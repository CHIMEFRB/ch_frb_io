#include <cstdlib>
#include <stdarg.h> // for va_start/va_end
#include <thread>
#include <pthread.h>
#include <functional>
#include <iostream>
#include <ostream>
#include <sstream>
#include <sys/time.h>
#include <unistd.h>
#include <zmq.hpp>
#include "chlog.hpp"

using namespace std;
using namespace ch_frb_io;

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


//////////////////////

std::mutex cout_mutex;

//////////////////////



static string 
vstringprintf(const char* format, va_list lst) {
    char temps[256];
    // truncates if length > size of 'temps'
    int n = vsnprintf(temps, sizeof(temps), format, lst);
    if (n < 0)
        throw runtime_error("vstringprintf failed: " + string(strerror(errno)));
    if (n < 256)
        return string(temps);
    // Try again with larger temp buffer
    char* temp2 = new char[n+1];
    n = vsnprintf(temp2, n+1, format, lst);
    if (n < 0)
        throw runtime_error("vstringprintf(2) failed: " + string(strerror(errno)));
    string s(temp2);
    delete[] temp2;
    return s;
}

static string 
__attribute__ ((format(printf,1,2)))
stringprintf(const char* format, ...) {
    va_list lst;
    va_start(lst, format);
    string s = vstringprintf(format, lst);
    va_end(lst);
    return s;
}


// Holds a ZeroMQ socket that is used to talk to multiple logging
// servers.
class chime_log_socket {
    typedef std::lock_guard<std::mutex> scoped_lock;

public:
    chime_log_socket() :
        _ctx(NULL),
        _socket(NULL),
        _mutex(),
        _stream(std::bind(&chime_log_socket::send,
                          this, std::placeholders::_1),
                std::bind(&chime_log_socket::unlock, this))
    {}

    ~chime_log_socket() {
        if (_socket)
            delete _socket;
        if (_ctx)
            delete _ctx;
    }

    void open_socket(zmq::context_t* ctx) {
        if (!ctx)
            ctx = _ctx = new zmq::context_t();
        _socket = new zmq::socket_t(*ctx, ZMQ_PUB);
    }

    void close_socket() {
        delete _socket;
        _socket = NULL;
    }
    
    void add_server(std::string port) {
        if (!_socket)
            throw runtime_error("chime_log_socket::add_server called but socket has not been initialized.");
        _socket->connect(port);
    }

    void remove_server(std::string port) {
        if (!_socket)
            throw runtime_error("chime_log_socket::remove_server called but socket has not been initialized.");
        _socket->disconnect(port);
    }

    void send(std::string msg) {
        if (_socket) {
            //cout << "send(" << msg << ")" << endl;
            //scoped_lock l(_mutex);
            _socket->send(static_cast<const void*>(msg.data()), msg.size());
        }
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::thread::id tid = std::this_thread::get_id();
            cout << "send(): pid/thread " << getpid() << "/" << tid << endl; //msg << endl;
        }
    }

    void unlock() {
        _mutex.unlock();
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::thread::id tid = std::this_thread::get_id();
            cout << "unlock(): pid/thread " << getpid() << "/" << tid << endl;
        }
    }


    chime_log_stream& get_stream() {
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::thread::id tid = std::this_thread::get_id();
            cout << "get_stream: pid/thread " << getpid() << "/" << tid << endl;
        }
        _mutex.lock();
        return _stream;
    }

protected:
    zmq::context_t* _ctx;
    zmq::socket_t* _socket;

    std::mutex _mutex;

    chime_log_stream _stream;
};


// GLOBALS for logging.
static chime_log_socket logsock;

void chime_log_quit() {
    logsock.close_socket();
}

void chime_log_init(zmq::context_t* ctx) {
    logsock.open_socket(ctx);
}

void chime_log_add_server(string port) {
    logsock.add_server(port);
}

void chime_log_remove_server(string port) {
    logsock.remove_server(port);
}

chime_log_stream& chime_log(log_level lev, const char* file, int line, const char* function) {
    struct timeval tv;
    string datestring;
    if (gettimeofday(&tv, NULL)) {
        cout << "Error calling gettimeofday(): " << strerror(errno) << endl;
    } else {
        struct tm cal;
        if (!gmtime_r(&tv.tv_sec, &cal)) {
            cout << "Error calling gmtime_r(): " << strerror(errno) << endl;
        } else {
            datestring = stringprintf("%04i-%02i-%02i-%02i:%02i:%02i.%03i",
                                      cal.tm_year + 1900, cal.tm_mon + 1, cal.tm_mday,
                                      cal.tm_hour, cal.tm_min, cal.tm_sec,
                                      int(tv.tv_usec / 1000));
        }
    }
    string levstring = (lev == log_level_debug ? "DEBUG" :
                        (lev == log_level_info ? "INFO" :
                         (lev == log_level_warn ? "WARN" :
                          (lev == log_level_err ? "ERROR" : "XXX"))));
    string header = stringprintf("%s %s:%i [%s] %s ", levstring.c_str(), file, line, function, datestring.c_str());

    chime_log_stream& s = logsock.get_stream();
    s << header;
    return s;
}

void
__attribute__ ((format(printf,5,6)))
chime_logf(enum log_level lev, const char* file, int line, const char* function, const char* pattern, ...) {
    va_list lst;
    va_start(lst, pattern);
    string s = vstringprintf(pattern, lst);
    va_end(lst);
    chime_log_stream& ll = chime_log(lev, file, line, function);
    ll << s;
    ll.done();
}


} // namespace

//////////////////// Part of the main() demo -- a log server. ////////////////////////

static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

//std::mutex cout_mutex;

static void* server_main(zmq::context_t* ctx, string port) {
    zmq::socket_t sock(*ctx, ZMQ_SUB);
    sock.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    sock.bind(port);
    cout << "server_main: bound " << port << endl;
    for (;;) {
        zmq::message_t msg;
        //cout << "server_main: waiting to receive message." << endl;
        try {
            if (!sock.recv(&msg)) {
                cout << "log server: failed to receive message" << endl;
                break;
            }
        } catch (const zmq::error_t& e) {
            cout << "log server: error receiving message: " << e.what() << endl;
            break;
        }

        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            cout << "Server received message: " << msg_string(msg) << endl;
        }
    }
    cout << "server_main done." << endl;
    return NULL;
}

void log_client(int num) {
    for (int i=0; i<10; i++) {
        chlog("Hello I am client " << num << ", message " << i);
    }
}

int main() {
    zmq::context_t ctx;

    chime_log_init(&ctx);
    //chime_log_init(NULL);

    string port = "tcp://127.0.0.1:6666";
    thread serverthread(std::bind(server_main, &ctx, port));
    serverthread.detach();
    chime_log_add_server(port);

    // string port2 = "tcp://127.0.0.1:6667";
    // thread serverthread2(std::bind(server_main, &ctx, port2));
    // chime_log_add_server(port2);

    usleep(1000000);

    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);

    chlogf("Gotta love %s style; %03i", "printf", 7);

    usleep(1000000);

    chime_log_quit();

    chime_log_init(NULL);
    chime_log_add_server(port);

    usleep(1000000);

    thread logger1(std::bind(log_client, 1));
    thread logger2(std::bind(log_client, 2));

    logger1.join();
    logger2.join();

    usleep(3000000);

    cout << "main() finished" << endl;
    return 0;
}



