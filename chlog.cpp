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
        _local(true)
    {}

    ~chime_log_socket() {
        if (_socket)
            delete _socket;
        if (_ctx)
            delete _ctx;
    }

    void chime_log_local(bool loc) {
        scoped_lock l(_mutex);
        _local = loc;
    }

    void open_socket(zmq::context_t* ctx) {
        scoped_lock l(_mutex);
        if (_socket)
            throw runtime_error("chime_log_socket::open_socket called but socket is already open.");
        if (!ctx)
            ctx = _ctx = new zmq::context_t();
        _socket = new zmq::socket_t(*ctx, ZMQ_PUB);

        // Send buffer: ZMQ_SNDHWM: default 1000 messages
    }

    void close_socket() {
        scoped_lock l(_mutex);
        if (_socket)
            delete _socket;
        _socket = NULL;
    }
    
    void add_server(std::string port) {
        scoped_lock l(_mutex);
        if (!_socket)
            throw runtime_error("chime_log_socket::add_server called but socket has not been initialized.");
        _socket->connect(port);
    }

    void remove_server(std::string port) {
        scoped_lock l(_mutex);
        if (!_socket)
            throw runtime_error("chime_log_socket::remove_server called but socket has not been initialized.");
        _socket->disconnect(port);
    }

    void send(std::string header, std::string msg) {
        msg = header + " " + msg;
        {
            scoped_lock l(_mutex);
            if (_local)
                cout << msg << endl;
            if (!_socket)
                return;
            _socket->send(static_cast<const void*>(msg.data()), msg.size());
        }
    }

protected:
    zmq::context_t* _ctx;
    zmq::socket_t* _socket;

    // Mutex used to protect the socket (and also _socket/_ctx state)
    std::mutex _mutex;

    bool _local;
};


// GLOBALS for logging.
static chime_log_socket logsock;

void chime_log_close_socket() {
    logsock.close_socket();
}

void chime_log_open_socket(zmq::context_t* ctx) {
    logsock.open_socket(ctx);
}

void chime_log_add_server(string port) {
    logsock.add_server(port);
}

void chime_log_remove_server(string port) {
    logsock.remove_server(port);
}

void chime_log(log_level lev, const char* file, int line, const char* function,
               string msg) {
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
    string header = levstring + stringprintf(" %s:%i [%s]", file, line, function) + datestring;

    logsock.send(header, msg);
}

void
__attribute__ ((format(printf,5,6)))
chime_logf(enum log_level lev, const char* file, int line, const char* function, const char* pattern, ...) {
    va_list lst;
    va_start(lst, pattern);
    string s = vstringprintf(pattern, lst);
    va_end(lst);
    chime_log(lev, file, line, function, s);
}


} // namespace

//////////////////// Part of the main() demo -- a log server. ////////////////////////

std::mutex cout_mutex;

static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

static void* server_main(zmq::context_t* ctx, string port, int sid) {
    zmq::socket_t sock(*ctx, ZMQ_SUB);
    sock.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    sock.bind(port);
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << "server_main " << sid << ": bound " << port << endl;
    }
    for (;;) {
        zmq::message_t msg;
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
            cout << "Server " << sid << ": " << msg_string(msg) << endl;
        }
    }
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << "server_main " << sid << " done." << endl;
    }
    return NULL;
}

void log_client(int num) {
    for (int i=0; i<10; i++) {
        chlog("Hello I am client " << num << ", message " << i);
    }
}

int main() {
    zmq::context_t ctx;

    // Try opening the socket with a given context.
    chime_log_open_socket(&ctx);

    string port = "tcp://127.0.0.1:6666";
    thread serverthread(std::bind(server_main, &ctx, port, 1));
    serverthread.detach();
    chime_log_add_server(port);

    string port2 = "tcp://127.0.0.1:6667";
    thread serverthread2(std::bind(server_main, &ctx, port2, 2));
    serverthread2.detach();
    chime_log_add_server(port2);


    // What if we connect to a server that doesn't exist?
    string port3 = "tcp://127.0.0.1:6668";
    chime_log_add_server(port3);


    usleep(1000000);

    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);

    chlogf("Gotta love %s style; %03i", "printf", 7);

    usleep(100000);

    chime_log_close_socket();

    // Now re-open with a new context
    chime_log_open_socket(NULL);

    chime_log_add_server(port);
    chime_log_add_server(port2);

    // (wait for connect)
    usleep(1000000);

    thread logger1(std::bind(log_client, 1));
    thread logger2(std::bind(log_client, 2));

    logger1.join();
    logger2.join();

    usleep(1000000);

    cout << "main() finished" << endl;
    return 0;
}
