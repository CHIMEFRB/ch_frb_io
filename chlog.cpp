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

class chime_log_socket {
    typedef std::lock_guard<std::mutex> scoped_lock;

public:
    chime_log_socket() :
        _socket(NULL)
    {}

    ~chime_log_socket() {
        //cout << "~chime_log_socket" << endl;
        if (_socket)
            delete _socket;
    }

    void open_socket(zmq::context_t &ctx) {
        _socket = new zmq::socket_t(ctx, ZMQ_PUB);
    }

    void close_socket() {
        delete _socket;
        _socket = NULL;
    }
    
    void add_server(std::string port) {
        //cout << "chime_log connecting to " << port << endl;
        if (!_socket)
            throw runtime_error("chime_log_socket::add_server called by socket has not been initialized.");
        _socket->connect(port);
    }

    void send(std::string msg) {
        //cout << "chime_log sending message: " << msg << endl;
        if (!_socket)
            return;
        {
            scoped_lock l(_mutex);
            size_t n = _socket->send(static_cast<const void*>(msg.data()), msg.size());
        }
        //cout << "sent " << n << endl;
    }

protected:
    zmq::socket_t* _socket;

    std::mutex _mutex;
};

// GLOBALS for logging.
static chime_log_socket logsock;
static zmq::context_t* logctx = NULL;

void chime_log_quit() {
    logsock.close_socket();
}

// an atexit() function.
static void delete_logctx() {
    chime_log_quit();
    //cout << "delete_logctx()" << endl;
    if (logctx)
        delete logctx;
    logctx = NULL;
}

/*
 Initializes the CHIME/FRB distributed logging socket for this node.

 If a ZeroMQ context is passed in, it will be used to create the
 socket.  Otherwise, a new context will be created.

 NOTE, if a context is passed in, you likely will need to call
 chime_log_quit() to guarantee that the socket is closed before the
 ZeroMQ context is deleted.

 (When a context is not passed in, an atexit() handler is registered
 to clean up.)
 */
void chime_log_init(zmq::context_t* ctx = NULL) {
    if (!ctx) {
        if (!logctx) {
            logctx = new zmq::context_t();
            std::atexit(delete_logctx);
        }
        ctx = logctx;
    }
    logsock.open_socket(*ctx);
}

void chime_log_add_server(string port) {
    logsock.add_server(port);
}

static void chime_send_log(const string &s) {
    logsock.send(s);
}

// inspired by http://stackoverflow.com/questions/2212776/overload-handling-of-stdendl
class chime_log_stream : public std::ostream {
public:
    chime_log_stream()
        : std::ostream(&buffer),
          buffer()
    {}
protected:
    // A little string buffer subclass that intercepts the sync() call.
    class chstringbuf : public std::stringbuf {
    public:
        chstringbuf() {}

        virtual int sync() {
            chime_send_log(str());
            str("");
            return 0;
        }
    };
    chstringbuf buffer;
};

static chime_log_stream chstream;


ostream& chime_log(log_level lev, const char* file, int line, const char* function) {
    struct timeval tv;

    if (gettimeofday(&tv, NULL)) {
        cout << "Error calling gettimeofday(): " << strerror(errno) << endl;
        return chstream;
    }

    struct tm cal;
    if (!gmtime_r(&tv.tv_sec, &cal)) {
        cout << "Error calling gmtime_r(): " << strerror(errno) << endl;
        return chstream;
    }

    string datestring = stringprintf("%04i-%02i-%02i-%02i:%02i:%02i.%03i",
                                     cal.tm_year + 1900, cal.tm_mon + 1, cal.tm_mday,
                                     cal.tm_hour, cal.tm_min, cal.tm_sec,
                                     int(tv.tv_usec / 1000));

    chstream << "[" << file << ":" << line << " [" << function << "] "
           << datestring << "] ";
    return chstream;
}

void
__attribute__ ((format(printf,5,6)))
chime_logf(enum log_level lev, const char* file, int line, const char* function, const char* pattern, ...) {
    va_list lst;
    va_start(lst, pattern);
    string s = vstringprintf(pattern, lst);
    va_end(lst);
    chime_log(lev, file, line, function) << s << endl;
}


static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

std::mutex cout_mutex;

static void* server_main(zmq::context_t* ctx, string port) {
    /*
    void* varg) {
    std::tuple<zmq::context_t*, string> *args = reinterpret_cast<std::tuple<zmq::context_t*, string>* >(varg);
    zmq::context_t* ctx = std::get<0>(*args);
    string port = std::get<1>(*args);
     */

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

int main() {
    zmq::context_t ctx;
    //chime_log_init(&ctx);

    chime_log_init(NULL);

    string port = "tcp://127.0.0.1:6666";
    thread serverthread(std::bind(server_main, &ctx, port));
    serverthread.detach();
    chime_log_add_server(port);
    /*
     pthread_t serverthread;
     std::tuple<zmq::context_t*, string> args(&ctx, port);
     if (pthread_create(&serverthread, NULL, server_main, &args)) {
     cout << "Failed to create server thread." << endl;
     return -1;
     }
     */

    // string port2 = "tcp://127.0.0.1:6667";
    // thread serverthread2(std::bind(server_main, &ctx, port2));
    // chime_log_add_server(port2);

    usleep(1000000);

    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);

    chlogf("Gotta love %s style; %03i", "printf", 7);

    usleep(1000000);
    //ctx.close();

    //chime_log_quit();

    cout << "main() finished" << endl;
    return 0;
}



