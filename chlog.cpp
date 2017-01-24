#include <cstdlib>
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
__attribute__ ((format(printf,1,2)))
stringprintf(const char* format, ...) {
    va_list lst;
    char temps[256];
    va_start(lst, format);
    // truncates if length > size of 'temps'
    if (vsnprintf(temps, sizeof(temps), format, lst) < 0)
        throw runtime_error("stringprintf failed: " + string(strerror(errno)));
    va_end(lst);
    return string(temps);
}

class chime_log_socket {
public:
     chime_log_socket(zmq::context_t &ctx) :
         _socket(ctx, ZMQ_PUB)
    {}

    ~chime_log_socket() {
        cout << "~chime_log_socket" << endl;
    }

    void add_server(std::string port) {
        cout << "chime_log connecting to " << port << endl;
        _socket.connect(port);
    }

    void send(std::string msg) {
        //cout << "chime_log sending message: " << msg << endl;
        size_t n = _socket.send(static_cast<const void*>(msg.data()), msg.size());
        //cout << "sent " << n << endl;
    }

protected:
    zmq::socket_t _socket;
    
};

static chime_log_socket* logsock = NULL;

void chime_log_quit() {
    cout << "chime_log_quit" << endl;
    if (logsock)
        delete logsock;
    logsock = NULL;
}

void chime_log_init(zmq::context_t &ctx) {
    if (logsock)
        delete logsock;
    logsock = new chime_log_socket(ctx);
    std::atexit(chime_log_quit);
}

void chime_log_add_server(string port) {
    if (!logsock)
        return;
    logsock->add_server(port);
}

static void chime_send_log(const string &s) {
    if (logsock) {
        logsock->send(s);
    } else {
        cout << s;
        cout.flush();
        cout << "SEND" << endl;
    }
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


static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

std::mutex cout_mutex;

static void* server_main(void* varg) {
    std::tuple<zmq::context_t*, string> *args = reinterpret_cast<std::tuple<zmq::context_t*, string>* >(varg);
    zmq::context_t* ctx = std::get<0>(*args);
    string port = std::get<1>(*args);

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

    chime_log_init(ctx);

    /*
    string port = "tcp://127.0.0.1:6666";
    //thread serverthread(std::bind(server_main, &ctx, port));
    pthread_t serverthread;
    std::tuple<zmq::context_t*, string> args(&ctx, port);
    if (pthread_create(&serverthread, NULL, server_main, &args)) {
        cout << "Failed to create server thread." << endl;
        return -1;
    }
    chime_log_add_server(port);
    //serverthread.detach();
     */

    // string port2 = "tcp://127.0.0.1:6667";
    // thread serverthread2(std::bind(server_main, &ctx, port2));
    // chime_log_add_server(port2);

    usleep(1000000);

    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);

    usleep(1000000);
    //ctx.close();

    chime_log_quit();

    cout << "main() finished" << endl;
    return 0;
}



