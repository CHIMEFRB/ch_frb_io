#include <cstdlib>
#include <map>
#include <stdarg.h> // for va_start/va_end
#include <thread>
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

std::string 
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
        _local(true),
        _name(),
        _threadnames()
    {}

    ~chime_log_socket() {
        if (_socket)
            delete _socket;
        if (_ctx)
            delete _ctx;
    }

    void set_local(bool loc) {
        scoped_lock l(_mutex);
        _local = loc;
    }

    void set_name(std::string name) {
        _name = name;
    }

    void set_thread_name(std::string name) {
        _threadnames[std::this_thread::get_id()] = name;
    }

    std::string get_thread_name() {
        std::thread::id tid = std::this_thread::get_id();
        auto it = _threadnames.find(tid);
        if (it == _threadnames.end()) {
            std::ostringstream ss;
            ss << "Thread-" << tid;
            return ss.str();
        }
        return it->second;
    }

    void open_socket(zmq::context_t* ctx) {
        scoped_lock l(_mutex);
        if (_socket)
            throw runtime_error("chime_log_socket::open_socket called but socket is already open.");
        if (!ctx)
            ctx = _ctx = new zmq::context_t();
        _socket = new zmq::socket_t(*ctx, ZMQ_PUB);

        if (!_name.length()) {
            char tmp[256];
            gethostname(tmp, 256);
            _name = tmp;
            cout << "Set hostname: " << _name << endl;
        }
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

            string tname = get_thread_name();
            if (_local)
                cout << tname << " " << msg << endl;
            if (!_socket)
                return;
            msg = _name + " " + tname + " " + msg;
            _socket->send(static_cast<const void*>(msg.data()), msg.size());
        }
    }

protected:
    zmq::context_t* _ctx;
    zmq::socket_t* _socket;

    // Mutex used to protect the socket (and also _socket/_ctx state)
    std::mutex _mutex;

    // Print to cout also?
    bool _local;

    // Client name to send
    std::string _name;

    // Thread name map
    std::map<std::thread::id, std::string> _threadnames;
};


// GLOBALS for logging.
static chime_log_socket logsock;

void chime_log_open_socket(zmq::context_t* ctx) {
    logsock.open_socket(ctx);
}

void chime_log_close_socket() {
    logsock.close_socket();
}

void chime_log_local(bool loc) {
    logsock.set_local(loc);
}

void chime_log_set_name(std::string name) {
    logsock.set_name(name);
}

void chime_log_set_thread_name(std::string name) {
    logsock.set_thread_name(name);
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
    string header = (levstring +
                     stringprintf(" %s:%i [%s]", file, line, function) +
                     datestring);
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

