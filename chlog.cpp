#include <cstdlib>
#include <map>
#include <stdarg.h> // for va_start/va_end
#include <thread>
#include <mutex>
#include <functional>
#include <iostream>
#include <ostream>
#include <sstream>
#include <sys/time.h>
#include <unistd.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <zmq.hpp>
#include "chlog.hpp"

using namespace std;
using namespace ch_frb_io;

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

// in misc.cpp
std::string ip_to_string(const sockaddr_in &addr);

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

    void set_name(const std::string &name) {
        scoped_lock l(_mutex);
        _name = name;
    }

    void set_thread_name(const std::string &name) {
        scoped_lock l(_mutex);
        _threadnames[std::this_thread::get_id()] = name;
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
    
    void add_server(const std::string &port) {
        scoped_lock l(_mutex);
        if (!_socket)
            throw runtime_error("chime_log_socket::add_server called but socket has not been initialized.");
        _socket->connect(port.c_str());
    }

    void remove_server(const std::string &port) {
        scoped_lock l(_mutex);
        if (!_socket)
            throw runtime_error("chime_log_socket::remove_server called but socket has not been initialized.");
        _socket->disconnect(port.c_str());
    }

    void send(std::string header, std::string msg, bool do_assert) {
        scoped_lock l(_mutex);
        string tname = get_thread_name();
        if (_local) {
            if (do_assert)
                // Assume we want to know where it was...
                cout << header << endl;
            cout << "[" << tname << "] " << msg << endl;
        }
        if (!_socket)
            return;
        msg = _name + " " + tname + " " + header + " " + msg;
        _socket->send(static_cast<const void*>(msg.data()), msg.size());
    }

protected:
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

void chime_log_set_name(const std::string &name) {
    logsock.set_name(name);
}

void chime_log_set_thread_name(const std::string &name) {
    logsock.set_thread_name(name);
}

void chime_log_add_server(const std::string &port) {
    logsock.add_server(port);
}

void chime_log_remove_server(const std::string &port) {
    logsock.remove_server(port);
}

void chime_log(log_level lev, const char* file, int line, const char* function,
               const std::string &msg, bool do_assert) {
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
                     stringprintf(" %s:%i [%s] ", file, line, function) +
                     datestring);
    logsock.send(header, msg, do_assert);
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



chime_log_server::chime_log_server(std::ostream& out,
                                   zmq::context_t* ctx,
                                   const std::string &hostname,
                                   int port) :
    _out(out),
    _quit(false)
{
    if (!ctx)
        ctx = _ctx = new zmq::context_t();
    _socket = new zmq::socket_t(*ctx, ZMQ_SUB);
    _socket->setsockopt(ZMQ_SUBSCRIBE, "", 0);

    string addr = "tcp://" + hostname + ":";
    if (port == -1)
        addr = addr + "*";
    else
        addr = addr + std::to_string(port);

    _socket->bind(addr.c_str());

    char addrx[256];
    size_t addrsz = 256;
    _socket->getsockopt(ZMQ_LAST_ENDPOINT,
                        reinterpret_cast<void*>(addrx), &addrsz);
    _address = string(addrx);
    //cout << "Bound to address " << _address << endl;

    // When binding "*", the returned address is like "tcp://0.0.0.0:6534".
    // Replace "0.0.0.0" by my IP address
    size_t i = _address.find("0.0.0.0");
    if (i != std::string::npos) {
        string ipaddress = "0.0.0.0";
        char name[256];
        if (gethostname(name, 256)) {
            cout << "Failed to gethostname(): " << strerror(errno) << endl;
        } else {
            //cout << "Hostname: " << name << endl;
            struct addrinfo* addrs = NULL;
            if (getaddrinfo(name, NULL, NULL, &addrs)) {
                cout << "Failed to getaddrinfo(): " << gai_strerror(errno) << endl;
            } else {
                struct addrinfo* a;
                for(a = addrs; a != NULL; a = a->ai_next) {
                    if (a->ai_family != AF_INET)
                        continue;
                    struct sockaddr_in* sin = reinterpret_cast<struct sockaddr_in*>(a->ai_addr);
                    string ipstr = ip_to_string(*sin);
                    if (ipstr.size()) {
                        ipaddress = ipstr;
                        break;
                    }
                }
                freeaddrinfo(addrs);
            }
        }
        _address.replace(i, 7, ipaddress);
        //cout << "Replaced address: " << _address << endl;
    }
}

chime_log_server::~chime_log_server() {
    delete _socket;
    if (_ctx) {
        delete _ctx;
    }
}

std::string chime_log_server::get_address() {
    return _address;
}

static string msg_string(zmq::message_t &msg) {
    return string(static_cast<const char*>(msg.data()), msg.size());
}

void chime_log_server::run() {

    void* p_sock = _socket->operator void*();
    zmq_pollitem_t pollitems[] = {
        { p_sock, 0, ZMQ_POLLIN, 0 },
    };

    for (;;) {
        zmq::message_t msg;
        try {

            int r = zmq::poll(pollitems, 1, 1000);
            if (r == -1) {
                cout << "log server: zmq::poll error: " << strerror(errno) << endl;
                break;
            }

            if (_quit) {
                cout << "log server: _quit!" << endl;
                break;
            }

            if (!(pollitems[0].revents & ZMQ_POLLIN)) {
                cout << "log server: no input ready" << endl;
                continue;
            }

            if (!_socket->recv(&msg)) {
                _out << "log server: failed to receive message" << endl;
                break;
            }
        } catch (const zmq::error_t& e) {
            _out << "log server: error receiving message: " << e.what() << endl;
            break;
        }
        _out << msg_string(msg) << endl;
    }
    cout << "log server: exiting" << endl;
}

// Starts a new thread to run this server.
std::thread chime_log_server::start() {
    thread t(std::bind(&chime_log_server::run, this));
    t.detach();
    return t;
}

void chime_log_server::stop() {
    _quit = true;
}


} // namespace

