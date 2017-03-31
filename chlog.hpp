#ifndef _CH_FRB_LOG_HPP
#define _CH_FRB_LOG_HPP

#include <ostream>
#include <sstream>
#include <zmq.hpp>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

// printf into a string...
std::string 
__attribute__ ((format(printf,1,2)))
stringprintf(const char* format, ...);

enum log_level {
    log_level_debug = 1,
    log_level_info = 2,
    log_level_warn = 3,
    log_level_err = 4,
};

// Not meant to be called directly; called by preprocessor macro expansion
void chime_log(enum log_level lev, const char* file, int line, const char* function, const std::string &logmsg, bool do_assert=false);

// Not meant to be called directly; called by preprocessor macro expansion
void
__attribute__ ((format(printf,5,6)))
chime_logf(enum log_level lev, const char* file, int line, const char* function, const char* pattern, ...);

// This is the main logging method to use in CHIME/FRB code for
// distributed logging:
//
//  chlog("Hello world: " << 1 << ", " << 2 << ", " << 3);
//
#define chlog(...)                                                      \
    do {                                                                \
        std::ostringstream ss;                                          \
        ss << __VA_ARGS__;                                              \
        chime_log(ch_frb_io::log_level_info, __FILE__, __LINE__, __PRETTY_FUNCTION__, ss.str()); \
    } while(0)

// Like chlog, but takes a printf-style format string and arguments.
//
//  chlogf("Hello world: I am %03i", 7);
//
#define chlogf(pat, ...)                                                \
    do {                                                                \
        chime_logf(log_level_info, __FILE__, __LINE__, __PRETTY_FUNCTION__, pat, __VA_ARGS__); \
    } while(0)

// If we wanted to compile-out log messages below a certain log level,
// we could do:
#define chdebug(...) {}

// Assert that THING is true; log THING and bail out.
// FIXME -- bail out w/ assert() rather than exit(-1) ??
#define chassert(THING)                                                 \
    do { if (!(THING)) {                                                \
    chime_log(log_level_info, __FILE__, __LINE__, __PRETTY_FUNCTION__, "Assertion failed: " #THING, true); \
    usleep(1000000);                                                    \
    exit(-1);                                                           \
    } } while(0)

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
void chime_log_open_socket(zmq::context_t* ctx = NULL);

// Cleans up the distributed logging system.
void chime_log_close_socket();

// Enable/disable logging to cout as well as to the network.
void chime_log_local(bool);

// Sets a mnemonic name for this client; default is hostname
void chime_log_set_name(const std::string &name);

// Sets a mnemonic name for the current thread
void chime_log_set_thread_name(const std::string &name);

/*
 Starts sending log messages to the given server address (ZeroMQ
 address string, like "tcp://127.0.0.1:6667").
 */
void chime_log_add_server(const std::string &port);

void chime_log_remove_server(const std::string &port);



////////
// A simple log server that prints messages to cout.
////////
class chime_log_server {
public:
    chime_log_server(std::ostream& out = std::cout,
                     zmq::context_t* ctx = NULL,
                     const std::string &hostname = "*",
                     int port = -1);

    ~chime_log_server();

    std::string get_address();

    // Runs main serving loop.  Never returns.
    void run();

    // Starts a new thread to run this server.
    std::thread start();

    // Signals the server thread to stop.
    void stop();

protected:
    std::ostream& _out;
    zmq::socket_t* _socket;
    std::string _address;
    zmq::context_t* _ctx;
    bool _quit;
};



}

#endif
