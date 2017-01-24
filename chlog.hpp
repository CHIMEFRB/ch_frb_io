#include <ostream>

enum log_level {
    log_level_debug = 1,
    log_level_info = 2,
    log_level_warn = 3,
    log_level_err = 4,
};

std::ostream& chime_log(enum log_level lev, const char* file, int line, const char* function);

#define chlog(...) do { chime_log(log_level_info, __FILE__, __LINE__, __PRETTY_FUNCTION__) << __VA_ARGS__ << endl;} while(0)

#define chdebug(...) {}

