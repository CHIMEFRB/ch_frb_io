#include <iostream>
#include <ostream>
#include <sstream>
#include <sys/time.h>
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

#if 0
class chime_log_stream : public std::ostringstream {

public:
    chime_log_stream() {
        cout << "chime_log_stream() constructor" << endl;
    }
    virtual ~chime_log_stream() {
        cout << "~chime_log_stream()" << endl;
        cout << "string: " << str() << endl;
    }

    /*
    virtual basic_ostream& flush() {
        cout << "chime_log_stream::flush()" << endl;
        return *this;
    }

    template<class T>
    chime_log_stream& operator<<(T val) {
        cout << "chime_log_stream:: operator<< " << val << endl;
        return *this;
    }
     */
};
#endif

// inspired by http://stackoverflow.com/questions/2212776/overload-handling-of-stdendl
class chime_log_stream : public std::ostream {
    class MyStreamBuf: public std::stringbuf {
        //std::ostream& output;
    public:
        MyStreamBuf()
        //std::ostream& str)
        // : output(str)
        {}

        virtual int sync() {
            cout << "[blah]" << str();
            str("");
            cout.flush();
            cout << "SYNC" << endl;
            return 0;
        }
    };
    MyStreamBuf buffer;
public:
    chime_log_stream() //std::ostream& str)
        : std::ostream(&buffer),
          buffer()
    {}
};



chime_log_stream chstream;


ostream& chime_log(log_level lev, const char* file, int line, const char* function) {
    struct timeval tv;

    if (gettimeofday(&tv, NULL)) {
        cout << "Error calling gettimeofday(): " << strerror(errno) << endl;
        //return cout;
    }

    struct tm cal;
    if (!gmtime_r(&tv.tv_sec, &cal)) {
        cout << "Error calling gmtime_r(): " << strerror(errno) << endl;
        //return cout;
    }

    /*
     string datestring = stringprintf("%04i-%02i-%02i-%02i:%02i:%02i.%03i",
     cal.tm_year + 1900, cal.tm_mon + 1, cal.tm_mday,
     cal.tm_hour, cal.tm_min, cal.tm_sec,
     int(tv.tv_usec / 1000));
     return cout << "[" << file << ":" << line << " [" << function << "] "
     << datestring << "] ";
     */

    string datestring = stringprintf("%04i-%02i-%02i-%02i:%02i:%02i.%03i",
                                     cal.tm_year + 1900, cal.tm_mon + 1, cal.tm_mday,
                                     cal.tm_hour, cal.tm_min, cal.tm_sec,
                                     int(tv.tv_usec / 1000));

    //chime_log_stream stream;
    //stringbuf stream;
    //ostringstream* stream = new chime_log_stream(); //ostringstream();
    chstream << "[" << file << ":" << line << " [" << function << "] "
           << datestring << "] ";
    return chstream;
}


int main() {
    //chlog();
    chlog("Hello world");
    chlog("Hello " << 1 << ", " << 2 << ", 3");
    chdebug("Debug" << 42+43);
}



