#ifndef TIMER_H_AUG_17_2005
#define TIMER_H_AUG_17_2005

#include <ostream>
using namespace std;
#include <sys/times.h>
#include <unistd.h>

class Timer {
public:
    Timer() {
        times(&starttime_);
        ticksPerSecond_ = (int) sysconf(_SC_CLK_TCK);
    }
    ~Timer() {}

    void reset() {
        times(&starttime_);
    }

    // Returns the number of seconds of user time that have elapsed
    double time() {
        struct tms endtime;
        times(&endtime);
        clock_t diff = endtime.tms_utime - starttime_.tms_utime;
        double seconds = ((double)diff) / ticksPerSecond_;
        return seconds;
    }

    static void printTime(ostream &out, const double &sec) {
        if (sec < 60) {
            out << sec << " secs";
        } else if (sec < 3600)  out << int(sec / 60) << " mins, " << sec - int(sec / 60) * 60
                                        << " secs";
        else if (sec < 86400) out <<  int(sec / 3600) << " hrs, "
                                      << (sec - int(sec / 3600) * 3600) / 60.0 << " mins";
        else             out << int(sec / 86400) << " days, "
                                 << (sec - int(sec / 86400) * 86400) / 3600.0 << " hrs";
    }

private:
    struct tms starttime_;
    int ticksPerSecond_;
};


#endif
