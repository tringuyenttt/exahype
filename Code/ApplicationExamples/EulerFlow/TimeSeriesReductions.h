#ifndef __TIMESERIES_EULER__
#define __TIMESERIES_EULER__

/* A class for writing ASCII Time Series, lazily in headers. */
#include <vector>
#include <array>
#include <stdio.h>
#include <errno.h>
#include <string.h>

using namespace std;

class TimeSeriesReductions {
    enum index { tidx=0, time=1, l1=2, l2=3, max=4, min=5, avg=6, LEN=7 };
    const char * const colnames[LEN]  = { "1:plotindex ", "2:time ", "3:l1norm ", "4:l2norm ", "5:max ", "6:min ", "7:avg " };
    const char * const colformat[LEN] = { "%.0f\t",        "%e\t",     "%e\t",       "%e\t",       "%e\t",    "%e\t",    "%e\t"  };
    double data[LEN];
    int avgcnt;
    FILE* asc;

public:
    TimeSeriesReductions(const char* filename) {
	asc = fopen(filename, "w");
	if(asc == NULL) {
		fprintf(stderr, "ASCII writer: Could not open output file '%s'\n", filename);//, strerror(errno));
		exit(-1);
	}
	
	// print the header
	fprintf(asc, "# ");
	for(int i=0; i<LEN; i++) {
            fprintf(asc, colnames[i]);
        }
        fprintf(asc, "\n");

        // zero everything
        for(int i=0; i<LEN; i++)
            data[i] = 0;
    }

    void initRow(double current_time) {
        data[time] = current_time;
        data[l1] = data[l2] = 0;
        data[max] = data[min] = 0;
        data[avg] = 0;
        avgcnt = 0;
    }

    void addValue(double val, double dx) {
        data[l1] += abs(val) * dx;
        data[l2] += val * val * dx;
        if(val*dx > data[max]) data[max] = val*dx;
        if(val*dx < data[min]) data[min] = val*dx;
        data[avg] += val;
        avgcnt++;
    }

    void finishRow() {
        data[tidx]++;
        data[l2] = sqrt(data[l2]);
        data[avg] = data[avg] / avgcnt;
    }

    void writeRow() {
        finishRow();

	for(int i=0; i<LEN; i++) {
            fprintf(asc, colformat[i], data[i]);
        }
        fprintf(asc, "\n");
    }
};


#endif /* __TIMESERIES_EULER__ */
