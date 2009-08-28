const double Undef = -999.123;
//const double DERR = -999.999;

#include <string>
using namespace std;

void reportError();

void checkNArgs(const int n, const int nrhs);

template<class A>
inline int getInt(A* mxhndl) {
    return int(mxGetScalar(mxhndl));
}

template<class A>
inline double getDouble(A* mxhndl) {
    return double(mxGetScalar(mxhndl));
}

inline char* getString(const mxArray* p) {
    char* input_buf = 0;
    int status;
    int m = mxGetM(p);
    int n = mxGetN(p);
    int buflen = m*n + 1;
    string msg;

    if (m == 1) {
        input_buf = (char*)mxCalloc(buflen, sizeof(char));
        status = mxGetString(p, input_buf, buflen);
        if(status != 0) {
            msg = string(input_buf) 
                  + "\nNot enough space. String is truncated.";
            mexWarnMsgTxt(msg.c_str());
        }
    }
    else {
        mexErrMsgTxt("string must be a row vector");
    }
    return input_buf;
}



