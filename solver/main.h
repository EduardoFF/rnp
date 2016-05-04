#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <assert.h>

#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <pthread.h>
#include <glog/logging.h>


using namespace std;


#define ALL(x) x.begin(), x.end()
#define PV(v,type) copy(ALL(v),ostream_iterator<type>(cout," "))
#define IT(c) __typeof((c).begin())
#define FOREACH(i,c) for(__typeof((c).begin()) i=(c).begin();i!=(c).end();++i)

#define VERBOSE(x) if(param.verbose >= x) 
//#define DEBUG(x) if(param.debug >= x) 

#define DBG(l, x) if (DEBUG(l)) cout << x;


#define CEIL(VARIABLE) ( (VARIABLE - (int)VARIABLE)==0 ? (int)VARIABLE : (int)VARIABLE+1 )
#define EPSILON 1e-4

#define GETDIR(x) (x.find_last_of("/")!=string::npos?x.substr(0,x.find_last_of("/")):".")

#define DEBUG(t,m) printf("%s: %s\n",__FUNCTION__,m)
#define CRITICAL 1

#define DEBUGP(t,m, ...) printf("%s: ",__FUNCTION__);printf(m,__VA_ARGS__);if(t == CRITICAL) exit(-1)



#endif
