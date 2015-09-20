/*
 * =====================================================================================
 *
 *       Filename:  main.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/15/2010 11:25:05 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef _MAIN
#define _MAIN

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
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sys/stat.h>
#include <sys/wait.h>
#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include "boost/graph/adjacency_list.hpp"
#include <boost/config.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>


#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
//#include <boost/geometry/geometries/adapted/tuple_cartesian.hpp>
//#include <boost/geometry/geometries/adapted/c_array_cartesian.hpp>
//#include <boost/geometry/geometries/register/box.hpp>
//#include <boost/geometry/geometries/adapted/std_as_linestring.hpp>

//#include <boost/geometry/geometries/cartesian2d.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <boost/xpressive/xpressive.hpp>



#include <boost/graph/graphviz.hpp>



using namespace std;
using namespace boost;


#define ALL(x) x.begin(), x.end()
#define PV(v,type) copy(ALL(v),ostream_iterator<type>(cout," "))
#define FOREACH(i,c) for(__typeof((c).begin()) i=(c).begin();i!=(c).end();++i)
#define RFOREACH(i,c) for(__typeof((c).rbegin()) i=(c).rbegin();i!=(c).rend();++i)
#define ITERATOR(c) __typeof((c).begin())

#define VERBOSE(x) if(prob.verbose >= x) 
#define DEBUG(x) if(prob.debug >= x) 

static ofstream cdebug;

#define CEIL(VARIABLE) ( (VARIABLE - (int)VARIABLE)==0 ? (int)VARIABLE : (int)VARIABLE+1 )
#define EPSILON 1e-4
#endif
