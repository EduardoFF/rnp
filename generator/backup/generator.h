/*
 * =====================================================================================
 *
 *       Filename:  generator.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/01/2010 01:08:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _GENERATOR
#define _GENERATOR

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

#include "boost/lexical_cast.hpp"
#include <boost/random.hpp>
#include "boost/graph/adjacency_list.hpp"
#include <boost/config.hpp>

#include <boost/graph/graphviz.hpp>

#include "GetPot.h"


using namespace std;
using namespace boost;

#define ALL(x) x.begin(), x.end()
#define PV(v,type) copy(ALL(v),ostream_iterator<type>(cout," "))

#define VERBOSE(x) if(prob.verbose > x) 
#define DEBUG(x) if(prob.debug > x) 

#define GRAPHVIZ_CMD(x) "neato -Tps -s " #x "> "


// problem parameters (default value)
typedef struct{
     int st_nodes; // Number of static nodes (10)
     int b_nodes; // Number of base station nodes (2)
     int K; //Maximum number of relays to use (+INF)
     double dimX, dimY; //Width and length of field (10 x 10)
     double tx_range; //Transmission range (1.0)

     unsigned int seed; //Random seed (time(0))

     double grid_res; //Grid resolution (1.0)
     string model_file; // LP Model filename (model.mod)
     string data_file; // LP Data filename (model.dat)

     int verbose; // Verbose level (1)

     int debug; // Debug level (0)
     string instance_id;  // id to identify problem instance (default)
     bool use_gurobi; // use gurobi solver instead of glpk (false)
     int glpk_out; //Display glpk output (0) 
     int grb_out; //Display gurobi output (0)


} param;


     

#endif
