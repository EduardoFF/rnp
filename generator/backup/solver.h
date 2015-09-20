/*
 * =====================================================================================
 *
 *       Filename:  solver.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 12:27:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef _SOLVER
#define _SOLVER

#include <glpk.h>
#define VAR_RELAY 'y'  //Char used in the model to define binary variables for relay X -- i.e y[rX]
#define VAR_EDGE 'x'   //Char used in the model to define flow variables in edges (X,Y) -- i.e x[X, Y]

class GLPSol{
	public:
		vector<string> relaysToUse;
		vector< pair< string, string> > edgesToUse;
		enum vartypes{VARTYPE_NULL, VARTYPE_RELAY, VARTYPE_EDGE};
		double tx_rg;
		GLPSol(double t) :  tx_rg(t){}
		Graph *solve(int K);
	//	 *getSolutionGraph();
		int get_var_type(const char *);
		string get_relay_id(const char *);
		pair< string, string> get_edge_from_var(const char*);

		void write_to_sol(const char *, double);
		int createDataFile(int);
};


#endif
