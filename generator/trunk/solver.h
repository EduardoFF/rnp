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
#include "metric.h"

#ifdef _COMPILE_GUROBI

#include "gurobi_c++.h"

#endif


#ifdef _COMPILE_CPLEX

#include <cplex.h>
#include "cplex/cplex_consts.h"

#endif

#define VAR_RELAY 'y'  //Char used in the model to define binary variables for relay X -- i.e y[X]
#define VAR_EDGE 'x'   //Char used in the model to define flow variables in edges (X,Y) -- i.e x[X, Y]


class lp_solution_t{
	public:
	bool solved;
	double objval;
	map<string,double> value;
	map<string,char>  type;
	int n_relays;
	int n_edges;
	lp_solution_t(){}
	lp_solution_t(string filename);
	void write_to_file(string filename);

} ;

#ifdef _COMPILE_GUROBI

class mycallback: public GRBCallback
{
  public:

    GRBVar* vars;
    double *values;
    int numvars;
    int lastmsg;
    bool done;
    mycallback(GRBVar* xvars, double *xvalues, int xnumvars) {
      vars    = xvars;
      values = xvalues;
      numvars = xnumvars;
      lastmsg = -100;
      done = false;
    }
  protected:
      void callback ();
};


#endif
class Heuristic;

class GLPSol{
	private:

		string mps_filename;
		Network *net;

		bool has_mps;


		vector<string> relaysToUse;
		vector< pair< pair< string, string>, double> > edgesToUse;
		enum vartypes{VARTYPE_NULL, VARTYPE_RELAY, VARTYPE_EDGE};
		int get_var_type(const char *);
		string get_relay_id(const char *);
		void fix_glp_with_lpsol(glp_prob *, lp_solution_t *);
		pair< string, string> get_edge_from_var(const char*);
#ifdef _COMPILE_GUROBI
		lp_solution_t* solve_grb(glp_prob *);
		lp_solution_t* solve_grb(string, lp_solution_t*, grb_params_t &);
		lp_solution_t *solve_grb(glp_prob*, lp_solution_t*, grb_params_t &);
#endif
#ifdef _COMPILE_CPLEX
		lp_solution_t* solve_cplex(glp_prob *);
		lp_solution_t* solve_cplex(string, lp_solution_t*, cplex_params_t *);
		lp_solution_t *solve_cplex(glp_prob*, lp_solution_t*, cplex_params_t *);
#endif

		void solve_glp(glp_prob*);
		void write_to_sol(const char *, double);

	public:
		GLPSol(Network *n) :  net(n), has_mps(false){}
		GLPSol(Network *n, string fn) : mps_filename(fn),net(n), has_mps(true){}
		lp_solution_t *solve(int K, Heuristic *);
		static string create_data_file(Network *, int, Metric *);
		static void generate_and_save_mps(Network *net, string filename);
		Graph *create_graph_from_sol(lp_solution_t *);
	
};


#endif
