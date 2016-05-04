///
///      @file  solver.h
///     @brief  MILP Solver for generic problems
///
/// Support GLPK, CPLEX and GUROBI
/// Based on Matteo Salani work on glpk_wrapper
///
///    @author  Eduardo Feo (), eduardo@idsia.ch
///
///  @internal
///    Created  06/10/2011
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  IDSIA (Dalle Molle Institute for Artificial Intelligence)
///  Copyright  Copyright (c) 2011, Eduardo Feo
///
/// This program is free software; you can redistribute it and/or modify it 
/// under the terms of the GNU General Public License as published by the 
/// Free Software Foundation; either version 2 of the License, or (at your option) 
/// any later version. This program is distributed in the hope that it will be useful, 
/// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
/// or FITNESS FOR A PARTICULAR PURPOSE. 
/// See the GNU General Public License for more details at  
/// http://www.gnu.org/copyleft/gpl.html

///=====================================================================================
///


#ifndef SOLVER_H
#define SOLVER_H


#include "CpuTime.h"
#include "solverparameter.h"

#ifdef __WITH_GLPK
#include <glpk.h>
#endif

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/assume_abstract.hpp>


#ifdef _COMPILE_GUROBI
#include "gurobi_c++.h"
#endif


#ifdef _COMPILE_CPLEX
#include <ilcplex/cplex.h>
#include "cplex_consts.h"
#endif






/**
 * @class LpSolution
 * @brief Represents the solution of a MILP problem
 *
 * MORE COMMENTS COME HERE
 * 
 */
class LpSolution {

	// Serialization
	friend class boost::serialization::access;
	template<class Archive>
		void serialize(Archive & ar, const unsigned int /* file_version */);

	double m_cpuTimeModelReading;
	double m_cpuTimeDataReading;
	double m_cpuTimeModelGeneration;
	double m_cpuTimeSolving;
	double m_cpuTimePostSolve;
	double m_cpuTimeTotal;
	bool m_is_optimal;
	public:
	bool solved;
	double objval;
	map<string,double> value;
	map<string,char>  type;
	//int n_relays;
	//int n_edges;
	LpSolution() { 
		m_cpuTimeModelReading=0;
		m_cpuTimeDataReading = 0;
		m_cpuTimeModelGeneration = 0;
		m_cpuTimeSolving = 0;
		m_cpuTimePostSolve = 0;
		m_cpuTimeTotal = 0;
		m_is_optimal = false;
		solved=false;
		objval = 1000000000;
	}
	LpSolution(string filename);
	void print();
	void write_to_file(string filename) const;
	double cpuTimeModelReading() const { return m_cpuTimeModelReading;}
	double cpuTimeDataReading() const { return m_cpuTimeDataReading;}
	double cpuTimeModelGeneration() const { return m_cpuTimeModelGeneration;}
	double cpuTimeSolving() const { return m_cpuTimeSolving;}
	double cpuTimePostSolve() const { return m_cpuTimePostSolve;}
	double cpuTimeTotal() const { return m_cpuTimeTotal;}

	bool isOptimal() const { return m_is_optimal;}

	void is_optimal(bool b){ m_is_optimal = b;}
	void cpuTimeModelReading(double t) { m_cpuTimeModelReading = t;}
	void cpuTimeDataReading(double t) { m_cpuTimeDataReading = t;}
	void cpuTimeModelGeneration(double t) { m_cpuTimeModelGeneration = t;}
	void cpuTimeSolving(double t) { m_cpuTimeSolving = t;}
	void cpuTimePostSolve(double t) { m_cpuTimePostSolve = t;}
	void cpuTimeTotal(double t) { m_cpuTimeTotal = t;}



} ;

void save_lpsolution(const LpSolution &p, string &str);
void save_lpsolution(const LpSolution &p, std::vector<char> & , int&);
void restore_lpsolution(LpSolution &p, string &str);
void restore_lpsolution(LpSolution &p, char *, int );


#ifdef __WITH_GLPK

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

typedef struct {
	double timenow;
	bool has_incumbent;
	double bestobj;
	double gap;
} CallbackInfo;



typedef void(*WriteToDataFileFunc)(ofstream &); 
typedef int(*GLPSolSimpleCallback)(void *);


#ifdef _COMPILE_CPLEX
typedef std::map<int, std::pair<string, int> > CPX_VARIABLESINFO;

struct timeliminfo {
    double timestart;
    double timelim;
    double acceptablegap;
    int    aborted;
    int numvars;
    bool output_progress;
    GLPSolSimpleCallback simplecallback;
    CPX_VARIABLESINFO *varsinfo;
    LpSolution *incumbent;
};
typedef struct timeliminfo CPX_TIMELIMINFO, *CPX_TIMELIMINFOptr;

struct loginfo {
   double lastincumbent;
   int    lastlog;
   int    numcols;
};
typedef struct loginfo CPX_LOGINFO, *CPX_LOGINFOptr;

static int CPXPUBLIC
   timelimcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                    void *cbhandle);


#endif



class Heuristic;



class GLPSol{
	private:
		SolverParams param;
		string mps_filename;
		CpuTime cpu_time;
		double m_cpuTimeModelReading;
		double m_cpuTimeDataReading;
		double m_cpuTimeModelGeneration;
		double m_cpuTimeSolving;
		double m_cpuTimePostSolve;
		double m_cpuTimeTotal;
		GLPSolSimpleCallback m_simpleCallback;
		bool m_useSimpleCallback;





		//vector<string> relaysToUse;
		//vector< pair< pair< string, string>, double> > edgesToUse;
		//enum vartypes{VARTYPE_NULL, VARTYPE_RELAY, VARTYPE_EDGE};
		//

		int get_var_type(const char *);
		//string get_relay_id(const char *);
		//
		//
		void fix_glp_with_lpsol(glp_prob *, LpSolution *);
		//pair< string, string> get_edge_from_var(const char*);
#ifdef _COMPILE_GUROBI
		LpSolution* solve_grb(glp_prob *);
		LpSolution* solve_grb(string, LpSolution*, GrbParams &);
		LpSolution *solve_grb(glp_prob*, LpSolution*, GrbParams &);
#endif
#ifdef _COMPILE_CPLEX
		LpSolution* solve_cplex(glp_prob *);
		LpSolution* solve_cplex(string, LpSolution*, CplexParams *);
		LpSolution *solve_cplex(glp_prob*, LpSolution*, CplexParams *);
#endif

		void solve_glp(glp_prob*);
		void write_to_sol(const char *, double);

	public:
		GLPSol(SolverParams &param_) : param(param_) { m_simpleCallback=NULL; m_useSimpleCallback = false;}
		GLPSol(SolverParams &param_, string fn) : param(param_) { m_simpleCallback = NULL; m_useSimpleCallback = false;}
		LpSolution *solve(WriteToDataFileFunc, Heuristic *h = NULL );
//		static string create_data_file(Problem *);
		string create_data_file(WriteToDataFileFunc f);
		void generate_and_save_mps(WriteToDataFileFunc f, string filename);
		~GLPSol() {}
		void   setSimpleCallback(GLPSolSimpleCallback f) {m_simpleCallback=f; m_useSimpleCallback=(f!=NULL);}
		double cpuTimeModelReading() { return m_cpuTimeModelReading;}
		double cpuTimeDataReading() { return m_cpuTimeDataReading;}
		double cpuTimeModelGeneration() { return m_cpuTimeModelGeneration;}
		double cpuTimeSolving() { return m_cpuTimeSolving;}
		double cpuTimePostSolve() { return m_cpuTimePostSolve;}
		double cpuTimeTotal() { return m_cpuTimeTotal;}

	
};


#endif

#endif
