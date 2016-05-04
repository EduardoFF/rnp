///
///      @file  solverparameter.h
///     @brief  Parameters for solver
///
/// Basic parameters for solver library
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


#ifndef SOLVERPARAMETER_H
#define SOLVERPARAMETER_H


/** 
 * @enum SOLVER
 * @brief Supported solvers
   * More detailed enum description.
 */
typedef enum{
	SOLVER_GLPK=1, ///< GLPK 
	SOLVER_CPLEX=2, ///< CPLEX 
	SOLVER_GUROBI=3 ///< GUROBI
}SOLVER;





/** @brief Parameters for GUROBI solver
 */
typedef struct GrbParams{
	int cuts; ///< Gurobi cuts
} GrbParams;


/** @brief Parameters for CPLEX solver
 */



typedef struct CplexParams{
	char * name; 
	int cpx_index;
	char type;
	double db_value;
	int int_value;
	struct CplexParams* next;
} CplexParams;

// problem parameters (default value)
class SolverParams{
public:
	SolverParams(){
		//Default values
		seed = 12345;
		model_file = "model.mod";
		data_file = "model.dat";
		has_mps = false;
		mps_filename = "tmp.mps";	
		use_cplex = false;
		use_gurobi = true;
		which_solver = SOLVER_CPLEX;
		output_path=".";
		cplex_par = NULL;
		time_limit = 1e10;
		keep_mps = false;
		get_suboptimal = false;
		output_progress = 0;
//which_solver = SOLVER_GUROBI;
		verbose = 0;
		do_postsolve = 1;
		allow_parallel = 1;
		noMIP = 0;
		strong_timelim = 1;

		// Default gurobi
		grb.cuts = -1;

	}

     unsigned int seed; //Random seed (time(0))
     bool has_mps;
     string mps_filename;
     bool keep_mps;
     string model_file; // LP Model filename (model.mod)
     string data_file; // LP Data filename (model.dat)
     int verbose; // Verbose level (1)
     int debug; // Debug level (0)
     bool use_gurobi; // use gurobi solver instead of glpk (false)
     bool use_cplex; // use cplex solver (false)
     SOLVER which_solver;
     int glpk_out; //Display glpk output (0) 
     int grb_out; //Display gurobi output (0)
     string output_path; // Directory where to put output files (. )
     bool output_progress;
     string output_progress_file;
     GrbParams grb;
     CplexParams *cplex_par;
     double time_limit;
     bool do_postsolve;
     bool allow_parallel;
     bool noMIP;
     bool strong_timelim;

     bool get_suboptimal;
     //int use_heuristic;
     //string heuristic;
     string input_mps_file;
     string output_mps_file;
};

#endif
