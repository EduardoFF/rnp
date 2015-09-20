/*
 * =====================================================================================
 *
 *       Filename:  solver.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 11:51:31 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eduardo Feo
 *        Company:  
 *
 * =====================================================================================
 */

#include "main.h"
#include "generator.h"
#include "Plotter.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "solver.h"
#include "CpuTime.h"
#include "heuristic.h"

#ifdef _COMPILE_GUROBI
extern "C" {
#include "../wrapper/wrapper.h"
}
#endif



extern param prob;


extern ofstream cdebug;

#define MAX_LINK_NEIGHBORS 10.0
#define CPX_INT 0
#define CPX_DB 1

void exit_on_error(string msg, glp_prob *_P, glp_tran *_tran)
{
	fprintf(stderr, "msg.c_str()");

	if(_tran != NULL)
		glp_mpl_free_wksp(_tran);

	glp_delete_prob(_P);

	exit(-1);



}

void freeGLPK(glp_prob *_P, glp_tran *_tran)
{
	if(_tran != NULL)
		glp_mpl_free_wksp(_tran);

	glp_delete_prob(_P);
}

string var_status_to_string(glp_prob *_P, int col)
{
	switch(glp_get_col_type(_P, col))
	{
	case(GLP_FR):
		return ("free (unbounded) variable");
	case(GLP_LO):
		return("variable with lower bound");
	case(GLP_UP):
		return("variable with upper bound");
	case(GLP_DB):
		return("double-bounded variable");
	case(GLP_FX):
		return("fixed variable");
	}
	return("Unknown variable status");
}


string var_kind_to_string(glp_prob *_P, int col)
{
	// page 70
	switch(glp_get_col_kind(_P, col))
	{
	case(GLP_CV):
		return ("continuous");
	case(GLP_IV):
		return("integer");
	case(GLP_BV):
		return("binary");
	}
	return("Unknown variable type");
}



string problem_type_to_string(glp_prob *_P)
{
	// page 34
	switch(glp_get_obj_dir(_P))
	{
	case(GLP_MIN):
		return ("minimization");
	case(GLP_MAX):
		return("maximization");
	}
	return("Unknown variable type");
}


void lp_solution_t::write_to_file(string filename){
	ofstream outf(filename.c_str());
	
	cout << "Writing solution to " << filename << endl;

	outf << "SOLVED\t" << solved << endl;
	if(solved){
		outf << "OBJ_VAL\t" << objval << endl;
		outf << "N_VARS\t" << value.size() << endl;
		outf << "MAX_RELAYS\t" << prob.K << endl;
		outf << "N_RELAYS\t" << n_relays << endl;
		outf << "N_EDGES\t" << n_edges << endl;
		FOREACH(it,value){
			string var = it->first;
			outf << var << "\t" << type[var] << "\t" << value[var] << endl;
		}	
	}
	outf.close();
}


lp_solution_t::lp_solution_t(string filename){
	int nelem;
	ifstream ifile(filename.c_str());

	cout << "Loading solution from " << filename << endl;

	ifile >> solved;
	ifile >> objval;
	ifile >> nelem;
	for(int i=0;i<nelem;i++){
		string varid;
		double val;
		char typ;
		ifile >> varid >> typ >> val;
		value[varid] = val;
		type[varid] = typ;
	}
}

string GLPSol::create_data_file(Network *net, int K, Metric *metric){


	VERBOSE(1)printf("Creating data file...\n");
	string data_filename = prob.output_path + "/" + prob.instance_id + "-model.dat";
	ofstream dataFile(data_filename.c_str());
	if( dataFile.is_open()){

		/*  Creating sets */
		dataFile << "set S :=";

		for(int i = 0; i < net->numStatic(); i++){
			dataFile << " s" << i;

		}
		dataFile << ";\n";
		dataFile << "set B :=";

		for(int i=0; i< net->numBases();i++){
			dataFile << " b" << i;
		}
		dataFile << ";\n";

		dataFile << "set R :=";
		for(int i=0;i<net->numRelays();i++){
			Node &n = net->getRelay(i);
			dataFile << " " << n.id;
		}
		dataFile << ";\n";

		/* Creating edges */
		dataFile << "set edges :=";

		for(int i=0;i < net->size(); i++){
			Node &ni = net->getNode(i);
			for( int j=i+1;j<net->size();j++){
				Node &nj = net->getNode(j);
				if( net->inRange(ni,nj)){
					dataFile << " ( " << ni.id << " , " << nj.id << ")";
					dataFile << " ( " << nj.id << " , " << ni.id << ")";
				}


			}	
		}
		dataFile << " ;\n";


		
		/* Setting costs for edges */
		
		/*
		dataFile << "param d_f := ";
		for(int i=0;i < net->size(); i++){
			Node &ni = net->getNode(i);
			for( int j=i+1;j<net->size();j++){
				Node &nj = net->getNode(j);
				if( net->inRange(ni,nj)){
					double link_cost;



					double link_d = net->distanceLink(ni,nj);
				//	double link_n = net->neighborsLink(ni,nj);

					double d_f = (max(link_d,2.5)/net->tx_range); 
				//	double n_f = (min(link_n,MAX_LINK_NEIGHBORS) / MAX_LINK_NEIGHBORS);

					//link_cost = pow(link_cost, 2.0);

					// Calculate cost based on edge length 
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << d_f << " "; 
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << d_f << " "; 


				}


			}	
		}
		dataFile << ";\n";
		
*/
		dataFile << "param c := ";
		for(int i=0;i < net->size(); i++){
			Node &ni = net->getNode(i);
			for( int j=i+1;j<net->size();j++){
				Node &nj = net->getNode(j);
				if( net->inRange(ni,nj)){
					double link_cost;
					link_cost = (*metric)(ni.id,nj.id);
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << link_cost << " ";
					link_cost = (*metric)(nj.id,ni.id);
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << link_cost << " ";
				}
			}
		}
		dataFile << ";\n";
		dataFile << "param AVERAGE_LINK_COST := " << metric->average();
		dataFile << ";\n";
		dataFile << "param MAX_LINK_COST := " << metric->max();
		// cout << "param MAX_LINK_COST := " << metric->max() << endl;


		/** Calculate expected hop-count */
		dataFile << ";\n";
		int expected_hop_count = 0;
		for(int i = 0; i< net->numStatic();i++){
			Node &ns = net->getStatic(i);
			int hop_count = 100000000;
			for(int j=0;j < net->numBases();j++){
				Node &nb = net->getBase(j);
				double d = net->distanceLink(ns,nb);
				int hc = (int)ceil(1.0*d/net->getRange());
				hop_count = min(hc,hop_count);
			}
			expected_hop_count += hop_count;
		}
		dataFile << "param EXPECTED_HOP_COUNT := " << expected_hop_count;

		// cout << "param EXPECTED_HOP_COUNT := " << expected_hop_count << endl;
		


		
		
	
/*		
		dataFile << "param n_f := ";
		for(int i=0;i < net->size(); i++){
			Node &ni = net->getNode(i);
			for( int j=i+1;j<net->size();j++){
				Node &nj = net->getNode(j);
				if( net->inRange(ni,nj)){
					double link_cost;



				//	double link_d = net->distanceLink(ni,nj);
					double link_n = net->neighborsLink(ni,nj);

				//	double d_f = (max(link_d,2.5)/net->tx_range); 
					double n_f = (min(link_n,MAX_LINK_NEIGHBORS) / MAX_LINK_NEIGHBORS);

					//link_cost = pow(link_cost, 2.0);

					// Calculate cost based on edge length
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << n_f << " "; 
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << n_f << " "; 


				}


			}	
		}
*/
		dataFile << ";\n";
	

		dataFile << "param K := " << K << ";\n";
		dataFile << "end;\n";






		return data_filename;

	}else
		return "";




}

int GLPSol::get_var_type(const char *v){
	if(v[0] == VAR_RELAY)
		return VARTYPE_RELAY;
	else if(v[0] == VAR_EDGE)
		return VARTYPE_EDGE;
	else
		return VARTYPE_NULL;
}

string GLPSol::get_relay_id(const char *v){
	string var(v);
	return var.substr(var.find('[')+1,var.find(']')-var.find('[')-1);
}

pair<string, string> GLPSol::get_edge_from_var(const char* varname){
	string v(varname);
	string r1 = v.substr(v.find('[')+1,v.find(',') - v.find('[') -1);
	string r2 = v.substr(v.find(',')+1, v.find(']') - v.find(',') - 1 );
	return make_pair(r1,r2);
}

void GLPSol::write_to_sol(const char *varname, double value){
	/*  Just add relays that belong to an edge */	
	/*   if(get_var_type(varname) == VARTYPE_RELAY){
		string rid = get_relay_id(varname);
		DEBUG(1)cout << "Relayid " << rid << endl;
		relaysToUse.push_back(rid);
	}else 
	*/
	if(get_var_type(varname) == VARTYPE_EDGE){
		pair<string,string> e = get_edge_from_var(varname);
		DEBUG(1)cdebug << "Edge " << e.first << ", " << e.second << endl;
		edgesToUse.push_back(make_pair(e,value));
		if(e.first[0] == 'r')
			relaysToUse.push_back(e.first);
		if(e.second[0] == 'r')
			relaysToUse.push_back(e.second);
	}
	/*else
		printf("GLPSol: Invalid variable type\n");*/



}
void GLPSol::fix_glp_with_lpsol(glp_prob *mip, lp_solution_t *sol){
	double tmp;
	char type;
	int col_index;
	string nameVar;

	int numvars;

	int verbose = 0;

	/** GUROBI: Get variable values **/
	FOREACH(it,sol->value){

		tmp = it->second;

		nameVar = it->first;

		type = sol->type[nameVar];

		/** GLPK search variable index by name **/
		col_index = glp_find_col(mip, nameVar.c_str());

		if (col_index != 0){
			/** GLPK set variable bounds **/
			if ((type == 'B') || (type == 'I')){
				//if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", 
				//		     nameVar.c_str(), type, tmp, round(tmp));
				glp_set_col_bnds(mip, col_index, GLP_FX, round(tmp), round(tmp));
			}
			else{
				//if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", 
				//		     nameVar.c_str(), type, tmp, tmp);
				glp_set_col_bnds(mip, col_index, GLP_FX, tmp, tmp);
			}
		}
	}



}

#ifdef _COMPILE_GUROBI

void freeGRB(GRBModel *model, GRBEnv  *env){
	/** GUROBI: free structures **/
	if (model) delete model;
	if (env) delete env;
}


void mycallback::callback(){


      try {
        if (where == GRB_CB_MIPNODE && !done ) {
		setSolution(vars,values,numvars);
		done = true;

		

	}else if (where == GRB_CB_MIPSOL) {
          double obj     = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
          int    nodecnt = (int) getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
          double* x = getSolution(vars, numvars);
          cout << "**** New solution at node " << nodecnt << ", obj "
                             << obj <<  "****" << endl;
          delete[] x;


	}
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }

}
	


lp_solution_t *GLPSol::solve_grb(string mps_file, lp_solution_t *prev_sol, grb_params_t &params){
	int GLPK_out = 1;
	int GRB_out = 1;
	double obj_val;
	double dblattr;

	GRBVar *sol_vars;
	double *sol_values;


	char type;
	int col_index;
	double tmp,bound;
	
	lp_solution_t *ret_sol = new lp_solution_t();


	int glp_status, grb_status, numvars, ret, GRB_IsMIP;
	string nameVar;
	int verbose = 1;


	/************/
	/** GUROBI **/
	/************/

	GRBEnv *env;
	try {
		env = new GRBEnv();



		env->set(GRB_IntParam_OutputFlag,1);

		/** GUROBI: Read model **/
		GRBModel model =  GRBModel(*env,mps_file.c_str());
		/** GUROBI: Get environment **/
		GRBEnv mipenv = model.getEnv();


		/** GUROBI: Set parameters **/

		/** GUROBI: Ask for more precision **/
		mipenv.set(GRB_DoubleParam_FeasibilityTol,10E-6);
		mipenv.set(GRB_DoubleParam_IntFeasTol,10E-5);
		mipenv.set(GRB_DoubleParam_MIPGap,10E-6);


		/** Other parameters */
		mipenv.set(GRB_IntParam_Cuts,params.cuts);
		mipenv.set(GRB_IntParam_RootMethod,1);
		mipenv.set(GRB_IntParam_Symmetry,-1);


		/** GUROBI: get numvars and numrows **/
		numvars = model.get(GRB_IntAttr_NumVars);


		/** GUROBI: get model type **/
		GRB_IsMIP = model.get(GRB_IntAttr_IsMIP);

		/** Load previous solution when avaliable 
		 *  To do that - use callbacks 
		 *  */

		mycallback *cb;

		if(prev_sol){
			DEBUG(1) cdebug << "Registering callback" << endl;

			sol_vars = model.getVars();
			sol_values = (double *)malloc(numvars*sizeof(double));
			for(int i=0;i<numvars;i++){
				nameVar = sol_vars[i].get(GRB_StringAttr_VarName);
				sol_values[i] = prev_sol->value[nameVar];
				if(sol_values[i] > EPSILON)
					DEBUG(1)cdebug << "setting " << nameVar << " = " << sol_values[i] << endl;
			}
			cb = new mycallback(sol_vars, sol_values, numvars);
			model.setCallback(cb);
		}

		/** GUROBI: Optimize model **/
		model.optimize();

		/** GUROBI: Retreive the optimization status **/
		grb_status = model.get(GRB_IntAttr_Status);
		switch(grb_status){
		case GRB_OPTIMAL:
			break;
		case GRB_INFEASIBLE :
			fprintf(stderr, "Error GRB optimization failed with code GRB_INFEASIBLE\n");
		case GRB_INF_OR_UNBD :
			fprintf(stderr, "Error GRB optimization failed with code GRB_INF_OR_UNBD \n");
		case GRB_UNBOUNDED :
			fprintf(stderr, "Error GRB optimization failed with code GRB_UNBOUNDED \n");
		case GRB_CUTOFF :
			fprintf(stderr, "Error GRB optimization failed with code GRB_CUTOFF \n");
		case GRB_ITERATION_LIMIT :
			fprintf(stderr, "Error GRB optimization failed with code GRB_ITERATION_LIMIT \n");
		case GRB_NODE_LIMIT :
			fprintf(stderr, "Error GRB optimization failed with code GRB_NODE_LIMIT \n");
		case GRB_TIME_LIMIT :
			fprintf(stderr, "Error GRB optimization failed with code GRB_TIME_LIMIT \n");
		case GRB_SOLUTION_LIMIT :
			fprintf(stderr, "Error GRB optimization failed with code GRB_SOLUTION_LIMIT \n");
		case GRB_INTERRUPTED :
			fprintf(stderr, "Error GRB optimization failed with code GRB_INTERRUPTED \n");
		case GRB_SUBOPTIMAL :
			fprintf(stderr, "Error GRB optimization failed with code GRB_SUBOPTIMAL \n");
		case GRB_NUMERIC :
			fprintf(stderr, "Error GRB optimization failed with code GRB_NUMERIC \n");

			/** GUROBI: Quit in any case non optimal **/
			//freeGRB(env);
		}

		/** GUROBI: Get obj function value **/
		tmp = model.get(GRB_DoubleAttr_IntVio);
		bound = model.get(GRB_DoubleAttr_ObjBound);
		tmp = model.get(GRB_DoubleAttr_ObjVal);

		/* ********************** */

		obj_val = tmp;


		/* ************ */
		if (verbose) printf ("Objective %lf\n", tmp);
		if (verbose) printf ("Best bound %lf\n", bound);
		if (verbose) printf ("Absolute gap %lf\n", fabs(tmp - bound));

		/** Write solution to object lp_solution_t */

		ret_sol->objval = obj_val;
		ret_sol->solved = true;
		GRBVar *vars = model.getVars();
		for(int i=0;i<numvars;i++){
			nameVar = vars[i].get(GRB_StringAttr_VarName);
			double tmpdbl = vars[i].get(GRB_DoubleAttr_X);
			if(tmpdbl > EPSILON){
				ret_sol->value[nameVar] = vars[i].get(GRB_DoubleAttr_X);
				ret_sol->type[nameVar] = vars[i].get(GRB_CharAttr_VType);
				//if(fabs(ret_sol->value[nameVar]) > 0.0)
				DEBUG(1) cdebug << "Writing sol: " << nameVar 
					<< " has value " << ret_sol->value[nameVar] << endl;
			}
		}
		delete vars;



	}catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Error during optimization" << endl;
	}


	/** TO DO: Clean everything */
	delete env;


	return ret_sol;



}


lp_solution_t *GLPSol::solve_grb(glp_prob *mip, lp_solution_t *prev_sol, grb_params_t &params){
		/** GLPK: Generate Variable indexing **/
	glp_create_index(mip);

	/** GLPK: Generate LP **/

	glp_write_mps(mip, GLP_MPS_FILE, NULL, "tmp.mps");

	lp_solution_t *retval = solve_grb("tmp.mps",prev_sol,params);
	remove("tmp.mps");
	return retval;




}

#endif


#ifdef _COMPILE_CPLEX

void free_cplex(char * context, int value, CPXENVptr cpx_env,CPXLPptr cpx_lp){
	//cpx_udef * tmp;


	/** GLPK: free user defined parameters **/
	/* 
	while (cpx_user_def_params != NULL){
		tmp = cpx_user_def_params;
		cpx_user_def_params = cpx_user_def_params->next;
		free(tmp);
	}
*/
	/** CPLEX: free structures **/
	if (cpx_lp) CPXfreeprob(cpx_env, &cpx_lp);
	if (cpx_env) CPXcloseCPLEX(&cpx_env);

	printf("ERRORS OCCURRED.\nContext: %s, Value: %d.\nGLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\n", context, value);

	
}

lp_solution_t *GLPSol::solve_cplex(string mps_file, lp_solution_t *prev_sol, cplex_params_t *params){
	/** CPLEX: data structures **/
	CPXENVptr cpx_env = NULL;
	CPXLPptr cpx_lp = NULL;
	/** CPLEX: variables */
	int      numvars, j, retCPX, CPX_IsMIP, CPX_noMIP, CPX_probType,cpx_probe;
	double   db_tmp, bound,obj_val;
	char nameCPX[80];
	char type;
	double cpx_feasibilityTol = 10E-6;
	double cpx_intFeasTol = 10E-4;
	double cpx_mipGap = 10E-5;
	/** execution flags **/
	char CPX_out, GLPK_out, verbose;
  	CPX_out = GLPK_out = CPX_noMIP = verbose = 0;
	double cpx_timelim =3600; // Time limit one hour
	
	CPX_out = 1;
	verbose = 1;
	
	string nameVar;
	lp_solution_t *ret_sol = new lp_solution_t();


	/** CPLEX: user defined parameters **/

	cpx_env = CPXopenCPLEX(&retCPX);
	if (retCPX || cpx_env == NULL)
	{
		fprintf(stderr, "Error: could not create cpx_environment\n");
		free_cplex("CPX create environment", retCPX, cpx_env, cpx_lp);
	}

	retCPX = CPXsetintparam(cpx_env, CPX_PARAM_SCRIND, CPX_out?CPX_ON:CPX_OFF);
	if (retCPX) free_cplex("CPX set output", retCPX, cpx_env, cpx_lp);
	
	retCPX = CPXsetintparam(cpx_env, CPX_PARAM_PROBE, 3);
	if (retCPX) free_cplex("CPX set probe", retCPX, cpx_env, cpx_lp);

	/** CPLEX: Create LP model **/
	cpx_lp = CPXcreateprob(cpx_env, &retCPX, "cpx_lp");
	if (retCPX) free_cplex("CPX create problem", retCPX, cpx_env, cpx_lp);

	/** CPLEX: Read model **/
	retCPX = CPXreadcopyprob (cpx_env, cpx_lp, mps_file.c_str(), "MPS");
	if (retCPX) free_cplex("CPX load model", retCPX, cpx_env, cpx_lp);

	/** Remove utility files from disk **/
//	if (!keep_tmp_mps) remove("model.mps");
	/** CPLEX: Set parameters **/

	/** CPLEX: Ask for more precision **/
	retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPRHS , cpx_feasibilityTol);
	if (retCPX) free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
	retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPINT, cpx_intFeasTol);
	if (retCPX) free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
	retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPGAP, cpx_mipGap);
	if (retCPX) free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
	retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_TILIM, cpx_timelim);
	if (retCPX) free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);


	/** Set other user-defined parameters for Cplex **/
	cplex_params *tmp = params;
	while (tmp != NULL){
		if (tmp->type == CPX_INT){
			retCPX = CPXsetintparam(cpx_env, tmp->cpx_index, tmp->int_value);
			if (retCPX) free_cplex("CPX Set udef int para", retCPX, cpx_env, cpx_lp);
		}
		else if (tmp->type == CPX_DB){
			retCPX = CPXsetdblparam(cpx_env, tmp->cpx_index, tmp->db_value);
			if (retCPX) free_cplex("CPX Set udef dbl para", retCPX, cpx_env, cpx_lp);
		}
		tmp = tmp->next;
	}

	/** CPLEX: get numvars **/
	numvars = CPXgetnumcols(cpx_env, cpx_lp);

	/** CPLEX: get model type **/
	CPX_probType = CPXgetprobtype(cpx_env, cpx_lp);
	if ((CPX_probType != CPXPROB_MILP) && (CPX_probType != CPXPROB_LP))
	{
		fprintf(stderr, "Error: unable to solve problem type :%d\n", CPX_probType);
		free_cplex("CPX problem type", 0, cpx_env, cpx_lp);
	}
	if (verbose) printf ("Problem type: %d\n", CPX_probType);
	CPX_IsMIP = (CPXPROB_MILP == CPX_probType);

	/** CPLEX: Solve the linear relaxation **/
	if ((CPX_noMIP) || (!CPX_IsMIP)){
		//TODO allow for other simplex optimizers
		retCPX = CPXdualopt(cpx_env, cpx_lp);
		if (retCPX) free_cplex("CPX LP Optimize", retCPX, cpx_env, cpx_lp);

		/** CPLEX: Retreive the optimization status **/
		retCPX = CPXgetstat(cpx_env, cpx_lp);
		if (retCPX != CPX_STAT_OPTIMAL){
			/** CPLEX: Quit in any case non optimal **/
			free_cplex("Error CPX Lp optimization failed", retCPX, cpx_env, cpx_lp);
			ret_sol->solved = false;
			return ret_sol;
		}
	}
	else{
		retCPX = CPXmipopt(cpx_env, cpx_lp);
		if (retCPX) free_cplex("CPX MIP Optimize", retCPX, cpx_env, cpx_lp);

		/** CPLEX: Retreive the optimization status **/
		retCPX = CPXgetstat(cpx_env, cpx_lp);
		if ((retCPX != CPXMIP_OPTIMAL) && (retCPX != CPXMIP_OPTIMAL_TOL)){
			/** CPLEX: Quit in any case non optimal **/
			free_cplex("Error CPX MIP optimization failed", retCPX, cpx_env, cpx_lp);
			ret_sol->solved = false;
			return ret_sol;
		}
	}
	/** CPLEX: Get obj function value **/
	if (!CPX_noMIP){
		retCPX = CPXgetbestobjval(cpx_env, cpx_lp, &bound);
		if (retCPX) free_cplex("CPX accessing bound", retCPX, cpx_env, cpx_lp);
	}

	retCPX = CPXgetobjval(cpx_env, cpx_lp, &db_tmp);
	if (retCPX) free_cplex("CPX obj value", retCPX, cpx_env, cpx_lp);

	if (verbose) printf ("Objective %lf\n", db_tmp);
	if (verbose) printf ("Best bound %lf\n", bound);
	if (verbose) printf ("Absolute gap %lf\n", fabs(db_tmp - bound));
	
	obj_val = db_tmp;
	ret_sol->objval = obj_val;
	ret_sol->solved = true;
	/** CPLEX: Get variable values **/
	for (j = 0; j < numvars; ++j){


		retCPX = CPXgetx (cpx_env, cpx_lp, &db_tmp, j, j);
		if (retCPX) free_cplex("CPX Get var value", retCPX, cpx_env, cpx_lp);


		int surplus, colnamespace; 
		char ** colname = NULL;
		char * colnamestore = NULL;

		retCPX = CPXgetcolname(cpx_env, cpx_lp, NULL, NULL, 0, &surplus, j, j);
		if (( retCPX != CPXERR_NEGATIVE_SURPLUS ) && ( retCPX != 0 ))  {
			free_cplex("CPX Get var names", retCPX, cpx_env, cpx_lp);
		}

		colnamespace = - surplus;
		if ( colnamespace > 0 ) {

			colname = (char **) malloc (sizeof(char *));
			colnamestore = (char *)  malloc(colnamespace);

			if ( colname == NULL || colnamestore == NULL ) {
				free_cplex("CPX Allocating memory for col name" , 0, cpx_env, cpx_lp);
			}

			retCPX = CPXgetcolname (cpx_env, cpx_lp, colname, colnamestore, colnamespace, &surplus, j, j);
			if ( retCPX ) {
				free_cplex("CPX Get final var names", retCPX, cpx_env, cpx_lp);
			}
		}
		else {
			free_cplex("CPX no name associated", 0, cpx_env, cpx_lp);
		}

//		if (verbose) printf ("Processed variable %d name %s\n", j, colname[0]);

		sprintf(nameCPX, "%s", colname[0]);

		free(colnamestore);
		free(colname);

		char type[1];
		retCPX = CPXgetctype (cpx_env, cpx_lp, type, j,j);
		if (retCPX) free_cplex("CPX Accessing variable type", retCPX, cpx_env, cpx_lp);

		nameVar = nameCPX;	
	

		double tmpdbl = db_tmp;
		if(tmpdbl > EPSILON){
			ret_sol->value[nameVar] = tmpdbl;
			ret_sol->type[nameVar] = type[0];
			//if(fabs(ret_sol->value[nameVar]) > 0.0)
			DEBUG(1) cdebug << "Writing sol: " << nameVar 
				<< " has value " << ret_sol->value[nameVar] << endl;
		}

	}

	/** TO DO: Clean everything
	 * */
	return ret_sol;




}

lp_solution_t *GLPSol::solve_cplex(glp_prob *mip, lp_solution_t *prev_sol, cplex_params_t *params){
		/** GLPK: Generate Variable indexing **/
	glp_create_index(mip);

	/** GLPK: Generate LP **/

	glp_write_mps(mip, GLP_MPS_FILE, NULL, "tmp.mps");

	lp_solution_t *retval = solve_cplex("tmp.mps",prev_sol,params);
	remove("tmp.mps");
	return retval;
}




#endif
void GLPSol::solve_glp(glp_prob *mip){

	glp_iocp *iparm;	
	double glpk_iparm_mip_gap = 10E-4;
	double glpk_iparm_tol_int = 10E-4;
	double glpk_iparm_tol_obj = 10E-4;


	int ret,glp_status;
	
	/** GLPK initialize parameters **/
	iparm = (glp_iocp*) malloc(sizeof(glp_iocp));
	glp_init_iocp(iparm);
	iparm->presolve = GLP_ON;
	iparm->mip_gap = glpk_iparm_mip_gap;
	iparm->tol_int = glpk_iparm_tol_int;
	iparm->tol_obj = glpk_iparm_tol_obj;

	/** GLPK get the optimal integer solution **/
	ret = glp_intopt(mip, iparm);
	if (ret){
		fprintf(stderr, "glp_intopt, Error on optimizing the model : %d \n", ret);
		//freeMem();
	}

	glp_status = glp_mip_status(mip);
	switch (glp_status){
	case GLP_OPT:
		break;
	case GLP_FEAS:
		fprintf(stderr, "Error GLPK simplex is not optimal, GLP_FEAS, code %d\n", ret);
	case GLP_NOFEAS:
		fprintf(stderr, "Error GLPK simplex is not optimal, GLP_NOFEAS, code %d\n", ret);
	case GLP_UNDEF:
		fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNDEF, code %d\n", ret);
	}



}


Graph *GLPSol::create_graph_from_sol(lp_solution_t *sol){


/**  Commented: Done at GLPSOl::solve 
	FOREACH(it,sol->value){
		if(fabs(it->second) > EPSILON)
		write_to_sol((it->first).c_str(),it->second);
	}
	


	std::sort(relaysToUse.begin(), relaysToUse.end());
	std::vector<std::string>::iterator new_end_pos;
	new_end_pos = std::unique( relaysToUse.begin(), relaysToUse.end() );

	// The elements between new_end_pos and str.end() hold
	// hold their old values, and must be erased to complete
	// the operation:

	relaysToUse.erase( new_end_pos, relaysToUse.end() );
*/
	VERBOSE(2){
		printf("Relays to use\n");
		PV(relaysToUse,string);
		printf("\n");
		printf("Edges to use\n");
		for(int i=0;i<edgesToUse.size();i++)
			printf(" (%s %s) %f",edgesToUse[i].first.first.c_str(), edgesToUse[i].first.second.c_str(),
			       edgesToUse[i].second);
		printf("\n");
	}

	printf("Building solution graph...\n");
	Graph *g = new Graph(net->numStatic()+net->numBases()+relaysToUse.size());
	map<string,int> idToIdx;
	for(int i=0;i<net->numStatic();i++){
		g->set_vertex_id(i,net->getStatic(i).id);
		g->set_vertex_loc(i,net->getStatic(i).x / net->dimX ,net->getStatic(i).y / net->dimY);
		g->set_vertex_type(i,net->getStatic(i).t);
		idToIdx[net->getStatic(i).id] = i;
	}
	for(int i=0;i<net->numBases();i++){
		int ix = i + net->numStatic();
		g->set_vertex_id(ix,net->getBase(i).id);
		g->set_vertex_loc(ix,net->getBase(i).x / net->dimX ,net->getBase(i).y / net->dimY);
		g->set_vertex_type(ix,net->getBase(i).t);
		idToIdx[net->getBase(i).id] = i + net->numStatic();
	}
	for(int i=0;i<relaysToUse.size();i++){
		int ix = i + net->numStatic()+net->numBases();
		Node &n = net->getNodeById(relaysToUse[i]);
		g->set_vertex_id(ix,n.id);
		g->set_vertex_loc(ix,n.x / net->dimX,n.y/net->dimY);
		g->set_vertex_type(ix,n.t);
		idToIdx[relaysToUse[i]] = ix;
	}

	for(int i=0;i<edgesToUse.size();i++){
		int u,v;
		u = idToIdx[edgesToUse[i].first.first];
		v = idToIdx[edgesToUse[i].first.second];

		g->addEdge(u,v);
		stringstream ss;

		ss << edgesToUse[i].second;
		g->set_edge_label(u,v,ss.str());
	}
	return g;

}
lp_solution_t *GLPSol::solve(int max_relay, Heuristic *H){


	glp_prob *lp;
	glp_tran *tran;
	glp_iocp params;
	int ret;

	lp_solution_t *ret_sol = new lp_solution_t();

	CpuTime cpu_time;
	string out_time_file; 
	double timer;
	out_time_file = prob.output_path + "/" + prob.instance_id + ".sol.time";
	FILE *log_fp = fopen(out_time_file.c_str(),"w");


	double obj_val;



	cpu_time.start();  // start timer
	
	lp = glp_create_prob();
	tran = glp_mpl_alloc_wksp();

	if(!has_mps){
		/** Generate problem from network structure */

		ret = glp_mpl_read_model(tran, prob.model_file.c_str(), 1);
		if (ret != 0) 
			printf("Error on translating model\n");
		else
			printf("Model read: OK\n");
		timer = cpu_time.cpu_time_elapsed(log_fp, "MODEL_READING");
		printf("time elapsed: MODEL_READING  %g\n",timer);

		Metric *metric;
		printf("metric is %s\n",prob.metric.c_str());
		if( prob.metric == "PRR")
			metric = new EstPRR_Metric(net);
		else if(prob.metric == "DISTANCE_SQ")
			metric = new DistSQ_Metric(net);
		
		/* * Time recording */
		timer = cpu_time.cpu_time_elapsed(log_fp, "METRIC_CALC");
		printf("time elapsed: METRIC_CALC  %g\n",timer);
		
		string dataFile = GLPSol::create_data_file(net,max_relay,metric);
		if(dataFile == "")
			printf("Error creating data file\n");

		ret = glp_mpl_read_data(tran, dataFile.c_str());
		if (ret != 0) 
			printf("Error on translating data\n");
		else
			printf("Data read: OK\n");


		/* * Time recording */
		timer = cpu_time.cpu_time_elapsed(log_fp, "DATA READING");
		printf("time elapsed: DATA READING  %g\n",timer);

		if (glp_mpl_generate(tran, NULL) != 0)
			exit_on_error("Error on generating model\n", lp, tran);
		glp_mpl_build_prob(tran, lp);
		/* * Time recording */
		timer = cpu_time.cpu_time_elapsed(log_fp,"GENERATION");
		printf("time elapsed: GENERATION  %g\n",timer);
	}else{
		/** Load mps file */
		printf("Loading mps file %s\n",mps_filename.c_str());
		glp_read_mps(lp,  GLP_MPS_FILE, NULL,mps_filename.c_str());
		glp_create_index(lp);
		/* * Time recording */
		timer = cpu_time.cpu_time_elapsed(log_fp,"LOAD_MPS");
		printf("time elapsed: LOAD_MPS  %g\n",timer);

	}




	


	if(prob.use_gurobi){
#ifdef _COMPILE_GUROBI
		wrapper_params wp;
		wp.glp_out = prob.glpk_out;


		wp.grb_out = prob.grb_out;

		lp_solution_t *prev_sol = NULL;
		if(H){
			prev_sol = (*H)();
			int cnt = 0;
			FOREACH(jt,prev_sol->value){
				if(jt->second > 0)
					cnt++;
			}
			DEBUG(1)cdebug << "prev sol has " << cnt << " variables > 0" << endl;
		}
		lp_solution_t *sol;
		if(has_mps)
			sol = solve_grb(mps_filename,prev_sol,prob.grb);
		else
			sol = solve_grb(lp,prev_sol,prob.grb);

		ret_sol = sol;
		
		fix_glp_with_lpsol(lp,sol);
		obj_val = sol->objval;
#else
		printf("GUROBI SUPPORT NOT COMPILED\n");
#endif
		//obj_val = solve_glp_grb(lp, &wp);
	}else if(prob.use_cplex){
#ifdef _COMPILE_CPLEX
		printf("********* USING CPLEX ***********\n");
		lp_solution_t *prev_sol = NULL;
		
		lp_solution_t *sol;
		if(has_mps)
			sol = solve_cplex(mps_filename,prev_sol,prob.cplex_par);
		else
			sol = solve_cplex(lp,prev_sol,prob.cplex_par);

		ret_sol = sol;
		if(sol->solved){
			printf("********************* CPLEX ENDED - FIXING GLPK *******************\n");
			fix_glp_with_lpsol(lp,sol);
			obj_val = sol->objval;
		}else{
			printf("********************* CPLEX ENDED - NO SOLUTION FOUND *******************\n");
		
		}
#else
		printf("CPLEX SUPPORT NOT COMPILED\n");
#endif
			
	}else{

		glp_simplex(lp, NULL);

		obj_val = glp_get_obj_val(lp);


		bool feasible;
		switch(glp_get_status(lp))
		{
		case GLP_OPT:
			printf("---> LP solution is optimal and has value: %g\n", obj_val);
			//glp_mpl_postsolve(tran, lp, GLP_MIP);
			//lpx_print_sol(lp, "solution.sol");
			break;
		case GLP_NOFEAS:
			printf("---> LP problem has no feasible solution (proven by the solver)\n");
			feasible = false;
			break;
		case GLP_UNBND:
			printf("---> LP solution is unbounded\n");
			feasible = false;
			break;
		case GLP_FEAS:
			printf("---> LP solution is feasible and has value %g, however, its optimality"
			       " (or non-optimality) has not been proven\n", obj_val);
			feasible = false;
			break;
		case GLP_UNDEF:
			printf("---> LP solution is undefined\n");
			feasible = false;
			break;
		}      
		printf("%d variables are continuos, %d are integer (of which %d are binary)\n", 
		       (glp_get_num_cols(lp) - glp_get_num_int(lp)), glp_get_num_int(lp), glp_get_num_bin(lp));

	}



	//	glp_init_iocp(&params);
	//		params.presolve = GLP_OFF;



	//		glp_intopt(lp, &params);
	/* * Time recording */
	
	if(ret_sol->solved){
		solve_glp(lp);
	
	
		timer = cpu_time.cpu_time_elapsed(log_fp,"SOLVING");
		printf("time elapsed: SOLVING  %g\n",timer);






		if(!has_mps)
			glp_mpl_postsolve(tran, lp, GLP_MIP);
		string sol_output;
		sol_output = prob.output_path + "/solution.sol";
		lpx_print_sol(lp, sol_output.c_str());
		//sol_output = prob.output_path + "/sensitivity.sol";
		// lpx_print_sens_bnds(lp, sol_output.c_str());

		int num_integer_variables = 0;
		int binary_vars = 0;
		int ncols = glp_get_num_cols(lp);

		for(int c=1; c <= ncols; c++)
		{
			const char * col_name = glp_get_col_name(lp,c);

			VERBOSE(3){
				printf("Var %d is a %s, %s, with bounds [%g, %g]", c, (var_status_to_string(lp,c)).c_str(), 
				       (var_kind_to_string(lp,c)).c_str(),  glp_get_col_lb(lp, c),  glp_get_col_ub(lp, c));
				if(col_name != NULL)
					printf(" name [%s]",col_name);
			}
			double s = glp_mip_col_val(lp,c);
			VERBOSE(3){	
				printf(" value [%f]",s);
				printf("\n");
			}

			if(s > EPSILON){
				string str_col_name(col_name);
				//ret_sol->value[str_col_name] = s;

			}

			if( glp_get_col_kind(lp, c) != GLP_CV )
			{
				//integer_var_index[num_integer_variables++] = c;
				if( glp_get_col_kind(lp, c) == GLP_BV )
				{
					binary_vars++;
					//fprintf(stderr, "Name col %d: %s\n", c, glp_get_col_name(P, c));
				}
			}
		}

		ret_sol->objval = obj_val;
		ret_sol->solved = true;
	
		/* * Finishing time recording */
		fprintf(stderr, "time elapsed: cpu  %g  total \n",cpu_time.cpu_time_elapsed(), 
			cpu_time.total_time_elapsed()); 

		cpu_time.end(log_fp, "SOLVING");

		timer = cpu_time.total_time_elapsed();
		printf("time elapsed: POST SOLVE  %g\n",timer); 
		fclose (log_fp);

		/** Process lp_sol to extract solution info
		 *  Before It was done in create_graph_from_sol
		 *  */
		FOREACH(it,ret_sol->value){
			if(fabs(it->second) > EPSILON)
			write_to_sol((it->first).c_str(),it->second);
		}

		std::sort(relaysToUse.begin(), relaysToUse.end());
		std::vector<std::string>::iterator new_end_pos;
		new_end_pos = std::unique( relaysToUse.begin(), relaysToUse.end() );

		// The elements between new_end_pos and str.end() hold
		// hold their old values, and must be erased to complete
		// the operation:

		relaysToUse.erase( new_end_pos, relaysToUse.end() );

		/** Write additional info to lp_sol  */
		ret_sol->n_relays = relaysToUse.size();
		ret_sol->n_edges = edgesToUse.size();

		/** Output simple solution information */
		printf("########### solution #############\n");
		printf("number of relays: %d\n",relaysToUse.size());
		printf("number of edges: %d\n",edgesToUse.size());
		printf("objective function value: %f\n",obj_val);
	}else{
		printf("NO SOLUTION FOUND!!!!\n");
	}


	freeGLPK(lp,tran);


	printf("*********************  SOLVE FINISHED ****************** \n");
	return ret_sol;
}	


void GLPSol::generate_and_save_mps(Network *net, string fname){
	
	glp_prob *lp;
	glp_tran *tran;
	glp_iocp params;
	CpuTime cpu_time;
	double timer;
	string out_time_file; 
	int ret;
	out_time_file = prob.output_path + prob.instance_id + ".gen.time";
	FILE *log_fp = fopen(out_time_file.c_str(),"w");

	cpu_time.start();  // start timer


	lp = glp_create_prob();
	tran = glp_mpl_alloc_wksp();
	ret = glp_mpl_read_model(tran, prob.model_file.c_str(), 1);
	if (ret != 0) 
		printf("Error on translating model\n");
	else
		printf("Model read: OK\n");
	timer = cpu_time.cpu_time_elapsed(log_fp, "MODEL_READING");
	printf("time elapsed: METRIC_CALC  %g\n",timer);

	Metric *metric;
	if( prob.metric == "PRR")
		metric = new EstPRR_Metric(net);
	else if(prob.metric == "DISTANCE_SQ")
		metric = new DistSQ_Metric(net);

	timer = cpu_time.cpu_time_elapsed(log_fp, "METRIC_CALC");
	printf("time elapsed: METRIC_CALC  %g\n",timer);
	string dataFile = GLPSol::create_data_file(net, prob.K,metric);
	if(dataFile == "")
		printf("Error creating data file\n");

	ret = glp_mpl_read_data(tran, dataFile.c_str());
	if (ret != 0) 
		printf("Error on translating data\n");
	else
		printf("Data read: OK\n");


	/* * Time recording */
	timer = cpu_time.cpu_time_elapsed(log_fp, "DATA READING");
	printf("time elapsed: DATA READING  %g\n",timer);

	if (glp_mpl_generate(tran, NULL) != 0)
		exit_on_error("Error on generating model\n", lp, tran);


	glp_mpl_build_prob(tran, lp);

	
	/** GLPK: Generate Variable indexing **/
	glp_create_index(lp);

	/** GLPK: Generate LP **/

	glp_write_mps(lp, GLP_MPS_FILE, NULL, fname.c_str());

	cpu_time.end(log_fp, "MPS_GENERATION");
	freeGLPK(lp, tran);




}
