
#include "main.h"
#include "solver.h"

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <ios>
#define CPX_INT 0
#define CPX_DB 1



ofstream o_progress;
FILE *log_fp;

#ifdef __WITH_GLPK
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

string GLPSol::create_data_file(WriteToDataFileFunc f){
	string data_filename = param.data_file;
	ofstream dataFile(data_filename.c_str());


	if(param.verbose) printf("Creating data file ..\n");
	if( dataFile.is_open()){
		//prob->writeToDataFile(dataFile);
		(*f)(dataFile);

	}


	return data_filename;


}

int GLPSol::get_var_type(const char *v){
/*
 * 
 	if(v[0] == VAR_RELAY)
		return VARTYPE_RELAY;
	else if(v[0] == VAR_EDGE)
		return VARTYPE_EDGE;
	else
		return VARTYPE_NULL;

		*/
}




void GLPSol::write_to_sol(const char *varname, double value){
	printf("var %s %f\n",varname,value);

}
void GLPSol::fix_glp_with_lpsol(glp_prob *mip, LpSolution *sol){
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


#endif

void LpSolution::write_to_file(string filename) const{


}


void LpSolution::print(){
//	FOREACH(it,values){
//		cout << (*it)->first << " " << (*it)->second << endl;
//	}
}

LpSolution::LpSolution(string filename){

}




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

	if(value)
	{
		char *errstr;
		char buffer[4096];
		errstr = (char *)CPXgeterrorstring(cpx_env, value, buffer);
		if ( errstr != NULL ) {
			fprintf (stderr, "%s\n", buffer);
		}
		else {
			fprintf (stderr, "CPLEX Error %5d: Unknown error code\n",
				value);
		}

		fprintf(stderr,"Context: %s, Value: %d.\nGLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\n", context, value);
	}
	

	
}
	static int CPXPUBLIC
timelimcallback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{
	int status = 0;
	double gap = 0.0;
	double timenow = 0.0;
	double bestobj = 0.0;


	CPX_TIMELIMINFOptr info = (CPX_TIMELIMINFOptr) cbhandle;
	int            hasincumbent = 0;


	status = CPXgetcallbackinfo (env, cbdata, wherefrom,
				     CPX_CALLBACK_INFO_MIP_FEAS, &hasincumbent);
	if ( status )  goto TERMINATE;



	if ( !info->aborted ) {
		status = CPXgetcallbackinfo (env, cbdata, wherefrom,
					     CPX_CALLBACK_INFO_MIP_REL_GAP, &gap);
		if ( status )  goto TERMINATE;

		/* Turn the gap into a percentage */
		gap *= 100.0;

		status = CPXgettime (env, &timenow);
		if ( status )  goto TERMINATE;

		if ( hasincumbent)
		{
			status = CPXgetcallbackinfo(env, cbdata, wherefrom, 
						    CPX_CALLBACK_INFO_BEST_INTEGER, &bestobj);
			if(status)
				goto TERMINATE;

			if ( timenow - info->timestart > info->timelim  &&
			     gap < info->acceptablegap                    ) {

				fprintf (stderr,
					 "Good enough solution at time %.2fsec, gap = %g%%, obj %.2f\n",
					 timenow - info->timestart, gap, bestobj);

				/* callback may be called again during the clean up phase after the
				   abort has been issued, so remember that abort processing has
				   already occurred. */

				info->aborted = 1;


				/* 
				 * Obtaining a pointer to CPXLPptr Does not work
				 CPXLPptr cpx_lp;
				 printf("Trying to get incumbent - if something happened, there was an error\n");
				 status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_USER_PROBLEM,
				 &cpx_lp);
				 if ( status ) {
				 printf("Error CPX_CALLBACK_INFO_USER_PROBLEM %d\n", status);
				 goto TERMINATE;
				 }


*/

				int numvars = info->numvars;
				double *vals = (double*)malloc(numvars*sizeof(double));
				status = CPXgetcallbackincumbent(env,cbdata,wherefrom, vals, 0, numvars-1); 
				if(status)
				{
					printf("Error CPXgetcallbackincumbent %d\n", status);
				}
				else
				{
					printf("Good CPXgetcallbackincumbent\n");
				}
				// Get incumbent 
				LpSolution *incumbent = new LpSolution();
				CPX_VARIABLESINFO &varsinfo = (*info->varsinfo);

				for(int i=0; i< numvars; i++)
				{
					std::pair<string, char>  varinfo = varsinfo[i];
					string nameVar = varinfo.first;
					char type = varinfo.second;
					incumbent->value[nameVar] = vals[i];
					incumbent->type[nameVar] = type;
				}

				incumbent->is_optimal(false);
				incumbent->solved = true;
				incumbent->objval = bestobj;
				info->incumbent = incumbent;
				printf("Incumbent solution set\n");

				free(vals);
			}
		}

		// Check if we should report solving progress
		// If that's the case, the ofstream should be 
		// initialized and ready
		if(info->output_progress)
		{
			o_progress << (timenow - info->timestart);
			if( hasincumbent)
			{
				o_progress << " " << bestobj << " " << gap;
			}
			else
			{
				o_progress << " NULL";
			}
			o_progress << endl;
		}


		// Call use-defined callback
		CallbackInfo cbinfo;
		cbinfo.has_incumbent = hasincumbent;
		cbinfo.timenow = timenow - info->timestart;
		cbinfo.gap = gap;
		cbinfo.bestobj = bestobj;

		if(info->simplecallback != NULL)
		{
			(*info->simplecallback)(&cbinfo);
		}







		
	}
	return info->aborted;

TERMINATE:

		return (status);

} /* END timelimcallback */


CPX_VARIABLESINFO *CPXgetVariablesInfo(CPXENVptr env, CPXLPptr cpx_lp)
{

	/** CPLEX: get numvars **/

	int numvars = CPXgetnumcols(env, cpx_lp);
	int j;
	int retCPX;
	double db_tmp;
	char nameCPX[80];
	string nameVar;
	bool errorCPX = false;
	CPX_VARIABLESINFO *ret_map = 
		new CPX_VARIABLESINFO();

	for (j = 0; j < numvars; ++j)
	{


		int surplus, colnamespace; 
		char ** colname = NULL;
		char * colnamestore = NULL;

		retCPX = CPXgetcolname(env, cpx_lp, NULL, NULL, 0, &surplus, j, j);
		if (( retCPX != CPXERR_NEGATIVE_SURPLUS ) && ( retCPX != 0 ))  
		{
			fprintf(stderr, "Error CPX Get var names %d\n", retCPX);
			errorCPX = true;
			break;
		}

		colnamespace = - surplus;
		if ( colnamespace > 0 ) {

			colname = (char **) malloc (sizeof(char *));
			colnamestore = (char *)  malloc(colnamespace);

			if ( colname == NULL || colnamestore == NULL ) 
			{

				fprintf(stderr, "Error CPX Allocating memory for col name \n");
				
				if(colname) free( colname );
				if(colnamestore) free( colnamestore );
				errorCPX = true;
				break;
			}

			retCPX = CPXgetcolname (env, cpx_lp, colname, colnamestore, colnamespace, &surplus, j, j);
			if ( retCPX ) 
			{
				fprintf(stderr, "Error CPX Get final names %d\n", retCPX);
				free( colname );
				free( colnamestore );
				errorCPX = true;
				break;
			}
		}
		else 
		{

			fprintf(stderr, "Error CPX no name associated\n");
			errorCPX = true;
			break;
		}

		sprintf(nameCPX, "%s", colname[0]);

		free(colnamestore);
		free(colname);

		char type[1];
		retCPX = CPXgetctype (env, cpx_lp, type, j,j);
		if (retCPX) {

			fprintf(stderr, "Error CPX Accessing variable type %d\n", retCPX);
			errorCPX = true;
			break;
		}

		nameVar = nameCPX;
		(*ret_map)[j] = make_pair(nameVar, type[0]);	



	}
	if( errorCPX)
	{
		delete ret_map;
		return NULL;
	}

	return ret_map;

}

LpSolution *GLPSol::solve_cplex(string mps_file, LpSolution *prev_sol, CplexParams *params){
	/** CPLEX: data structures **/
	CPXENVptr cpx_env = NULL;
	CPXLPptr cpx_lp = NULL;
	/** CPLEX: variables */
	int      numvars, j, retCPX, CPX_IsMIP, CPX_noMIP, CPX_probType,cpx_probe;
	double   db_tmp, bound,obj_val;
	char nameCPX[80];
	char type;
	double cpx_feasibilityTol = 10E-8;
	double cpx_intFeasTol = 10E-8;
	double cpx_mipGap = 10E-8;
	double cpx_trelim = 2048; // MEMORY LIMIT 2GB
	/** execution flags **/
	char CPX_out, GLPK_out, verbose;
  	CPX_out = GLPK_out = CPX_noMIP = verbose = 0;
	double cpx_timelim = param.time_limit; // Time limit
	
	int usetimelimcallback = 
		(param.get_suboptimal || m_useSimpleCallback || param.output_progress);
	CPX_TIMELIMINFO mytimeliminfo;


	CPX_out = param.verbose;
	GLPK_out = param.verbose;

	CPX_noMIP = param.noMIP;

	verbose = param.verbose;
	

	/// Time recording
	m_cpuTimeModelGeneration = cpu_time.cpu_time_elapsed(log_fp,"GENERATION");


	string nameVar;
	LpSolution *ret_sol = new LpSolution();


	/** CPLEX: user defined parameters **/

	cpx_env = CPXopenCPLEX(&retCPX);
	if (retCPX || cpx_env == NULL)
	{
		fprintf(stderr, "Error: could not create cpx_environment\n");
		free_cplex("CPX create environment", retCPX, cpx_env, cpx_lp);
	}

	retCPX = CPXsetintparam(cpx_env, CPX_PARAM_SCRIND, CPX_out?CPX_ON:CPX_OFF);
	if (retCPX) free_cplex("CPX set output", retCPX, cpx_env, cpx_lp);
	
	retCPX = CPXsetintparam(cpx_env, CPX_PARAM_PROBE, 0);
	if (retCPX) free_cplex("CPX set probe", retCPX, cpx_env, cpx_lp);

	if(!param.allow_parallel)
	{
		retCPX = CPXsetintparam(cpx_env, CPX_PARAM_THREADS, 1);
		if (retCPX) free_cplex("CPX set threads", retCPX, cpx_env, cpx_lp);
	}

	retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_TRELIM , cpx_trelim);
	if (retCPX) free_cplex("CPX set trelim", retCPX, cpx_env, cpx_lp);

	/** CPLEX: Create LP model **/
	cpx_lp = CPXcreateprob(cpx_env, &retCPX, "cpx_lp");
	if (retCPX) free_cplex("CPX create problem", retCPX, cpx_env, cpx_lp);

	/** CPLEX: Read model **/
	retCPX = CPXreadcopyprob (cpx_env, cpx_lp, mps_file.c_str(), "MPS");
	if (retCPX)
	{
		free_cplex("CPX load model", retCPX, cpx_env, cpx_lp);
		fprintf(stderr, "Error reading mps %s\n",mps_file.c_str());
	}

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
	if(param.strong_timelim)
	{
		retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_TILIM, cpx_timelim);
		if (retCPX) free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
	}


	/** Set other user-defined parameters for Cplex **/
	CplexParams *tmp = params;
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

	/** CPLEX: Callback setting **/

	if ( usetimelimcallback ) {
		double t;
		retCPX = CPXgettime (cpx_env, &t);
		if ( retCPX ) {
			fprintf (stderr, "Failed to initialize timer.\n");
			free_cplex("CPX to initialize timer\n", retCPX, cpx_env, cpx_lp);
		}
		mytimeliminfo.acceptablegap = 100.0;
		mytimeliminfo.aborted       = 0;
		mytimeliminfo.timestart     = t;
		mytimeliminfo.timelim       = param.time_limit;
		mytimeliminfo.numvars       = numvars;
		mytimeliminfo.varsinfo      = CPXgetVariablesInfo(cpx_env, cpx_lp);
		mytimeliminfo.incumbent = NULL;
		mytimeliminfo.simplecallback = (m_useSimpleCallback?m_simpleCallback:NULL);
		mytimeliminfo.output_progress = param.output_progress;

		retCPX = CPXsetinfocallbackfunc (cpx_env, timelimcallback, &mytimeliminfo);
		if ( retCPX ) {
			fprintf (stderr, "Failed to set time limit callback function.\n");
			free_cplex("CPX Failed callback set", retCPX, cpx_env, cpx_lp);
			
		}
	}



	/** CPLEX: get model type **/
	CPX_probType = CPXgetprobtype(cpx_env, cpx_lp);
	if ((CPX_probType != CPXPROB_MILP) && (CPX_probType != CPXPROB_LP))
	{
		fprintf(stderr, "Error: unable to solve problem type :%d\n", CPX_probType);
		free_cplex("CPX problem type", 1, cpx_env, cpx_lp);
	}
	if (verbose) printf ("Problem type: %d\n", CPX_probType);
	CPX_IsMIP = (CPXPROB_MILP == CPX_probType);

	/** CPLEX: Solve the linear relaxation **/
	if ((!CPX_IsMIP)){
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
		if(CPX_noMIP)
		{

			retCPX = CPXchgprobtype(cpx_env, cpx_lp,CPXPROB_LP);
			if (retCPX) free_cplex("CPX change problem type", retCPX, cpx_env, cpx_lp);

			retCPX = CPXlpopt(cpx_env, cpx_lp);
			if (retCPX) free_cplex("CPX LP Optimize", retCPX, cpx_env, cpx_lp);

			/** CPLEX: Retreive the optimization status **/
			retCPX = CPXgetstat(cpx_env, cpx_lp);
			if ((retCPX == CPX_STAT_UNBOUNDED) || (retCPX == CPX_STAT_INFEASIBLE)
			    || (retCPX == CPX_STAT_INForUNBD)){
				printf("Model is infeasible or unbounded\n");

				/** CPLEX: Quit in any case non optimal **/
				free_cplex("Error CPX LP optimization failed", retCPX, cpx_env, cpx_lp);
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
				printf("Optimal NOT FOUND checking for SUBOPTIMAL\n");
				if( param.get_suboptimal)
				{
					if(mytimeliminfo.incumbent != NULL)
					{

						delete ret_sol;
						ret_sol = mytimeliminfo.incumbent;
						if(verbose)
						{
							printf("Obtained suboptimal\n");
						}

						free_cplex("Finished", 0, cpx_env, cpx_lp);
						return ret_sol;
					}
					else
					{
						if(verbose)
						{
							printf("Could not obtain sub-optimal\n");
						}
					}
				}

				/** CPLEX: Quit in any case non optimal **/
				free_cplex("Error CPX MIP optimization failed", retCPX, cpx_env, cpx_lp);
				ret_sol->solved = false;
				return ret_sol;
			}
		}
	}
	/** CPLEX: Get obj function value **/
	retCPX = CPXgetobjval(cpx_env, cpx_lp, &db_tmp);
	if (retCPX) free_cplex("CPX obj value", retCPX, cpx_env, cpx_lp);
	if (verbose) printf ("Objective %lf\n", db_tmp);

	if (!CPX_noMIP){
		retCPX = CPXgetbestobjval(cpx_env, cpx_lp, &bound);
		if (retCPX) free_cplex("CPX accessing bound", retCPX, cpx_env, cpx_lp);
		if (verbose) printf ("Best bound %lf\n", bound);
		if (verbose) printf ("Absolute gap %lf\n", fabs(db_tmp - bound));
	}


	
	obj_val = db_tmp;
	ret_sol->objval = obj_val;
	ret_sol->is_optimal(true);
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
				free_cplex("CPX Allocating memory for col name" , 1, cpx_env, cpx_lp);
			}

			retCPX = CPXgetcolname (cpx_env, cpx_lp, colname, colnamestore, colnamespace, &surplus, j, j);
			if ( retCPX ) {
				free_cplex("CPX Get final var names", retCPX, cpx_env, cpx_lp);
			}
		}
		else {
			free_cplex("CPX no name associated", 1, cpx_env, cpx_lp);
		}

		if (verbose) printf ("Processed variable %d name %s\n", j, colname[0]);

		sprintf(nameCPX, "%s", colname[0]);

		free(colnamestore);
		free(colname);

		char type[1];
		if( !CPX_noMIP)
		{
			retCPX = CPXgetctype (cpx_env, cpx_lp, type, j,j);
			if (retCPX) free_cplex("CPX Accessing variable type", retCPX, cpx_env, cpx_lp);
		}
		else
		{
			type[0] = 'C'; // continuous variable
		}


		nameVar = nameCPX;	
	

		double tmpdbl = db_tmp;
		if(tmpdbl > EPSILON){
			ret_sol->value[nameVar] = tmpdbl;
			ret_sol->type[nameVar] = type[0];
			//if(fabs(ret_sol->value[nameVar]) > 0.0){
			//cout << "Writing sol: " << nameVar 
			//	<< " has value " << ret_sol->value[nameVar] << endl;
			//}
		//	DEBUG(1) cdebug << "Writing sol: " << nameVar 
		//		<< " has value " << ret_sol->value[nameVar] << endl;
		}

	}

	/** TO DO: Clean everything
	 * */
	
	free_cplex("Finished", 0, cpx_env, cpx_lp);
	return ret_sol;




}

LpSolution *GLPSol::solve_cplex(glp_prob *mip, LpSolution *prev_sol, CplexParams *params){
		/** GLPK: Generate Variable indexing **/
	glp_create_index(mip);

	/** GLPK: Generate LP **/

	glp_write_mps(mip, GLP_MPS_FILE, NULL, param.mps_filename.c_str());

	LpSolution *retval = solve_cplex(param.mps_filename,prev_sol,params);
	if(!param.keep_mps)
		remove(param.mps_filename.c_str());
	return retval;
}




#endif



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
	

//TODO SET TIMELIMIT

LpSolution *GLPSol::solve_grb(string mps_file, LpSolution *prev_sol, GrbParams &params){
	int GLPK_out = 1;
	int GRB_out = 1;
	double obj_val;
	double dblattr;

	GRBVar *vars = 0;
	double *sol_values;


	char type;
	int col_index;
	double tmp,bound;
	
	LpSolution *ret_sol = new LpSolution();

	/// Time recording
	m_cpuTimeModelGeneration = cpu_time.cpu_time_elapsed(log_fp,"GENERATION");


	int glp_status, grb_status, numvars, ret, GRB_IsMIP;
	string nameVar;
	int verbose = param.verbose;


	/************/
	/** GUROBI **/
	/************/

	GRBEnv *env;
	try {
		env = new GRBEnv();



		env->set(GRB_IntParam_OutputFlag,verbose);

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
		mipenv.set(GRB_IntParam_Method,1);
		mipenv.set(GRB_IntParam_Symmetry,-1);


		/** GUROBI: get numvars and numrows **/
		numvars = model.get(GRB_IntAttr_NumVars);
		vars = model.getVars();

		/** GUROBI: get model type **/
		GRB_IsMIP = model.get(GRB_IntAttr_IsMIP);


		/** Load previous solution when avaliable 
		 *  To do that - use callbacks 
		 *  */

		mycallback *cb;
		if(prev_sol)
		{
			DEBUG(0,"Registering callback\n");

			sol_values = (double *)malloc(numvars*sizeof(double));
			for(int i=0;i<numvars;i++){
				nameVar = vars[i].get(GRB_StringAttr_VarName);
				sol_values[i] = prev_sol->value[nameVar];
				if(sol_values[i] > EPSILON){
					DEBUGP(0,"setting %s = %f\n",nameVar.c_str(),sol_values[i]);
				}
			}
			cb = new mycallback(vars, sol_values, numvars);
			model.setCallback(cb);
			free(sol_values);
		}

		/** GUROBI: Optimize model **/
		model.optimize();

		/** GUROBI: Retreive the optimization status **/
		grb_status = model.get(GRB_IntAttr_Status);

		// Print status message
		if( verbose)
		{
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


			}
		}

		if( grb_status != GRB_OPTIMAL)
		{
			ret_sol->solved = false;
		}
		else
		{



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

			/** Write solution to object LpSolution */

			ret_sol->objval = obj_val;
			ret_sol->solved = true;
			ret_sol->is_optimal(true);
			for(int i=0;i<numvars;i++){
				nameVar = vars[i].get(GRB_StringAttr_VarName);
				double tmpdbl = vars[i].get(GRB_DoubleAttr_X);
				if(tmpdbl > EPSILON){
					ret_sol->value[nameVar] = vars[i].get(GRB_DoubleAttr_X);
					ret_sol->type[nameVar] = vars[i].get(GRB_CharAttr_VType);
				}
			}
		}



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
	delete [] vars;
	delete env;


	return ret_sol;



}


LpSolution *GLPSol::solve_grb(glp_prob *mip, LpSolution *prev_sol, GrbParams &params){
		/** GLPK: Generate Variable indexing **/
	glp_create_index(mip);

	/** GLPK: Generate LP **/

	glp_write_mps(mip, GLP_MPS_FILE, NULL, param.mps_filename.c_str());

	LpSolution *retval = solve_grb(param.mps_filename,prev_sol,params);
	if(!param.keep_mps)
		remove(param.mps_filename.c_str());
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
	free(iparm);



}
LpSolution *GLPSol::solve( WriteToDataFileFunc f, Heuristic *H)
{


	glp_prob *lp;
	glp_tran *tran;
	glp_iocp params;
	int ret;

	LpSolution *ret_sol;

	string out_time_file; 
	double timer;
	//out_time_file = param.output_path + "/" + param.instance_id + ".sol.time";
	//
	out_time_file = "sol.time";
	log_fp = fopen(out_time_file.c_str(),"w");

	bool verbose = param.verbose;

	if(param.output_progress)
	{
		o_progress.open(param.output_progress_file.c_str());
	}


	double obj_val;



	cpu_time.start();  // start timer
	
	lp = glp_create_prob();
	tran = glp_mpl_alloc_wksp();

	glp_term_out(param.verbose?GLP_ON:GLP_OFF);

	if(!param.has_mps)
	{
		/** Generate problem from network structure */

		ret = glp_mpl_read_model(tran, param.model_file.c_str(), 1);
		if (ret != 0) 
			printf("Error on translating model\n");
		else
		{
			if(verbose) printf("Model read: OK\n");
		}
			
		timer = cpu_time.cpu_time_elapsed(log_fp, "MODEL_READING");
		m_cpuTimeModelReading = timer;
		if( verbose) printf("time elapsed: MODEL_READING  %g\n",timer);

		
		string dataFile = create_data_file(f);
		if(dataFile == "")
			printf("Error creating data file\n");

		ret = glp_mpl_read_data(tran, dataFile.c_str());
		if (ret != 0) 
			printf("Error on translating data\n");
		else
		{
			if(verbose) printf("Data read: OK\n");
			
		}


		/* * Time recording */
		timer = cpu_time.cpu_time_elapsed(log_fp, "DATA READING");
		m_cpuTimeDataReading = timer;

		if(verbose) printf("time elapsed: DATA READING  %g\n",timer);

		if (glp_mpl_generate(tran, NULL) != 0)
			exit_on_error("Error on generating model\n", lp, tran);
		glp_mpl_build_prob(tran, lp);


		if(verbose) printf("time elapsed: GENERATION  %g\n",timer);
		
	}
	else
	{
		/** Load mps file */
		printf("Loading mps file %s\n",mps_filename.c_str());
		glp_read_mps(lp,  GLP_MPS_FILE, NULL,mps_filename.c_str());
		glp_create_index(lp);

		timer = cpu_time.cpu_time_elapsed(log_fp,"LOAD_MPS");
		if(param.verbose) printf("time elapsed: LOAD_MPS  %g\n",timer);

	}


	if(param.which_solver == SOLVER_GUROBI)
	{
#ifdef _COMPILE_GUROBI

		if(param.verbose)
		{
			printf("********* USING GUROBI ***********\n");
		}


		LpSolution *prev_sol = NULL;
		/*
		if(H){
			prev_sol = (*H)();
			int cnt = 0;
			FOREACH(jt,prev_sol->value){
				if(jt->second > 0)
					cnt++;
			}
			DEBUG(1)cdebug << "prev sol has " << cnt << " variables > 0" << endl;
		}
		*/
		LpSolution *sol;
		if(param.has_mps)
		{
			sol = solve_grb(mps_filename,prev_sol,param.grb);
		}
		else{
			sol = solve_grb(lp,prev_sol,param.grb);
		}

		ret_sol = sol;
		if(param.verbose)
			printf("********************* GUROBI ENDED - FIXING GLPK *******************\n");
		fix_glp_with_lpsol(lp,sol);
		obj_val = sol->objval;
#else
		printf("GUROBI SUPPORT NOT COMPILED\n");
#endif
		//obj_val = solve_glp_grb(lp, &wp);
	}
	else if(param.which_solver == SOLVER_CPLEX)
	{
#ifdef _COMPILE_CPLEX
		if(param.verbose)
		{
			printf("********* USING CPLEX ***********\n");
		}
		LpSolution *prev_sol = NULL;
		
		LpSolution *sol;
		if(param.has_mps)
			sol = solve_cplex(param.mps_filename,prev_sol,param.cplex_par);
		else
			sol = solve_cplex(lp,prev_sol,param.cplex_par);

		ret_sol = sol;
		if(sol->solved){
			if(param.verbose)
			{
				printf("********************* CPLEX ENDED - FIXING GLPK *******************\n");
			}
			if(param.do_postsolve)
				fix_glp_with_lpsol(lp,sol);
			obj_val = sol->objval;
		}else{
			if(param.verbose)
			{
				printf("********************* CPLEX ENDED - NO SOLUTION FOUND *******************\n");
			}
		
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



	
	if(ret_sol->solved)
	{
		
		if(param.do_postsolve)
			solve_glp(lp);
		
	
		timer = cpu_time.cpu_time_elapsed(log_fp,"SOLVING");
		m_cpuTimeSolving = timer;

		if(param.verbose) printf("time elapsed: SOLVING  %g\n",timer);






		if(!param.has_mps && param.do_postsolve)
			glp_mpl_postsolve(tran, lp, GLP_MIP);
		string sol_output;

		//sol_output = param.output_path + "/solution.sol";
		//lpx_print_sol(lp, sol_output.c_str());
		//sol_output = param.output_path + "/sensitivity.sol";
		// lpx_print_sens_bnds(lp, sol_output.c_str());

		int num_integer_variables = 0;
		int binary_vars = 0;
		int ncols = glp_get_num_cols(lp);

		if( param.do_postsolve)
		{
			for(int c=1; c <= ncols; c++)
			{
				const char * col_name = glp_get_col_name(lp,c);

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
		}
/*
		ret_sol->objval = obj_val;
		ret_sol->is_optimal(true);
		ret_sol->solved = true;
*/
		timer = cpu_time.cpu_time_elapsed(log_fp,"POSTSOLVE");
		m_cpuTimePostSolve = timer;

		if(param.verbose) printf("time elapsed: POSTSOLVE  %g\n",timer);

	
		/* * Finishing time recording */
		if(param.verbose)
		{
			fprintf(stderr, 
				"time elapsed: cpu  %g  total %g\n",
				cpu_time.cpu_time_elapsed(),
				cpu_time.total_time_elapsed()); 
		}
	
		cpu_time.end(log_fp, "SOLVING");

		timer = cpu_time.cpu_time_elapsed();
		m_cpuTimeTotal = timer;

		fclose (log_fp);

	}else
	{
		timer = cpu_time.cpu_time_elapsed(log_fp,"SOLVING");
		m_cpuTimeSolving = timer;

		timer = cpu_time.cpu_time_elapsed(log_fp,"POSTSOLVE");
		m_cpuTimePostSolve = timer;

		timer = cpu_time.cpu_time_elapsed();
		m_cpuTimeTotal = timer;


		if(param.verbose) printf("time elapsed: SOLVING  %g\n",timer);


		if(param.verbose)
		{
			printf("Solver: NO SOLUTION FOUND\n");
		}
	}

	glp_delete_index(lp);
	freeGLPK(lp,tran);

	if(param.output_progress)
	{
		o_progress.close();
	}


	if(param.verbose) {
		printf("*********************  SOLVE FINISHED ****************** \n");
		if(ret_sol->solved)
		{
			printf("Problem solved with obj_val %.2f\n", ret_sol->objval);
		}
	}

	/// Recording cpu times
	ret_sol->cpuTimeModelReading(m_cpuTimeModelReading);
	ret_sol->cpuTimeDataReading(m_cpuTimeDataReading - m_cpuTimeModelReading);
	ret_sol->cpuTimeModelGeneration(m_cpuTimeModelGeneration - m_cpuTimeDataReading);
	ret_sol->cpuTimeSolving(m_cpuTimeSolving - m_cpuTimeModelGeneration);
	ret_sol->cpuTimePostSolve(m_cpuTimePostSolve - m_cpuTimeSolving);
	ret_sol->cpuTimeTotal(m_cpuTimeTotal);
	return ret_sol;
}	


void GLPSol::generate_and_save_mps( WriteToDataFileFunc f, string fname){
	
	glp_prob *lp;
	glp_tran *tran;
	glp_iocp params;
	double timer;
	string out_time_file; 
	int ret;
	//out_time_file = param.output_path + param.instance_id + ".gen.time";
	FILE *log_fp = fopen(out_time_file.c_str(),"w");

	cpu_time.start();  // start timer


	lp = glp_create_prob();
	tran = glp_mpl_alloc_wksp();
	ret = glp_mpl_read_model(tran, param.model_file.c_str(), 1);
	if (ret != 0) 
		printf("Error on translating model\n");
	else
		printf("Model read: OK\n");
	timer = cpu_time.cpu_time_elapsed(log_fp, "MODEL_READING");
	m_cpuTimeModelReading = timer;
	if(param.verbose) printf("time elapsed: METRIC_CALC  %g\n",timer);

	//EstPRR_Metric *metric = new EstPRR_Metric(net);

	//timer = cpu_time.cpu_time_elapsed(log_fp, "METRIC_CALC");
	//printf("time elapsed: METRIC_CALC  %g\n",timer);
	//
	//string dataFile = GLPSol::create_data_file(net, param.K,metric);
	string dataFile = create_data_file(f);
	if(dataFile == "")
		printf("Error creating data file\n");

	ret = glp_mpl_read_data(tran, dataFile.c_str());


	if (ret != 0) 
		printf("Error on translating data\n");
	else
		printf("Data read: OK\n");


	/* * Time recording */
	timer = cpu_time.cpu_time_elapsed(log_fp, "DATA READING");
	m_cpuTimeDataReading = timer;
	if(param.verbose) printf("time elapsed: DATA READING  %g\n",timer);

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
template<class Archive>
void 
LpSolution::serialize(Archive & ar, const unsigned int /* file_version */)
{
	ar & m_cpuTimeModelReading;
	ar & m_cpuTimeDataReading;
	ar & m_cpuTimeModelGeneration;
	ar & m_cpuTimeSolving;
	ar & m_cpuTimePostSolve;
	ar & m_cpuTimeTotal;
	ar & m_is_optimal;
	ar & solved;
	ar & objval;
	ar & value;
	//	printf("LpSolution value len %d\n", value.size());
	//printf("LpSolution value type %d\n", type.size());
//	ar & type;
}

#define MYLIB_ENABLE_COMPRESSION 1
void restore_lpsolution(LpSolution &p, string &str)
{

  //	printf("restore TEXT LpSolution len %\n", str.size());
	namespace bio=boost::iostreams;

	bool enable_compression = MYLIB_ENABLE_COMPRESSION;
	try {
	
		// open the archive
		//std::istringstream mystream(str.c_str());
		std::istringstream mystream;
		mystream.str(str);
		bio::filtering_istream f;
		if(enable_compression)
		{
			f.push(bio::gzip_decompressor());
		}
		f.push(mystream);
	
		//mystream.flush();
		boost::archive::text_iarchive ia(f);
		ia >> p;
	} catch(const std::exception & ex)
	{
		fprintf(stderr, "restore_lpsolution (text) catch expection %s\n", ex.what());
	}

}


void restore_lpsolution(LpSolution &p, char * data, int len)
{

  //	printf("restore BINARY LpSolution len %d\n", len);
	namespace bio=boost::iostreams;
	bool enable_compression = MYLIB_ENABLE_COMPRESSION;
	try {
		typedef boost::iostreams::basic_array_source<char> Device;
		boost::iostreams::stream_buffer<Device> my_stream(data, len);
		bio::filtering_istream f;
		if(enable_compression)
		{
			f.push(bio::gzip_decompressor());
		}
		f.push(my_stream);
		boost::archive::binary_iarchive ia(f);

		// open the archive
		//std::istringstream mystream(str.c_str());
		//mystream.flush();
		//boost::archive::text_iarchive ia(mystream);
		ia >> p;
	} catch(const std::exception & ex)
	{
		fprintf(stderr, "restore_lpsolution catch expection %s\n", ex.what());
	}

}

void 
save_lpsolution(const LpSolution &p, std::vector<char> &my_data, int &buffer_size)
{

	bool enable_compression = MYLIB_ENABLE_COMPRESSION;
//	memset(buffer, 0x0, buffer_size);
	my_data.clear();

	std::cout << "my data has size " 
		<< my_data.size() << std::endl;
	try {

		namespace io=boost::iostreams;
		io::back_insert_device<std::vector<char> > inserter(my_data);
		io::stream<io::back_insert_device<std::vector<char> > > my_stream(inserter);
//		io::stream<io::basic_array<char> > source(buffer);
		{
			io::filtering_stream<io::output> f;
			if(enable_compression)
			{
				f.push(io::gzip_compressor());
			}
			f.push(my_stream);

			//		boost::archive::text_oarchive oa(mystream);
			boost::archive::binary_oarchive oa(f);
			//
			printf("binary serialize lpsolution start\n");
			oa << p;
			printf("before flushing\n");
			f.flush();
		}

		//buffer_size =  io::seek(source, 0, std::ios_base::end);
		buffer_size = my_data.size();

		std::cout << "my data has size " 
			<<buffer_size << std::endl;

	} catch(const std::exception & ex)
	{
		fprintf(stderr, "save_lpsolution (binary) catch expection %s\n", ex.what());
	}

}

void save_lpsolution(const LpSolution &p, string &str)
{

	bool enable_compression = MYLIB_ENABLE_COMPRESSION;
	try {
		// open the archive
		std::ostringstream mystream;

		{
			namespace bio=boost::iostreams;
			bio::filtering_stream<bio::output> f;

			if(enable_compression)
			{
				f.push(bio::gzip_compressor());
			}
			f.push(mystream);

			//		boost::archive::text_oarchive oa(mystream);
			boost::archive::text_oarchive oa(f);
			//
			printf("text serialize lpsolution start\n");
			oa << p;
			f.flush();
		}
		str = mystream.str();

		//		printf("text serialize lpsolution len %d\n",
		//     str.size());

	} catch(const std::exception & ex)
	{
		fprintf(stderr, "restore_lpsolution catch expection %s\n", ex.what());
	}
}

