/***************************************************************************
  wrapper.c  -  description
  -------------------
begin                : May 2010
copyright            : Matteo Salani
email                : matteo.salani@idsia.ch
description          : This is a wrapper to process a MILP model
written in GNU - MathProg (http://www.gnu.org/software/glpk/)
solve it with Gurobi (http://gurobi.com/)
post-process the output with GLPK
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/** please, have GLPK and GUROBI installed **/
#include "glpk.h"
#include "gurobi_c.h"


#include "wrapper.h"

/** GLPK: data structures **/
glp_prob *mip = NULL;
glp_tran *tran = NULL;
glp_smcp *parm;
glp_iocp *iparm;
int ret, col_index;
double glpk_iparm_mip_gap = 10E-4;
double glpk_iparm_tol_int = 10E-4;
double glpk_iparm_tol_obj = 10E-4;

/** GUROBI: data structures **/
GRBenv   *env   = NULL;
GRBenv   *mipenv = NULL;
GRBmodel *model = NULL;
int      numvars, j, retGRB, GRB_IsMIP;
double   tmp, bound;
char * nameGRB;
char type;

/** execution flags **/
char GRB_out, GLPK_out, verbose;

void freeMem(){

	/** GLPK: free structures **/
	if (tran) glp_mpl_free_wksp(tran);
	if (mip) glp_delete_prob(mip);

	if (retGRB) printf("ERROR: %s\n", GRBgeterrormsg(env));

	/** GUROBI: free structures **/
	if (model) GRBfreemodel(model);
	if (env) GRBfreeenv(env);

	printf("ERRORS OCCURRED.\nGLPK -> GUROBI -> GLPK wrapper v0.1 (2010)\n");

	exit(1);
}

void usage(){
	fprintf(stderr, "GLPK -> GUROBI -> GLPK wrapper v0.1 (2010)\nby Matteo Salani (matteo.salani@idsia.ch)\n\n"
		"Usage: \n"
		"wrapper_GLPK_GUROBI [OPTIONS] -m <model.mod> -d <data.dat>\n"
		"with OPTIONS:\n"
		"--glpk_out (enable GLPK output, default disabled)\n"
		"--grb_out (enable GRB output, default disabled)\n"
		"--glpk_mip_gap <value> (default 10E-4)\n"
		"--glpk_tol_int <value> (defaule 10E-4)\n"
		"--glpk_tol_obj <value> (defaule 10E-4)\n");
}


//* Solve GLPK prob using GUROBI


void gurobi_set_basis(GRBmodel *mod){
	int ncons;

	int error;

	error = GRBgetintattr(mod, "NumConstrs", &ncons);
	if (error) freeMem();
	

	printf("Number Constraints %d\n",ncons);



  
	for(int i=0; i < ncons; i++){
		error = GRBsetintattrelement(mod, "CBasis", i, 0);

		if(error) freeMem();
	}

	int nvars;

	error = GRBgetintattr(mod, "NumVars", &nvars);
        printf("Number Variables %d\n",ncons);
	for(int i =0;i<nvars;i++){
		error = GRBsetintattrelement(mod, "VBasis", i, 0);
		if(error) freeMem();
	}



}

GRBmodel *fixed_model(GRBmodel *mdl0)
{
	GRBenv *env;
	GRBmodel *mdl;
	double f, *y;
	int i;
	static char *statusname[] = {
		"infeasible",
		"infeasible or unbounded",
		"unbounded",
		"cutoff",
		"iteration limit",
		"node limit",
		"time limit",
		"solution limit",
		"interrupted",
		"numeric difficulty"
		};

	if (!(mdl = GRBfixedmodel(mdl0)))
		return 0;
	if (!(env = GRBgetenv(mdl))) {
		//dpf(d, "\nGRBgetenv failed in fixed_model().");
 badret:
		GRBfreemodel(mdl);
		return 0;
		}
	if (GRBsetintparam(env, "Presolve", 0)) {
		//intbasis_fail(d, "setintparam(\"Presolve\")");
		goto badret;
		}

	gurobi_set_basis(mdl);
	if (GRBoptimize(mdl)) {
		//intbasis_fail(d, "optimize()");
		goto badret;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i)) {
		//intbasis_fail(d, "getintattr()");
		goto badret;
		}
	if (i != GRB_OPTIMAL) {
//		if (i >= GRB_INFEASIBLE && i <= GRB_NUMERIC)
			//dpf(d, "\nGRBoptimize of fixed model: %s.",
			//	statusname[i-GRB_INFEASIBLE]);
//		else
			//dpf(d, "\nSurprise status %d after GRBoptimize of fixed model.",
				//i);
		goto badret;
		}
/*  	if (d->missing & 2 && (y = d->y0)
	 && !GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, n_con, y)) {
		d->y = y;
		d->missing &= ~2;
		}
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f)) {
		if (f > 0.)
//			dpf(d, "\nplus %.0f simplex iteration%s for intbasis",
//				f, "s" + (f == 1.));
		}
		*/
	return mdl;
	}


void gurobi_sens_output(GRBmodel *mod, char *filename){
#define DBL_INF 10e10
	int ncons;

	int error;

	error = GRBgetintattr(mod, "NumConstrs", &ncons);
	if (error) freeMem();
	
	FILE *outfile = fopen(filename, "w");	

	fprintf(outfile, "Sensitivity Analysis\n");

	fprintf(outfile,"Name\tSlack\tSARHSLow\tSARHSUp\n");
	printf("Number of constraints: %d\n",ncons);

	for(int i=0; i < ncons; i++){

		char *consname;
  		error = GRBgetstrattrelement(mod, "ConstrName", i, &consname);
		if (error) freeMem();
		double slack=0, sarhslow=0, sarhsup=0;


	
  		error = GRBgetdblattrelement(mod, "Slack", i, &slack);
		if (error) freeMem();
  		error = GRBgetdblattrelement(mod, "SARHSLow", i, &sarhslow);
		if (error) freeMem();
  		error = GRBgetdblattrelement(mod, "SARHSUp", i, &sarhsup);
		if (error) freeMem();

 
		
		//printf("%s\t%lf\t%lf\t%lf\n",consname,slack,sarhslow,sarhsup);
	
		fprintf(outfile,"%s\t%f\t",consname,slack);
		if(fabs(sarhslow) < DBL_INF)
			fprintf(outfile, "%f",sarhslow);
		else{
			fprintf(outfile, "%c",(sarhslow>0?'+':'-'));
			fprintf(outfile, "inf");
		}
		fprintf(outfile,"\t");

		if(fabs(sarhsup) < DBL_INF)
			fprintf(outfile, "%f",sarhsup);
		else{
			fprintf(outfile, "%c",(sarhsup>0?'+':'-'));
			fprintf(outfile, "inf");
		}
		fprintf(outfile,"\n");
		       	




	}
	printf("\n");
	fclose(outfile);


}

double solve_glp_grb(glp_prob *mip, wrapper_params *par){


	GLPK_out = par->glp_out;
	GRB_out = par->grb_out;
	double obj_val;



	/** GLPK: Generate Variable indexing **/
	glp_create_index(mip);

	/** GLPK: Generate LP **/
	glp_write_mps(mip, GLP_MPS_FILE, NULL, "tmp.mps");


	/************/
	/** GUROBI **/
	/************/

	retGRB = GRBloadenv(&env, NULL);
	if (retGRB || env == NULL)
	{
		fprintf(stderr, "Error: could not create environment\n");
		exit(1);
	}

	retGRB = GRBsetintparam(env, "OutputFlag", GRB_out?1:0);
	if (retGRB) freeMem();

	//retGRB = GRBsetintparam(env, "Sensitivity", 1);
	//if (retGRB) freeMem();

	/** GUROBI: Read model **/
	retGRB = GRBreadmodel(env, "tmp.mps", &model);
	if (retGRB) freeMem();

	/** Remove utility files from disk **/
	//remove("tmp.mps");

	/** GUROBI: Get environment **/
	mipenv = GRBgetenv(model);
	if (!mipenv) freeMem();

	/** GUROBI: Set parameters **/

	/** GUROBI: Ask for more precision **/
	retGRB = GRBsetdblparam(mipenv, "FeasibilityTol", 10E-6);
	if (retGRB) freeMem();
	retGRB = GRBsetdblparam(mipenv, "IntFeasTol", 10E-5);
	if (retGRB) freeMem();
	retGRB = GRBsetdblparam(mipenv, "MIPgap", 10E-6);
	if (retGRB) freeMem();

	/* * Playing with gurobi parameters and attr*/

	//gurobi_set_basis();
	retGRB = GRBsetintparam(mipenv, "Cuts", 3);
	if (retGRB) freeMem();

	retGRB = GRBsetintparam(mipenv, "RootMethod", 1);
	if (retGRB) freeMem();

	retGRB = GRBsetintparam(mipenv, "Symmetry", -1);
	if (retGRB) freeMem();

	

	/** GUROBI: get numvars and numrows **/
	retGRB = GRBgetintattr(model, "NumVars", &numvars);
	if (retGRB) freeMem();


	/** Test variable names */
	for(int j=0;j<numvars;j++){	
		retGRB = GRBgetstrattrelement(model, "VarName", j, &nameGRB);
		printf("GRB Var %d Name %s\n",j,nameGRB); 
	}
	/** GUROBI: get model type **/
	retGRB = GRBgetintattr(model, "IsMIP", &GRB_IsMIP);
	if (retGRB) freeMem();

	/** GUROBI: Optimize model **/
	retGRB = GRBoptimize(model);
	if (retGRB) freeMem();

	
	
	/** GUROBI: Retreive the optimization status **/
	GRBgetintattr(model, "Status", &retGRB);
	switch(retGRB){
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
		freeMem();
	}

	/** GUROBI: Get obj function value **/
	retGRB = GRBgetdblattr(model, "IntVio", &tmp);
	if (retGRB) freeMem();


	retGRB = GRBgetdblattr(model, "ObjBound", &bound);
	if (retGRB) freeMem();

	retGRB = GRBgetdblattr(model, "ObjVal", &tmp);
	if (retGRB) freeMem();

	/* ********************** */

	obj_val = tmp;


	/* ************ */
	if (verbose) printf ("Objective %lf\n", tmp);
	if (verbose) printf ("Best bound %lf\n", bound);
	if (verbose) printf ("Absolute gap %lf\n", fabs(tmp - bound));

	/** GUROBI: Get variable values **/
	for (j = 0; j < numvars; ++j){

		retGRB = GRBgetdblattrelement(model, "X", j, &tmp);
		if (retGRB) freeMem();

		retGRB = GRBgetstrattrelement(model, "VarName", j, &nameGRB);
		printf("GRB Var %d Name %s\n",j,nameGRB); 
		if (retGRB) freeMem();

		retGRB = GRBgetcharattrelement(model, "VType", j, &type);
		if (retGRB) freeMem();

		/** GLPK search variable index by name **/
		col_index = glp_find_col(mip, nameGRB);

		if (col_index != 0){
			/** GLPK set variable bounds **/
			if ((type == 'B') || (type == 'I')){
				if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", nameGRB, type, tmp, round(tmp));
				glp_set_col_bnds(mip, col_index, GLP_FX, round(tmp), round(tmp));
			}
			else{
				if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", nameGRB, type, tmp, tmp);
				glp_set_col_bnds(mip, col_index, GLP_FX, tmp, tmp);
			}
		}
	}

	if (GRB_IsMIP){

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
			freeMem();
		}

		ret = glp_mip_status(mip);
		switch (ret){
		case GLP_OPT:
			break;
		case GLP_FEAS:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_FEAS, code %d\n", ret);
			freeMem();
		case GLP_NOFEAS:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_NOFEAS, code %d\n", ret);
			freeMem();
		case GLP_UNDEF:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNDEF, code %d\n", ret);
			freeMem();
		}
	}
	else{

		/*GLPK initialize parameters */
		parm = (glp_smcp*) malloc(sizeof(glp_smcp));
		glp_init_smcp(parm);
		parm->meth = GLP_DUALP;
		parm->tol_bnd = 10E-4;
		parm->tol_dj = 10E-4;

		/* GLPK get the optimal basis */
		//ret = glp_simplex(mip, parm);
		if (ret){
			fprintf(stderr, "glp_simplex, Error on optimizing the model : %d \n", ret);
			freeMem();
		}
		ret = glp_get_status(mip);
		switch (ret){
		case GLP_OPT:
			break;
		case GLP_FEAS:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_FEAS, code %d\n", ret);
			freeMem();
		case GLP_INFEAS:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_INFEAS, code %d\n", ret);
			freeMem();
		case GLP_NOFEAS:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_NOFEAS, code %d\n", ret);
			freeMem();
		case GLP_UNBND:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNBND, code %d\n", ret);
			freeMem();
		case GLP_UNDEF:
			fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNDEF, code %d\n", ret);
			freeMem();
		}


	}

	//GRBmodel *fmod = fixed_model(model);
	//gurobi_sens_output(fmod, "/tmp/sens.sol");
        GRBwrite(model, "/tmp/model.sol");






	/** GUROBI: free structures **/
	if (model) GRBfreemodel(model);
	if (env) GRBfreeenv(env);

	return obj_val;
}


int main2(int argc, char *argv[]){
	char * file_model;
	char * file_data;

	if (argc < 5){
		usage();
		exit(1);
	}

	file_model = file_data = NULL;
	GRB_out = GLPK_out = verbose = 1;

	for(int i=1 ; i<=argc-1 ; i++){
		if (strcmp(argv[i],"-m")==0) file_model = argv[i+1];
		if (strcmp(argv[i],"-d")==0) file_data = argv[i+1];
		if (strcmp(argv[i],"-v")==0) verbose = 1;
		if (strcmp(argv[i],"--glpk_out")==0) GLPK_out = 1;
		if (strcmp(argv[i],"--grb_out")==0) GRB_out = 1;
		if (strcmp(argv[i],"--glpk_mip_gap")==0) glpk_iparm_mip_gap = atof(argv[i+1]);
		if (strcmp(argv[i],"--glpk_tol_int")==0) glpk_iparm_tol_int = atof(argv[i+1]);
		if (strcmp(argv[i],"--glpk_tol_obj")==0) glpk_iparm_tol_obj = atof(argv[i+1]);
	}

	if ((file_model==NULL) || (file_data == NULL)){
		usage();
		fprintf(stderr, "Error no model or data files provided\n");
		freeMem();
	}

	/** GLPK: Open environment **/
	mip = glp_create_prob();
	tran = glp_mpl_alloc_wksp();

	glp_term_out(GLPK_out?GLP_ON:GLP_OFF);

	/** GLPK: Read model written in MathProg **/
	ret = glp_mpl_read_model(tran, file_model, 1);

	if (ret){
		fprintf(stderr, "Error on translating model\n");
		freeMem();
	}

	/** GLPK: Read data for MathProg **/
	ret = glp_mpl_read_data(tran, file_data);
	if (ret){
		fprintf(stderr, "Error on translating data\n");
		freeMem();
	}

	/** GLPK: Generate model (merge data an model) **/
	ret = glp_mpl_generate(tran, NULL);
	if (ret){
		fprintf(stderr, "Error on generating model\n");
		freeMem();
	}

	/** GLPK: Generate Build Model **/
	glp_mpl_build_prob(tran, mip);

	wrapper_params wpar;
	wpar.grb_out = GRB_out;
	wpar.glp_out = GLPK_out;
	solve_glp_grb(mip,&wpar);
	/** GLPK: Perform postprocessing **/
	ret = glp_mpl_postsolve(tran, mip, GLP_MIP);
	if (ret != 0) fprintf(stderr, "Error on postsolving model\n");

	/** GLPK: free structures **/
	if (tran) glp_mpl_free_wksp(tran);
	if (mip) glp_delete_prob(mip);

	if (retGRB) printf("ERROR: %s\n", GRBgeterrormsg(env));

	/** GUROBI: free structures **/
	if (model) GRBfreemodel(model);
	if (env) GRBfreeenv(env);

	printf("Done.\nGLPK -> GUROBI -> GLPK wrapper v0.1 (2010)\n");

	exit(0);
}

/*int main(int argc, char *argv[]){
	main2(argc,argv);
}
*/
