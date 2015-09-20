/***************************************************************************
                          wrapper_cplex.c  -  description
                             -------------------
    begin                : October 2010
    copyright            : Matteo Salani
    email                : matteo.salani@idsia.ch
    description          : This is a wrapper to process a MILP model 
                           written in GNU - MathProg (http://www.gnu.org/software/glpk/)
                           solve it with CPLEX
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

/** please, have GLPK and CPLEX installed **/
#include "glpk.h"
#include "cplex.h"

#include "cplex_consts.h"

#define CPX_INT 0
#define CPX_DB 1

/** GLPK: data structures **/
glp_prob *mip = NULL;
glp_tran *tran = NULL;
glp_smcp *parm;
glp_iocp *iparm;
int ret, col_index;
double glpk_iparm_mip_gap = 10E-4;
double glpk_iparm_tol_int = 10E-4;
double glpk_iparm_tol_obj = 10E-4;

/** CPLEX: data structures **/
CPXENVptr cpx_env = NULL;
CPXLPptr cpx_lp = NULL;

int      numvars, j, retCPX, CPX_IsMIP, CPX_noMIP, CPX_probType;
double   db_tmp, bound;
char nameCPX[80];
char type;
double cpx_feasibilityTol = 10E-6;
double cpx_intFeasTol = 10E-5;
double cpx_mipGap = 10E-4;

/** CPLEX: user defined parameters **/
typedef struct _cpx_udef {
  char * name;
  int cpx_index;
  char type;
  double db_value;
  int int_value;
  struct _cpx_udef * next;
} cpx_udef;

cpx_udef * cpx_user_def_params;

/** execution flags **/
char CPX_out, GLPK_out, verbose;

int keep_tmp_mps;

void freeMem(char * context, int value){
cpx_udef * tmp;

  /** GLPK: free user defined parameters **/
  while (cpx_user_def_params != NULL){
    tmp = cpx_user_def_params;
    cpx_user_def_params = cpx_user_def_params->next;
    free(tmp);
  }

  /** GLPK: free structures **/
  if (tran) glp_mpl_free_wksp(tran);
  if (mip) glp_delete_prob(mip);

  /** CPLEX: free structures **/
  if (cpx_lp) CPXfreeprob(cpx_env, &cpx_lp);
  if (cpx_env) CPXcloseCPLEX(&cpx_env);

  printf("ERRORS OCCURRED.\nContext: %s, Value: %d.\nGLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\n", context, value);

  exit(1);
}

void usage(){
    fprintf(stderr, "GLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\nby Matteo Salani (matteo.salani@idsia.ch)\n\n"
                    "Usage: \n"
                    "wrapper_GLPK_CPLEX [OPTIONS] -m <model.mod> -d <data.dat>\n"
                    "with OPTIONS:\n"
                    "-v verbose\n"
                    "--glpk_out (enable GLPK output, default disabled)\n"
                    "--cpx_out (enable CPX output, default disabled)\n"
                    "--glpk_mip_gap <value> GLPK mip gap (default 10E-4)\n"
                    "--glpk_tol_int <value> GLPK integer tolerance (default 10E-4)\n"
                    "--glpk_tol_obj <value> GLPK objective tolerance (default 10E-4)\n"
                    "--cpx_feas_tol <value> CPLEX feasibility tolerance (default 10E-6)\n"
                    "--cpx_int_tol <value> CPLEX integer tolerance (default 10E-5)\n"
                    "--cpx_set_db_para <parameter name> <value> Set CPLEX double paramenter (example --cpx_set_db_para IntFeasTol 10E-5)\n"
                    "--cpx_set_int_para <parameter name> <value> Set CPLEX integer paramenter (example --cpx_set_int_para ModelSense -1)\n"
                    "--cpx_nomip tells Cplex to solve the linear relaxation only \n"
                    "--keep_tmp_files Keep the temporary mps file\n");
}

int main(int argc, char *argv[]){
char * file_model;
char * file_data;
cpx_udef * tmp;

  if (argc < 5){
    usage();
    exit(1);
  }

  printf("GLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\n\n");

  //pupulate cplex lookup table
  populate_cpx_const();

  cpx_user_def_params = NULL;
  file_model = file_data = NULL;
  CPX_out = GLPK_out = CPX_noMIP = verbose = keep_tmp_mps = 0;

  for(int i=1 ; i<=argc-1 ; i++){
    if (strcmp(argv[i],"-m")==0) file_model = argv[i+1];
    if (strcmp(argv[i],"-d")==0) file_data = argv[i+1];
    if (strcmp(argv[i],"-v")==0) verbose = 1;
    if (strcmp(argv[i],"--keep_tmp_files")==0) keep_tmp_mps = 1;
    if (strcmp(argv[i],"--glpk_out")==0) GLPK_out = 1;
    if (strcmp(argv[i],"--cpx_out")==0) CPX_out = 1;
    if (strcmp(argv[i],"--glpk_mip_gap")==0) glpk_iparm_mip_gap = atof(argv[i+1]);
    if (strcmp(argv[i],"--glpk_tol_int")==0) glpk_iparm_tol_int = atof(argv[i+1]);
    if (strcmp(argv[i],"--glpk_tol_obj")==0) glpk_iparm_tol_obj = atof(argv[i+1]);
    if (strcmp(argv[i],"--cpx_feas_tol")==0) cpx_feasibilityTol =  atof(argv[i+1]);
    if (strcmp(argv[i],"--cpx_int_tol")==0) cpx_intFeasTol =  atof(argv[i+1]);
    if (strcmp(argv[i],"--cpx_mip_gap")==0) cpx_mipGap = atof(argv[i+1]) ;
    if (strcmp(argv[i],"--cpx_set_db_para")==0) {
      cpx_udef * tmp = (cpx_udef *) malloc(sizeof(cpx_udef));
      tmp->name = argv[i+1];
      tmp->cpx_index = cplex_parameters[tmp->name];
      tmp->type = CPX_DB;
      tmp->db_value = atof(argv[i+2]);
      tmp->next = cpx_user_def_params; 
      cpx_user_def_params = tmp;
    }
    if (strcmp(argv[i],"--cpx_set_int_para")==0) {
      cpx_udef * tmp = (cpx_udef *) malloc(sizeof(cpx_udef));
      tmp->name = argv[i+1];
      tmp->cpx_index = cplex_parameters[tmp->name];
      tmp->type = CPX_INT;
      tmp->int_value = atoi(argv[i+2]);
      tmp->next = cpx_user_def_params; 
      cpx_user_def_params = tmp;
    }
    if (strcmp(argv[i],"--cpx_nomip")==0) CPX_noMIP = 1;
  } 

  if ((file_model==NULL) || (file_data == NULL)){ 
    usage();
    fprintf(stderr, "Error no model or data files provided\n");
    freeMem("Input data", 0);
  }

  printf("Wrapper parameters:\n");
  printf("  model: %s\n", file_model);
  printf("  data: %s\n", file_data);
  printf("  verbosity level: %s\n", verbose ? "quiet" : "verbose");
  printf("  keep temporary files: %s\n", keep_tmp_mps ? "on" : "off");
  printf("  GLPK parameters:\n");
  printf("    GLPK output: %s\n", GLPK_out ? "on" : "off");
  printf("    GLPK Mip Gap: %lf\n", glpk_iparm_mip_gap);
  printf("    GLPK Int Tolerance: %lf\n", glpk_iparm_tol_int);
  printf("    GLPK Objective Tolerance: %lf\n", glpk_iparm_tol_obj);
  printf("  CPLEX parameters:\n");
  printf("    CPLEX output: %s\n", CPX_out ? "on" : "off");
  printf("    CPLEX Feasibility tolerance: %lf\n", cpx_feasibilityTol);
  printf("    CPLEX Integer Feasibility tolerance: %lf\n", cpx_intFeasTol);
  if (cpx_user_def_params != NULL){
    printf("  USER DEF CPLEX parameters:\n");
    tmp = cpx_user_def_params;
    while (tmp != NULL){
      if (tmp->type == CPX_INT){
        printf("    CPLEX %s[%d] : %d\n", tmp->name, tmp->cpx_index, tmp->int_value);
      }
      else if (tmp->type == CPX_DB){
        printf("    CPLEX %s[%d] : %lf\n", tmp->name, tmp->cpx_index, tmp->db_value);
      }
      tmp = tmp->next;
    }
  }
  printf("\n\n");

  /** GLPK: Open cpx_environment **/
  mip = glp_create_prob();
  tran = glp_mpl_alloc_wksp();

  glp_term_out(GLPK_out?GLP_ON:GLP_OFF);

  /** GLPK: Read model written in MathProg **/
  ret = glp_mpl_read_model(tran, file_model, 1);

  if (ret){ 
    fprintf(stderr, "Error on translating model\n");
    freeMem("GLPK read model", ret);
  }
  
  /** GLPK: Read data for MathProg **/
  ret = glp_mpl_read_data(tran, file_data);
  if (ret){ 
    fprintf(stderr, "Error on translating data\n");
    freeMem("GLPK read data", ret);
  }

  /** GLPK: Generate model (merge data an model) **/
  ret = glp_mpl_generate(tran, NULL);
  if (ret){ 
    fprintf(stderr, "Error on generating model\n");
    freeMem("GLPK mpl generate", ret);
  }

  /** GLPK: Generate Build Model **/
  glp_mpl_build_prob(tran, mip);
          
  /** GLPK: Generate Variable indexing **/
  glp_create_index(mip);

  /** GLPK: Generate LP **/
  glp_write_mps(mip, GLP_MPS_FILE, NULL, "model.mps");


  /***********/
  /** CPLEX **/
  /***********/

  cpx_env = CPXopenCPLEX(&retCPX);
  if (retCPX || cpx_env == NULL)
  {
    fprintf(stderr, "Error: could not create cpx_environment\n");
    freeMem("CPX create environment", retCPX);
  }

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_SCRIND, CPX_out?CPX_ON:CPX_OFF);
  if (retCPX) freeMem("CPX set output", retCPX);

  /** CPLEX: Create LP model **/
  cpx_lp = CPXcreateprob(cpx_env, &retCPX, "cpx_lp");
  if (retCPX) freeMem("CPX create problem", retCPX);

  /** CPLEX: Read model **/
  retCPX = CPXreadcopyprob (cpx_env, cpx_lp, "model.mps", "MPS");
  if (retCPX) freeMem("CPX load model", retCPX);

  /** Remove utility files from disk **/
  if (!keep_tmp_mps) remove("model.mps");

  /** CPLEX: Set parameters **/

  /** CPLEX: Ask for more precision **/
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPRHS , cpx_feasibilityTol);
  if (retCPX) freeMem("CPX set dbl param", retCPX);
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPINT, cpx_intFeasTol);
  if (retCPX) freeMem("CPX set dbl param", retCPX);
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPGAP, cpx_mipGap);
  if (retCPX) freeMem("CPX set dbl param", retCPX);

  /** Set other user-defined parameters for Cplex **/
  tmp = cpx_user_def_params;
  while (tmp != NULL){
    if (tmp->type == CPX_INT){
      retCPX = CPXsetintparam(cpx_env, tmp->cpx_index, tmp->int_value);
      if (retCPX) freeMem("CPX Set udef int para", retCPX);
    }
    else if (tmp->type == CPX_DB){
      retCPX = CPXsetdblparam(cpx_env, tmp->cpx_index, tmp->db_value);
      if (retCPX) freeMem("CPX Set udef dbl para", retCPX);
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
    freeMem("CPX problem type", 0);
  }
  if (verbose) printf ("Problem type: %d\n", CPX_probType);
  CPX_IsMIP = (CPXPROB_MILP == CPX_probType);

  /** CPLEX: Solve the linear relaxation **/
  if ((CPX_noMIP) || (!CPX_IsMIP)){
    //TODO allow for other simplex optimizers
    retCPX = CPXdualopt(cpx_env, cpx_lp);
    if (retCPX) freeMem("CPX LP Optimize", retCPX);

    /** CPLEX: Retreive the optimization status **/
    retCPX = CPXgetstat(cpx_env, cpx_lp);
    if (retCPX != CPX_STAT_OPTIMAL){
      /** CPLEX: Quit in any case non optimal **/
      freeMem("Error CPX Lp optimization failed", retCPX);
    }
  }
  else{
    retCPX = CPXmipopt(cpx_env, cpx_lp);
    if (retCPX) freeMem("CPX MIP Optimize", retCPX);

    /** CPLEX: Retreive the optimization status **/
    retCPX = CPXgetstat(cpx_env, cpx_lp);
    if ((retCPX != CPXMIP_OPTIMAL) && (retCPX != CPXMIP_OPTIMAL_TOL)){
      /** CPLEX: Quit in any case non optimal **/
      freeMem("Error CPX MIP optimization failed", retCPX);
    }
  }

  /** CPLEX: Get obj function value **/
  if (!CPX_noMIP){
    retCPX = CPXgetbestobjval(cpx_env, cpx_lp, &bound);
    if (retCPX) freeMem("CPX accessing bound", retCPX);
  }

  retCPX = CPXgetobjval(cpx_env, cpx_lp, &db_tmp);
  if (retCPX) freeMem("CPX obj value", retCPX);
  
  if (verbose) printf ("Objective %lf\n", db_tmp);
  if (verbose) printf ("Best bound %lf\n", bound);
  if (verbose) printf ("Absolute gap %lf\n", fabs(db_tmp - bound));

  /** CPLEX: Get variable values **/
  for (j = 0; j < numvars; ++j){


      retCPX = CPXgetx (cpx_env, cpx_lp, &db_tmp, j, j);
      if (retCPX) freeMem("CPX Get var value", retCPX);


      int surplus, colnamespace; 
      char ** colname = NULL;
      char * colnamestore = NULL;
      
      retCPX = CPXgetcolname(cpx_env, cpx_lp, NULL, NULL, 0, &surplus, j, j);
      if (( retCPX != CPXERR_NEGATIVE_SURPLUS ) && ( retCPX != 0 ))  {
        freeMem("CPX Get var names", retCPX);
      }
      
      colnamespace = - surplus;
      if ( colnamespace > 0 ) {
           
           colname = (char **) malloc (sizeof(char *));
           colnamestore = (char *)  malloc(colnamespace);
           
           if ( colname == NULL || colnamestore == NULL ) {
            freeMem("CPX Allocating memory for col name" , 0);
           }
           
           retCPX = CPXgetcolname (cpx_env, cpx_lp, colname, colnamestore, colnamespace, &surplus, j, j);
           if ( retCPX ) {
              freeMem("CPX Get final var names", retCPX);
           }
      }
      else {
        freeMem("CPX no name associated", 0);
      }
  
      if (verbose) printf ("Processed variable %d name %s\n", j, colname[0]);
      
      sprintf(nameCPX, "%s", colname[0]);
	  
	    free(colnamestore);
      free(colname);

    char type[1];
    retCPX = CPXgetctype (cpx_env, cpx_lp, type, j,j);
    if (retCPX) freeMem("CPX Accessing variable type", retCPX);

    /** GLPK search variable index by name **/
    col_index = glp_find_col(mip, nameCPX);

    if (col_index != 0){ 
      /** GLPK set variable bounds **/
      if ((type[0] == 'B') || (type[0] == 'I')){
        if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", nameCPX, type[0], db_tmp, round(db_tmp));
        glp_set_col_bnds(mip, col_index, GLP_FX, round(db_tmp), round(db_tmp));
      }
      else{
        if (verbose) printf ("Variable %s is of type %c value %lf fixed to %lf\n", nameCPX, type[0], db_tmp, db_tmp);
        glp_set_col_bnds(mip, col_index, GLP_FX, db_tmp, db_tmp);
      }
    }
  }

  if ((CPX_IsMIP) && (!CPX_noMIP)){
          
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
      freeMem("GLPK Error re-optimization", ret);
    }

    ret = glp_mip_status(mip);
    switch (ret){
            case GLP_OPT:
            break;
            case GLP_FEAS:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_FEAS, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_NOFEAS:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_NOFEAS, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_UNDEF:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNDEF, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
    }
  }
  else{

    /** GLPK initialize parameters **/
    parm = (glp_smcp*) malloc(sizeof(glp_smcp));
    glp_init_smcp(parm);
    parm->meth = GLP_DUALP;
    parm->tol_bnd = 10E-4;
    parm->tol_dj = 10E-4;
  
    /** GLPK get the optimal basis **/
    //ret = glp_simplex(mip, parm);
    if (ret){ 
     fprintf(stderr, "glp_simplex, Error on optimizing the model : %d \n", ret);
     freeMem("GLPK Error re-optimization", ret);
    }
    ret = glp_get_status(mip);
    switch (ret){
            case GLP_OPT:
            break;
              case GLP_FEAS:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_FEAS, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_INFEAS:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_INFEAS, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_NOFEAS:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_NOFEAS, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_UNBND:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNBND, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
            case GLP_UNDEF:
              fprintf(stderr, "Error GLPK simplex is not optimal, GLP_UNDEF, code %d\n", ret);
              freeMem("GLPK Error re-optimization", ret);
    }
  }

  /** GLPK: Perform postprocessing **/
  ret = glp_mpl_postsolve(tran, mip, GLP_MIP);
  if (ret != 0) fprintf(stderr, "Error on postsolving model\n");

  /** GLPK: free structures **/
  if (tran) glp_mpl_free_wksp(tran);
  if (mip) glp_delete_prob(mip);

  /** CPLEX: free structures **/
  if (cpx_lp) CPXfreeprob(cpx_env, &cpx_lp);
  if (cpx_env) CPXcloseCPLEX(&cpx_env);

  printf("Done.\nGLPK -> CPLEX -> GLPK wrapper v0.2 (2010)\n");

  exit(0);
}
