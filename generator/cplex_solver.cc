#include "main.h"
#include "generator.h"
#include "solver.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "matheuristic.h"
#include "cplex_solver.h"
#include "rnpsolution.h"

extern param prob;

ofstream of_progress;

typedef struct timeliminfo CPX_TIMELIMINFO, *CPX_TIMELIMINFOptr;
struct timeliminfo {
  double 		timestart;
  double 		timelim;
  double 		acceptablegap;
  int    		aborted;
  int 			numvars;
  bool 			output_progress;
  double 		last_progress_report;
  double 		progress_report_interval;
  CPXsimpleCallback 	infocallbackfunc;
  CPXsimpleCallback 	heuristiccallbackfunc;
  CPXsimpleCallback 	incumbentcallbackfunc;

  CPX_VARIABLESINFO* 	varsinfo;
  bool 			has_incumbent;
  double 		bestobj;
  LpSolutionPtr 	incumbent;
};


typedef struct {
  double 		timenow;
  bool 			has_incumbent;
  double 		bestobj;
  double 		gap;
  LpSolutionPtr 	incumbent;
} InfoCallbackData;


typedef struct {
  double 		timenow;
  bool 			has_incumbent;
  double 		bestobj;
  double 		gap;
  LpSolutionPtr 	incumbent;
  LpSolutionPtr 	heuristic;
} HeuristicCallbackData;

typedef struct {
  double 		timenow;
  bool 			has_incumbent;
  double 		bestobj;
  double 		gap;
  LpSolutionPtr 	incumbent;
} IncumbentCallbackData;

Matheuristic *MyCplexSolver::m_mathh = NULL;


static LpSolutionPtr getIncumbent(double *vals, int numvars, CPX_VARIABLESINFO &varsinfo, double objval);

  static int CPXPUBLIC
timelimcallback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{
  int status = 0;
  double gap = 0.0;
  double timenow = 0.0;
  int solpooln = 0;
  bool TRY_SOL_POOL = false;


  CPX_TIMELIMINFOptr info = (CPX_TIMELIMINFOptr) cbhandle;
  int hasincumbent = 0;


  status = CPXgetcallbackinfo (env, cbdata, wherefrom,
			       CPX_CALLBACK_INFO_MIP_FEAS, &hasincumbent);
  if ( status )  goto TERMINATE;

  //printf("timelimcallback hasincumbent %d\n", hasincumbent);
  if( TRY_SOL_POOL)
  {
    CPXCLPptr lp_ptr;
    status = CPXgetcallbackinfo( env, cbdata, wherefrom, CPX_CALLBACK_INFO_USER_PROBLEM,
				 &lp_ptr);
    if ( status )  {
      fprintf(stderr, "getcallbackinfo LP failed\n");
      goto TERMINATE;
    }
    solpooln = CPXgetsolnpoolnumsolns(env, lp_ptr);
    printf("MIP has %d solutions in SOLUTION POOL\n", solpooln);
  }
  if( hasincumbent)
  {
    double bestobj=0.0;
    status = CPXgetcallbackinfo(env, cbdata, wherefrom, 
				CPX_CALLBACK_INFO_BEST_INTEGER, 
				&bestobj);
    if ( status )  goto TERMINATE;
    /// Determine if we need to load incumbent
    bool load_incumbent = false;
    if( !info->has_incumbent )
      load_incumbent = true;
    if( info->has_incumbent)
    {
      //printf("timelimcallback best_obj %f info->bestobj %f\n",
      //     bestobj, info->bestobj);
      if( fabs(bestobj - info->bestobj) > EPSILON)
	load_incumbent = true;
    }
    /// If incumbent callback is enabled, this is never executed
    /// RECALL: incumbent callback is called before accepting incumbent
    if( load_incumbent)
    {
      //printf("timelimcallback loading incumbent\n");
      int numvars = info->numvars;
      double *vals = (double*)malloc(numvars*sizeof(double));
      status = CPXgetcallbackincumbent(env,cbdata,wherefrom, vals, 0, numvars-1); 
      if(status)
      { 
	fprintf(stderr, "CPXgetcallback incumbent failed %d\n", status);
	goto TERMINATE;
      }

      LpSolutionPtr incumbent = getIncumbent(vals, numvars,(*info->varsinfo), bestobj);
      info->bestobj = incumbent->objval;
      info->incumbent = incumbent;
      info->has_incumbent = true;
      free(vals);
    }
  }

  if ( !info->aborted ) {
    status = CPXgetcallbackinfo (env, cbdata, wherefrom,
				 CPX_CALLBACK_INFO_MIP_REL_GAP, &gap);
    if ( status )  goto TERMINATE;

    /* Turn the gap into a percentage */
    gap *= 100.0;

    status = CPXgettime (env, &timenow);
    if ( status )  goto TERMINATE;

    if ( info->has_incumbent)
    {
      if( info->timelim > 0)
      {

	if ( timenow - info->timestart > info->timelim  &&
	     gap < info->acceptablegap                    ) {

	  fprintf (stderr,
		   "Good enough solution at time %.2fsec, gap = %g%%, obj %.2f\n",
		   timenow - info->timestart, gap, info->bestobj);

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

	}
      }
    }

    // Check if we should report solving progress
    // If that's the case, the ofstream should be 
    // initialized and ready
    if(info->output_progress)
    {
      if( (timenow - info->last_progress_report) 
	  > info->progress_report_interval)
      {
	info->last_progress_report = timenow;
	of_progress << (timenow - info->timestart);
	if( hasincumbent)
	{
	  of_progress << " " << info->bestobj << " " << gap;
	}
	else
	{
	  of_progress << " NULL";
	}
	of_progress << endl;
      }
    }



  }
  return info->aborted;

TERMINATE:
  fprintf(stderr,"timelimcallback failed with status %d\n",status);
  return (status);

} /* END timelimcallback */

  static int CPXPUBLIC 
heuristiccallback(CPXCENVptr env, void *cbdata, int wherefrom,
		  void *cbhandle, double *objval_p, double *x,
		  int *checkfeas_p, int *useraction_p)
{

  int status = 0;
  double bestobj = 0.0;

  double gap = 0.0;
  double timenow = 0.0;
  printf("Heuristic callback called\n");

  int hasincumbent;
  // Call use-defined callback

  CPX_TIMELIMINFOptr info = (CPX_TIMELIMINFOptr) cbhandle;
  HeuristicCallbackData cbinfo;

  *useraction_p = CPX_CALLBACK_DEFAULT;

  if ( !info->aborted ) {
    status = CPXgetcallbackinfo (env, cbdata, wherefrom,
				 CPX_CALLBACK_INFO_MIP_FEAS, &hasincumbent);
    if ( status )  goto TERMINATE;
    info->has_incumbent = hasincumbent;

    status = CPXgetcallbackinfo (env, cbdata, wherefrom,
				 CPX_CALLBACK_INFO_MIP_REL_GAP, &gap);
    if ( status )  goto TERMINATE;
    /* Turn the gap into a percentage */
    gap *= 100.0;
    //info->gap = gap;


    status = CPXgettime (env, &timenow);
    if ( status )  goto TERMINATE;


    if( info->heuristiccallbackfunc)
    {

      cbinfo.has_incumbent = hasincumbent;
      cbinfo.timenow = timenow - info->timestart;
      cbinfo.gap = gap;
      cbinfo.bestobj = info->bestobj;
      cbinfo.incumbent = info->incumbent;

      bool gotone = (*info->heuristiccallbackfunc)(&cbinfo);
      if(gotone)
      {
	printf("Heuristic callback pulled sol %.2f\n",
	       cbinfo.heuristic->objval);
	// Filling in pulled solution
	*objval_p = cbinfo.heuristic->objval;
	int numvars = info->numvars;
	CPX_VARIABLESINFO &varsinfo = (*info->varsinfo);
	for(int i=0; i< numvars; i++)
	{
	  std::pair<string, char>  varinfo = varsinfo[i];
	  string nameVar = varinfo.first;
	  char type = varinfo.second;
	  double val = cbinfo.heuristic->value[nameVar];
	  x[i] = val;
	}
	// Check feasibility, just in case
	*checkfeas_p = 0;
	*useraction_p = CPX_CALLBACK_SET;
      }
    }
  }
  return (status);
TERMINATE:
  return status;
}

  static LpSolutionPtr 
getIncumbent(double *vals, int numvars, CPX_VARIABLESINFO &varsinfo, double objval)
{

  //double *vals = (double*)malloc(numvars*sizeof(double));
  //status = CPXgetcallbackincumbent(env,cbdata,wherefrom, vals, 0, numvars-1); 
  //if(status)
  //{
  //	printf("Error CPXgetcallbackincumbent %d\n", status);
  //}
  //else
  //{
  //	printf("Good CPXgetcallbackincumbent\n");
  //}
  // Get incumbent 
  LpSolutionPtr incumbent(new LpSolution());

  int validvars = 0;
  for(int i=0; i< numvars; i++)
  {
    std::pair<string, char>  varinfo = varsinfo[i];
    string nameVar = varinfo.first;
    char type = varinfo.second;
    if( fabs(vals[i]) > EPSILON)
    {
      validvars++;
      incumbent->value[nameVar] = vals[i];
      incumbent->type[nameVar] = type;
    }
  }

  incumbent->is_optimal(false);
  incumbent->solved = true;
  incumbent->objval = objval;
  printf("Incumbent solution set. valid vars %d\n", validvars);
  return incumbent;
}


  static int CPXPUBLIC 
incumbentcallback(CPXCENVptr env, void *cbdata, int wherefrom,
		  void *cbhandle, double objval_p, double *x,
		  int *checkfeas_p, int *useraction_p)
{

  int status = 0;

  double gap = 0.0;
  double timenow = 0.0;

  printf("incumbent callback called. obj_val %.2f\n", objval_p);
  int hasincumbent;
  // Call use-defined callback

  CPX_TIMELIMINFOptr info = (CPX_TIMELIMINFOptr) cbhandle;

  //	status = CPXgetcallbackinfo(env, cbdata, wherefrom, 
  //				    CPX_CALLBACK_INFO_BEST_INTEGER, &bestobj
  //

  IncumbentCallbackData cbinfo;
  int numvars = info->numvars;
  double *vals = x;
  LpSolutionPtr incumbent = getIncumbent(vals, numvars,(*info->varsinfo), objval_p);
  //double *vals = (double*)malloc(numvars*sizeof(double));
  //status = CPXgetcallbackincumbent(env,cbdata,wherefrom, vals, 0, numvars-1); 
  //if(status)
  //{
  //	printf("Error CPXgetcallbackincumbent %d\n", status);
  //}
  //else
  //{
  //	printf("Good CPXgetcallbackincumbent\n");
  //}
  // Get incumbent 
  /*
    LpSolutionPtr incumbent(new LpSolution());
    CPX_VARIABLESINFO &varsinfo = (*info->varsinfo);

    int validvars = 0;
    for(int i=0; i< numvars; i++)
    {
    std::pair<string, char>  varinfo = varsinfo[i];
    string nameVar = varinfo.first;
    char type = varinfo.second;
    if( fabs(vals[i]) > EPSILON)
    {
    validvars++;
    incumbent->value[nameVar] = vals[i];
    incumbent->type[nameVar] = type;
    }
    }

    incumbent->is_optimal(false);
    incumbent->solved = true;
    incumbent->objval = objval_p;
    info->incumbent = incumbent;
    info->has_incumbent = true;
    printf("Incumbent solution set. valid vars %d\n", validvars);
    */

  // Here we set the incumbent objval				    
  info->bestobj = incumbent->objval;
  info->incumbent = incumbent;
  info->has_incumbent = true;


  if(info->incumbentcallbackfunc != NULL)
  {
    printf("Preparing incumbemtcallbackfunc\n");
    cbinfo.has_incumbent = true;
    cbinfo.timenow = timenow - info->timestart;
    cbinfo.gap = gap;
    cbinfo.bestobj = objval_p;
    cbinfo.incumbent = incumbent;
    (*info->incumbentcallbackfunc)(&cbinfo);

  }



  return status;
}
CPX_VARIABLESINFO *MyCPXgetVariablesInfo(CPXENVptr env, CPXLPptr cpx_lp)
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



void my_free_cplex(char * context, int value, CPXENVptr cpx_env,CPXLPptr cpx_lp){
  //cpx_udef * tmp;


  /** GLPK: free user defined parameters **/
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

    fprintf(stderr,"Context: %s, Value: %d.\n", context, value);
  }



}

void MyCplexSolver::addEdgeVar(string i, string j)
{

  int retCPX;
  int ccnt;
  char ctype;
  char **colname;
  colname = (char **)malloc(1*sizeof(char *));
  colname[0] = (char *)malloc(16*sizeof(char));	
  //	char colname[1][16];

  int ix = CPXgetnumcols(m_cpxenv,m_cpxlp);
  if(prob.ilp_model.integer_flow)
    ctype = 'I';
  else
    ctype = 'C';
  ccnt = 1;
  sprintf(colname[0], "x[%s,%s]",i.c_str(),j.c_str());

  //	printf("Adding variable %s\n",colname[0]);
  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
  if (retCPX) my_free_cplex("CPX create problem", retCPX, m_cpxenv, m_cpxlp);
  m_varX[make_pair(i,j)]=ix;

  ctype = 'B';
  ccnt = 1;
  sprintf(colname[0], "bx[%s,%s]",i.c_str(),j.c_str());

  //printf("Adding variable %s\n",colname[0]);
  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
  if (retCPX) my_free_cplex("CPX create problem", retCPX, m_cpxenv, m_cpxlp);
  m_varBX[make_pair(i,j)]=ix+1;

  free(colname[0]);
  free(colname);

}



void 
MyCplexSolver::addRelayVar(string rid)
{

  int retCPX;
  int ccnt;
  char ctype;
  char **colname;
  colname = (char **)malloc(1*sizeof(char *));
  colname[0] = (char *)malloc(16*sizeof(char));	

  int ix = CPXgetnumcols(m_cpxenv,m_cpxlp);

  ctype = 'B';
  ccnt = 1;
  sprintf(colname[0], "y[%s]",rid.c_str());
  //printf("Adding variable %s\n",colname[0]);
  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
  if (retCPX) my_free_cplex("CPX create problem", retCPX, m_cpxenv, m_cpxlp);
  m_varY[rid] = ix;

  free(colname);

}

void 
MyCplexSolver::addOtherVars()
{
  int retCPX;
  int ccnt;
  char ctype;
  char **colname;
  colname = (char **)malloc(1*sizeof(char *));
  colname[0] = (char *)malloc(16*sizeof(char));	

  int ix = CPXgetnumcols(m_cpxenv,m_cpxlp);

  if( prob.ilp_model.allow_incomplete_delivery )
  {
    /// we add a flow variable from every SN/RN to a fake sink, 
    /// thus, all data that can not be delivered is going to that sink
    ///BUG ALERT! Do not add a link to DUMP for RNs because it might
    /// cause some problems - in any case, it it enough for SNs
    ctype = 'C';
    for(int i=0; i< m_netWithRelays->size();i++)
    {	
      Node &ni = m_netWithRelays->getNode(i);
      if( ni.t == BASE || ni.t == RELAY)
	continue;
      ccnt = 1;
      sprintf(colname[0], "x[%s,DUMP]",ni.id.c_str());

      //	printf("Adding variable);
      retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
      if (retCPX) my_free_cplex("CPX create problem", retCPX, m_cpxenv, m_cpxlp);
      m_varXD[make_pair(ni.id,"DUMP")]=ix++;
    }
/*
    sprintf(colname[0], "D");
    //printf("Adding variable %s\n",colname[0]);
    retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
    if (retCPX) my_free_cplex("CPX new cols", retCPX, m_cpxenv, m_cpxlp);
    m_varD = ix;
    printf("added D\n");
    fflush(stdout);
    */
  }

  if( prob.ilp_model.use_flow_neighbor_sink )
    {
      for(int i=0; i< m_netWithRelays->size();i++)
	{	
	  Node &ni = m_netWithRelays->getNode(i);
	  if( ni.t == BASE )
	    continue;
	  ccnt = 1;
	  sprintf(colname[0], "x[%s,FDUMP]",ni.id.c_str());

	  //	printf("Adding variable);
	  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
	  if (retCPX) my_free_cplex("CPX adding variable FDUMP", retCPX, m_cpxenv, m_cpxlp);
	  m_varXF[make_pair(ni.id,"FDUMP")]=ix++;
	}     
    }

  free(colname);
}

void 
MyCplexSolver::addLFVar(string id)
{
  int retCPX;
  int ccnt;
  char ctype;
  char **colname;
  colname = (char **)malloc(1*sizeof(char *));
  colname[0] = (char *)malloc(16*sizeof(char));	

  int ix = CPXgetnumcols(m_cpxenv,m_cpxlp);

  ctype = 'C';
  ccnt = 1;
  sprintf(colname[0], "lf[%s]",id.c_str());
  //printf("Adding variable %s\n",colname[0]);
  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
  if (retCPX) my_free_cplex("CPX new cols", retCPX, m_cpxenv, m_cpxlp);
  m_varLF[id] = ix;

  free(colname);
}

int 
MyCplexSolver::getLFVarIndex(string i)
{
  ITERATOR(m_varLF) it = m_varLF.find(i);
  assert(it != m_varLF.end());
  return it->second;
}

void 
MyCplexSolver::addFXVar(string id)
{
  int retCPX;
  int ccnt;
  char ctype;
  char **colname;
  colname = (char **)malloc(1*sizeof(char *));
  colname[0] = (char *)malloc(16*sizeof(char));	

  int ix = CPXgetnumcols(m_cpxenv,m_cpxlp);

  ctype = 'B';
  ccnt = 1;
  sprintf(colname[0], "fx[%s]",id.c_str());
  //printf("Adding variable %s\n",colname[0]);
  retCPX = CPXnewcols (m_cpxenv, m_cpxlp, ccnt, NULL, NULL, NULL,&ctype, colname);
  if (retCPX) my_free_cplex("CPX new cols", retCPX, m_cpxenv, m_cpxlp);
  m_varFX[id] = ix;

  free(colname);
}

int 
MyCplexSolver::getFXVarIndex(string i)
{
  ITERATOR(m_varFX) it = m_varFX.find(i);
  assert(it != m_varFX.end());
  return it->second;
}

int 
MyCplexSolver::getYVarIndex(string i)
{
  ITERATOR(m_varY) it = m_varY.find(i);
  assert(it != m_varY.end());
  return it->second;
}

int 
MyCplexSolver::getXVarIndex(string i, string j)
{
  if( j == "DUMP")
  {
    ITERATOR(m_varXD) it = m_varXD.find( make_pair(i, "DUMP"));
    assert(it != m_varXD.end());
    return it->second;
  }
  if( j == "FDUMP" )
    {
      ITERATOR(m_varXF) it = m_varXF.find( make_pair(i, "FDUMP"));
      assert(it != m_varXF.end());
      return it->second;      
    }
  //  printf("Retrieving Var X[%s,%s] index\n",i.c_str(),j.c_str());
  ITERATOR(m_varX) it = m_varX.find(make_pair(i,j));
  assert(it != m_varX.end());
  return it->second;
}


int 
MyCplexSolver::getBXVarIndex(string i, string j)
{
  ITERATOR(m_varBX) it = m_varBX.find(make_pair(i,j));
  assert(it != m_varBX.end());
  return it->second;
}

/* 
 * int CPXaddrows(CPXCENVptr env, CPXLPptr lp, int ccnt, int rcnt, int nzcnt, const
 * double * rhs, const char * sense, const int * rmatbeg, const int * rmatind, const
 * double * rmatval, char ** colname, char ** rowname)
 * */

/*
 * s.t. flow{k in R}: sum{(k,j) in edges} x[k,j] = sum{(i,k) in edges} x[i,k];
 * */
void 
MyCplexSolver::addFlowCons()
{
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(16*sizeof(char));
  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc((m_varX.size()+m_varXD.size())*sizeof(int));
  double *rmatval = (double *)malloc((m_varX.size()+m_varXD.size())*sizeof(double));

  printf("Adding %d flow constraints\n", m_netWithRelays->numRelays());
  fflush(stdout);
  for(int i=0; i< m_netWithRelays->numRelays();i++)
  {
    Node &nr = m_netWithRelays->getRelay(i);
    rmatbeg[0]=0;
    nzcnt=0;
    for(int j=0; j< m_netWithRelays->size(); j++)
    {
      Node &nn = m_netWithRelays->getNode(j);
      if( m_netWithRelays->inRange(nr,nn) && (nr.id != nn.id))
      {
	rmatind[nzcnt] = getXVarIndex(nr.id, nn.id);
	rmatval[nzcnt++] = 1.0;
	/*
	  rmatind[nzcnt] = getXVarIndex(nn.id, nr.id);
	  rmatval[nzcnt++] = -1.0;
	  */
      }
    }

    for(int j=0; j< m_netWithRelays->size(); j++)
    {
      Node &nn = m_netWithRelays->getNode(j);
      if( m_netWithRelays->inRange(nn,nr) && (nr.id != nn.id))
      {
	rmatind[nzcnt] = getXVarIndex(nn.id, nr.id);
	rmatval[nzcnt++] = -1.0;
      }
    }
    if( prob.ilp_model.use_flow_neighbor_sink)
      {
	//! 
	rmatind[nzcnt] = getXVarIndex(nr.id, "FDUMP");
	rmatval[nzcnt++] = 1.0;
	//!
      }
 

//    if( prob.ilp_model.allow_incomplete_delivery )
//    {
//      rmatind[nzcnt] = getXVarIndex(nr.id, "DUMP");
//      rmatval[nzcnt++] = 1.0;
//    }
    ccnt = 0;
    rcnt = 1;
    rhs = 0.0;
    sense = 'E';
    colname = NULL;
    sprintf(rowname[0],"flow[%s]",nr.id.c_str());
    printf("Adding row %s\n", rowname[0]);
    fflush(stdout);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}


void 
MyCplexSolver::addFlowEdgeCons()
{
  /* s.t. flowedge{(i,j) in edges} : x[i,j] <= MAX_CAPACITY * bx[i,j];
   * */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  double *rhs = (double *)malloc(m_varX.size()*sizeof(double));;
  char *sense = (char *)malloc(m_varX.size()*sizeof(char));
  char ** colname;
  char **rowname;
  rowname = (char **)malloc(m_varX.size()*sizeof(char *));
  for( int i=0; i<m_varX.size();i++)
    rowname[i] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(m_varX.size()*sizeof(int));
  int *rmatind = (int *)malloc(2*m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(2*m_varX.size()*sizeof(double));


  printf("Adding %d flowedge constraints\n", m_varX.size());
  fflush(stdout);
  rcnt=0;
  nzcnt=0;
  FOREACH(it, m_varX)
  {
    rmatbeg[rcnt] = nzcnt;
    const pair<string,string> &p = it->first;
    rmatind[nzcnt] = it->second;
    rmatval[nzcnt] = 1.0;
    nzcnt++;

    rmatind[nzcnt] = getBXVarIndex(p.first, p.second);
    rmatval[nzcnt] = -1.0 * prob.ilp_model.max_capacity;
    nzcnt++;
    rhs[rcnt] = 0.0;
    sense[rcnt] = 'L';
    sprintf(rowname[rcnt],"flowedge[(%s,%s)]",p.first.c_str(),p.second.c_str());
    rcnt++;

  }
  printf("adding ...\n");
  fflush(stdout);
  colname = NULL;
  ccnt = 0;
  retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
		      rhs,sense, rmatbeg,rmatind,rmatval,colname,rowname);
  if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);

  free(rhs);
  free(sense);
  for(int i=0;i<m_varX.size();i++)
    free(rowname[i]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);


}

void MyCplexSolver::addNodeDegreeCons()
{
  /*
   * s.t. nodedegree{k in R union S}: (sum{(j,k) in edges} bx[j,k]) <= D;
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));

  //printf("Adding %d node degree constraints\n", 
  //       m_netWithRelays->numStatic() + m_netWithRelays->numBases());

  for(int i=0; i< nnodes; i++)
  {
    Node &ni = m_netWithRelays->getNode(i);
    rmatbeg[0]=0;
    nzcnt=0;
    if( ni.t != BASE)
    {
      for(int j=0; j< nnodes; j++)
      {
	Node &nj = m_netWithRelays->getNode(j);

	if( m_netWithRelays->inRange(nj,ni) && i!=j)
	{
	  rmatind[nzcnt] = getBXVarIndex(nj.id, ni.id);
	  rmatval[nzcnt++] = 1.0;

	}
      }
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_degree;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"nodedeg[%s]",ni.id.c_str());
      //printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);

    }
  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}


void 
MyCplexSolver::addFlowNeighborCons()
{
  /*
   * s.t. flowneighbour{i in  S}: 
   * (sum{ (i,j) in edges diff edgesR} sum{(j,k) in edges diff edgesR} x[j,k]) - MAX_LOC_FLOW <= fx[i] * MAX_FLOW;
   *
   * s.t. flowneighbor2{i in  S}: 
   * ((sum{ (i,j) in edges diff edgesR} sum{(j,k) in edges diff edgesR} x[j,k]) - MAX_LOC_FLOW) >= (fx[i] -1)*MAX_FLOW;
   */

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  int nstatic = m_netWithRelays->numStatic();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));

  printf("Adding %d flowneighbor constraints\n", nstatic*2);
  printf("max_flow %f  max_local_flow %f\n", 
	 prob.ilp_model.max_flow, 
	 prob.ilp_model.max_local_flow);
  fflush(stdout);

  for(int i=0; i< nstatic; i++)
    {
      //      printf("flown[%d] ",i);
      Node &ni = m_netWithRelays->getStatic(i);
      rmatbeg[0]=0;
      nzcnt=0;
      for(int j=0; j< nnodes; j++)
	{
	  Node &nj = m_netWithRelays->getNode(j);

	  ///  BUG ALERT!: why considering only the static for flow neighbor constraints ? 
	  //      if(nj.t != RELAY && ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	  if(ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	    {
	      rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	      rmatval[nzcnt++] = 1.0;
		      
	      for(int k=0; k<nnodes; k++)
		{
		  Node &nk = m_netWithRelays->getNode(k);
		  ///  BUG ALERT!: why considering only the static for flow neighbor constraints ?
		  //	  if( nk.t != RELAY && nj.id != nk.id && m_netWithRelays->inRange(nj,nk))
		  if( nj.id != nk.id && m_netWithRelays->inRange(nj,nk))	  
		    {

		      rmatind[nzcnt] = getXVarIndex(nj.id, nk.id);
		      rmatval[nzcnt++] = 1.0;

		      //          WE FOUND A BUG HERE! 12.02.2015
		      //	    rmatind[nzcnt] = getFXVarIndex(ni.id);
		      //	    rmatval[nzcnt++] = -prob.ilp_model.max_flow;
		    }
		}
	    }
	}

      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -prob.ilp_model.max_flow;
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"flown[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);


      rmatbeg[0]=0;
      nzcnt=0;
      for(int j=0; j< nnodes; j++)
	{
	  Node &nj = m_netWithRelays->getNode(j);
	  /// BUG ALERT: see above
	  //      if(nj.t != RELAY && ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	  if(ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	    {
	      rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	      rmatval[nzcnt++] = 1.0;

	      for(int k=0; k<nnodes; k++)
		{
		  Node &nk = m_netWithRelays->getNode(k);
		  /// BUG ALERT: see above 
		  //	  if( nk.t != RELAY && nj.id != nk.id && m_netWithRelays->inRange(nj,nk))
		    if( nj.id != nk.id && m_netWithRelays->inRange(nj,nk))
		    {

		      rmatind[nzcnt] = getXVarIndex(nj.id, nk.id);
		      rmatval[nzcnt++] = 1.0;

		      //WE FOUND A VERY OLD BUG HERE :(
		      //rmatind[nzcnt] = getFXVarIndex(ni.id);
		      //rmatval[nzcnt++] = -prob.ilp_model.max_flow;
		    }
		}
	    }
	}
      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -prob.ilp_model.max_flow;
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow - prob.ilp_model.max_flow;
      sense = 'G';
      colname = NULL;
      sprintf(rowname[0],"flown2[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
    }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);

}

void
MyCplexSolver::addStrictFlowNeighborCons()
{
  /*
   *   fx[i] <= 0 
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  int nstatic = m_netWithRelays->numStatic();
  int nrelays = m_netWithRelays->numRelays();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(2*sizeof(int));
  double *rmatval = (double *)malloc(2*sizeof(double));
      
  FOREACH(it, m_varLF)
    {
      rmatbeg[0]=0;
      nzcnt=0;
      string nid = it->first;
      int ix = it->second;
      
      rmatind[nzcnt] = ix;
      rmatval[nzcnt++] = 1.0;
      //
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow+1;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"strictflown[%s]", nid.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
    }
 
}

void
MyCplexSolver::addImFlowNeighborCons()
{
  /*
   * s.t. localflow{i in  S U R}: 
   * lf[i] = sum{j in (i,j)} sum{ k in (j,k), k!= j} x[j,k] + sum{j in (i,j)} x[i,j]
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  int nstatic = m_netWithRelays->numStatic();
  int nrelays = m_netWithRelays->numRelays();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));

  printf("Adding %d improved flow neighbo constraints\n", nstatic*2);

  for(int i=0; i< nstatic; i++)
    {
      Node &ni = m_netWithRelays->getStatic(i);
      ////////
      //! if LF[i] >= max_local_flow then FX[i] = 1
      rmatbeg[0]=0;
      nzcnt=0;
      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = 1.0;
      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -prob.ilp_model.max_flow;
      //! 
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"flown1[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);

      ///////
      //! if FX[i] == 1 then LF[i] >= max_local_flow
      rmatbeg[0]=0;
      nzcnt=0;
      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = 1.0;
      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -prob.ilp_model.max_flow;
      //! 
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow - prob.ilp_model.max_flow;
      sense = 'G';
      colname = NULL;
      sprintf(rowname[0],"flown2[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
      
    }

  ///! now for relays
  for(int i=0; i< nrelays; i++)
    {
      Node &ni = m_netWithRelays->getRelay(i);
      ////////
      //! if LF[i] >= max_local_flow + (1-y[i])*max_flow then FX[i] = 1
      //! LF[i] - max_local_flow + y[i]*max_flow - max_flow <= FX[i]*max_flow
      //! LF[i] + y[i]*max_flow - FX[i]*max_flow <= max_flow + max_local_flow 
      rmatbeg[0]=0;
      nzcnt=0;
      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = 1.0;
      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -prob.ilp_model.max_flow;
      rmatind[nzcnt] = getYVarIndex(ni.id);
      rmatval[nzcnt++] = prob.ilp_model.max_flow;
      //! 
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow + prob.ilp_model.max_flow;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"flown1[%s]",ni.id.c_str());
      printf("Adding row %s: lf[%s] + %ffx[%s] + %fy[%s] < %f \n",
	     rowname[0], ni.id.c_str(), -prob.ilp_model.max_flow,
	     ni.id.c_str(), prob.ilp_model.max_flow,ni.id.c_str(),
	     rhs);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);

      ///////
      //! if FX[i] == 1 then LF[i] >= max_local_flow + (1-y[i])*max_flow
      //! note: we must use 2*max_flow for big-m here, o.w. it could be infeasible in the free case
      //! LF[i] - max_local_flow - (1-y[i])*max_flow >= (FX[i]-1)*2*max_flow
      //! LF[i] + y[i]*max_flow - FX[i]*2*max_flow >= -max_flow + max_local_flow
      rmatbeg[0]=0;
      nzcnt=0;
      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = 1.0;
      rmatind[nzcnt] = getFXVarIndex(ni.id);
      rmatval[nzcnt++] = -2*prob.ilp_model.max_flow;
      rmatind[nzcnt] = getYVarIndex(ni.id);
      rmatval[nzcnt++] = prob.ilp_model.max_flow;
      //! 
      ccnt = 0;
      rcnt = 1;
      rhs = prob.ilp_model.max_local_flow - prob.ilp_model.max_flow;
      sense = 'G';
      colname = NULL;
      sprintf(rowname[0],"flown2[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
      
    }
  

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);  
}

void
MyCplexSolver::addSinkFlowCons()
{

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));
  ///
  double allD = 0.0;
  for(int i=0; i< m_netWithRelays->numStatic();i++)
    {
      allD += m_netWithRelays->getDemand(i);
    }
  
  FOREACH(it, m_varXF)
    {
      string nid = (it->first).first;
      rmatbeg[0]=0;
      nzcnt=0;
      rmatind[nzcnt] = it->second;
      rmatval[nzcnt++] = 1.0;
      rmatind[nzcnt] = getFXVarIndex(nid);
      rmatval[nzcnt++] = -allD;
      //
      ccnt = 0;
      rcnt = 1;
      rhs = 0;
      sense = 'L';
      colname = NULL;
      sprintf(rowname[0],"sinkflow[%s]",nid.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
      
    }
}

void
MyCplexSolver::addLocalFlowCons()
{
  /*
   * s.t. localflow{i in  S U R}: 
   * lf[i] = sum{j in (i,j)} sum{ k in (j,k), k!= j} x[j,k] + sum{j in (i,j)} x[i,j]
   */

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  int nstatic = m_netWithRelays->numStatic();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc((m_varX.size()+1)*sizeof(int));
  double *rmatval = (double *)malloc((m_varX.size()+1)*sizeof(double));
  
  //bool verbose = false;
  
  printf("Adding %d local flow constraints\n", nnodes);
  fflush(stdout);

  for(int i=0; i< nnodes; i++)
    {
      Node &ni = m_netWithRelays->getNode(i);
      //TODO shall we include BSs here?
      if( ni.t == BASE )
	continue;
      rmatbeg[0]=0;
      nzcnt=0;
      //      printf("ni = %s\n", ni.id.c_str());
      //fflush(stdout);

      /// first,compute neighborhod
      std::map<int, int > neigh;
      neigh[i] = 0;
      if( prob.ilp_model.im_flow_neighbor_use_interference_range )
	{
	  for( int j=0; j< nnodes; j++)
	    {
	      if( j==i )
		continue;
	      Node &nj = m_netWithRelays->getNode(j);
	      if (nj.t == BASE )
		continue;
	      double dist = m_netWithRelays->distanceLink(ni,nj);
	      if( dist > prob.ilp_model.im_flow_neighbor_interference_range )
		continue;
	      neigh[j]=1;
	    }
	}
      else
	{
	  int nh = prob.ilp_model.im_flow_neighbor_nh+1;
	  int cnh = 1;

	  while( nh--)
	    {
	      //printf("nh %d\n",nh);
	      std::set<int> newones;
	      FOREACH(it, neigh )
		{
		  int j= it->first;
		  Node &nj = m_netWithRelays->getNode(j);
		  if( nj.t == BASE )
		    continue;
		  //printf("checkj %s\n", nj.id.c_str());
		  //fflush(stdout);
	  
		  for( int k=0; k< nnodes; k++)
		    {
		      if( k==j )
			continue;
		      Node &nk = m_netWithRelays->getNode(k);
		      if( !m_netWithRelays->inRange(nj, nk) )
			continue;
		      if( neigh.find(k) == neigh.end() )
			{
			  //printf("new neighbor %d %s\n",
			  //     k, nk.id.c_str());
			  //fflush(stdout);
			  newones.insert(k);
			  //neigh[k] = cnh+1;
			}
		    }
		}
	      FOREACH(it, newones)
		{
		  neigh[*it] = cnh+1;
		}
	      cnh++;
	    }
	}
      printf("neighborhod size %d\n", neigh.size());
      fflush(stdout);
      /// then, add the oflows
      FOREACH(it, neigh)
	{
	  int j= it->first;
	  Node &nj = m_netWithRelays->getNode(j);
	  for( int k=0; k< nnodes; k++)
	    {
	      if( k==j )
		continue;
	      Node &nk = m_netWithRelays->getNode(k);
	      if( !m_netWithRelays->inRange(nj, nk) )
		continue;
		
	      rmatind[nzcnt] = getXVarIndex(nj.id, nk.id);
	      rmatval[nzcnt++] = 1.0;
	    }
	}
      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = -1.0;
      ccnt = 0;
      rcnt = 1;
      rhs = 0;
      sense = 'E';
      colname = NULL;
      sprintf(rowname[0],"lflow[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
    }


#if 0
      for(int j=0; j< nnodes; j++)
	{
	  Node &nj = m_netWithRelays->getNode(j);

	  ///  BUG ALERT!: why considering only the static for flow neighbor constraints ? 
	  //      if(nj.t != RELAY && ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	  if(ni.id != nj.id && m_netWithRelays->inRange(ni,nj))
	    {
	      rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	      rmatval[nzcnt++] = 1.0;
		      
	      for(int k=0; k<nnodes; k++)
		{
		  Node &nk = m_netWithRelays->getNode(k);
		  ///  BUG ALERT!: why considering only the static for flow neighbor constraints ?
		  //	  if( nk.t != RELAY && nj.id != nk.id && m_netWithRelays->inRange(nj,nk))
		  if( nj.id != nk.id && m_netWithRelays->inRange(nj,nk))	  
		    {

		      rmatind[nzcnt] = getXVarIndex(nj.id, nk.id);
		      rmatval[nzcnt++] = 1.0;

		      //          WE FOUND A BUG HERE! 12.02.2015
		      //	    rmatind[nzcnt] = getFXVarIndex(ni.id);
		      //	    rmatval[nzcnt++] = -prob.ilp_model.max_flow;
		    }
		}
	    }
	}

      rmatind[nzcnt] = getLFVarIndex(ni.id);
      rmatval[nzcnt++] = -1.0;
      ccnt = 0;
      rcnt = 1;
      rhs = 0;
      sense = 'E';
      colname = NULL;
      sprintf(rowname[0],"lflow[%s]",ni.id.c_str());
      printf("Adding row %s\n", rowname[0]);
      retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			  &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
      if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);

    }
#endif

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
  
}

void MyCplexSolver::addRelayBoundCons()
{
  /*
   * s.t. bound: sum{k in R} y[k] <= K;
   * s.t. bound_low: sum{k in R} y[k] >= minK;
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));

  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_netWithRelays->numRelays()*sizeof(int));
  double *rmatval = (double *)malloc(m_netWithRelays->numRelays()*sizeof(double));

  int row_index = CPXgetnumrows(m_cpxenv, m_cpxlp);

  //printf("Adding relay bounds constraints\n");
  rmatbeg[0]=0;
  nzcnt=0;
  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    Node &ni = m_netWithRelays->getRelay(i);
    rmatind[nzcnt] = getYVarIndex(ni.id);
    rmatval[nzcnt++] = 1.0;

  }

  ccnt = 0;
  rcnt = 1;
  rhs = prob.K;
  sense = 'L';
  colname = NULL;
  sprintf(rowname[0],"bound");
  //printf("Adding row %s\n", rowname[0]);
  retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
		      &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
  if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
  m_relayUpperBndRow = row_index++;


  rmatbeg[0]=0;
  nzcnt=0;
  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    Node &ni = m_netWithRelays->getRelay(i);
    rmatind[nzcnt] = getYVarIndex(ni.id);
    rmatval[nzcnt++] = 1.0;

  }

  ccnt = 0;
  rcnt = 1;
  rhs = prob.minK;
  sense = 'G';
  colname = NULL;
  sprintf(rowname[0],"bound_low");
  //printf("Adding row %s\n", rowname[0]);
  retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
		      &rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
  if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
  m_relayLowerBndRow = row_index;

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}


void 
MyCplexSolver::addBaseCons()
{
  /*
   * s.t. base: sum {k in B} sum{(i,k) in edges} x[i,k] = sum{i in S} d[i] - D;
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  /// plus one, maybe D var
  int allxs = m_varX.size() + m_varXD.size()+ m_varXF.size()+1;
  int *rmatind = (int *)malloc(allxs*sizeof(int));
  double *rmatval = (double *)malloc(allxs*sizeof(double));

  double allD = 0.0;
  for(int i=0; i< m_netWithRelays->numStatic();i++)
  {
    allD += m_netWithRelays->getDemand(i);
  }

  //printf("Adding base constraint - sum d[i] = %f\n", allD);
  rmatbeg[0]=0;
  nzcnt=0;

  for(int i=0; i< m_netWithRelays->numBases(); i++)
  {
    Node &ni = m_netWithRelays->getBase(i);
    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(nj,ni) && ni.id != nj.id)
      {
	rmatind[nzcnt] = getXVarIndex(nj.id, ni.id);
	rmatval[nzcnt++] = 1.0;
      }
    }
  }

  if( prob.ilp_model.allow_incomplete_delivery )
  {
    /// deficit of delivered demand
    for(int i=0; i< m_netWithRelays->size();i++)
    {	
      Node &ni = m_netWithRelays->getNode(i);
      if( ni.t == BASE || ni.t == RELAY)
	continue;
      rmatind[nzcnt] = getXVarIndex(ni.id, "DUMP");
      rmatval[nzcnt++] = 1.0;
    }
  }

  if( prob.ilp_model.use_flow_neighbor_sink )
    {
      FOREACH(it, m_varXF )
	{
	  rmatind[nzcnt] = it->second;
	  rmatval[nzcnt++] = 1.0;
	}
    }

  if( nzcnt )
  {
    printf("BaseCons nzcnt %d\n", nzcnt);
    ccnt = 0;
    rcnt = 1;
    rhs = allD;
    sense = 'E';
    colname = NULL;
    sprintf(rowname[0],"base");
    printf("Adding row %s\n", rowname[0]);
    fflush(stdout);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
  }


  printf("going to free\n");
  fflush(stdout);
  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
  printf("freed\n");
  fflush(stdout);
}


void MyCplexSolver::addDemandCons()
{
  /*
   * s.t. innodes{k in S}: sum{(k,j) in edges} x[k,j] - sum{(i,k) in edges} x[i,k] = d[k];
   */

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc((m_varX.size()+m_varXD.size())*sizeof(int));
  double *rmatval = (double *)malloc( (m_varX.size()+m_varXD.size())*sizeof(double));

  //printf("Adding %d innode constraints\n", 
  //       m_netWithRelays->numStatic());

  for(int i=0; i< m_netWithRelays->numStatic(); i++)
  {
    Node &ni = m_netWithRelays->getStatic(i);
    rmatbeg[0]=0;
    nzcnt=0;
    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(ni,nj) && ni.id != nj.id)
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	rmatval[nzcnt++] = 1.0;
	/*
	  rmatind[nzcnt] = getXVarIndex(nj.id, ni.id);
	  rmatval[nzcnt++] = -1.0;
	  */
      }
    }

    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(nj,ni) && ni.id != nj.id)
      {
	/*
	  rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	  rmatval[nzcnt++] = 1.0;
	  */

	rmatind[nzcnt] = getXVarIndex(nj.id, ni.id);
	rmatval[nzcnt++] = -1.0;

      }
    }
    if( prob.ilp_model.allow_incomplete_delivery )
    {
      rmatind[nzcnt] = getXVarIndex(ni.id, "DUMP");
      rmatval[nzcnt++] = 1.0;
    }

    if( prob.ilp_model.use_flow_neighbor_sink )
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, "FDUMP");
	rmatval[nzcnt++] = 1.0;
      }

    ccnt = 0;
    rcnt = 1;
    rhs = m_netWithRelays->getDemand(i);
    printf("Demand %d = %.2f\n",i,rhs);
    sense = 'E';
    colname = NULL;
    sprintf(rowname[0],"innode[%s]",ni.id.c_str());
    //printf("Adding row %s\n", rowname[0]);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);


  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}

void MyCplexSolver::addRelaySelCons()
{
  /*  y[k] => some flow passed through k
   * s.t. relayselected{k in R} : FLOW_RELAY_CONSTRAINTS * ((sum{(i,k) in edges} x[i,k]) - y[k]) >= 0;
   */

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));


  //printf("Adding %d relaysel constraints\n", 
  //      m_netWithRelays->numRelays());

  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    Node &ni = m_netWithRelays->getRelay(i);
    rmatbeg[0]=0;
    nzcnt=0;
    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(ni,nj) && ni.id != nj.id)
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	rmatval[nzcnt++] = 1.0;
      }
    }
    if( prob.ilp_model.use_flow_neighbor_sink )
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, "FDUMP");
	rmatval[nzcnt++] = 1.0;
      }
    rmatind[nzcnt] = getYVarIndex(ni.id);
    rmatval[nzcnt++] = -1.0;

    ccnt = 0;
    rcnt = 1;
    rhs = 0.0;
    sense = 'G';
    colname = NULL;
    sprintf(rowname[0],"relaysel[%s]",ni.id.c_str());
    printf("Adding row %s\n", rowname[0]);
    fflush(stdout);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);


  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}

void 
MyCplexSolver::addFlowSelCons()
{
  /*
   * some flow passed through k => y[k]
   * s.t. flowselected{k in R}: FLOW_RELAY_CONSTRAINTS * (sum{(i,k) in edges} x[i,k]) <= sum{i in S} d[i] * y[k];
   */

  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;
  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));

  double allD = 0.0;
  for(int i=0; i< m_netWithRelays->numStatic();i++)
  {
    allD += m_netWithRelays->getDemand(i);
  }

  //printf("Adding %d flowsel constraints\n", 
  //       m_netWithRelays->numRelays());

  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    Node &ni = m_netWithRelays->getRelay(i);
    rmatbeg[0]=0;
    nzcnt=0;
    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(ni,nj) && ni.id != nj.id)
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, nj.id);
	rmatval[nzcnt++] = 1.0;


      }
    }
    if( prob.ilp_model.use_flow_neighbor_sink )
      {
	rmatind[nzcnt] = getXVarIndex(ni.id, "FDUMP");
	rmatval[nzcnt++] = 1.0;
      }
    rmatind[nzcnt] = getYVarIndex(ni.id);
    rmatval[nzcnt++] = -allD;

    ccnt = 0;
    rcnt = 1;
    rhs = 0.0;
    sense = 'L';
    colname = NULL;
    sprintf(rowname[0],"flowsel[%s]",ni.id.c_str());
    printf("Adding row %s\n", rowname[0]);
    fflush(stdout);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);


  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}

void 
MyCplexSolver::setObjFunc()
{
  /*
   * int CPXchgobj(CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const
   * double * values)
   */

  /*
   * (sum{(i,j) in edges} c[i,j] * x[i,j])
   * + sum{i in S union R} fx[i] * FLOW_NEIGHBOUR_PENALTY + sum{k in R} y[k] * PENALTY_RELAY;
   */

  int retCPX;
  CPXchgobjsen(m_cpxenv, m_cpxlp, CPX_MIN);

  printf("setting objecive\n");
  fflush(stdout);

  int allxs = m_varX.size() + m_varFX.size()+ m_varY.size() + m_varXD.size() + m_varXF.size();
  int *indeces = (int *)malloc(allxs*sizeof(int));
  double *values = (double *)malloc(allxs*sizeof(double));



  double allD = 0.0;
  for(int i=0; i< m_netWithRelays->numStatic();i++)
  {
    allD += m_netWithRelays->getDemand(i);
  }

  int cnt=0;

  if(  prob.ilp_model.allow_incomplete_delivery )
  {
    double incomplete_delivery_penalty =
      m_netWithRelays->maxPathLength()
      * m_metric->max()
      * prob.ilp_model.incomplete_delivery_penalty;
    printf("using incomplete_delivery_penalty %.2f\n",
	   incomplete_delivery_penalty);
    printf("maxPathLength %d metricMax %f penalty_w %.2f\n",
	   m_netWithRelays->maxPathLength(),
	   m_metric->max(),
	   prob.ilp_model.incomplete_delivery_penalty);

    for(int i=0; i< m_netWithRelays->size();i++)
    {	
      Node &ni = m_netWithRelays->getNode(i);
      if( ni.t == BASE || ni.t == RELAY)
	continue;
      indeces[cnt] = getXVarIndex(ni.id, "DUMP");
      values[cnt++] = incomplete_delivery_penalty;
    }
  }

  if( prob.ilp_model.use_flow_neighbor_sink )
    {
      FOREACH(it, m_varXF)
	{
	  indeces[cnt] = it->second; 
	  values[cnt++] = allD*m_metric->max()*prob.ilp_model.lost_flow_penalty;
	}
    }

  /// minimize cost of flow
  FOREACH(it, m_varX)
  {
    const pair<string, string> &p = it->first;
    int ix = it->second;

    double link_cost;
    link_cost = (*m_metric)(p.first,p.second);
    //		printf("link_cost %s %s %.2f\n",p.first.c_str(),p.second.c_str(),link_cost);

    indeces[cnt] = ix;
    values[cnt++] = link_cost;
  }

  /// + flow neighbor penalties
  double flow_neighbor_penalty = (m_netWithRelays->expectedHopCount() * m_metric->max() * allD);
  printf("flow_neighbor_penalty %f w %f\n",flow_neighbor_penalty,
	 prob.ilp_model.flow_neighbor_penalty_w);
  flow_neighbor_penalty *= prob.ilp_model.flow_neighbor_penalty_w;

  //flow_neighbor_penalty /= m_netWithRelays->numStatic();
  FOREACH(it, m_varFX)
  {
    int ix = it->second;

    indeces[cnt] = ix;
    values[cnt++] = flow_neighbor_penalty;
  }

  FOREACH(it, m_varY)
  {
    int ix = it->second;
    indeces[cnt] = ix;
    values[cnt++] = prob.ilp_model.penalty_relay;
  }
  retCPX = CPXchgobj(m_cpxenv, m_cpxlp, cnt,indeces,values);

  if (retCPX) my_free_cplex("CPX chgobj", retCPX, m_cpxenv, m_cpxlp);
  printf("objecive SET\n");
  fflush(stdout);

}

  void 
MyCplexSolver::addBWCapacityCons()
{
  /*
   * s.t. bandwidth_capacity{i in nodes}: (sum{(i,j) in edges} x[i,j] + sum{(k,i) in edges} x[k,i]) <= MAX_CAPACITY;
   */
  int ccnt;
  int rcnt;
  int nzcnt;
  int retCPX;
  int nnodes = m_netWithRelays->size();
  double rhs;
  char sense;
  char **colname;
  char **rowname;

  /** CPLEX: get numvars **/
  int var_index;
  //int row_index = CPXgetnumrows(cpx_env, cpx_lp);

  rowname = (char **)malloc(1*sizeof(char *));
  rowname[0] = (char *)malloc(25*sizeof(char));


  int *rmatbeg = (int *)malloc(1*sizeof(int));
  int *rmatind = (int *)malloc(m_varX.size()*sizeof(int));
  double *rmatval = (double *)malloc(m_varX.size()*sizeof(double));

  printf("Adding %d bandwidth capacity constraints\n", nnodes);
  fflush(stdout);

  for(int i=0; i< nnodes; i++)
  {
    Node &ni = m_netWithRelays->getNode(i);
    rmatbeg[0]=0;
    nzcnt=0;
    for(int j=0; j< nnodes; j++)
    {
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(ni,nj) && i!=j)
      {
	var_index = getXVarIndex(ni.id, nj.id);

	rmatind[nzcnt] = var_index; 				
	rmatval[nzcnt++] = 1.0;
      }

      if( m_netWithRelays->inRange(nj,ni) && i!=j)
      {
	//m_coefByVar[var_index].insert(row_index);

	var_index = getXVarIndex(nj.id, ni.id);
	rmatind[nzcnt] = var_index; 				
	rmatval[nzcnt++] = 1.0;
	//m_coefByVar[var_index].insert(row_index);
      }
    }
    ccnt = 0;
    rcnt = 1;
    rhs = prob.ilp_model.max_capacity;
    sense = 'L';
    colname = NULL;
    sprintf(rowname[0],"bwcap[%s]",ni.id.c_str());
    printf("Adding row %s\n", rowname[0]);
    fflush(stdout);
    retCPX = CPXaddrows(m_cpxenv, m_cpxlp,ccnt,rcnt,nzcnt,
			&rhs,&sense, rmatbeg,rmatind,rmatval,colname,rowname);
    if (retCPX) my_free_cplex("CPX add rows", retCPX, m_cpxenv, m_cpxlp);
  }

  free(rowname[0]);
  free(rowname);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
}

void 
MyCplexSolver::setRelayBnds(int mink, int maxk)
{
  /*
   * int CPXchgrhs(CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const
   * double * values)
   */
  int retCPX;
  int indices[2];

  double values[2];

  indices[0] = m_relayUpperBndRow;
  indices[1] = m_relayLowerBndRow;

  values[0] = maxk;
  values[1] = mink;

  retCPX = CPXchgrhs(m_cpxenv, m_cpxlp, 2, indices, values);
  if( retCPX)
    my_free_cplex("CPX chgrhs", retCPX, m_cpxenv, m_cpxlp);

}

void 
MyCplexSolver::setFlowRelayCons()
{
  addRelaySelCons();
}

void MyCplexSolver::initialize()
{
  CpuTime cpu_time;
  cpu_time.start();
  for(int i=0; i< m_netWithRelays->size();i++)
  {	
    Node &ni = m_netWithRelays->getNode(i);
    for(int j=0; j<m_netWithRelays->size();j++)
    {
      if( i== j)
	continue;
      Node &nj = m_netWithRelays->getNode(j);

      if( m_netWithRelays->inRange(ni,nj))
      {
	addEdgeVar(ni.id,nj.id);
      }

    }
    if(ni.t == RELAY)
    {
      addRelayVar(ni.id);
    }
    else
    {
      if( prob.ilp_model.use_flow_neighbor_constraints )
	addFXVar(ni.id);
    }
    if( prob.ilp_model.use_flow_neighbor_constraints &&
	prob.ilp_model.use_im_flow_neighbor_constraints)
      {
	printf("ERROR: can not use old and improved flow neighbor constraints at the same time\n");
	exit(-1);
      }
    ///! local flow variables
    ///! 
    if( prob.ilp_model.use_im_flow_neighbor_constraints )
      {
	if( ni.t != BASE )
	  {
	  addFXVar(ni.id);
	  addLFVar(ni.id);
	  }
      }

  }
  addOtherVars();

  printf("added %d X Variables\n", m_varX.size());
  fflush(stdout);
  
  addFlowCons();
  addFlowEdgeCons();
  addBWCapacityCons();
  //    if( !m_isDisabled[FLOW_NEIGHBOR_CONST])
  if( prob.ilp_model.use_flow_neighbor_constraints)
    addFlowNeighborCons();
  //    if( !m_isDisabled[NODE_DEGREE_CONST])

  if( prob.ilp_model.use_im_flow_neighbor_constraints )
    {
      addLocalFlowCons();
      addImFlowNeighborCons();
    }
  if( prob.ilp_model.strict_flow_neighbor_constraints )
    {
      addStrictFlowNeighborCons();
    }
  if( prob.ilp_model.use_flow_neighbor_sink )
    {
      addSinkFlowCons();
    }
  if( prob.ilp_model.use_node_degree_constraints)
    addNodeDegreeCons();
  addDemandCons();

  addFlowSelCons();
  if( prob.ilp_model.flow_relay_constraints)
  {

    addRelaySelCons();
  }
  printf("added Relay Sel Constraints\n");
  fflush(stdout);
  ///------------------------
  addBaseCons();
  printf("added Base Constraints\n");
  fflush(stdout);
  ///------------------------
  addRelayBoundCons();
  printf("added Base Constraints\n");
  fflush(stdout);
  ///------------------------
  
  printf("Using metric: %s\n", prob.metric.c_str());
  fflush(stdout);
  if( prob.metric == "PRR")
  {
    m_metric = new EstPRR_Metric(m_netWithRelays);
  }
  else if(prob.metric == "DISTANCE_SQ")
  {
    m_metric = new DistSQ_Metric(m_netWithRelays);
    //cout << *m_metric << endl;
  }
  else if(prob.metric == "LINK_WEIGHT")
  {
    m_metric = new LinkWeight_Metric(m_netWithRelays);
  }
  else if(prob.metric == "CONSTANT")
  {
    m_metric = new Constant_Metric(m_netWithRelays, 1.0);
  }
  else
  {
    printf("Invalid metric %s\n", prob.metric.c_str());
    exit(-1);
  }



  setObjFunc();
  //enableAllRelays();

  m_cpuTimeModelGeneration = cpu_time.cpu_time_elapsed();
  m_keepVarsInfo = false;
  m_varsinfo = NULL;

}

  void
MyCplexSolver::keepVarsInfo(bool keep)
{
  m_keepVarsInfo = keep;
}
  void 
MyCplexSolver::enableRelay(int rid)
{
  Node &nr = m_netWithRelays->getRelay(rid);
  enableRelay(nr.id);
}
void MyCplexSolver::enableRelay(string rid)
{
  /*
   * CPXtightenbds(CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const char * lu, const double * bd)
   */
  int retCPX;
  int cnt;
  int var_ix = getYVarIndex(rid);
  char bnd;
  double bnd_val;

  bnd = 'L';
  bnd_val = 0.0;
  cnt = 1;
  retCPX = CPXtightenbds(m_cpxenv,m_cpxlp,cnt,&var_ix, &bnd,&bnd_val);
  if( retCPX)
    my_free_cplex("CPX tightenbds", retCPX, m_cpxenv, m_cpxlp);


  bnd = 'U';
  bnd_val = 1.0;
  retCPX = CPXtightenbds(m_cpxenv,m_cpxlp,cnt,&var_ix,&bnd,&bnd_val);
  if( retCPX)
    my_free_cplex("CPX tightenbds", retCPX, m_cpxenv, m_cpxlp);
}

void MyCplexSolver::fixRelay(string rid)
{
  /*
   * CPXtightenbds(CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const char * lu, const double * bd)
   */
  int retCPX;
  int cnt;
  int var_ix = getYVarIndex(rid);
  char bnd;
  double bnd_val;
  cnt = 1;
  bnd = 'B';
  bnd_val = 1.0;
  retCPX = CPXtightenbds(m_cpxenv,m_cpxlp,cnt,&var_ix,&bnd,&bnd_val);
  if( retCPX)
    my_free_cplex("CPX tightenbds", retCPX, m_cpxenv, m_cpxlp);
}

void MyCplexSolver::fixRelay(int rid)
{

  Node &nr= m_netWithRelays->getRelay(rid);
  fixRelay(nr.id);
}

void MyCplexSolver::disableRelay(int rid)
{
  /*
   * CPXtightenbds(CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const char * lu, const double * bd)
   */
  int retCPX;
  int cnt;
  int var_ix = getYVarIndex(m_netWithRelays->getRelay(rid).id);
  char bnd;
  double bnd_val;
  cnt = 1;
  bnd = 'B';
  bnd_val = 0.0;
  retCPX = CPXtightenbds(m_cpxenv,m_cpxlp,cnt,&var_ix,&bnd,&bnd_val);
  if( retCPX)
    my_free_cplex("CPX tightenbds", retCPX, m_cpxenv, m_cpxlp);


}
void MyCplexSolver::enableAllRelays()
{
  //printf("Enabling all relays\n");
  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    enableRelay(i);
  }
}

void MyCplexSolver::disableAllRelays()
{
  for(int i=0; i< m_netWithRelays->numRelays(); i++)
  {
    disableRelay(i);
  }
}

MyCplexSolver::MyCplexSolver(Network *static_net, Network *net_with_relays):  
  m_staticNet(static_net),m_netWithRelays(net_with_relays)
{
  m_metric = NULL;
  m_simpleCallback = NULL;
  m_heuristicCallback = NULL;
  m_incumbentCallback = NULL;
  m_useInitial = false;

  int retCPX;
  m_cpxlp = NULL;

  m_cpxenv = CPXopenCPLEX(&retCPX);
  if (retCPX || m_cpxenv == NULL)
  {
    fprintf(stderr, "Error: could not create cpx_environment\n");
    my_free_cplex("CPX create environment", retCPX, m_cpxenv, m_cpxlp);
  }


  /// CPLEX: Create LP model 
  m_cpxlp = CPXcreateprob(m_cpxenv, &retCPX, "cpx_lp");
  if (retCPX) my_free_cplex("CPX create problem", retCPX, m_cpxenv, m_cpxlp);

  initialize();
}

double MyCplexSolver::getObjVal()
{
  assert(m_sol.get() != NULL);
  return m_sol->getObjVal();
}

RNPSolutionPtr MyCplexSolver::solution()
{
  return m_sol;
}

  void
MyCplexSolver::useInitial(RNPSolutionPtr initsol)
{
  m_useInitial = true;
  m_initSol = initsol;
}

bool MyCplexSolver::solve()
{
  /// CPLEX: data structures 
  CPXENVptr cpx_env = m_cpxenv;
  /// CPLEX: variables
  int      numvars, j, retCPX, CPX_IsMIP, CPX_noMIP, CPX_probType;
  double   db_tmp, bound,obj_val;
  char nameCPX[80];
  char type;
  double cpx_feasibilityTol = 10E-6;
  double cpx_intFeasTol = 10E-6;
  double cpx_mipGap = prob.mycplex_params.mipGap;

  /// Memory limit
  double cpx_trelim = prob.mycplex_params.trelim; 
  double cpx_workmem = prob.mycplex_params.workmem; 
  int cpx_nodefileind = prob.mycplex_params.nodefileind; 

  ///
  int cpx_mipemphasis = prob.mycplex_params.mipemphasis;
  int cpx_threads = prob.mycplex_params.threads; 
  int cpx_preind = prob.mycplex_params.preind;
  int cpx_lpmethod = prob.mycplex_params.lpmethod;
  int cpx_parallelmode = prob.mycplex_params.parallelmode;
  int cpx_mipdisplay = prob.mycplex_params.mipdisplay;
  int cpx_probe = prob.mycplex_params.probe;
  int cpx_writemps = prob.lp_params.save_mps;
  int cpx_writeparam = prob.mycplex_params.save_cpx_params;

  /// execution flags 
  char CPX_out, GLPK_out, verbose;
  CPX_out = GLPK_out = CPX_noMIP = verbose = 0;
  double cpx_strongtimelim = prob.lp_params.strong_timelim; // Strong Time limit
  double cpx_weaktimelim = prob.lp_params.weak_timelim; // weak time limit
  CpuTime cpu_time;

  int usetimelimcallback = 
    ((m_simpleCallback != NULL) 
     || prob.lp_params.log_progress 
     || cpx_weaktimelim > 0 
     || prob.lp_params.get_suboptimal);
  if( prob.lp_params.log_progress)
  {
    string progress_file=prob.output_path  + prob.instance_id + ".cpxlog";
    of_progress.open(progress_file.c_str(), ios_base::out | ios_base::trunc );
  }


  int useheuristiccallback = (m_heuristicCallback != NULL);
  int useincumbentcallback = (m_incumbentCallback != NULL);
  CPX_TIMELIMINFO mytimeliminfo;


  //	CPX_out = param.verbose;
  //	GLPK_out = param.verbose;

  CPX_noMIP = prob.lp_params.noMIP;
  verbose = prob.lp_params.verbose;
  CPX_out = verbose;

  string cpxout_file;
  CPXFILEptr fp_cpxout;
  if( prob.mycplex_params.logToFile)
  {
    cpxout_file = prob.output_path + 
      prob.instance_id + ".cpxout";
    fp_cpxout = CPXfopen(cpxout_file.c_str(), "w");
    retCPX =  CPXsetlogfile(m_cpxenv, fp_cpxout);
    if(retCPX )
      fprintf(stderr, "Cannot open cpxout file\n");
    else
      fprintf(stdout, "Logging CPLEX in %s\n", cpxout_file.c_str());
  }


  /// Time recording
  cpu_time.start();

  CPXLPptr cpx_lp = CPXcloneprob(m_cpxenv,m_cpxlp,&retCPX);
  //CPXLPptr cpx_lp= m_cpxlp;
  if (retCPX) my_free_cplex("CPX clone prob", retCPX, cpx_env, cpx_lp);

  string nameVar;

  LpSolutionPtr ret_sol(new LpSolution());


  //uncomment
  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_SCRIND, CPX_out?CPX_ON:CPX_OFF);
  if (retCPX) my_free_cplex("CPX set output", retCPX, cpx_env, cpx_lp);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_PROBE, cpx_probe);
  if (retCPX) my_free_cplex("CPX set probe", retCPX, cpx_env, cpx_lp);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_MIPDISPLAY, cpx_mipdisplay);
  if (retCPX) my_free_cplex("CPX set mip display", retCPX, cpx_env, cpx_lp);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
  if (retCPX) my_free_cplex("CPX set mipcbredlp", retCPX, cpx_env, cpx_lp);

  if(prob.lp_params.force_branch_and_cut)
  {
    printf("CPX forcing branch and cut\n");
    retCPX = CPXsetintparam(cpx_env, CPX_PARAM_MIPSEARCH, 
			    CPX_MIPSEARCH_TRADITIONAL);
    if (retCPX) my_free_cplex("CPX set search strategy", retCPX, cpx_env, cpx_lp);
  }

  //uncomment
  if(!prob.lp_params.allow_parallel)
  {
    retCPX = CPXsetintparam(cpx_env, CPX_PARAM_THREADS, 1);
    if (retCPX) my_free_cplex("CPX set threads", retCPX, cpx_env, cpx_lp);
  }
  else
  {
    if( cpx_threads > 0 )
    {
      retCPX = CPXsetintparam(cpx_env, CPX_PARAM_THREADS, cpx_threads);
      if (retCPX) my_free_cplex("CPX set threads", retCPX, cpx_env, cpx_lp);
    }
  }

  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_TRELIM , cpx_trelim);
  if (retCPX) my_free_cplex("CPX set trelim", retCPX, cpx_env, cpx_lp);
  printf("CPX MEMORY: set trelim %f\n", cpx_trelim);

  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_WORKMEM , cpx_workmem);
  if (retCPX) my_free_cplex("CPX set workmem", retCPX, cpx_env, cpx_lp);
  printf("CPX MEMORY: set workmem %f\n", cpx_workmem);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_NODEFILEIND , cpx_nodefileind);
  if (retCPX) my_free_cplex("CPX set nodefileind", retCPX, cpx_env, cpx_lp);
  printf("CPX MEMORY: set nodefileind %d\n", cpx_nodefileind);


  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_MIPEMPHASIS , cpx_mipemphasis);
  if (retCPX) my_free_cplex("CPX set mipemphasis", retCPX, cpx_env, cpx_lp);

  
  /// CPLEX: Ask for more precision 
   //uncomment
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPRHS , cpx_feasibilityTol);
  if (retCPX) my_free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPINT, cpx_intFeasTol);
  if (retCPX) my_free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
  retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_EPGAP, cpx_mipGap);
  if (retCPX) my_free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);

  if(cpx_strongtimelim > 0)
  {
    printf("Using CPX_PARAM_TILIM %f\n", cpx_strongtimelim);
    retCPX = CPXsetdblparam(cpx_env, CPX_PARAM_TILIM, cpx_strongtimelim);
    if (retCPX) my_free_cplex("CPX set dbl param", retCPX, cpx_env, cpx_lp);
  }

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_LPMETHOD, cpx_lpmethod);
  printf("CPLEX: Seting lpmethod to %d\n", cpx_lpmethod);
  if (retCPX) my_free_cplex("CPX set lpmethod", retCPX, cpx_env, cpx_lp);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, cpx_preind);
  if (retCPX) my_free_cplex("CPX set preind", retCPX, cpx_env, cpx_lp);

  retCPX = CPXsetintparam(cpx_env, CPX_PARAM_PARALLELMODE, cpx_parallelmode);
  if (retCPX) my_free_cplex("CPX set parallelmode", retCPX, cpx_env, cpx_lp);

  /// CPLEX: get numvars
  numvars = CPXgetnumcols(cpx_env, cpx_lp);


  double timestart;
  retCPX = CPXgettime (cpx_env, &timestart);
  if ( retCPX )
  {
    my_free_cplex("CPX to initialize timer\n", retCPX, cpx_env, cpx_lp);
  }

  /// CPLEX: Callback setting

  if( usetimelimcallback || useheuristiccallback || useincumbentcallback)
  {

    mytimeliminfo.acceptablegap = 100.0;
    mytimeliminfo.aborted       = 0;
    mytimeliminfo.timestart     = timestart;
    mytimeliminfo.timelim       = cpx_weaktimelim;
    mytimeliminfo.numvars       = numvars;
    if( m_keepVarsInfo )
    {
      if( m_varsinfo == NULL )
	m_varsinfo = MyCPXgetVariablesInfo(cpx_env, cpx_lp);
      mytimeliminfo.varsinfo = m_varsinfo;
    }
    else
    {
      mytimeliminfo.varsinfo = MyCPXgetVariablesInfo(cpx_env, cpx_lp);
    }
    mytimeliminfo.infocallbackfunc = m_simpleCallback;
    mytimeliminfo.heuristiccallbackfunc = m_heuristicCallback;
    mytimeliminfo.incumbentcallbackfunc = m_incumbentCallback;
    mytimeliminfo.output_progress = prob.lp_params.log_progress;
    mytimeliminfo.last_progress_report = 0; /// in seconds
    mytimeliminfo.progress_report_interval = 2; /// in seconds
    mytimeliminfo.has_incumbent = false;
    mytimeliminfo.bestobj = 100000000;
  }
  if ( usetimelimcallback ) 
  {
    if (verbose) printf ("setting-up info callback\n");
    retCPX = CPXsetinfocallbackfunc (cpx_env, timelimcallback, &mytimeliminfo);
    if ( retCPX ) {
      fprintf (stderr, "Failed to set time limit callback function.\n");
      my_free_cplex("CPX Failed callback set", retCPX, cpx_env, cpx_lp);

    }
  }

  if ( useheuristiccallback )
  {
    if (verbose) printf ("setting-up heuristic callback\n");
    retCPX = CPXsetheuristiccallbackfunc(cpx_env, heuristiccallback, 
					 &mytimeliminfo);
    if ( retCPX ) {
      fprintf (stderr, "Failed to set heuristic callback function.\n");
      my_free_cplex("CPX Failed callback set", retCPX, cpx_env, cpx_lp);

    }
  }

  if ( useincumbentcallback )
  {

    if (verbose) printf ("setting-up incumbent callback\n");
    retCPX = CPXsetincumbentcallbackfunc(cpx_env, incumbentcallback, 
					 &mytimeliminfo);
    if ( retCPX ) {
      fprintf (stderr, "Failed to set incumbent callback function.\n");
      my_free_cplex("CPX Failed callback set", retCPX, cpx_env, cpx_lp);

    }
  }

  /** CPLEX: get model type **/
  CPX_probType = CPXgetprobtype(cpx_env, cpx_lp);
  if ((CPX_probType != CPXPROB_MILP) && (CPX_probType != CPXPROB_LP))
  {
    fprintf(stderr, "Error: unable to solve problem type :%d\n", CPX_probType);
    my_free_cplex("CPX problem type", 1, cpx_env, cpx_lp);
  }
  if (verbose) printf ("MyCPLEX: Problem type: %d\n", CPX_probType);
  CPX_IsMIP = (CPXPROB_MILP == CPX_probType);


  /** Assume it is going to be solved */
  bool got_solved = true;

  /** Assume we get a (sub-) optimal solution */
  bool has_solution = true;

  /** CPLEX: Solve the linear relaxation **/
  if ((!CPX_IsMIP))
  {
    //TODO allow for other simplex optimizers
    retCPX = CPXdualopt(cpx_env, cpx_lp);
    if (retCPX) my_free_cplex("CPX LP Optimize", retCPX, cpx_env, cpx_lp);

    /** CPLEX: Retreive the optimization status **/
    retCPX = CPXgetstat(cpx_env, cpx_lp);
    if (retCPX != CPX_STAT_OPTIMAL)
    {
      /** This is not an error **/
      //my_free_cplex("Error CPX Lp optimization failed", retCPX, cpx_env, cpx_lp);
      got_solved = false;
      has_solution = false;
      //return ret_sol;
    }
  }
  else
  {
    if(CPX_noMIP)
    {
      if (verbose) printf ("changing MIP to LP\n" );
      /// It is a MIP Problem, but we change it to LP
      retCPX = CPXchgprobtype(cpx_env, cpx_lp,CPXPROB_LP);
      if (retCPX) my_free_cplex("CPX change problem type", retCPX, cpx_env, cpx_lp);

      retCPX = CPXlpopt(cpx_env, cpx_lp);
      if (retCPX) my_free_cplex("CPX LP Optimize", retCPX, cpx_env, cpx_lp);

      /// CPLEX: Retreive the optimization status
      retCPX = CPXgetstat(cpx_env, cpx_lp);
      if ((retCPX == CPX_STAT_UNBOUNDED) || (retCPX == CPX_STAT_INFEASIBLE)
	  || (retCPX == CPX_STAT_INForUNBD)){
	printf("Model is infeasible or unbounded\n");

	/// This is not an error 
	//my_free_cplex("Error CPX LP optimization failed", retCPX, cpx_env, cpx_lp);
	got_solved = false;
	has_solution = false;
	//return ret_sol;
      }
    }
    else
    {
#define DEBUGMIPSTART 0
      if( m_useInitial )
      {
	LpSolutionPtr initsol_lp = m_initSol->getLpSolution();
	int mcnt = 1;
	int nzcnt = initsol_lp->value.size();
	printf("Adding MIP Start with %d values\n",
	       nzcnt);
	if( nzcnt == 0)
	{
	  my_free_cplex("CPX MIPStart empty", retCPX, cpx_env, cpx_lp);
	}
	int *beg = (int *) malloc( sizeof(int) );
	beg[0] = 0;

	int *effortlevel = (int *) malloc( sizeof(int) );
	//effortlevel[0] = CPX_MIPSTART_SOLVEMIP ;
	effortlevel[0] = CPX_MIPSTART_AUTO ;
#if DEBUGMIPSTART
	int varindices[1];
	double values[1];
	nzcnt = 1;
	printf("DEBUGMIPSTART only one variable\n");
#else
	int *varindices = (int *) malloc( nzcnt * sizeof(int) );
	double *values = (double *) malloc( nzcnt * sizeof(double) );
#endif
	int ix=0;
	FOREACH(it, initsol_lp->value)
	{
	  string vname = it->first;
	  int vix;
	  retCPX = CPXgetcolindex( cpx_env, cpx_lp, vname.c_str(), &vix );
	  if (retCPX) my_free_cplex("CPX getcolindex", retCPX, cpx_env, cpx_lp);
	  //printf("index of %s is %d\n",
	  //vname.c_str(), vix);
	  varindices[ix] = vix;
	  values[ix] = it->second;
	  ix++;
#if DEBUGMIPSTART
	  printf("DEBUGMIPSTART variable set %d %f\n",
		 varindices[0], values[0]);
	  break;
#endif
	}
	//char ** mipstartname = NULL;
	printf("about to add mip start\n");
	retCPX = CPXaddmipstarts( cpx_env, cpx_lp, mcnt, nzcnt, beg, 
				  varindices, values, effortlevel, NULL );
	if (retCPX) my_free_cplex("CPX addmipstart", retCPX, cpx_env, cpx_lp);

	printf("added mip start\n");
	
#if !DEBUGMIPSTART
	free( varindices );
	free( values );
#endif
	free( beg );
	free( effortlevel );
	
      }

      if( cpx_writemps)
      {
	string mps_file=prob.output_path  + prob.instance_id + ".mps";
	//string mps_file = prob.lp_params.output_mps_file;
	if (verbose) printf ("writing mps to %s\n", mps_file.c_str());
	retCPX = CPXwriteprob (cpx_env, cpx_lp, mps_file.c_str(), NULL);
	if (retCPX) my_free_cplex("CPX writeprob", retCPX, cpx_env, cpx_lp);
      }


      if( cpx_writeparam)
      {
	string cpxpar_file=prob.output_path  + prob.instance_id + ".cpxpar";
	if (verbose) printf ("writing params to %s\n", cpxpar_file.c_str());
	retCPX = CPXwriteparam(cpx_env, cpxpar_file.c_str());
	if (retCPX) my_free_cplex("CPX writeparam", retCPX, cpx_env, cpx_lp);
      }
      if (verbose ) printf("CPLEX: about to optimize\n");
      retCPX = CPXmipopt(cpx_env, cpx_lp);
      if (verbose )printf("CPLEX: mip opt DONE!\n");
      if (retCPX) my_free_cplex("CPX MIP Optimize", retCPX, cpx_env, cpx_lp);

      /// CPLEX: retrieve the optimization status/
      retCPX = CPXgetstat(cpx_env, cpx_lp);
      if ((retCPX != CPXMIP_OPTIMAL) && (retCPX != CPXMIP_OPTIMAL_TOL))
      {
	got_solved = false;
	has_solution = false;
	printf("CPLEX: optimal NOT FOUND checking for SUBOPTIMAL\n");

	if( prob.lp_params.get_suboptimal)
	{
	  if(mytimeliminfo.incumbent != NULL)
	  {
	    // old ret_sol should die (shared_ptr)
	    ret_sol = mytimeliminfo.incumbent;
	    has_solution = true;
	    if(verbose)
	    {
	      printf("CPLEX: obtained suboptimal\n");
	    }

	    //my_free_cplex("Finished", 0, cpx_env, cpx_lp);
	    //return ret_sol;
	  }
	  else
	  {
	    if(verbose)
	    {
	      printf("CPLEX: could not obtain sub-optimal\n");
	    }
	  }
	}


	//my_free_cplex("Error CPX MIP optimization failed", retCPX, cpx_env, cpx_lp);
	//et_sol->solved = false;
	//return ret_sol;
      }
    }
  }
  // Uncomment to disable code below
  //got_solved = false;

  /// If got solve, retrieve INFO and fill LpSol
  if( got_solved)
  {
    /// CPLEX: Get obj function value 
    retCPX = CPXgetobjval(cpx_env, cpx_lp, &db_tmp);
    if (retCPX) my_free_cplex("CPX obj value", retCPX, cpx_env, cpx_lp);
    if (verbose) printf ("Objective %lf\n", db_tmp);

    if (!CPX_noMIP)
    {
      retCPX = CPXgetbestobjval(cpx_env, cpx_lp, &bound);
      if (retCPX) my_free_cplex("CPX accessing bound", retCPX, cpx_env, cpx_lp);
      if (verbose) printf ("Best bound %lf\n", bound);
      if (verbose) printf ("Absolute gap %lf\n", fabs(db_tmp - bound));
    }

    obj_val = db_tmp;

    double gap;
    retCPX = CPXgetmiprelgap(cpx_env, cpx_lp, &gap);
    if (retCPX) my_free_cplex("CPX mip rel gap", retCPX, cpx_env, cpx_lp);

    /// Turn the gap into a percentage
    gap *= 100.0;

    double timeend;
    retCPX = CPXgettime(cpx_env, &timeend);
    if (retCPX) my_free_cplex("CPX get time", retCPX, cpx_env, cpx_lp);

    if( prob.lp_params.log_progress)
    {
      /// log final
      of_progress << (timeend-timestart) << " " << obj_val << " " << gap << endl;
    }


    ret_sol->objval = obj_val;
    ret_sol->is_optimal(true);
    ret_sol->solved = true;
    /// CPLEX: Get variable values
    for (j = 0; j < numvars; ++j)
    {
      retCPX = CPXgetx (cpx_env, cpx_lp, &db_tmp, j, j);
      if (retCPX) my_free_cplex("CPX Get var value", retCPX, cpx_env, cpx_lp);


      int surplus, colnamespace; 
      char ** colname = NULL;
      char * colnamestore = NULL;

      retCPX = CPXgetcolname(cpx_env, cpx_lp, NULL, NULL, 0, &surplus, j, j);
      if (( retCPX != CPXERR_NEGATIVE_SURPLUS ) && ( retCPX != 0 ))  {
	my_free_cplex("CPX Get var names", retCPX, cpx_env, cpx_lp);
      }

      colnamespace = - surplus;
      if ( colnamespace > 0 ) 
      {
	colname = (char **) malloc (sizeof(char *));
	colnamestore = (char *)  malloc(colnamespace);

	if ( colname == NULL || colnamestore == NULL ) 
	{
	  my_free_cplex("CPX Allocating memory for col name" , 1, cpx_env, cpx_lp);
	}

	retCPX = 
	  CPXgetcolname (cpx_env, cpx_lp, colname, colnamestore, colnamespace, &surplus, j, j);
	if ( retCPX ) 
	{
	  my_free_cplex("CPX Get final var names", retCPX, cpx_env, cpx_lp);
	}
      }
      else 
      {
	my_free_cplex("CPX no name associated", 1, cpx_env, cpx_lp);
      }

      //if (verbose) printf ("Processed variable %d name %s\n", j, colname[0]);

      sprintf(nameCPX, "%s", colname[0]);

      free(colnamestore);
      free(colname);

      char type[1];
      if( !CPX_noMIP)
      {
	retCPX = CPXgetctype (cpx_env, cpx_lp, type, j,j);
	if (retCPX) my_free_cplex("CPX Accessing variable type", retCPX, cpx_env, cpx_lp);
      }
      else
      {
	type[0] = 'C'; // continuous variable
      }
      nameVar = nameCPX;	

      double tmpdbl = db_tmp;
      if(tmpdbl > EPSILON)
      {
	ret_sol->value[nameVar] = tmpdbl;
	ret_sol->type[nameVar] = type[0];
	//if(fabs(ret_sol->value[nameVar]) > 0.0){
	//cout << "Writing sol: " << nameVar 
	//	<< " has value " << ret_sol->value[nameVar] << endl;
	//	}
	//	DEBUG(1) cdebug << "Writing sol: " << nameVar 
	//		<< " has value " << ret_sol->value[nameVar] << endl;
      }

    }
  }

  ret_sol->cpuTimeSolving(cpu_time.cpu_time_elapsed());
  ret_sol->cpuTimeTotal(cpu_time.cpu_time_elapsed());

  if( has_solution && m_mathh != NULL )
  {
    /// push last solution
    int error = m_mathh->pushSolution(ret_sol);
    if(error)
    {
      printf("Error pushing last solution\n");
    }
    else
    {
      printf("MyCplexSolver: Pushed last solution with obj %f\n", ret_sol->objval);
    }
  }
  /** Here we clean all **/
  if (cpx_lp) CPXfreeprob(cpx_env, &cpx_lp);

  if( usetimelimcallback || useheuristiccallback || useincumbentcallback)
  {

    if( !m_keepVarsInfo )
      /// Clean timelim struct 
      delete mytimeliminfo.varsinfo;
  }
  if( of_progress.is_open())
  {
    of_progress.close();
  }

  if( prob.mycplex_params.logToFile)
  {
    retCPX =  CPXfclose(fp_cpxout);
    if(retCPX )fprintf(stderr, "CPLEX: cannot close cpxout file\n");
  }

  RNPSolutionPtr sol( new RNPSolution(ret_sol));
  m_sol = sol;
  return m_sol->solved();
}

MyCplexSolver::~MyCplexSolver()
{
  delete m_metric;
  my_free_cplex("Finished", 0, m_cpxenv, m_cpxlp);

  if( m_keepVarsInfo )
  {
    /// Clean timelim struct 
    if( m_varsinfo != NULL )
      delete m_varsinfo;
  }
}

  void
MyCplexSolver::enableConstraint( ConstraintType cons, bool isEnabled)
{
  m_isDisabled[cons] = !isEnabled;

}

  void 
MyCplexSolver::setMatheuristic(Matheuristic *math) 
{
  m_mathh = math;
  m_heuristicCallback = &mathPullCallback;
  m_incumbentCallback = &mathPushCallback;

  //	m_simpleCallback = &mathPushCallback;


}

  int 
MyCplexSolver::mathPullCallback(void *data)
{

  HeuristicCallbackData *cbinfo = (HeuristicCallbackData*)(data);

  printf("CPLEX: %.2f\n", cbinfo->timenow);
  printf("Pulling solution\n");
  bool valid;
  LpSolutionPtr lp_sol = m_mathh->pullSolution(&valid);
  if( valid)
  {
    printf("Got valid solution\n");
    cbinfo->heuristic = lp_sol;
  } else
  {
    printf("Not valid\n");
  }
  return valid;
}

  int 
MyCplexSolver::mathPushCallback(void *data)
{

  InfoCallbackData *cbinfo = (InfoCallbackData*)(data);
  //	m_mathh->testPush(cbinfo->bestobj);
  if( cbinfo->bestobj < m_mathh->lastPull())
  {
    int error = m_mathh->pushSolution(cbinfo->incumbent);
    if(error)
    {
      printf("Error pushing solution\n");
    }
    else
      printf("MyCplexSolver: Pushed %f\n", cbinfo->bestobj);
    return error;
  } else {
    printf("Attempting to push %f aborted - lastPull %f\n", 
	   cbinfo->bestobj, m_mathh->lastPull());
    return 1;
  }
}
