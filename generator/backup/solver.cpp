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

#include "generator.h"
#include "graph.h"
#include "solver.h"
#include "network.h"
extern "C" {
#include "../wrapper/wrapper.h"
}
//extern void solve_glp_grb(glp_prob *mip);
extern Network net;
extern param prob;

extern Graph *G;

void exit_on_error(string msg, glp_prob *_P, glp_tran *_tran)
{
  fprintf(stderr, "msg.c_str()");

  if(_tran != NULL)
    glp_mpl_free_wksp(_tran);
 
  glp_delete_prob(_P);

  exit(-1);
}
void free_space(glp_prob *_P, glp_tran *_tran)
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

int GLPSol::createDataFile(int K){


	VERBOSE(1)printf("Creating data file...\n");
	ofstream dataFile ("model.dat");
	if( dataFile.is_open()){

		/*  Creating sets */
		dataFile << "set S :=";

		for(int i = 0; i < net.numStatic(); i++){
			dataFile << " s" << i;

		}
		dataFile << ";\n";
		dataFile << "set B :=";
	
		for(int i=0; i< net.numBases();i++){
			dataFile << " b" << i;
		}
		dataFile << ";\n";

		dataFile << "set R :=";
		for(int i=0;i<net.numRelays();i++){
			dataFile << " r" << i;
		}
		dataFile << ";\n";
		
		/* Creating edges */
		dataFile << "set edges :=";
		graph_traits<Graph>::edge_iterator e_first,e_end,e;
		tie(e_first,e_end) = edges(*G);
		for(e = e_first; e != e_end; e++){
			graph_traits < Graph >::vertex_descriptor
				u = source(*e, *G), v = target(*e, *G);
			dataFile << " ( " << net.getNode(u).id << " , " << net.getNode(v).id << ")";
			dataFile << " ( " << net.getNode(v).id << " , " << net.getNode(u).id << ")";

		}
		dataFile << " ;\n";
		dataFile << "param K := " << K << ";\n";
		dataFile << "end;\n";
		





		return 0;
		
	}else
		return -1;


	

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
	if(get_var_type(varname) == VARTYPE_RELAY){
		string rid = get_relay_id(varname);
		DEBUG(1)cout << "Relayid " << rid << endl;
		relaysToUse.push_back(rid);
	}else if(get_var_type(varname) == VARTYPE_EDGE){
		pair<string,string> e = get_edge_from_var(varname);
		DEBUG(1)cout << "Edge " << e.first << ", " << e.second << endl;
		edgesToUse.push_back(e);
	}else
		printf("GLPSol: Invalid variable type\n");



}

Graph *GLPSol::solve(int K){


	glp_prob *lp;
	glp_tran *tran;
	glp_iocp params;
	int ret;
	lp = glp_create_prob();
	tran = glp_mpl_alloc_wksp();
	ret = glp_mpl_read_model(tran, "model.mod", 1);
	if (ret != 0) 
		printf("Error on translating model\n");
	else
		printf("Model read: OK\n");

	if(createDataFile(K))
		printf("Error creating data file\n");

	ret = glp_mpl_read_data(tran, "model.dat");
	if (ret != 0) 
		printf("Error on translating data\n");
	else
		printf("Data read: OK\n");

	if (glp_mpl_generate(tran, NULL) != 0)
		exit_on_error("Error on generating model\n", lp, tran);


	glp_mpl_build_prob(tran, lp);

	if(prob.use_gurobi){

		solve_glp_grb(lp, prob.glpk_out, prob.grb_out);
	}else{



		glp_simplex(lp, NULL);

		double obj_val = glp_get_obj_val(lp);
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
		glp_init_iocp(&params);
		params.presolve = GLP_OFF;

		//glp_simplex(P, NULL);


		glp_intopt(lp, &params);
	}




	glp_mpl_postsolve(tran, lp, GLP_MIP);
	lpx_print_sol(lp, "solution.sol");

	int num_integer_variables = 0;
	int binary_vars = 0;
	int ncols = glp_get_num_cols(lp);

	for(int c=1; c <= ncols; c++)
	{
		const char * col_name = glp_get_col_name(lp,c);

		VERBOSE(2){
			printf("Var %d is a %s, %s, with bounds [%g, %g]", c, (var_status_to_string(lp,c)).c_str(), 
			       (var_kind_to_string(lp,c)).c_str(),  glp_get_col_lb(lp, c),  glp_get_col_ub(lp, c));
			if(col_name != NULL)
				printf(" name [%s]",col_name);
		}
		double s = glp_mip_col_val(lp,c);
		VERBOSE(2){	
			printf(" value [%f]",s);
			printf("\n");
		}

		if(s > 0)
			write_to_sol(col_name, s);

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


  VERBOSE(2){
	  printf("Relays to use\n");
	  PV(relaysToUse,string);
	  printf("\n");
	  printf("Edges to use\n");
	  for(int i=0;i<edgesToUse.size();i++)
		  printf(" (%s %s)",edgesToUse[i].first.c_str(), edgesToUse[i].second.c_str());
	  printf("\n");
  }
   
  printf("Building solution graph...");
  Graph *g = new Graph(net.numStatic()+net.numBases()+relaysToUse.size());
  map<string,int> idToIdx;
  for(int i=0;i<net.numStatic();i++){
	  g->set_vertex_id(i,net.getStatic(i).id);
	  idToIdx[net.getStatic(i).id] = i;
  }
  for(int i=0;i<net.numBases();i++){
         g->set_vertex_id(i + net.numStatic(),net.getBase(i).id);
         idToIdx[net.getBase(i).id] = i + net.numStatic();
  }
  for(int i=0;i<relaysToUse.size();i++){
	  g->set_vertex_id(i+net.numStatic()+net.numBases(),relaysToUse[i]);
	  idToIdx[relaysToUse[i]] = i+net.numStatic()+net.numBases();
  }

  for(int i=0;i<edgesToUse.size();i++)
	g->addEdge(idToIdx[edgesToUse[i].first],idToIdx[edgesToUse[i].second]);


  /** Output simple solution information */
  printf("########### solution #############\n");
  printf("number of relays: %d\n",relaysToUse.size());
  printf("number of edges: %d\n",edgesToUse.size());
  printf("objective function value: %f\n",glp_get_obj_val(lp));
  


  return g;


 


  
 

	


	/*  
	LPX *lp;  // LP Problem
	
	lp = lpx_create_prob();

	int ret = lpx_read_model(lp, NULL, "model.mod");
        if (ret != 0){
		cout << "Solver: Error reading model\n";
	}


	//lpx_set_prob_name(lp, "prob_name");

*/

}	



