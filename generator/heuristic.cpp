/*
 * =====================================================================================
 *
 *       Filename:  heuristic.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/08/2010 06:29:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "main.h"
#include "generator.h"
//#include "Plotter.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "lpsolver.h"
#include "solver.h"
#include "rnpsolution.h"

#include "heuristic.h"

extern param prob;
void generate_grid(Network *net, double dimX, double dimY, double r );
void output_graphs(Graph *G, string filename);
Network * convex_hull_restriction(Network *);

LpSolution *Low_res_heuristic::update_solution(LpSolution *cnt_sol, map<int,int> index_map){

	using namespace boost::xpressive;
	int cnt = 0;

 	LpSolution *ret = new LpSolution();

	/*printf("Relay ids substitutions\n");
	FOREACH(it,index_map){
		cout << it->first << " -> " << it->second << endl;
	}*/
        sregex expression, variable, listnodeid, nodeid, nodeindex, nodetype;

/*  	nodetype = alpha;
	nodeindex = +digit;
	nodeid = nodetype >> nodeindex;
	listnodeid = (nodeid) | (nodeid) >> ',' >> listnodeid;
	variable = +alpha;
	expression = variable >> '[' >> listnodeid >> ']';
*/


	sregex separators;
	separators = as_xpr(',') | '[' | ']';
	
	FOREACH(it,cnt_sol->value){
		// match regex
		string var(it->first);

		if(it->second > 0)
			DEBUG(1)
				cdebug << var << ":" << it->second << endl;

		// for each match, the token iterator should first take the value of
		// the first marked sub-expression followed by the value of the second
		// marked sub-expression
		int const subs[] = { 1, 2 };

		sregex_token_iterator cur( var.begin(), var.end(), separators,-1);
		sregex_token_iterator end;
		std::ostringstream os;
 
		os << *cur << "[";
		cur++;

		bool fi = false;
		for( ; cur != end; ++cur )
		{
			if(fi)
				os << ",";
			else
				fi = true;
			string toparse = *cur;
			if(toparse[0]=='r'){
				int ix = atoi(toparse.substr(1,toparse.size()-1).c_str());
				int nix = index_map[ix];
				os << "r" << nix;

			}else
				os << *cur;
		}
		os << "]";
		//cout << os.str() << endl;
		//

		if(it->second > 0){
			DEBUG(1)
				cdebug << "replacing " << var << " with " << os.str() << endl;
			cnt++;
		}

		ret->value[os.str()] = it->second;




		//std::ostream_iterator< std::string > out_iter( std::cout, "\n" );
		//std::copy( cur, end, out_iter );


	}
	DEBUG(1)
		cdebug << "Replaced " << cnt << " variables " << endl;
	return ret;


	


}

Low_res_heuristic::Low_res_heuristic(LpSolution *prev_sol, double actual_res, int factor, bool silent){
	/* Create map for relay indexes */
	int rx,ry; // number of relays per column, row in simplified problem
	double iter_res = actual_res*factor;
	rx = (int)floor(prob.dimX/iter_res) + 1;
	ry = (int)floor(prob.dimY/iter_res) + 1;

	int nrx, nry; // number of relays per column, row in the actual problem
	nrx = (int)floor(prob.dimX/(actual_res)) + 1;
	nry = (int)floor(prob.dimY/(actual_res)) + 1;

	map<int,int> relay_index;

	for(int idx=0;idx<rx*ry;idx++){
		int new_idx;
		int rel_col = factor * (idx/ry);
		int rel_row = factor * (idx%ry);
		new_idx = rel_col*nry + rel_row;
		relay_index[idx] = new_idx;
	}


	hvalue = update_solution(prev_sol,relay_index);
	int cnt=0;

	FOREACH(jt,hvalue->value){
		if(jt->second > 0)
			cnt++;
	}
	DEBUG(1)
		cdebug << "Heuristic sol has " << cnt << " variables > 0" << endl;



}

Low_res_heuristic::Low_res_heuristic(Network *net, double actual_res, int factor, bool silent){
	Network *iter_net = new Network(*net);

	// Generate grid

	double iter_res = actual_res*factor;
	generate_grid(iter_net, net->dimX, net->dimY, iter_res);
	VERBOSE(1) printf("Low_res_heuristic (%lf) number of relays %d\n",iter_res, iter_net->numRelays());
	if(prob.use_convex_hull){
		Network *it_net_hull = convex_hull_restriction(iter_net);
		delete iter_net;
		iter_net = it_net_hull;
		VERBOSE(1) printf("Low_res_heuristic (%lf) number of relays (hull) %d\n",
				  iter_res,iter_net->numRelays());
	}
	/*  Display resulting graph */

	std::ostringstream os;
	Graph *G;
	if(!silent){
		G = iter_net->createGraph();

		os << "Low_res_heu-" << iter_res << "-static+grid";

		output_graphs(G,os.str());
	}

	MyLpSolver lp(iter_net, prob.K);

	//LpSolution *sol = lp.solve(prob.K,NULL);
	lp.solve();
	LpSolutionPtr sol = lp.solution()->getLpSolution();

	if(!silent){
		//G = lp.create_graph_from_sol();
		//os.str("");
		//os << "Low_res_heu-" << iter_res << "-sol";
		//output_graphs(G,os.str());
	}

	/* Create map for relay indexes */
	int rx,ry; // number of relays per column, row in simplified problem
	rx = (int)floor(net->dimX/iter_res) + 1;
	ry = (int)floor(net->dimY/iter_res) + 1;

	int nrx, nry; // number of relays per column, row in the actual problem
	nrx = (int)floor(prob.dimX/(actual_res)) + 1;
	nry = (int)floor(prob.dimY/(actual_res)) + 1;

	map<int,int> relay_index;

	for(int idx=0;idx<rx*ry;idx++){
		int new_idx;
		int rel_col = factor * (idx/ry);
		int rel_row = factor * (idx%ry);
		new_idx = rel_col*nry + rel_row;
		relay_index[idx] = new_idx;
	}


	hvalue = update_solution(sol.get(),relay_index);
	int cnt=0;

	FOREACH(jt,hvalue->value){
		if(jt->second > 0)
			cnt++;
	}
	DEBUG(1)
		cdebug << "Heuristic sol has " << cnt << " variables > 0" << endl;


}

LpSolution *Low_res_heuristic::operator()(){
	return hvalue;
}
