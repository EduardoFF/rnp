
#include "main.h"
#include "generator.h"
#include "solver.h"
#include "graph.h"
#include "network.h"
#include "rnpsolution.h"
extern param prob;

enum vartypes{VARTYPE_NULL, VARTYPE_RELAY, VARTYPE_EDGE};

#define DRAW_DUMP 0

int get_var_type(const char *v){
	if(v[0] == VAR_RELAY)
		return VARTYPE_RELAY;
	else if(v[0] == VAR_EDGE)
		return VARTYPE_EDGE;
	else
		return VARTYPE_NULL;
}


string get_relay_id(const char *v){
	string var(v);
	return var.substr(var.find('[')+1,var.find(']')-var.find('[')-1);
}

pair<string, string> get_edge_from_var(const char* varname){
	string v(varname);
	string r1 = v.substr(v.find('[')+1,v.find(',') - v.find('[') -1);
	string r2 = v.substr(v.find(',')+1, v.find(']') - v.find(',') - 1 );
	return make_pair(r1,r2);
}



bool RNPSolution::isOptimal()
{
	return m_lpsol->isOptimal();
}

double RNPSolution::getObjVal()
{
	return m_lpsol->objval;
}

bool RNPSolution::solved()
{
	return m_lpsol->solved;
}

double RNPSolution::cpuTimeTotal()
{
	return m_lpsol->cpuTimeTotal();
}
double RNPSolution::cpuTimeSolving()
{
	return m_lpsol->cpuTimeSolving();
}
void RNPSolution::writeTimeLog(string filename)
{

	ofstream outf(filename.c_str());
	outf << "cpuTimeModelReading " << m_lpsol->cpuTimeModelReading() << endl;
	outf << "cpuTimeDataReading " << m_lpsol->cpuTimeDataReading() << endl;
	outf << "cpuTimeModelGeneration " << m_lpsol->cpuTimeModelGeneration() << endl;
	outf << "cpuTimeSolving " << m_lpsol->cpuTimeSolving() << endl;
	outf << "cpuTimeTotal " << m_lpsol->cpuTimeTotal() << endl;

	outf.close();
}

void
RNPSolution::writeToFile(string filename)
{

  ofstream outf(filename.c_str());
	
  if(prob.verbose)
    {
      cout << "Writing solution to " << filename << endl;
    }
  outf << "SOLVED\t" << m_lpsol->solved << endl;
  if(m_lpsol->solved)
    {
      outf << "OPTIMAL\t" << m_lpsol->isOptimal() << endl;
      outf << "OBJ_VAL\t" << m_lpsol->objval << endl;
      //outf << "N_VARS\t" << m_lpsol->value.size() << endl;
      outf << "MAX_RELAYS\t" << prob.K << endl;
      outf << "N_RELAYS\t" << numberRelays() << endl;
      outf << "N_EDGES\t" << edgesToUse.size() << endl;
      FOREACH(it,m_lpsol->value)
	{
	  string var = it->first;
	  outf << var << "\t" 
	       << m_lpsol->type[var] 
	       << "\t" << m_lpsol->value[var] 
	       << endl;
	}	
    }
  outf.close();
}

RNPSolution::RNPSolution(LpSolutionPtr lp_sol):m_lpsol(lp_sol)
{

	/** Process lp_sol to extract solution info
	 *  Before It was done in create_graph_from_sol
	 *  */
	m_validVariables = 0;
	m_hasDump = false;
	FOREACH(it,lp_sol->value){
		if(fabs(it->second) > EPSILON)
		{
			write_to_sol((it->first).c_str(),it->second);
			m_validVariables++;
		}
	}

	std::sort(relaysToUse.begin(), relaysToUse.end());
	std::vector<std::string>::iterator new_end_pos;
	new_end_pos = std::unique( relaysToUse.begin(), relaysToUse.end() );

	// The elements between new_end_pos and str.end() hold
	// hold their old values, and must be erased to complete
	// the operation:

	relaysToUse.erase( new_end_pos, relaysToUse.end() );

	/** Write additional info to lp_sol  */
	//lp_sol->n_relays = relaysToUse.size();
	//lp_sol->n_edges = edgesToUse.size();

	/** Output simple solution information */
	if(prob.verbose)
	{
		printf("########### solution #############\n");
		printf("number of relays: %d\n",relaysToUse.size());
		printf("number of edges: %d\n",edgesToUse.size());
		printf("objective function value: %f\n",lp_sol->objval);
	}
}

void RNPSolution::write_to_sol(const char *varname, double value){
	if(get_var_type(varname) == VARTYPE_EDGE){
		pair<string,string> e = get_edge_from_var(varname);
		DEBUG(1)cdebug << "Edge " << e.first << ", " << e.second << endl;
		if( e.second == "DUMP" )
		{
		  m_hasDump = true;
		  if( DRAW_DUMP )
		    edgesToUse.push_back(make_pair(e,value));
		} else
		    edgesToUse.push_back(make_pair(e,value));
		if(e.first[0] == 'r')
			relaysToUse.push_back(e.first);
		if(e.second[0] == 'r')
			relaysToUse.push_back(e.second);
		if(e.second[0] == 'r' && e.first[0] == 'r')
		{
			edgeSet.insert(make_pair(e.first, e.second));
		}
		if(value > 0)
		  outgoing_flow[e.first] += value;
	}
}


Graph* RNPSolution::create_graph_from_sol(Network *net)
{

	if(prob.verbose) printf("Building solution graph...\n");
	Graph *g = new Graph(net->numStatic()+net->numBases()+relaysToUse.size()+(m_hasDump&&DRAW_DUMP?1:0));
	map<string,int> idToIdx;
	int ix=0;
	for(int i=0;i<net->numStatic();i++){
		g->set_vertex_id(i,net->getStatic(i).id);
		g->set_vertex_loc(i,net->getStatic(i).x / net->dimX ,net->getStatic(i).y / net->dimY);
		g->set_vertex_type(i,net->getStatic(i).t);
		idToIdx[net->getStatic(i).id] = ix++;
	}
	for(int i=0;i<net->numBases();i++){
		g->set_vertex_id(ix,net->getBase(i).id);
		g->set_vertex_loc(ix,net->getBase(i).x / net->dimX ,net->getBase(i).y / net->dimY);
		g->set_vertex_type(ix,net->getBase(i).t);
		idToIdx[net->getBase(i).id] = ix++;
	}
	
	for(int i=0;i<relaysToUse.size();i++){
		Node &n = net->getNodeById(relaysToUse[i]);
		g->set_vertex_id(ix,n.id);
		g->set_vertex_loc(ix,n.x / net->dimX,n.y/net->dimY);
		g->set_vertex_type(ix,n.t);
		idToIdx[relaysToUse[i]] = ix++;
	}

	if( DRAW_DUMP )
	{
	  g->set_vertex_id(ix,"DUMP");
	  g->set_vertex_loc(ix,1,1);
	  g->set_vertex_type(ix,BASE);
	  idToIdx["DUMP"] = ix++;
	}
	for(int i=0;i<edgesToUse.size();i++){
		int u,v;
		if( !DRAW_DUMP )
		  if( edgesToUse[i].first.first == "DUMP" ||
		      edgesToUse[i].first.second == "DUMP" )
		    continue;
		u = idToIdx[edgesToUse[i].first.first];
		v = idToIdx[edgesToUse[i].first.second];

		g->addEdge(u,v);
		stringstream ss;

		ss << edgesToUse[i].second;
		g->set_edge_label(u,v,ss.str());
	}
	
	return g;

}

RNPSolution::~RNPSolution()
{
	//delete m_lpsol;
}


