#ifndef _RNPSOLUTION_H
#define _RNPSOLUTION_H
#include "main.h"

class LpSolution;

#define VAR_RELAY 'y'  //Char used in the model to define binary variables for relay X -- i.e y[X]
#define VAR_EDGE 'x'   //Char used in the model to define flow variables in edges (X,Y) -- i.e x[X, Y]
typedef boost::shared_ptr<LpSolution> LpSolutionPtr;
class RNPSolution
{
	LpSolutionPtr m_lpsol;
	vector<string> relaysToUse;
	vector< pair< pair< string, string>, double> > edgesToUse;
	int m_validVariables;
	set< pair<string, string> > edgeSet;
	std::map<string, double> outgoing_flow;
	bool m_hasDump;
	void write_to_sol(const char *varname, double value);
	void process_solution();

	public:
	RNPSolution(LpSolutionPtr);
	
	bool isOptimal();
	bool solved();
	int nValidVariables(){ return m_validVariables;}
	void writeTimeLog(string fname);
	void writeToFile(string filename);
	int numberRelays(){ return relaysToUse.size();}

	const vector<string> &getRelaysToUse(){ return relaysToUse;}
	const vector< pair< pair< string, string> , double > > getEdgesToUse(){
		return edgesToUse;
	}

	Graph* create_graph_from_sol(Network *net);
	double cpuTimeTotal();
	double cpuTimeSolving();
	LpSolutionPtr getLpSolution(){ return m_lpsol;}

	double getOutgoingFlow(string nodeid){ return outgoing_flow[nodeid];}
	bool isEdge(string nid1, string nid2){ return edgeSet.find(make_pair(nid1,nid2)) != edgeSet.end();}
	double getObjVal();
	~RNPSolution();

};
typedef boost::shared_ptr<RNPSolution> RNPSolutionPtr;

#endif
