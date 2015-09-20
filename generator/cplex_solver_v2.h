#ifndef _CPLEXSOLVERV2_H
#define _CPLEXSOLVERV2_H

class Network;
class Metric;
class Graph;
class Matheuristic;
//struct CallbackInfo;

//#include <ilcplex/cplex.h>
#include "ilcplex/ilocplex.h"
#include "rnpsolution.h"

typedef int(*CPXsimpleCallback)(void *);


typedef std::pair<string, string> VariableX;

class MyCplexSolverV2
{
	Network *m_staticNet;
	Network *m_netWithRelays;

	int m_varCnt;
	IloNumVarArray m_vars;

	std::map< VariableX, int > m_XVars;
	std::map< int, VariableX > m_rXVars;
//	IloEnv   env


//	CPXsimpleCallback m_simpleCallback;
//	Matheuristic *m_mathh;


	 
	double m_cpuTimeModelGeneration;



	int m_k;
	int m_minK;
	Metric *m_metric;
	RNPSolutionPtr m_sol;
	void initialize();
	void addEdgeVar(string i, string j);
	void addRelayVar(string i);
	void addFXVar(string i);

	void addFlowCons();
	void addFlowEdgeCons();
	void addBWCapacityCons();
	void addFlowNeighborCons();
	void addNodeDegreeCons();
	void addDemandCons();
	void addFlowSelCons();
	void addRelaySelCons();
	void addBaseCons();
	void addRelayBoundCons();
	void addRangedCons();


	void setObjFunc();


//	void setCallback(CPXsimpleCallback cb){ m_simpleCallback = cb;}

	int getXVarIndex(string i, string j);
	int getBXVarIndex(string i, string j);
	int getFXVarIndex(string i);
	int getYVarIndex(string i);

	int m_relayUpperBndRow;
	int m_relayLowerBndRow;

	std::map<string, double> outgoing_flow;
//	void write_to_sol(const char *varname, double value);
//	void process_soution();

	public:
	MyCplexSolverV2(Network *static_net, Network *net_with_relays);
	bool solve();
	void enableRelay(int rid);
	void enableRelay(string rid);

	void disableRelay(int rid);
	void enableAllRelays();
	void disableAllRelays();
	void fixRelay(int);
	void fixRelay(string);
	void setFlowRelayCons();

	~MyCplexSolverV2();

	void setRelayBnds(int low, int max);

	double getObjVal();
	RNPSolutionPtr solution();
};





#endif
