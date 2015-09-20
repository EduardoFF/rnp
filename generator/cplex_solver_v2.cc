
#include "main.h"
#include "generator.h"
#include "solver.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "matheuristic.h"
#include "cplex_solver_v2.h"
#include "rnpsolution.h"

extern param prob;

namespace Cplexv2
{
	ofstream o_progress;
}

IloEnv g_env;
inline VariableX VarX(string s1, string s2)
{

	VariableX varx = make_pair(s1,s2);
	return varx;
}


void MyCplexSolverV2::addEdgeVar(string si, string sj)
{
	char varname[50];
	sprintf(varname,"X[%s,%s]",si.c_str(), sj.c_str());
	cout << "Added variable X["<< si << "," 
		<< sj << "]" << endl;
	m_vars.add(IloNumVar(g_env, 0, 1, IloNumVar::Bool, varname));
	VariableX varx = VarX(si, sj); 
	m_XVars[ varx ] = m_varCnt;
	m_rXVars[m_varCnt] = varx;
	m_varCnt++;
}



void MyCplexSolverV2::addRelayVar(string rid)
{

}

void MyCplexSolverV2::addFXVar(string id)
{
}

int MyCplexSolverV2::getFXVarIndex(string i)
{
}


int MyCplexSolverV2::getYVarIndex(string i)
{
}
int MyCplexSolverV2::getXVarIndex(string i, string j)
{
}

int MyCplexSolverV2::getBXVarIndex(string i, string j)
{
}

/*
 * s.t. flow{k in R}: sum{(k,j) in edges} x[k,j] = sum{(i,k) in edges} x[i,k];
 * */
void MyCplexSolverV2::addFlowCons()
{
}

/* s.t. flowedge{(i,j) in edges} : x[i,j] <= MAX_CAPACITY * bx[i,j];
 * */

void MyCplexSolverV2::addFlowEdgeCons()
{

}

/*
 * s.t. flowneighbour{i in  S}: 
 * (sum{ (i,j) in edges diff edgesR} sum{(j,k) in edges diff edgesR} x[j,k]) - MAX_LOC_FLOW <= fx[i] * MAX_FLOW;
 *
 * s.t. flowneighbor2{i in  S}: 
 * ((sum{ (i,j) in edges diff edgesR} sum{(j,k) in edges diff edgesR} x[j,k]) - MAX_LOC_FLOW) >= (fx[i] -1)*MAX_FLOW;
 */

void MyCplexSolverV2::addFlowNeighborCons()
{

}

/*
 * s.t. nodedegree{k in R union S}: (sum{(j,k) in edges} bx[j,k]) <= D;
 */
void MyCplexSolverV2::addNodeDegreeCons()
{
}

/*
 * s.t. bound: sum{k in R} y[k] <= K;
 * s.t. bound_low: sum{k in R} y[k] >= minK;
 */
void MyCplexSolverV2::addRelayBoundCons()
{
}


/*
 * s.t. base: sum {k in B} sum{(i,k) in edges} x[i,k] = sum{i in S} d[i];
 */
void MyCplexSolverV2::addBaseCons()
{
}

/*
 * s.t. innodes{k in S}: sum{(k,j) in edges} x[k,j] - sum{(i,k) in edges} x[i,k] = d[k];
 */
void MyCplexSolverV2::addDemandCons()
{
}


/*
 * s.t. relayselected{k in R} : FLOW_RELAY_CONSTRAINTS * ((sum{(i,k) in edges} x[i,k]) - y[k]) >= 0;
 */
void MyCplexSolverV2::addRelaySelCons()
{
}

/*
 * s.t. flowselected{k in R}: FLOW_RELAY_CONSTRAINTS * (sum{(i,k) in edges} x[i,k]) <= sum{i in S} d[i] * y[k];
 */
void MyCplexSolverV2::addFlowSelCons()
{

}

/*
 * (sum{(i,j) in edges} c[i,j] * x[i,j])
 * + sum{i in S union R} fx[i] * FLOW_NEIGHBOUR_PENALTY + sum{k in R} y[k] * PENALTY_RELAY;
 */
void MyCplexSolverV2::setObjFunc()
{

}

/*
 * s.t. bandwidth_capacity{i in nodes}: (sum{(i,j) in edges} x[i,j] + sum{(k,i) in edges} x[k,i]) <= MAX_CAPACITY;
 */

void MyCplexSolverV2::addBWCapacityCons()
{
}

void MyCplexSolverV2::setRelayBnds(int mink, int maxk)
{
}

void MyCplexSolverV2::setFlowRelayCons()
{
	addRelaySelCons();
}

void MyCplexSolverV2::initialize()
{
	CpuTime cpu_time;
	cpu_time.start();
	for(int i=0; i< m_netWithRelays->size();i++)
	{	
		Node &ni = m_netWithRelays->getNode(i);
		for(int j=i+1; j<m_netWithRelays->size();j++)
		{
			Node &nj = m_netWithRelays->getNode(j);

			if( m_netWithRelays->inRange(ni,nj))
			{
				addEdgeVar(ni.id,nj.id);
				addEdgeVar(nj.id,ni.id);
			}

		}
		if(ni.t == RELAY)
		{
			addRelayVar(ni.id);
		}
		else
		{
			addFXVar(ni.id);
		}
		
	}

	IloEnv   env;
	try {
		IloModel model(env);
		IloCplex cplex(env);


		IloObjective   obj;
		IloNumVarArray var(env);
		IloRangeArray  rng(env);
		IloSOS1Array   sos1(env);
		IloSOS2Array   sos2(env);

		//cplex.use(SOSbranch(env, sos1));
		cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);

		cplex.extract(model);
	//	cplex.solve();

	//	IloNumArray vals(env);
	//	cplex.getValues(vals, var);
	//	env.out() << "Solution status = " << cplex.getStatus() << endl;
	//	env.out() << "Solution value  = " << cplex.getObjValue() << endl;
	//	env.out() << "Values          = " << vals << endl;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();


	addFlowCons();
	addFlowEdgeCons();
	addBWCapacityCons();
	addFlowNeighborCons();
	addNodeDegreeCons();
	addDemandCons();

	addFlowSelCons();
	if( prob.ilp_model.flow_relay_constraints)
	{

		addRelaySelCons();
	}
	addBaseCons();
	addRelayBoundCons();

	
	if( prob.metric == "PRR")
		m_metric = new EstPRR_Metric(m_netWithRelays);
	else if(prob.metric == "DISTANCE_SQ")
		m_metric = new DistSQ_Metric(m_netWithRelays);
	

	setObjFunc();
	//enableAllRelays();

	m_cpuTimeModelGeneration = cpu_time.cpu_time_elapsed();

}

void MyCplexSolverV2::enableRelay(int rid)
{
	Node &nr = m_netWithRelays->getRelay(rid);
	enableRelay(nr.id);
}
void MyCplexSolverV2::enableRelay(string rid)
{
	}

void MyCplexSolverV2::fixRelay(string rid)
{
	}

void MyCplexSolverV2::fixRelay(int rid)
{

	Node &nr= m_netWithRelays->getRelay(rid);
	fixRelay(nr.id);
}

void MyCplexSolverV2::disableRelay(int rid)
{
}
void MyCplexSolverV2::enableAllRelays()
{
	//printf("Enabling all relays\n");
	for(int i=0; i< m_netWithRelays->numRelays(); i++)
	{
		enableRelay(i);
	}
}

void MyCplexSolverV2::disableAllRelays()
{
	for(int i=0; i< m_netWithRelays->numRelays(); i++)
	{
		disableRelay(i);
	}
}

MyCplexSolverV2::MyCplexSolverV2(Network *static_net, Network *net_with_relays):  
	m_staticNet(static_net),m_netWithRelays(net_with_relays)
{
	m_metric = NULL;
	//m_simpleCallback = NULL;


	initialize();
}

double MyCplexSolverV2::getObjVal()
{
	assert(m_sol.get() != NULL);
	return m_sol->getObjVal();
}

RNPSolutionPtr MyCplexSolverV2::solution()
{
	return m_sol;
}

bool MyCplexSolverV2::solve()
{
	//RNPSolutionPtr sol( new RNPSolution(ret_sol));
	//m_sol = sol;
	//return m_sol->solved();
}

MyCplexSolverV2::~MyCplexSolverV2()
{
	delete m_metric;
}

