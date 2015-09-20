#include "main.h"
#include "generator.h"
#include "solver.h"
#include "graph.h"
#include "network.h"
#include "metric.h"
#include "matheuristic.h"
#include "lpsolver.h"
#include "rnpsolution.h"
extern param prob;



class LpWrapper
{
	public:
		static MyLpSolver *lpsolver;
		LpWrapper(){}
		static void writeDataFile(ofstream& os) {
			LpWrapper::lpsolver->writeDataFile(os);
		}
		static int simpleCallback(void *cbinfo) {
			return LpWrapper::lpsolver->simpleCallback(cbinfo);
		}
};

MyLpSolver *LpWrapper::lpsolver = NULL;

/*
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

void MyLpSolver::write_to_sol(const char *varname, double value){
	if(get_var_type(varname) == VARTYPE_EDGE){
		pair<string,string> e = get_edge_from_var(varname);
		DEBUG(1)cdebug << "Edge " << e.first << ", " << e.second << endl;
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
*/

void MyLpSolver::setMetric(Metric *metric)
{
	m_metric = metric;
}

MyLpSolver::~MyLpSolver()
{
	if(m_metric)
	{
		delete m_metric;
	}
}
int MyLpSolver::simpleCallback(void *info)
{

	CallbackInfo *cbinfo = (CallbackInfo *)info;
	printf("simpleCallback called at time %.2f\n", cbinfo->timenow);
	//while(m_matheuristic->Synchronize(cbinfo->timenow)!=0){sleep(10);}
	return 0;
}
bool MyLpSolver::solve()
{
	SolverParams params;

	params.seed = prob.seed;
	if(prob.lp_params.solver == "GUROBI")
	{
		params.which_solver = SOLVER_GUROBI;
	}
	else if(prob.lp_params.solver == "CPLEX")
	{
		params.which_solver = SOLVER_CPLEX;
	}else
	{
		/// In other case, use CPLEX
		params.which_solver = SOLVER_CPLEX;
	}
	params.model_file = prob.lp_params.model_file;
	//cout << "Model " << params.model_file << endl;
	params.data_file = prob.lp_params.data_file;

	if(params.data_file == "default")
	{
		params.data_file = prob.lp_params.output_path + 
			prob.instance_id + ".lpdat";
	}
	params.output_path = prob.lp_params.output_path;
	params.mps_filename = prob.lp_params.output_path + 
		prob.instance_id + ".mps";
	params.get_suboptimal = prob.lp_params.get_suboptimal;
	params.time_limit = prob.lp_params.strong_timelim;
	params.verbose = prob.lp_params.verbose;
	params.output_progress = prob.lp_params.log_progress;
	params.do_postsolve = prob.lp_params.do_postsolve;
	params.allow_parallel = prob.lp_params.allow_parallel;
	params.noMIP = prob.lp_params.noMIP;
	params.strong_timelim = prob.lp_params.strong_timelim;


	if(params.output_progress)
	{
		params.output_progress_file = prob.lp_params.output_path +
			prob.instance_id + ".lp_sol";
	}
	

	if( m_metric == NULL)
	{
		if( prob.metric == "PRR")
			m_metric = new EstPRR_Metric(m_net);
		else if(prob.metric == "DISTANCE_SQ")
			m_metric = new DistSQ_Metric(m_net);
	}

	GLPSol *lp = new GLPSol(params);
	if( m_matheuristic != NULL)
	{
		printf("Setting callback\n");

		lp->setSimpleCallback(&LpWrapper::simpleCallback);
	}

	LpWrapper::lpsolver = this;

	LpSolution *lp_sol = lp->solve(&LpWrapper::writeDataFile);
	
	LpSolutionPtr lp_sol_ptr(lp_sol);
	RNPSolutionPtr sol(new RNPSolution(lp_sol_ptr));
	m_sol = sol;
	

	
	if( prob.verbose)
	{
		printf("Finished Solving\n");
	}
	if(!prob.lp_params.keep_datafile)
	{
		remove(params.data_file.c_str());
	}
	delete m_metric;
	m_metric = NULL;
	delete lp;
	return m_sol->solved();
}

void MyLpSolver::write_time_to_file(string filename)
{
	m_sol->writeTimeLog(filename);
}




void MyLpSolver::write_solution_to_file(string filename)
{
	m_sol->writeToFile(filename);

}





void MyLpSolver::writeDataFile(ofstream &dataFile)
{
	VERBOSE(1)
	{
		printf("MyLpSolver Creating data file...\n");
	}

	//string data_filename = prob.output_path + "/" + prob.instance_id + "-model.dat";
	//ofstream dataFile(data_filename.c_str());
	
	if( dataFile.is_open()){

		/*  Creating sets */
		dataFile << "set S :=";

		for(int i = 0; i < m_net->numStatic(); i++){
			dataFile << " s" << i;

		}
		dataFile << ";\n";
		dataFile << "set B :=";

		for(int i=0; i< m_net->numBases();i++){
			dataFile << " b" << i;
		}
		dataFile << ";\n";

		dataFile << "set R :=";
		for(int i=0;i<m_net->numRelays();i++){
			Node &n = m_net->getRelay(i);
			dataFile << " " << n.id;
		}
		dataFile << ";\n";

		/* Creating edges */
		dataFile << "set edges :=";

		for(int i=0;i < m_net->size(); i++){
			Node &ni = m_net->getNode(i);
			for( int j=i+1;j<m_net->size();j++){
				Node &nj = m_net->getNode(j);
				if( m_net->inRange(ni,nj)){
					dataFile << " ( " << ni.id << " , " << nj.id << ")";
					dataFile << " ( " << nj.id << " , " << ni.id << ")";
				}


			}	
		}
		dataFile << " ;\n";


		
		/* Setting costs for edges */
		
		/*
		dataFile << "param d_f := ";
		for(int i=0;i < m_net->size(); i++){
			Node &ni = m_net->getNode(i);
			for( int j=i+1;j<m_net->size();j++){
				Node &nj = m_net->getNode(j);
				if( m_net->inRange(ni,nj)){
					double link_cost;



					double link_d = m_net->distanceLink(ni,nj);
				//	double link_n = m_net->neighborsLink(ni,nj);

					double d_f = (max(link_d,2.5)/m_net->tx_range); 
				//	double n_f = (min(link_n,MAX_LINK_NEIGHBORS) / MAX_LINK_NEIGHBORS);

					//link_cost = pow(link_cost, 2.0);

					// Calculate cost based on edge length 
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << d_f << " "; 
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << d_f << " "; 


				}


			}	
		}
		dataFile << ";\n";
		
*/
		dataFile << "param c := ";
		for(int i=0;i < m_net->size(); i++){
			Node &ni = m_net->getNode(i);
			for( int j=i+1;j<m_net->size();j++){
				Node &nj = m_net->getNode(j);
				if( m_net->inRange(ni,nj)){
					double link_cost;
					link_cost = (*m_metric)(ni.id,nj.id);
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << link_cost << " ";
					link_cost = (*m_metric)(nj.id,ni.id);
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << link_cost << " ";
				}
			}
		}
		dataFile << ";\n";
		dataFile << "param AVERAGE_LINK_COST := " << m_metric->average();
		dataFile << ";\n";
		dataFile << "param MAX_LINK_COST := " << m_metric->max();
		// cout << "param MAX_LINK_COST := " << m_metric->max() << endl;


		/** Calculate expected hop-count */
		dataFile << ";\n";
		int expected_hop_count = 0;
		for(int i = 0; i< m_net->numStatic();i++){
			Node &ns = m_net->getStatic(i);
			int hop_count = 100000000;
			for(int j=0;j < m_net->numBases();j++){
				Node &nb = m_net->getBase(j);
				double d = m_net->distanceLink(ns,nb);
				int hc = (int)ceil(1.0*d/m_net->getRange());
				hop_count = min(hc,hop_count);
			}
			expected_hop_count += hop_count;
		}
		dataFile << "param EXPECTED_HOP_COUNT := " << expected_hop_count << endl;

		dataFile << ";\n";
		// cout << "param EXPECTED_HOP_COUNT := " << expected_hop_count << endl;
		
		
		dataFile << "param FLOW_RELAY_CONSTRAINTS := " << prob.ilp_model.flow_relay_constraints;

		
		
	
/*		
		dataFile << "param n_f := ";
		for(int i=0;i < m_net->size(); i++){
			Node &ni = m_net->getNode(i);
			for( int j=i+1;j<m_net->size();j++){
				Node &nj = m_net->getNode(j);
				if( m_net->inRange(ni,nj)){
					double link_cost;



				//	double link_d = m_net->distanceLink(ni,nj);
					double link_n = m_net->neighborsLink(ni,nj);

				//	double d_f = (max(link_d,2.5)/m_net->tx_range); 
					double n_f = (min(link_n,MAX_LINK_NEIGHBORS) / MAX_LINK_NEIGHBORS);

					//link_cost = pow(link_cost, 2.0);

					// Calculate cost based on edge length
					dataFile << "[ " << ni.id << ", " << nj.id << "] ";
					dataFile << n_f << " "; 
					dataFile << "[ " << nj.id << ", " << ni.id << "] ";
					dataFile << n_f << " "; 


				}


			}	
		}
*/
		dataFile << ";\n";
	

		dataFile << "param K := " << m_k << ";\n";
		dataFile << "param minK := " << m_minK << ";\n";
		dataFile << "end;\n";





	}



}


