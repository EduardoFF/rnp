#ifndef _LPSOLVER_H
#define _LPSOLVER_H

#include "rnpsolution.h"
class Network;
class Metric;
class Graph;
class Matheuristic;
//struct CallbackInfo;

class MyLpSolver
{
	Network *m_net;
	int m_k;
	int m_minK;
	Metric *m_metric;
	RNPSolutionPtr m_sol;
	Matheuristic *m_matheuristic;

	public:


	MyLpSolver(Network *net,  
		   int maxK, 
		   int minK = 0) : 
		m_net(net),m_k(maxK),m_minK(minK) {
		m_metric = NULL;
		m_matheuristic = NULL;
	}
	void setMetric(Metric *);
	void setMatheuristic(Matheuristic *mth) { m_matheuristic = mth;}
	void writeDataFile(ofstream &);
	int simpleCallback(void *);
	bool solve();
	RNPSolutionPtr solution() {return m_sol;}

	void write_time_to_file(string filename);
	void write_solution_to_file(string filename);
	~MyLpSolver();
};




#endif

