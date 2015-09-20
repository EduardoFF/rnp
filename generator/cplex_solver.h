#ifndef _CPLEXSOLVER_H
#define _CPLEXSOLVER_H

class Network;
class Metric;
class Graph;
class Matheuristic;
//struct CallbackInfo;

#include <ilcplex/cplex.h>
#include "rnpsolution.h"

typedef int(*CPXsimpleCallback)(void *);

typedef std::map<int, std::pair<string, int> > CPX_VARIABLESINFO;

class MyCplexSolver
{
  Network *m_staticNet;
  Network *m_netWithRelays;
  typedef enum {
    FLOW_NEIGHBOR_CONST,
    NODE_DEGREE_CONST,
  } ConstraintType;

  std::map< ConstraintType, bool> m_isDisabled;

  CPXENVptr m_cpxenv;

  CPXLPptr m_cpxlp;

  CPXsimpleCallback m_simpleCallback;
  CPXsimpleCallback m_heuristicCallback;
  CPXsimpleCallback m_incumbentCallback;
  static Matheuristic *m_mathh;



  double m_cpuTimeModelGeneration;


  std::map< pair<string, string>, int> m_varX;
  std::map< string, int> m_varFX;
  std::map< pair<string, string>, int> m_varBX;
  std::map< string, int> m_varY;

  std::map< pair<string, string>, int > m_varXD;
  std::map< pair<string, string>, int > m_varXF;
//  int m_varD;

  //! LF is an auxiliary variable (continuous) that
  //! indicates the amount of flow surrounding a node
  //! -> i.e., "occupancy of the channel"
  std::map< std::string, int > m_varLF;
  
  int m_k;
  int m_minK;
  Metric *m_metric;
  bool m_useInitial;

  bool m_keepVarsInfo; /// useful for GA
  CPX_VARIABLESINFO *m_varsinfo;

  RNPSolutionPtr m_initSol;

  RNPSolutionPtr m_sol;
  void initialize();
  void addEdgeVar(string i, string j);
  void addRelayVar(string i);
  void addFXVar(string i);
  void addLFVar(std::string i);
  /// misc / auxiliary variables
  void addOtherVars();
  
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

  void addLocalFlowCons();
  void addImFlowNeighborCons();
  void addStrictFlowNeighborCons();
  void addSinkFlowCons();

  static int mathPushCallback(void *);
  static int mathPullCallback(void *);
  void setObjFunc();


  void setCallback(CPXsimpleCallback cb){ m_simpleCallback = cb;}

  int getXVarIndex(string i, string j);
  int getBXVarIndex(string i, string j);
  int getFXVarIndex(string i);
  int getLFVarIndex(string i);
  int getYVarIndex(string i);

  int m_relayUpperBndRow;
  int m_relayLowerBndRow;

  std::map<string, double> outgoing_flow;
  //	void write_to_sol(const char *varname, double value);
  //	void process_soution();

  public:
  MyCplexSolver(Network *static_net, Network *net_with_relays);
  bool solve();
  void enableRelay(int rid);
  void enableRelay(string rid);

  void enableConstraint( ConstraintType, bool);

  void disableRelay(int rid);
  void enableAllRelays();
  void disableAllRelays();
  void fixRelay(int);
  void fixRelay(string);
  void setFlowRelayCons();

  ~MyCplexSolver();

  void setRelayBnds(int low, int max);
  void setMatheuristic(Matheuristic *);

  double getObjVal();
  RNPSolutionPtr solution();
  void useInitial(RNPSolutionPtr);

  void keepVarsInfo(bool keep);
};

#endif
