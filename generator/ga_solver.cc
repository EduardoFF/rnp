#include "main.h"
#include "generator.h"
#include "Plotter.h"
#include "graph.h"
#include "network.h"
#include "metric.h"

#include "solver.h"
#include "lpsolver.h"
#include "CpuTime.h"
#include "heuristic.h"
#include "ga_solver.h"
#include "Random.h"
#include "matheuristic.h"
#include "cplex_solver.h"
#include "rnpsolution.h"

#include <ga/GASimpleGA.h>	// we're going to use the simple GA
#include <ga/GA2DBinStrGenome.h> // and the 2D binary string genome

/// for detecting chains
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

float Objective(GAGenome &);	// This is the declaration of our obj function.
// The definition comes later in the file.

extern param prob;
Network *ga_network_with_relays;
Network *ga_network_static;


MyCplexSolver *g_cplexSolver;

/// This maps serves to compute
/// the mutations for the relay locations
std::map<int, vector<int> > adyacent_relays; 

/// Stores the evaluated genomes to prevent 
/// re-evaluations
/// O(log(N)) to retrieve, which for N < 10000 is OK
std::map<string, double> evaluated_genomes;


std::map<int, set<int> > partner_relays;
UniformRandom *u_rand;
static double MUTATION_DIST_TH = 2 * prob.tx_range;

// For initializations and random selection of relays
vector<int> initializeRelaySet;

/// For the 3 way crossover
set<pair<double, int> > candidateRelaySet;
set<int> candidateRelays;

/// seem to be DEPRECATED
set<int> g_redList;
set<int> g_greenList;

std::map<int, double> max_flow_for_relay;
double min_candidate_flow=0;
int min_candidate = -1;


std::deque<MyGenome *> expat_set;

#define GA_UNDEFINED_VAL 10e7

typedef struct {
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_relay_flow;
  double max_relay_flow;
  int max_candidates;
  set<int> relays;
  set<int> candidates;
  int min_candidate;
} CandidateRelayRegion;

std::vector< CandidateRelayRegion > g_candidateRegions;
std::map< int, vector<int> > g_relayToRegions;

/// maps relay to set of sets of relays which have appear in a chain
std::map<int, set< set<int> > > g_relayChains;

void generate_simulation(Graph *G, Network *net, string filename);
void output_graphs(Graph *G, string filename);
int MyGenome::repairs = 0;

/// This amount of time corresponds to time to be deducted from 
/// total cpu time of GA algorithm
CpuTime g_cpuTime;
double g_cpuTimeToDeduct = 0;

/// Time spent solving MILP model for evaluations
double g_cpuTimeSolvingMILP = 0;

double g_cpuTimeNextFlush = prob.ga_params.report_score_interval;
double g_bestScore = GA_UNDEFINED_VAL;
RNPSolutionPtr g_bestGARNPSol;
int g_nEvaluations = 0;
int g_nUnfeasible = 0;
int g_nPopulations = 0;
int g_nUsedPartners = 0;
int g_nDetectedPartners = 0;
int g_nDetectedChains = 0;
int g_nUsedCandidates = 0;
int g_nUsedChains = 0;
int g_nConflictsInCrossover=0;

ostream *g_olog;
ostream *g_odebug;
ofstream o_population_log;
ofstream o_generation_log;
Matheuristic *g_matheuristic = 0;
LpSolutionPtr g_lastPulledSol;
double g_lastPulledSolValue = GA_UNDEFINED_VAL;

double g_lastPushedSolValue = GA_UNDEFINED_VAL;

ofstream o_objective_timeline;

std::map< int, set<int> > g_conflictMap;

RNPSolutionPtr EvaluateUsingCplexSolver(MyGenome &);
RNPSolutionPtr EvaluateUsingWrapper(MyGenome &);
RNPSolutionPtr EvaluateUsingDummy(MyGenome &);
double fixRNPSolution(MyGenome &, RNPSolutionPtr sol);


  int 
myCallback()
{

  double solvingTime = g_cpuTime.cpu_time_elapsed();
  double realSolvingTime = solvingTime - g_cpuTimeToDeduct;
  if(g_matheuristic)
  {
    printf("GA: %.2f\n", realSolvingTime);
    printf("Pulling solution\n");
    bool valid;
    LpSolutionPtr lp_sol = g_matheuristic->pullSolution(&valid);
    if( valid)
    {
      g_lastPulledSol = lp_sol;
      printf("Got valid solution\n");
      double obj_val = lp_sol->objval;
      g_lastPulledSolValue = obj_val;
      printf("GA Solver: Pulled %f\n", obj_val);
      if( obj_val < g_bestScore || g_bestScore == GA_UNDEFINED_VAL)
      {
	// Introduce new guy
	MyGenome *mygen = new MyGenome(prob.ga_params.genome_size);
	mygen->fromLpSolution(lp_sol);

	expat_set.push_back(mygen);
	mygen->evaluate();
	printf("New expat. Expat size %d\n", expat_set.size());
      }
      else
      {
	printf("GA has better score %f (%f)\n", g_bestScore, obj_val);
	if(prob.ga_params.always_accept_expat)
	{
	  // Introduce new guy
	  MyGenome *mygen = new MyGenome(prob.ga_params.genome_size);
	  mygen->fromLpSolution(lp_sol);
	  expat_set.push_back(mygen);
	  printf("New expat. Expat size %d\n", expat_set.size());
	}

      }


    } else
      printf("Not valid solution\n");
#if 0
    double testvalue;

    if(g_matheuristic->testPull(&testvalue))
    {
      printf("GA_SOLVER: TestValue %f\n", testvalue);
    }else
    {
      printf("GA_SOLVER: Pull failed\n");
    }


    //	while(g_matheuristic->Synchronize((int)realSolvingTime) != 0){sleep(10);}
#endif
    return 0;
  }
}

  void 
GABuildInitializeSet()
{
  initializeRelaySet.clear();
  int nrelays = ga_network_with_relays->numRelays();
  for(int i=0; i< nrelays;i++)
  {
    initializeRelaySet.push_back(i);
  }
  /// Shuffle 3 times

  std::random_shuffle( initializeRelaySet.begin(), initializeRelaySet.end());
  std::random_shuffle( initializeRelaySet.begin(), initializeRelaySet.end());
  std::random_shuffle( initializeRelaySet.begin(), initializeRelaySet.end());
  // Done
}


  void 
initializeCandidateRegions()
{
  g_candidateRegions.clear();
  int c_x = prob.ga_params.candidate_region_x;
  int c_y = prob.ga_params.candidate_region_y;
  int n_regions = c_x * c_y;
  printf("Initializing %d regions\n", n_regions);

  int max_candidates = prob.ga_params.candidate_set_size / n_regions;

  printf("candidates for each region %d\n", max_candidates);

  double d_x = 1.0*prob.dimX / c_x;
  double d_y = 1.0*prob.dimY / c_y;
  double x = 0.0;
  for(int x=0; x < c_x; x++)
  {
    for(int y=0; y < c_y; y++)
    {
      double xx = x*d_x;
      double yy = y*d_y;
      CandidateRelayRegion cR;
      cR.min_x = xx;
      cR.max_x = xx + d_x;
      cR.min_y = yy;
      cR.max_y = yy + d_y;
      cR.min_relay_flow = GA_UNDEFINED_VAL;
      cR.max_relay_flow = 0.0;
      cR.max_candidates = max_candidates;
      g_candidateRegions.push_back(cR);


      if(prob.ga_params.debug)
      {
	(*g_odebug) << "CandidateRegion " << endl;
	(*g_odebug) << "min_x " << cR.min_x << endl;
	(*g_odebug) << "max_x " << cR.max_x << endl;
	(*g_odebug) << "min_y " << cR.min_y << endl;
	(*g_odebug) << "max_y " << cR.max_y << endl;
	(*g_odebug) << "max_candidates " << cR.max_candidates << endl;
      }
    }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) << g_candidateRegions.size() << " CandidateRegions " << endl;
  }


  for(int i=0; i< ga_network_with_relays->numRelays(); i++)
  {
    Node &nr = ga_network_with_relays->getRelay(i);
    bool associated = false;
    for(int k=0; k< g_candidateRegions.size(); k++)
    {
      CandidateRelayRegion &cR = g_candidateRegions[k];
      /// If relay is inside
      if( cR.min_x <= nr.x && nr.x <= cR.max_x && cR.min_y <= nr.y && nr.y <= cR.max_y)
      {
	associated = true;
	g_relayToRegions[i].push_back(k);
	cR.relays.insert(i);
      }
    }
    if( !associated)
    {
      fprintf(stderr, "Relay %s not associated to any region\n", nr.id.c_str());
      exit(1);
    }
  }

  for(int k=0; k< g_candidateRegions.size(); k++)
  {
    CandidateRelayRegion &cR = g_candidateRegions[k];
    if(prob.ga_params.debug)
    {
      (*g_odebug) << "CandidateRegion " << k << " has " 
	<< cR.relays.size() << " relays" << endl;
    }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "CandidateRegions initialized" << endl;
  }
}


  string 
GAGetRelayId(int r)
{
  Node &nr = ga_network_with_relays->getRelay(r);
  return nr.id;
}

  GABoolean
MyGATerminator(GAGeneticAlgorithm & ga)
{
  double solvingTime = g_cpuTime.cpu_time_elapsed();
  double realSolvingTime = solvingTime - g_cpuTimeToDeduct;


  MyGenome& best = (MyGenome&)(ga.statistics().bestIndividual());

  o_generation_log << "GAPOP " << g_nPopulations 
    << " " << realSolvingTime << " " << best.score() << " ";
  o_generation_log << best.toString(" ") << endl;

  printf("flushFreq %d\n", ga.flushFrequency());
  printf("scoreFreq %d\n", ga.scoreFrequency());
  ga.flushScores();
  o_generation_log << ga.statistics() << endl;

  o_generation_log.flush();
  g_nPopulations++;

  //if(realSolvingTime > g_cpuTimeNextFlush)
  //
  /* Done at Evaluate
     if(solvingTime > g_cpuTimeNextFlush)
     {
     double bestScore = ga.statistics().bestIndividual().score();
     o_objective_timeline << realSolvingTime << " " << bestScore << endl;	
     printf("(%f) %f: %f\n",solvingTime,realSolvingTime,bestScore);
     g_cpuTimeNextFlush += prob.ga_params.report_score_interval;
     }
     */

  //printf("(%f) %f\n",solvingTime,realSolvingTime);
  if(realSolvingTime > prob.ga_params.time_limit)
  {
    if(prob.ga_params.debug)
    {
      (*g_odebug) << "Time limit reached -> terminating" << endl << std::flush;
    }
    return gaTrue;
  }
  else
  {
    return gaFalse;
  }
}

void MyDefaultInitializer(GAPopulation & p){
  /** DO ABSOLUTELY NOTHING!!! */
}
GASol::GASol(Network *net_s, Network *net_r)
{
  ga_network_with_relays = net_r;
  ga_network_static = net_s;
  for(int i=0; i< net_r->numRelays(); i++)
  {
    Node &r1 = net_r->getRelay(i);
    for(int j=i+1; j<net_r->numRelays(); j++)
    {
      Node &r2 = net_r->getRelay(j);
      double dist = net_r->distanceLink(r1,r2);
      if( dist <= prob.ga_params.adjacent_txrange * prob.tx_range)
      {
	adyacent_relays[i].push_back(j);
	adyacent_relays[j].push_back(i);
      }
    }
  }
  u_rand=new UniformRandom(prob.ga_params.seed);

  m_cpuTimeSolving = 0;
  m_matheuristic = NULL;
  m_hasSolution = false;
  GARandomSeed(prob.ga_params.seed);


  //
  //	for(int i=0; i< net_r->numRelays(); i++)
  //	{
  //		printf("relay %d has %d adjacents\n", i, adyacent_relays[i].size());
  //	}



}

int 
getRandomRelay() {

  /// If empty, create initializeRelaySet
  if(!initializeRelaySet.size())
  {
    GABuildInitializeSet();

    if(prob.ga_params.debug)
    {
      (*g_odebug) << "RandomRelaySet initialized" << endl;
    }
  }

  int g = initializeRelaySet.back();
  initializeRelaySet.pop_back();
  return g;
}


void 
randomGenome(MyGenome &mg, int k)
{
  mg[0] = k;

  for(int i=1; i<=k; i++)
  {
    int g = getRandomRelay();
    mg.gene(i, g);
  }
}


double GASol::solve()
{

  /// Start measuring cpu time
  g_cpuTime.start();

  GAParameterList params;
  GAGeneticAlgorithm *ga;
  GAScalingScheme *ga_scaling_scheme;
  GASelectionScheme *ga_selection_scheme;
  if( prob.ga_params.genome_size <= 0 )
  {
    prob.ga_params.genome_size = prob.K;
  }
  MyGenome genome(prob.ga_params.genome_size);

  MyCplexSolver *cplex_solver_k0 = NULL;

  ofstream flog;
  ofstream fdebug;

  if( prob.ga_params.log_to_file)	
  {
    string ga_log_file = prob.output_path + prob.ga_params.log_file;

    flog.open(ga_log_file.c_str());
    cout << "Logging to file " << ga_log_file << endl;
    g_olog = &flog;
    assert(flog.is_open());
  }
  else
  {
    g_olog = &cout;
  }
  (*g_olog) << "Logging initialized" << endl;


  if( prob.ga_params.debug)	
  {
    string ga_dbg_file = prob.output_path + prob.ga_params.debug_file;

    fdebug.open(ga_dbg_file.c_str());
    cout << "Debugging info to file " << ga_dbg_file << endl;
    g_odebug = &fdebug;
    assert(fdebug.is_open());
  }
  else
  {
    g_odebug = &cout;
  }



  if( prob.ga_params.use_candidate_relays)
    initializeCandidateRegions();

  if(m_matheuristic)
  {
    // set global variable
    // Not nicest way to do it,
    // but it should work
    g_matheuristic = m_matheuristic;
  }

  /// Disable relay flow constraints in ILP model
  //TODO: Restore after GA?
  prob.ilp_model.flow_relay_constraints = 0;

  // We dont need postsolve
  prob.lp_params.do_postsolve = false;


  // Limit the amount of time to optimize during evaluation
  prob.lp_params.get_suboptimal = true;
  prob.lp_params.weak_timelim = prob.ga_params.max_evaluation_time; //
  prob.lp_params.strong_timelim = prob.ga_params.max_evaluation_time*2; // just in case

  if( prob.ga_params.which_solver == "cplexsolver")
  {
    g_cplexSolver = new MyCplexSolver(ga_network_static, 
				      ga_network_with_relays);

    if( prob.ga_params.evaluate_use_initial_k0 )
    {
      RNPSolutionPtr init_sol;
      cplex_solver_k0 = new MyCplexSolver(ga_network_static, ga_network_static);
      bool solved_k0  = cplex_solver_k0->solve();
      if( solved_k0 )
      {
	init_sol = cplex_solver_k0->solution();
      }
      else
      {
	printf("Can't find initial solution - ending\n");
	exit(-1);
      }
      g_cplexSolver->useInitial(init_sol);
      //delete cplex_solver_k0;
    }
    /// saves time
    g_cplexSolver->keepVarsInfo(true);
  }else
  {

  }



  string ga_output_file = prob.output_path + prob.ga_params.output_file;
  o_objective_timeline.open(ga_output_file.c_str());

  if(prob.ga_params.log_population)
  {

    string ga_population_file = prob.output_path + 
      prob.ga_params.log_population_file;
    o_population_log.open(ga_population_file.c_str());
  }

  if(prob.ga_params.log_generation)
  {
    string ga_generation_file = prob.output_path +
      prob.ga_params.log_generation_file;
    o_generation_log.open(ga_generation_file.c_str());
  }

  /// Genome random initialization (shall we do it here?)


  randomGenome(genome, prob.ga_params.genome_size);
  GAPopulation myPopulation(genome);

  //myPopulation.initializer(&MyDefaultInitializer);

  if( prob.ga_params.algorithm == "GASimple")
  {
    GASimpleGA::registerDefaultParameters(params);
    ga = new GASimpleGA(myPopulation);

  }
  else if(prob.ga_params.algorithm == "GASteadyState")
  {
    GASteadyStateGA::registerDefaultParameters(params);
    ga = new GASteadyStateGA(myPopulation);
  }


  if( prob.ga_params.ga_scaling_scheme == "GANoScaling")
  {
    ga_scaling_scheme = new GANoScaling();
  }
  else if(prob.ga_params.ga_scaling_scheme == "GALinearScaling")
  {
    ga_scaling_scheme = new GALinearScaling();
  }

  if(prob.ga_params.ga_selection_scheme == "GARankSelector")
  {
    ga_selection_scheme = new GARankSelector();
  }
  else if(prob.ga_params.ga_selection_scheme == "GATournamentSelector")
  {
    ga_selection_scheme = new GATournamentSelector();
  } 
  else if(prob.ga_params.ga_selection_scheme == "GARouletteWheelSelector")
  {
    ga_selection_scheme = new GARouletteWheelSelector();
  } 
  else
  {
    exit(1);
  }

  /// Setting parameters
  params.set(gaNpopulationSize, prob.ga_params.populationSize);    // population size
  params.set(gaNpCrossover, prob.ga_params.pCrossover);       // probability of crossover 
  params.set(gaNpMutation, prob.ga_params.pMutation);       // probability of mutation
  params.set(gaNnGenerations, prob.ga_params.nGenerations);    // number of generations
  params.set(gaNminimaxi,-1);           // minimize obj val

  /// always print the parameters
  if( true )
  {
    (*g_olog) << "Solving GA-LIB with:\n";
    (*g_olog) << "genome_size " << prob.ga_params.genome_size << endl;
    (*g_olog) << "time_limit " << prob.ga_params.time_limit << endl;
    (*g_olog) << "Algorithm " << prob.ga_params.algorithm << endl;
    (*g_olog) << "Selection Scheme " << prob.ga_params.ga_selection_scheme << endl;
    (*g_olog) << "Scaling Scheme " << prob.ga_params.ga_scaling_scheme << endl;
    (*g_olog) << "gaNpopulationSize " << prob.ga_params.populationSize << endl;
    (*g_olog) << "gaNnGenerations " << prob.ga_params.nGenerations << endl;
    (*g_olog) << "gaNpCrossover " << prob.ga_params.pCrossover << endl;
    (*g_olog) << "Reporting interval " <<  prob.ga_params.report_score_interval << endl;
    (*g_olog) << "Candidate set size " <<  prob.ga_params.candidate_set_size << endl;
    (*g_olog) << "pUniformCross " <<  prob.ga_params.pUniformCross << endl;
    (*g_olog) << "pUniformMutate " <<  prob.ga_params.pUniformMutate << endl;
    (*g_olog) << "pMutSizeChange " <<  prob.ga_params.pMutSizeChange << endl;
    (*g_olog) << "pMutSizeIncrease " <<  prob.ga_params.pMutSizeIncrease << endl;
    (*g_olog) << "adjacent_txrange " <<  prob.ga_params.adjacent_txrange << endl;
    (*g_olog) << "candidate_set_size " <<  prob.ga_params.candidate_set_size << endl;
    (*g_olog) << "pCandidateRelay " <<  prob.ga_params.pCandidateRelay << endl;
    (*g_olog) << "pUseChain " <<  prob.ga_params.pUseChain << endl;
    (*g_olog) << "candidate_neighbor_check " <<  prob.ga_params.candidate_neighbor_check << endl;
    (*g_olog) << "pPartnerRelay " <<  prob.ga_params.pPartnerRelay << endl;
    (*g_olog) << "pAdoption " <<  prob.ga_params.pAdoption << endl;
    (*g_olog) << "pUseExpat " <<  prob.ga_params.pUseExpat << endl;
    (*g_olog) << "max_evaluation_time " <<  prob.ga_params.max_evaluation_time << endl;
    (*g_olog) << "use_conflicts " <<  prob.ga_params.use_conflicts << endl;
    (*g_olog) << "use_candidate_relays " <<  prob.ga_params.use_candidate_relays << endl;
    (*g_olog) << "use_partners " <<  prob.ga_params.use_partners << endl;
    (*g_olog) << "use_chains " <<  prob.ga_params.use_chains << endl;

    (*g_olog) << "do_update_conflitmap " <<  prob.ga_params.do_update_conflictmap << endl;
    (*g_olog) << "do_update_candidate " <<  prob.ga_params.do_update_candidate << endl;
    (*g_olog) << "do_update_partner " <<  prob.ga_params.do_update_partner << endl;
    (*g_olog) << "do_update_chain " <<  prob.ga_params.do_update_chain << endl;


    (*g_olog) << "evaluate_use_initial_k0 " <<  prob.ga_params.evaluate_use_initial_k0 << endl;

  }
  g_olog->flush();

  printf("Stating GA\n");

  ga->parameters(params);
  ga->scoreFrequency(1);
  ga->scaling(*ga_scaling_scheme);
  ga->selector(*ga_selection_scheme);
  ga->terminator(&MyGATerminator);




  ga->evolve();


  if(prob.ga_params.debug)
  {
    (*g_odebug) << "Evolve ended" << endl << std::flush;
  }

  MyGenome& best = (MyGenome&)(ga->statistics().bestIndividual());
  if(prob.ga_params.verbose)
  {
    double solvingTime = g_cpuTime.cpu_time_elapsed();
    double realSolvingTime = solvingTime - g_cpuTimeToDeduct;

    (*g_olog) << ga->statistics() << endl;
    (*g_olog) << "best guy: " << best << ":\n";
    (*g_olog) << "score: " << best.score() << endl;
    (*g_olog) << "undefined childs: " << MyGenome::repairs << endl;
    (*g_olog) << "solving time " << realSolvingTime << endl;
    (*g_olog) << "total time " << solvingTime << endl;
  }

  // Now, get LpSolution for best guy and generate output

  double ret_val = best.score();
  if(prob.ga_params.do_local_search)
  {
    if(prob.ga_params.verbose)
    {
      (*g_olog) << "Doing local search" << endl;
    }

    //Graph *G = eval_net->createGraph();
    //output_graphs(G,best.toString("_"));
    //delete G;

    prob.lp_params.get_suboptimal = true;
    prob.lp_params.weak_timelim = prob.ga_params.local_search_time;
    prob.lp_params.do_postsolve = false;
    prob.lp_params.noMIP = 0; // solve MIP
    prob.lp_params.strong_timelim = prob.ga_params.local_search_time*2; 
    prob.lp_params.verbose = 1;

    /// Restore relay flow constraints in ILP model
    //TODO: Do we really need to do this?
    prob.ilp_model.flow_relay_constraints = 1;

    RNPSolutionPtr sol;
    bool solved;

    Network *eval_net = new Network(*ga_network_static);
    int k = best[0];
    for(int i= k; i>0; i--)
    {
      Node &nr = ga_network_with_relays->getRelay(best[i]);
      eval_net->addNode(nr.x, nr.y,RELAY,nr.id);
    }
    // Add also adyacent relays
    // TODO: Can be a more elaborated process?

    for(int i=k; i>0; i--)
    {
      int rid = best[i];
      std::vector<int> &adyacent = adyacent_relays[rid];
      if(!adyacent.size())
	continue;
      for(int j=0; j<adyacent.size(); j++)
      {
	Node &nr = ga_network_with_relays->getRelay(adyacent[j]);
	eval_net->addNode(nr.x,nr.y,RELAY,nr.id);
      }
    }

    printf("GA Done - doing local search\n");

    // Note:  Number of relays in Genome might be 
    // less than K - given that there might be repetitions
    // therefore, define new value of K for model data file

    if( prob.ga_params.which_solver == "cplexsolver")
    {
      g_cplexSolver->disableAllRelays();
      g_cplexSolver->setFlowRelayCons();
      for(int k=0; k< eval_net->numRelays(); k++)
      {
	Node &nr = eval_net->getRelay(k);
	g_cplexSolver->enableRelay(nr.id);
      }
      g_cplexSolver->setRelayBnds(0, prob.K);

      solved = g_cplexSolver->solve();
      m_hasSolution = solved;
      if( solved)
      {
	sol = g_cplexSolver->solution();
	m_sol = sol;
      }

    }
    else
    {
#if 0
      // Now take best guy and do optimal local search with time-limit
      // Set ilp solver parameters (reset after?)


      int genome_k = eval_net->numRelays();


      /* Replaced by MyCplexSolver*/
      MyLpSolver *lp = new MyLpSolver(eval_net,genome_k);
      Metric *metric;
      if( prob.metric == "PRR")
	metric = new EstPRR_Metric(ga_network_with_relays);
      else if(prob.metric == "DISTANCE_SQ")
	metric = new DistSQ_Metric(ga_network_with_relays);
      lp->setMetric(metric);

      //MyCplexSolver *lp = new MyCplexSolver(eval_net, genome_k);

      solved = lp->solve(); // solve lp

      sol = (lp->solution());
      delete lp;
#endif
    }
    delete eval_net;


    if(solved)
    {
      ret_val = sol->getObjVal();
      //if(prob.ga_params.scale_obj)
      //{
      //	ret_val = ret_val / prob.ga_params.obj_scale_factor;
      //}

      //G = lp->create_graph_from_sol();
      //output_graphs(G,best.toString()+"-solution");
      //delete G;

    }
    else
    {
      printf("This should not happen! :( \n");
      ret_val = GA_UNDEFINED_VAL;
      exit(1);
    }
    double solvingTime = g_cpuTime.cpu_time_elapsed();
    double realSolvingTime = solvingTime - g_cpuTimeToDeduct;

    o_objective_timeline << realSolvingTime << " " << ret_val
      << " " << g_nEvaluations << " " << g_nUnfeasible 
      << " " << g_cpuTimeSolvingMILP 
      << " " << g_nConflictsInCrossover 
      << " " << g_nUsedCandidates
      << " " << g_nUsedPartners
      << " " << g_nDetectedPartners
      << " " << g_nDetectedChains
      << " " << g_nUsedChains
      << endl;
    o_objective_timeline.flush();


    // Clean mess

  } else
  {
    m_hasSolution = true;
    m_sol = g_bestGARNPSol;
  }



  o_objective_timeline.close();


  if( prob.ga_params.do_update_partner )
  {
    (*g_olog) << "Partner relays: " << endl;
    FOREACH(it, partner_relays) 
    {
      (*g_olog) << "[" << it->first << "]";
      FOREACH(jt, it->second)
      {
	(*g_olog) << " " << *jt;
      }
      (*g_olog) << endl;
    }
    (*g_olog) << " ====================== " << endl;
  }

  if( prob.ga_params.do_update_chain )
  {
    (*g_olog) << "Relay Chains: " << endl;
    FOREACH(it, g_relayChains) 
    {
      (*g_olog) << "[" << it->first << "]";
      FOREACH(jt, it->second)
      {
	(*g_olog) << " <";
	FOREACH(kt, *jt)
	  (*g_olog) << " " << *kt;
	(*g_olog) << " >";
      }
      (*g_olog) << endl;
    }
    (*g_olog) << " ====================== " << endl;
  }

  printf("Ready to clean all\n");

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "Ready to clean" << endl << std::flush;
  }

  /// Clean the mess
  ///NOTE: delete cplex causes SEGFAULT in ARES, WHY?>????????????
  //TODO Find out why
  //delete g_cplexSolver;

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "deleted CPX" << endl << std::flush;
  }
  delete ga_scaling_scheme;
  delete ga_selection_scheme;
  delete ga;

  /// right now, let it be
  /// it was causing segfault in ares when deleting after 
  /// solving
  //if( prob.ga_params.evaluate_use_initial_k0 )
  //{
    //if( cplex_solver_k0 != NULL )
      //delete cplex_solver_k0;
  //}

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "All clean, good bye" << endl << std::flush;
  }

  if( prob.ga_params.log_to_file)	
    flog.close();


  printf("GA: Solve terminated\n");
  return ret_val;

}
double GASol::test()
{

  // Disable logging
  prob.ga_params.log_to_file = 0;
  prob.ga_params.log_population = 0;
  prob.ga_params.verbose = 0;

  /// Disable relay flow constraints in ILP model
  //TODO: Restore after GA?
  prob.ilp_model.flow_relay_constraints = 0;

  MyGenome genome(prob.ga_params.genome_size);      // test default constructor (if we have one)
  cout << "genome after creation:\n" << genome << endl;

  for(int i=prob.ga_params.genome_size;i>0;i--)
  {
    int ri = GARandomInt(0,ga_network_with_relays->numRelays()-1);
    genome[i] = ri;
  }
  genome[0] = prob.ga_params.genome_size;

  genome.initialize();
  cout << "genome after manual initialization:\n" << genome << endl;
  genome.repair();

  cout << "genome after repair:\n" << genome << endl;
  float score = genome.score();

  cout << "Score is " << score << endl;


  genome.mutate(0.5);      // test the mutator
  cout << "genome after mutation:\n" << genome << endl;

  score = genome.evaluate();

  cout << "Score after mutation is " << score << endl;

  MyGenome* a = new MyGenome(genome);   // test copy constructor
  MyGenome* b = new MyGenome(genome);
  score = a->score();
  cout << "score after copy a" << score << endl;

  score = b->evaluate();

  cout << "score after copy b" << score << endl;

  b->initialize();


  cout << "score after initialize b" << score << endl;


  for(int i=(*b)[0];i>0;i--)
  {

    int ri = GARandomInt(0,ga_network_with_relays->numRelays()-1);
    (*b)[i] = ri;
  }
  (*b)[0] = prob.ga_params.genome_size / 2;
  score = b->score();
  cout << "Score after changing b " << score << endl;
  b->evaluate(gaTrue);
  score = b->score();
  cout << "Score after re-evaluating b " << score << endl;



  MyGenome* c = (MyGenome *)genome.clone(GAGenome::CONTENTS);
  cout << "clone contents:\n" << *c << "\n";
  MyGenome* d = (MyGenome *)genome.clone(GAGenome::ATTRIBUTES);
  cout << "clone attributes:\n" << *d << "\n";

  a->initialize();
  b->initialize();
  cout << "parents:\n" << *a << "\n" << *b << "\n";

  MyGenome::Cross(*a, *b, c, d);   // test two child crossover
  cout << "children of crossover:\n" << *c << "\n" << *d << "\n";
  double comp = c->compare(*d);       // test the comparator
  cout << "Comparing:\n" << *c << "\n" << *d << "\n" << comp << "\n";


  MyGenome::Cross(*a, *b, c, 0);   // test single child crossover
  cout << "child of crossover:\n" << *c << "\n";

  double eval = MyGenome::Evaluate(*c);
  cout << "Evaluating\n" << *c << "\n" << eval << "\n";


  delete a;
  delete b;
  delete c;
  delete d;

  return 0;
}


int findGeneMatch(int r, GA1DArrayGenome<int>  &mg, set<int> &exclusionSet)
{
  double max_dist = 0;
  int max_dist_r = r;
  int k = mg[0];

  Node &r1 = ga_network_with_relays->getRelay(r);
  for(int i=k;i>0;i--)
  {
    if(exclusionSet.find(mg[i]) == exclusionSet.end())
    {

      Node &r2 = ga_network_with_relays->getRelay(mg[i]);
      double d = ga_network_with_relays->distanceLink(r1,r2);
      if( d > max_dist)
      {
	max_dist = d;
	max_dist_r = mg[i];
      }
    }
  }
  return max_dist_r;

}
void 
MyGenome::Init(GAGenome &a)
{
  // your initializer here
  MyGenome &mg = (MyGenome&)a;

  int randomSize = GARandomInt(1,prob.ga_params.genome_size);
  randomGenome(mg, randomSize);
  mg.repair();
  //mg.resizeBehaviour(0,mg.length());
}
int 
MyGenome::equal(const GAGenome& a) const
{
  MyGenome &mg = (MyGenome&)a;
  if( mg[0] != (*this)[0])
    return 0;
  for(int i=mg[0]; i>0; i--)
    if(mg[i] != (*this)[0])
      return 0;
  return 1;
}


bool 
MutateSize( GAGenome& c, set<int> &mygenomeset, float pmut) 
{
  /// Number of attempts to find a valid mutation
  const int mutation_trials = 5;
  bool size_changed = false;
  MyGenome &child=(MyGenome &)c;
  if(GAFlipCoin(prob.ga_params.pMutSizeIncrease))
  {
    // Is there room for one more relay?
    if(child[0] < prob.ga_params.genome_size)
    {
      int mylen = child[0];
      /// Limited trials to get a relay which is not in the genome
      for(int i=0; i< mutation_trials; i++)
      {
	int ri = getRandomRelay();
	if( mygenomeset.find(ri) == mygenomeset.end())
	{
	  mylen++;
	  size_changed = true;
	  child[mylen] = ri;
	  mygenomeset.insert(ri);
	  break;
	}
      }
      /// Update size
      child[0] = mygenomeset.size();
    }
  }
  else
  {
    /// Can we remove a relay?
    /// Choose a random pos and swap to end
    if(child[0] > 1)
    {
      int ri = GARandomInt(1,mygenomeset.size());
      int rid = child[ri];
      int last_relay = child[mygenomeset.size()];
      /// swap last and dead relay
      child[ri] = last_relay;

      /// delete from set
      mygenomeset.erase(rid);
      size_changed = true;

      /// change size
      child[0] = mygenomeset.size();

    }
  }
  return size_changed;
}

int 
MyUniformMutation(GAGenome& c, std::set<int>& mygenomeset, float pmut) 
{
  /// Number of attempts to find a valid mutation
  const int mutation_trials = 5;
  MyGenome &child=(MyGenome &)c;
  /// First, decide if genome changes its size
  bool size_changed = false;
  if(GAFlipCoin(prob.ga_params.pMutSizeChange))
  {
    size_changed = MutateSize(c,mygenomeset,pmut);
  }

  if( child[0] == 0)
  {
    /// nothing to do
    return 0;
  }
  /// Then do normal UNIFORM mutation
  int nMut=0;
  for(int i=child[0]; i>=1; i--){
    if(GAFlipCoin(pmut)){
      int rid = child[i];
      /// Limited trials to find mutation relay
      for(int j=0; j<mutation_trials;j++)
      {
	int n_rid = getRandomRelay();
	if( mygenomeset.find(n_rid) == mygenomeset.end())
	{
	  child[i]=n_rid;
	  mygenomeset.erase(rid);
	  mygenomeset.insert(n_rid);
	  nMut++;
	  break;
	}
      }
    }
  }
  if(nMut || size_changed)
  {
    //c._evaluated=false;
    child.evaluate(gaTrue);
    //child.repair();
  }
  //printf("Mutated\n");
  if(!nMut && size_changed)
    nMut=1;
  return nMut;
}

int 
EnhancedMutation(GAGenome& c, std::set<int> &mygenomeset, float pmut) 
{
  /// Number of attempts to find a valid mutation
  const int mutation_trials = 5;

  MyGenome &child=(MyGenome &)c;
  /// First, decide if genome changes its size
  bool size_changed = false;
  if(GAFlipCoin(prob.ga_params.pMutSizeChange))
  {
    size_changed = MutateSize(c,mygenomeset,pmut);
  }
  if( child[0] == 0)
  {
    return 0;
  }
  /// Then do normal ADYACENT mutation ...
  int nMut=0;
  for(int i=child[0]; i>=1; i--){
    if(GAFlipCoin(pmut)) {
      int rid = child[i];

      if(GAFlipCoin(prob.ga_params.pEnhancedMutationRandomRelay))
      {
	/// Change relay by random relay

	/// Limited trials to find mutation relay
	for(int j=0; j<mutation_trials;j++)
	{
	  int n_rid = getRandomRelay();
	  if( mygenomeset.find(n_rid) == mygenomeset.end())
	  {
	    child[i]=n_rid;
	    mygenomeset.erase(rid);
	    mygenomeset.insert(n_rid);
	    nMut++;
	    break;
	  }
	}
      }
      else 
      {

	std::vector<int> &adyacent = adyacent_relays[rid];
	if(!adyacent.size())
	  continue;

	/// Limited trials to find mutation relay
	for(int j=0; j<mutation_trials;j++)
	{
	  int ri = GARandomInt(0,(int)(adyacent.size()-1));
	  int n_rid = adyacent[ri];
	  if( mygenomeset.find(n_rid) == mygenomeset.end())
	  {
	    child[i]=n_rid;
	    mygenomeset.erase(rid);
	    mygenomeset.insert(n_rid);
	    nMut++;
	    break;
	  }
	}
      }
    }
  }
  if(nMut || size_changed)
  {
    //c._evaluated=false;
    child.evaluate(gaTrue);
    //child.repair();
  }
  //printf("Mutated\n");
  if(!nMut && size_changed)
    nMut=1;
  return nMut;
}

int 
MyGenome::Mutate(GAGenome& c, float pmut) 
{
  MyGenome &child=(MyGenome &)c;
  if(pmut <= 0.0) return 0;

  std::set<int> mygenomeset;
  /// Copy Genome contents to SET O(Log(n)) to determine membership!
  for(int i=child[0]; i>0; i--)
    mygenomeset.insert(child[i]);
  /// Lets keep consistency, sort genome
  child[0] = mygenomeset.size();
  int k=1;
  FOREACH(it, mygenomeset)
  {
    child[k++] = *it;
  }


  if(GAFlipCoin(prob.ga_params.pUniformMutate))
  {
    //printf("Doing Uniform Mutation\n");
    return MyUniformMutation(c, mygenomeset, pmut);
  } 
  else 
  {
    return EnhancedMutation(c, mygenomeset, pmut);
  }

}

void MyGenome::repair()
{
  /// Sort and update size
  std::set<int> mygenomeset;

  /// Copy Genome contents to SET O(Log(n)) to determine membership!
  for(int i=(*this)[0]; i>0; i--)
    mygenomeset.insert((*this)[i]);

  (*this)[0] = mygenomeset.size();
  int k=1;
  FOREACH(it, mygenomeset)
  {
    (*this)[k++] = *it;
  }

  // DO NOTHING 
  /*
    float my_val = score();
    if(my_val == GA_UNDEFINED_VAL)
    {
    MyGenome::repairs+=1;

    }
    */


  /*
    sort(GAArray<int>::a, GAArray<int>::a + length());

    int *ix = unique(GAArray<int>::a, GAArray<int>::a + length());
    if(ix != (GAArray<int>::a+length()))
    {
    cout << "trying to resize " << (*this) << endl;
    printf("%d\n",ix-GAArray<int>::a);
  //resize(ix-GAArray<int>::a);
  }
  */


}

float 
MyGenome::Compare(const GAGenome& a, const GAGenome& b )
{
  // your comparison here
  MyGenome &sis=(MyGenome &)a;
  MyGenome &bro=(MyGenome &)b;
  if(sis[0] != bro[0]) return -1;
  float count = 0.0;
  for(int i=sis[0]; i>0; i--)
    count += ((sis[i] == bro[i]) ? 0 : 1);
  return count/sis[0];

}

string MyGenome::toString(string separator) {
  stringstream ss;
  int k = (*this)[0];
  for(int i=k;i>1; i--)
  {
    int rid = (*this)[i];
    Node &nr = ga_network_with_relays->getRelay(rid);
    ss << nr.id << separator;
  }

  int rid = (*this)[1];
  Node &nr = ga_network_with_relays->getRelay(rid);
  ss  << nr.id;
  return ss.str();
}
double 
fixRNPSolution(MyGenome &mg, RNPSolutionPtr sol)
{

  int k = mg[0];
  if(prob.ga_params.debug)
  {
    (*g_odebug) << "Fixing RNPSol " << endl;
  }

  LpSolutionPtr lpsol = sol->getLpSolution();
  double objval = lpsol->objval;
  int not_used = 0;
  for(int i= k; i>0; i--)
  {
    // Do the candidateRelaySet
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);

    if( fabs(oflow) <= EPSILON )
    {

      if(prob.ga_params.debug)
      {
	(*g_odebug) << "FIXRNP Relay " << nr.id << " not used " << endl;
      }
      printf("fixRNP found relay not used... fixing\n");
      stringstream ss; 
      ss << VAR_RELAY << "[" << nr.id << "]";
      string varname = ss.str();
      printf("setting %s to zero\n",varname.c_str());

      if(prob.ga_params.debug)
      {
	(*g_odebug) << "FIXRNP changing " << varname << " to 0 " << endl;
      }
      ITERATOR(lpsol->value) it = (lpsol->value).find(varname);
      if( it == (lpsol->value).end())
      {
	/// this is not more an error because we are not fixing the relays
	/// if genome_size > k, then, we assume that it is not used (var=0)
	//fprintf(stderr, "Fatal error fixing RNPSol. Variable %s not found\n",
		//varname.c_str());
	//exit(1);
      }
      else
      {
	it->second = 0;	
	//lpsol->value[varname] = 0;
	not_used += 1;
      }
    }
  }
  if(not_used)
  {

    if(prob.ga_params.debug)
    {
      (*g_odebug) << "FIXRNP fixing objval from  " 
	<< objval << " to ";
    }
    printf("Fixing objval from %.2f to ", objval);
    objval = objval - not_used*prob.ilp_model.penalty_relay;
    lpsol->objval = objval;
    printf("%.2f\n", objval);
    if(prob.ga_params.debug)
    {
      (*g_odebug) << objval << endl;
    }
  } else
  {
    printf("Nothing to fix\n");
  }
  return objval;
#if 0
  printf("Attempting to fix RNPSol\n");
  const vector<string> &relays = sol->getRelaysToUse();


  if(prob.ga_params.debug)
  {
    (*g_odebug) << "RNPSol has " << relays.size() << " relays" << endl;
  }
  LpSolutionPtr lpsol = sol->getLpSolution();
  double objval = lpsol->objval;
  int not_used = 0;
  FOREACH(rid, relays)
  {
    double oflow = sol->getOutgoingFlow(*rid);

    if(prob.ga_params.debug)
    {
      (*g_odebug) << "FIXRNP Relay " << *rid << " has oflow " << oflow << endl;
    }
    if( fabs(oflow) <= EPSILON )
    {

      if(prob.ga_params.debug)
      {
	(*g_odebug) << "FIXRNP Relay " << *rid << " not used " << endl;
      }
      printf("fixRNP found relay not used... fixing\n");
      stringstream ss; 
      ss << VAR_RELAY << "[" << *rid << "]";
      string varname = ss.str();
      printf("setting %s to zero\n",varname.c_str());

      if(prob.ga_params.debug)
      {
	(*g_odebug) << "FIXRNP changing " << varname << " to 0 " << endl;
      }
      ITERATOR(lpsol->value) it = (lpsol->value).find(varname);
      if( it == (lpsol->value).end())
      {
	fprintf(stderr, "Fatal error fixing RNPSol. Variable %s not found\n",
		varname.c_str());
	exit(1);
      }
      it->second = 0;	
      //lpsol->value[varname] = 0;
      not_used += 1;
    }
  }
  if(not_used)
  {
    printf("Fixing objval from %.2f to ", objval);
    objval = objval - not_used*prob.ilp_model.penalty_relay;
    printf("%.2f\n", objval);
  } else
  {
    printf("Nothing to fix\n");
  }
#endif
}



bool 
updateConflictMap(MyGenome &mg, RNPSolutionPtr sol)
{

  int k = mg[0];
  set<int> used, notused;
  for(int i= k; i>0; i--)
  {
    // Do the candidateRelaySet
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);
    if( fabs(oflow) < EPSILON)
    {
      notused.insert(mrelay);
    } else {
      used.insert(mrelay);
    }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "ConflictMap used " << used.size() 
      << " notused  " << notused.size() << endl;

  }
  /// Add conflicts
  FOREACH(ru, used)
  {
    set<int> &csetu = g_conflictMap[*ru];
    FOREACH(rn, notused)
    {
      csetu.insert(*rn);
      set<int> &csetn = g_conflictMap[*rn];
      csetn.insert(*ru);
    }
  }

  /// Here we found an ugly bug

  /// Remove conflicts

  //FOREACH(ru1, used)
  //{
    //set<int> &cset1 = g_conflictMap[*ru1];
    //for(ITERATOR(used) ru2=ru1; ru2 != used.end(); ru2++)
    //{
      //cset1.erase(*ru1);
      //set<int> &cset2 = g_conflictMap[*ru2];
      //cset2.erase(*ru2);
    //}
  //}

  FOREACH(ru1, used)
  {
    set<int> &cset1 = g_conflictMap[*ru1];
    for(ITERATOR(used) ru2=ru1; ru2 != used.end(); ru2++)
    {
      cset1.erase(*ru2);
      set<int> &cset2 = g_conflictMap[*ru2];
      cset2.erase(*ru1);
    }
  }

  if(prob.ga_params.debug)
  {
    for(int i= k; i>0; i--)
    {
      int mrelay = mg[i];

      set<int> &cset = g_conflictMap[mrelay];
      (*g_odebug) << "CM[ " << mrelay << "] "
	<< cset.size() << endl;
    }
  }
  return (notused.size() > 0);
}


/// In this version, we partition the area into regions
/// each relay is assigned to one region and we have a
/// candidate set for each region
bool updateCandidateRelaySetV2(MyGenome &mg, RNPSolutionPtr sol)
{
  int k = mg[0];
  int candidate_set_size = prob.ga_params.candidate_set_size;
  int candidate_regions = g_candidateRegions.size();
  bool update_candidates = false;

  bool update_lists = false;
  for(int i= k; i>0; i--)
  {
    // Do the candidateRelaySet
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);
    if(prob.ga_params.debug)
    {
      (*g_odebug) << "Relay " << nr.id << " has out_flow "
	<< oflow << endl;

    }
    if( update_lists)
    {
      if( fabs(oflow) < EPSILON)
      {
	ITERATOR(g_greenList) it = g_greenList.find(mrelay);
	/// If relay not in green list
	/// i.e never seen with positive flow
	if( it == g_greenList.end())
	{
	  /// flag it as red (useless)
	  g_redList.insert(mrelay);
	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "Relay " << nr.id 
	      << " in red_list "<< endl;

	  }
	}
      } 
      else 
      {
	/// insert in green list
	g_greenList.insert(mrelay);
	/// if it was in redlist
	ITERATOR(g_redList) it = g_redList.find(mrelay);
	if( it != g_redList.end())
	{
	  /// remove it
	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "Relay " << nr.id 
	      << " removed from red_list "<< endl;

	  }
	  g_redList.erase(it);
	}
      }
    }

    double moflow = max_flow_for_relay[mrelay];
    if( oflow > moflow)
    {
      if(prob.ga_params.debug)
      {
	(*g_odebug) << "Relay " << nr.id << " has new max_flow "
	  << oflow << " (" << moflow << ")" << endl;
      }
      max_flow_for_relay[mrelay]=oflow;
    } else 
    {
      continue;
    }
    vector<int> &regions = g_relayToRegions[mrelay];

    if(prob.ga_params.debug)
    {
      (*g_odebug) << "Relay " << nr.id << " is associated to " 
	<< regions.size() << " regions: " <<  endl;
    }
    FOREACH(cr, regions)
    {
      CandidateRelayRegion &cR = g_candidateRegions[*cr];
      bool entered=false;

      if( cR.candidates.size() < cR.max_candidates)
      {
	/// Candidate enters directly
	cR.candidates.insert(mrelay);
	entered = true;

	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "Relay " << nr.id << " enters directly cs " 
	    << *cr << " size: "
	    << cR.candidates.size() << endl;
	}
      } 
      else
      {
	if( oflow >= cR.min_relay_flow )
	{
	  /// Relay should enter
	  /// We must make space for new candidate
	  /// so lets remove minimum guy
	  cR.candidates.erase(cR.min_candidate);
	  cR.candidates.insert(mrelay);
	  entered = true;

	  if(prob.ga_params.debug)
	  {

	    Node &nn = ga_network_with_relays->getRelay(cR.min_candidate);
	    (*g_odebug) << "Relay " << nr.id << " replaces minimum "
	      << nn.id << " in cs " << *cr << endl;
	  }
	}
      }
      if( entered )
	update_candidates = true;
      /// Only update if we are full
      if( entered && cR.candidates.size() == cR.max_candidates)
      {
	/// Must update minimum
	/// This could be more efficient,
	//TODO: Improve
	double minimum_flow = GA_UNDEFINED_VAL;
	int minimum_guy;
	FOREACH(crelay, cR.candidates)
	{
	  double cflow = max_flow_for_relay[*crelay];
	  if( cflow < minimum_flow)
	  {
	    minimum_flow = cflow;
	    minimum_guy = (*crelay);
	  }
	}
	cR.min_candidate = minimum_guy;
	cR.min_relay_flow =  minimum_flow;

	if(prob.ga_params.debug)
	{
	  Node &nn = ga_network_with_relays->getRelay(cR.min_candidate);
	  (*g_odebug) << "Updated region candidates. new minimum " << nn.id << " with flow "
	    << cR.min_relay_flow << endl;
	}
      }
    }
  }
  if( update_candidates)
  {
    candidateRelays.clear();
    int crx = 0;
    FOREACH(cr, g_candidateRegions)
    {
      CandidateRelayRegion &cR = (*cr);
      candidateRelays.insert( cR.candidates.begin(), cR.candidates.end());

      if(prob.ga_params.debug)
      {
	(*g_odebug) << "CS " << crx << " (" << cR.candidates.size() << ")";
	FOREACH(it, cR.candidates)
	{
	  double flow = max_flow_for_relay[*it];
	  (*g_odebug) << "[ " << flow
	    << ", " 
	    << GAGetRelayId(*it) 
	    << " ]";
	}
	(*g_odebug) << endl;
      }

      crx++;
    }
  }

  if(prob.ga_params.debug)
  {

    (*g_odebug) << "RL " << " (" << g_redList.size() << ")";
    FOREACH(it, g_redList)
    {
      (*g_odebug) << " " << GAGetRelayId(*it) ;
    }
    (*g_odebug) << endl;

    (*g_odebug) << "GL " << " (" << g_greenList.size() << ")";
    FOREACH(it, g_greenList)
    {
      (*g_odebug) << " " << GAGetRelayId(*it) ;
    }
    (*g_odebug) << endl;
  }

  return update_candidates;
}


/// Different from partners -> extract the 'chains' (e.g., 2 or more interconnected relays)
bool 
updateChains(MyGenome &mg, RNPSolutionPtr sol)
{
  bool update_chain = false;
  int csize = 0; /// size of largest chain
  int k = mg[0];
  /// Check for partner relays

  std::set< std::pair<int, int> > edges;
  for(int i= k; i>0; i--)
  {
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);

    if( fabs(oflow) <= EPSILON )
      continue;

    for(int j=i-1;j>0;j--)
    {
      Node &nr2 = ga_network_with_relays->getRelay(mg[j]);
      if( sol->isEdge(nr.id, nr2.id) || sol->isEdge(nr2.id, nr.id))
      {
	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "connected relays " <<  nr.id
	    << " " << nr2.id << endl;
	}
	int r_a = min(mrelay, mg[j]);
	int r_b = max(mrelay, mg[j]);
	edges.insert( make_pair( r_a, r_b) );
      }
    }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "detected  " <<  edges.size()
      << " edges involving relays" << endl;
  }
  if( !edges.size() )
    return false;

  /// connected components
  using namespace boost;
  {
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    /// first, create a (bi)map for indeces
    std::map< int, int> rix, r_rix;
    int ri=0;
    FOREACH(it, edges)
    {
      if( rix.find(it->first) == rix.end())
      {
	rix[it->first]=ri;
	r_rix[ri]=it->first;
	ri++;
      }
      if( rix.find(it->second) == rix.end())
      {
	rix[it->second]=ri;
	r_rix[ri]=it->second;
	ri++;
      }
    }
    /// then, create graph

    Graph G;
    FOREACH(it, edges)
    {
      add_edge(rix[it->first], rix[it->second], G);
    }

    std::vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);
    g_nDetectedChains+=num;

    std::vector<int>::size_type i;

    if(prob.ga_params.debug)
      (*g_odebug) << "Total number of components: " << num << endl;

    std::vector< std::set<int> > C(num);
    csize=0;
    /// 2 passes
    /// 1 - find connected components C
    for (i = 0; i != component.size(); ++i)
    {
      if( component[i] >= C.size())
      {
	printf("WTF? component >\n");
	exit(-1);
      }
      C[component[i]].insert(r_rix[i]);
      if( C[component[i]].size() > csize)
	csize = C[component[i]].size();
    }

    ///2 - fill out the g_relayChains
    for (i = 0; i != component.size(); ++i)
    {
      int rid = r_rix[i];
      g_relayChains[rid].insert( C[component[i]] );
    }


    if(prob.ga_params.debug)
    {
      for (i = 0; i < C.size();i++)
      {
	(*g_odebug) << "Component [";
	FOREACH(it, C[i])
	{
	  (*g_odebug) << " " << *it;
	}
	(*g_odebug) << "]" << endl;
      }
    }
  }

  if(prob.ga_params.debug)
    (*g_odebug) << "Largest component: " << csize << endl;

  return csize>2;
}

bool 
updatePartners(MyGenome &mg, RNPSolutionPtr sol)
{

  bool update_partner = false;
  int k = mg[0];
  /// Check for partner relays

  for(int i= k; i>0; i--)
  {
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);
    for(int j=i-1;j>0;j--)
    {
      Node &nr2 = ga_network_with_relays->getRelay(mg[j]);
      if( sol->isEdge(nr.id, nr2.id) || sol->isEdge(nr2.id, nr.id))
      {
	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "partners " <<  nr.id
	    << " " << nr2.id << endl;
	}
	if( partner_relays[mrelay].count(mg[j]) )
	  g_nDetectedPartners++;
	partner_relays[mrelay].insert(mg[j]);
	partner_relays[mg[j]].insert(mrelay);
	update_partner = true;
      }
    }
  }
  return update_partner;
}
void 
updateCandidateRelaySet(MyGenome &mg, RNPSolutionPtr sol)
{
  /// Update candidate relay set

  int k = mg[0];
  int candidate_set_size = prob.ga_params.candidate_set_size;
  for(int i= k; i>0; i--)
  {
    // Do the candidateRelaySet
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    int mrelay = mg[i];
    double oflow = sol->getOutgoingFlow(nr.id);

    if(prob.ga_params.debug)
    {
      (*g_odebug) << "Relay " << nr.id << " has out_flow "
	<< oflow << endl;
    }
    double moflow = max_flow_for_relay[mrelay];
    if( oflow > moflow)
    {
      if(prob.ga_params.debug)
      {
	(*g_odebug) << "Relay " << nr.id << " has new max_flow "
	  << oflow << " (" << moflow << ")" << endl;
      }
      max_flow_for_relay[mrelay]=oflow;
      // New maximum
      // Check candidateSet
      if(oflow > min_candidate_flow || 
	 candidateRelaySet.size() < candidate_set_size)
      {

	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "Relay " << nr.id << " qualifies for candidate" << endl;
	}
	bool candidateOK = true;
	if(prob.ga_params.candidate_neighbor_check)
	{

	  /// Here we check if some adyacent guy is inside already
	  vector<int> &ady = adyacent_relays[mrelay];
	  set<int> myadjacents;
	  myadjacents.insert( ady.begin(), ady.end());
	  set<int> myintersect;
	  set_intersection(myadjacents.begin(), myadjacents.end(),
			   candidateRelays.begin(), candidateRelays.end(),
			   insert_iterator<set<int> >(myintersect,myintersect.begin()));
	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "Relay " << nr.id << " has " << myadjacents.size() << " adjacent relays and"
	      << myintersect.size() << " are candidates" << endl;
	  }

	  /// If there are too much adyacent relays 
	  /// inside candidate set, what should we do?
	  if( myintersect.size() > prob.ga_params.max_candidate_neighbor)
	  {

	    if(prob.ga_params.debug)
	    {
	      (*g_odebug) << "Relay " << nr.id << " has been discarded" << endl;
	    }
	    /// We could check those adyacents and, if the current relay
	    /// is a good one, replace an adjacent
	    /// but too much work
	    //TODO 
	    /// For the moment, discard candidate
	    candidateOK = false;
	  }
	}


	if( candidateOK)
	{
	  /// We need to check if same relay is already
	  /// in the set, so try to remove it

	  candidateRelaySet.erase(make_pair(moflow, mrelay));
	  candidateRelays.erase(mrelay);

	  if(candidateRelaySet.size() < candidate_set_size)
	  {
	    /// add new candidate
	    //printf("adding new candidate\n");
	    candidateRelaySet.insert(make_pair(oflow, mrelay));
	    candidateRelays.insert(mrelay);
	  }
	  else if(candidateRelaySet.size() == candidate_set_size)
	  {
	    /// Else, remove minimum and insert

	    candidateRelaySet.erase(make_pair(min_candidate_flow, min_candidate));
	    candidateRelays.erase(min_candidate);

	    candidateRelaySet.insert(make_pair(oflow, mrelay));
	    candidateRelays.insert(mrelay);
	  }
	  else
	  {
	    printf("Candidate set oversized\n");

	    // This should not happen
	    // Just in case, resize set
	    int surplus = 
	      candidateRelaySet.size() - candidate_set_size;
	    ITERATOR(candidateRelaySet) it = candidateRelaySet.begin();
	    while(surplus--)
	    {
	      candidateRelaySet.erase(it);
	      candidateRelays.erase(it->second);
	      it++;
	    }
	  }
	  ITERATOR(candidateRelaySet) it = candidateRelaySet.begin();
	  min_candidate_flow = (*it).first;
	  min_candidate = (*it).second;

	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "CS: ";
	    FOREACH(it, candidateRelaySet)
	    {
	      (*g_odebug) << "[ " 
		<<  (*it).first 
		<< ", " 
		<< GAGetRelayId((*it).second) 
		<< " ]";
	    }
	    (*g_odebug) << endl;
	  }

	}
      }
    }
  }
}


  RNPSolutionPtr
EvaluateUsingDummy(MyGenome& mg)
{
  LpSolutionPtr lpsol( new LpSolution());
  lpsol->solved = false;
  RNPSolutionPtr sol( new RNPSolution(lpsol));
  return sol;
}


  RNPSolutionPtr
EvaluateUsingCplexSolver(MyGenome& mg)
{
  printf("Evaluate with CPLEX solver\n");


  g_cplexSolver->disableAllRelays();
  printf("Disabled all relays\n");

  int k = mg[0];
  for(int i= k; i>0; i--)
  {
    //Node &nr = ga_network_with_relays->getRelay(mg[i]);

    /// if genome_size > prob.K, we can not fix relay, but only enable it
    if( k > prob.K )
    {
      g_cplexSolver->enableRelay(mg[i]);
    }
    else
    {
      g_cplexSolver->fixRelay(mg[i]);
    }
  }

  if( k > prob.K )
  {
    printf("Enabled relays\n");
    g_cplexSolver->setRelayBnds(prob.K,prob.K);
    printf("Relay bounds set\n");
  }
  else
  {
    printf("Fixed relays\n");
    g_cplexSolver->setRelayBnds(k,k);
    printf("Relay bounds set\n");
  }

  bool solved = g_cplexSolver->solve(); // solve lp
  RNPSolutionPtr sol = g_cplexSolver->solution();
  double solvingTime = sol->cpuTimeTotal();
  printf("Done with evaluation. MILP Solving time %.2f\n", solvingTime);
  g_cpuTimeSolvingMILP += solvingTime;
  // Before returning, call callback
  myCallback();

  return sol;


}
  string 
myGenomeToString(MyGenome &mg)
{
  stringstream ss;
  int k = mg[0];
  for(int i= k; i>0; i--)
  {
    if(mg[i] != Network::INVALID_NODE_INDEX)
    {
      if(ga_network_with_relays->isRelay(mg[i]))
      {
	Node &nr = ga_network_with_relays->getRelay(mg[i]);
	ss << " " << nr.id;
      } else
	ss << " CRAPPYNULL (" << mg[i] << ")";
    } else
      ss << " NULL";
  }
  return ss.str();
}

  void
MyGenome::fromLpSolution(LpSolutionPtr lpsol)
{
  printf("Filling MyGenome with LpSolution\n");
  RNPSolution rnpsol(lpsol);
  int k = rnpsol.numberRelays();
  printf("RNPSol has %d relays\n", k);
  //(*this)[0]=0;

  (*this)[0] = k;
  const std::vector<string> &relaysToUse = rnpsol.getRelaysToUse();
  if(relaysToUse.size() != rnpsol.numberRelays())
  {
    printf("ERROR here\n");
  }
  for(int i= k; i>0; i--)
  {
    int ix = ga_network_with_relays->getRelayIndex(relaysToUse.at(i-1));
    printf("fromLpSolution, checking relay index  %s (%d)\n", relaysToUse.at(i-1).c_str(), ix);
    if( ix == Network::INVALID_NODE_INDEX)
    {
      fprintf(stderr, "ERROR: Network invalid index\n");
    }
    (*this)[i] = ix;
  }

}






  float 
MyGenome::Evaluate(GAGenome& g)
{
  float ret_val;
  MyGenome &mg = (MyGenome&)g;
  string sg = myGenomeToString(mg);



  ITERATOR(evaluated_genomes) it = evaluated_genomes.find(sg);
  if( it != evaluated_genomes.end())
  {
    ret_val = it->second;
    return ret_val;
  }
  if(prob.ga_params.log_population)
  {
    o_population_log << sg;
  }

  RNPSolutionPtr sol;
  if( prob.ga_params.which_solver == "cplexsolver" 
      || prob.ga_params.which_solver == "default")
  {
    sol = EvaluateUsingCplexSolver(mg);
  }
  else if(prob.ga_params.which_solver == "dummy")
  {
    sol = EvaluateUsingDummy(mg);
  } 
  else
  {
#if 0
    sol = EvaluateUsingWrapper(mg);
#endif
    fprintf(stderr, "Evaluate with wrapper not supported anymore\n");
    exit(1);
  }

  if(sol->solved())
  {
    ret_val = sol->getObjVal();
    if(prob.ga_params.scale_obj)
    {
      ret_val = ret_val / prob.ga_params.obj_scale_factor;
    }


    Graph *G = sol->create_graph_from_sol(ga_network_with_relays);
    // Fix obj value cost considering non-used relays
    if( prob.ga_params.fix_rnpsol)
    {
      double old_ret_val = ret_val;
      ret_val = fixRNPSolution(mg, sol);
      if( ret_val == old_ret_val)
      {
	if(prob.ga_params.log_population)
	{
	  /// print fixed
	  o_population_log << " ( ok )";
	  //generate_simulation(G, ga_network_with_relays, sg.c_str());
	}
      }

    }

    if( g_matheuristic)
    {
      if( ret_val < g_lastPushedSolValue)
      {
	if( ret_val < g_lastPulledSolValue)
	{

	  printf("GA: Pushing solution with objval %.2f and valid variables %d\n",
		 ret_val, sol->nValidVariables());
	  if( !g_matheuristic->pushSolution(sol->getLpSolution()))
	    g_lastPushedSolValue = ret_val;
	}
      }
    }


    if(prob.ga_params.debug)
    {
      (*g_odebug) <<  "Checking " << mg.toString("_") << endl;
    }
    bool uCandidate, uPartner, uConflict, uChain;
    uCandidate=uPartner=uConflict=uChain=false;
    if( prob.ga_params.do_update_candidate)
    {
      printf("Updating candidate relay set\n");
      uCandidate = updateCandidateRelaySetV2(mg, sol);
    }

    if( prob.ga_params.do_update_partner)
    {
      printf("Updating partners relay set\n");
      uPartner = updatePartners(mg, sol);
    }

    if( prob.ga_params.do_update_chain)
    {
      printf("Updating chains\n");
      uChain = updateChains(mg, sol);
      if(prob.ga_params.debug && uChain)
      {
	(*g_odebug) <<  "GENERATING " << mg.toString("_") << endl;
	output_graphs(G,mg.toString("_"));
      }
    }

    if( prob.ga_params.do_update_conflictmap)
    {
      printf("Upating conflict map\n");
      uConflict = updateConflictMap(mg, sol);
    }

    if( uCandidate && uPartner && uConflict)
    {
      if(prob.ga_params.debug)
      {
	(*g_odebug) <<  "GENERATING " << mg.toString("_") << endl;
	output_graphs(G,mg.toString("_"));
      }
    }
  }
  else
  {
    ret_val = GA_UNDEFINED_VAL;
    g_nUnfeasible++;
  }

  if( ret_val < g_bestScore)
  {
    g_bestScore = ret_val;
    g_bestGARNPSol = sol;
  }

  //g.score(ret_val);
  g_nEvaluations++;

  /// Store score to prevent re-evaluations
  evaluated_genomes[sg] = ret_val;

  if(prob.ga_params.log_population)
  {
    o_population_log << " " << ret_val << endl;
  }
  double gaTime = g_cpuTime.cpu_time_elapsed();
  double realGATime = gaTime - g_cpuTimeToDeduct;


  if(prob.ga_params.verbose)
  {
    (*g_olog) << "CpuTimeTotal " << sol->cpuTimeTotal()
      << " CpuTimeSolving " << sol->cpuTimeSolving()
      << " CpuTimeToDeduct " << g_cpuTimeToDeduct
      << " TotalTime " << gaTime << endl;
  }

  if(gaTime > g_cpuTimeNextFlush)
  {
    o_objective_timeline << realGATime << " " << g_bestScore
      << " " << g_nEvaluations << " " << g_nUnfeasible
      << " " << g_cpuTimeSolvingMILP 
      << " " << g_nConflictsInCrossover 
      << " " << g_nUsedCandidates
      << " " << g_nUsedPartners
      << " " << g_nDetectedPartners 
      << " " << g_nDetectedChains
      << " " << g_nUsedChains
      << endl;
    o_objective_timeline.flush();
    //printf("(%f) %f: %f\n",solvingTime,realSolvingTime,g_bestScore);
    g_cpuTimeNextFlush += prob.ga_params.report_score_interval;
  }
  return ret_val;
}

  RNPSolutionPtr
EvaluateUsingWrapper(MyGenome& mg)
{
  float ret_val;
  // your evaluation here
  //cout << "Evaluating:\n" << mg << endl;

  Network *eval_net = new Network(*ga_network_static);
  int k = mg[0];
  for(int i= k; i>0; i--)
  {
    Node &nr = ga_network_with_relays->getRelay(mg[i]);
    eval_net->addNode(nr.x, nr.y,RELAY,nr.id);
  }
  //Graph *G = eval_net->createGraph();
  //output_graphs(G,mg.toString("_"));
  //delete G;



  // Note:  Number of relays in Genome might be 
  // less than K - given that there might be repetitions
  // therefore, define new value of K for model data file
  int genome_k = eval_net->numRelays();

  // Ask LP solver to find solution using all relays in genome
  MyLpSolver *lp = new MyLpSolver(eval_net,genome_k, genome_k);

  Metric *metric;
  if( prob.metric == "PRR")
    metric = new EstPRR_Metric(ga_network_with_relays);
  else if(prob.metric == "DISTANCE_SQ")
    metric = new DistSQ_Metric(ga_network_with_relays);
  lp->setMetric(metric);

  bool solved = lp->solve(); // solve lp
  RNPSolutionPtr sol = lp->solution();


  g_cpuTimeToDeduct += 
    (lp->solution()->cpuTimeTotal() - lp->solution()->cpuTimeSolving());

  delete lp;
  delete eval_net;

  // Before returning, call callback
  myCallback();

  return sol;




}


/// Current working version
/// Version 3 of enhanced crossover operator
  void 
EnhancedCrossoverV3(GA1DArrayGenome<int> &mom, 
		    GA1DArrayGenome<int> &dad, 
		    GA1DArrayGenome<int> *child)
{

  if(prob.ga_params.debug)
    (*g_odebug) <<  "Entering crossover " << endl; 
  set<int> new_genes;
  int max_size = max(mom[0],dad[0]);
  int min_size = min(mom[0],dad[0]);
  double score_dad, score_mom;

  score_dad = dad.score();
  score_mom = mom.score();

  GA1DArrayGenome<int> &best_parent = (score_dad>score_mom?mom:dad);
  GA1DArrayGenome<int> &bad_parent = (score_dad>score_mom?dad:mom);
  /// If parent are the same, create random child
  if( mom.equal(dad))
  {
    (*g_olog) << "Equal parents, creating random child" << endl;
    cout << "Equal parents, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  } 
  else if( score_dad == score_mom)
  {
    (*g_olog) << "Equal scores, creating random child" << endl;
    cout << "Equal scores, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  }

  int my_size = GARandomInt(min_size, max_size);

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "my_size " << my_size << " min_size "
      << min_size << " max_size " << max_size << endl;
  }
  set<int> parent_pool;
  for(int i=mom[0];i>0;i--)
    parent_pool.insert(mom[i]);

  for(int i=dad[0];i>0;i--)
    parent_pool.insert(dad[i]);

  set<int> gene_pool;
  gene_pool.insert( parent_pool.begin(), parent_pool.end());

  /// Put then inside vector
  vector<int> gene_pool_vec(gene_pool.size());
  int k=0;
  FOREACH(it, gene_pool)
  {
    gene_pool_vec[k++] = *it;
  }

  /// ... and shuffle
  std::random_shuffle( gene_pool_vec.begin(), gene_pool_vec.end());
  std::random_shuffle( gene_pool_vec.begin(), gene_pool_vec.end());

  /// number of genes to select from
  int n_uniform_genes = my_size;

  if(prob.ga_params.debug)
    (*g_odebug) <<  "my_size: " << my_size << endl; 
  if(prob.ga_params.debug)
  {
    (*g_odebug) << "UniformGenes: " << n_uniform_genes  
      << " from genePool of size " << gene_pool.size() << endl;
  }
  int c=0;
  bool check_conflicts = prob.ga_params.use_conflicts;
  int conflicts_found=0;
  while(new_genes.size() < n_uniform_genes && c < gene_pool_vec.size())
  {
    int mgene = gene_pool_vec[c++];
    bool conflict = false;
    if( check_conflicts)
    {
      set<int> &cset = g_conflictMap[mgene];
      FOREACH(rr, new_genes)
      {
	if( cset.find(*rr) != cset.end())
	{
	  /// We have a conflict
	  conflict = true;
	  break;
	}
      }
    }
    if( conflict)
    {
      g_nConflictsInCrossover++;
      conflicts_found++;
      /// Remove it from gene pool
      gene_pool.erase(mgene);
      continue;
    }

    new_genes.insert(mgene);
    /// Remove it from gene pool
    gene_pool.erase(mgene);
    /// If we have space for another gene, try partner relays
    if(new_genes.size() < n_uniform_genes && prob.ga_params.use_partners) 
    {
      if(GAFlipCoin(prob.ga_params.pPartnerRelay))
      {
	/// Choose a partner
	set<int> &partner = partner_relays[mgene];
	if(partner.size())
	{
	  set<int>::const_iterator it(partner.begin());
	  int rp = GARandomInt(0,partner.size()-1);
	  std::advance(it,rp);
	  /// Partner found
	  new_genes.insert(*it);
	  /// Just in case, try to remove it from gene_pool
	  /// if it was there
	  gene_pool.erase(*it);
	  /// log used partner
	  g_nUsedPartners++;
	}
      }
    }
  }

  int remaining = n_uniform_genes - (int)new_genes.size();

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "conflicts_found " << conflicts_found 
      << " remaining " << remaining << endl;
  }

  int chain_used = 0;
  if( prob.ga_params.use_chains && remaining > 0 )
  {
    if( GAFlipCoin(prob.ga_params.pUseChain ) )
    {
      /// SO WE ARE GOING TO USE A CHAIN... 
      vector<int> new_genes_v;
      new_genes_v.insert( new_genes_v.end(), new_genes.begin(), new_genes.end());

      /// ... and shuffle
      std::random_shuffle( new_genes_v.begin(), new_genes_v.end());
      std::random_shuffle( new_genes_v.begin(), new_genes_v.end());
      FOREACH(it, new_genes_v)
      {
	ITERATOR(g_relayChains) jt = g_relayChains.find(*it);
	if(jt != g_relayChains.end() )
	{

	  set< set<int> > &chains = jt->second;
	  ITERATOR(chains) kt(chains.begin());
	  int rp = GARandomInt(0, chains.size()-1);
	  std::advance(kt,rp);
	  if( (*kt).size() <= remaining )
	  {
	    /// Chain found
	    if(prob.ga_params.debug)
	    {
	      (*g_odebug) <<  "Using chain ";
	      FOREACH(lt, (*kt) )
	      {
		(*g_odebug) <<  " " << *lt;
	      }
	      (*g_odebug) <<  endl;
	    }
	    
	    new_genes.insert( (*kt).begin(), (*kt).end() );
	    chain_used = (*kt).size();
	    g_nUsedChains++;
	    break;
	  }
	}
      }
    }
  }


  remaining = n_uniform_genes - (int)new_genes.size();
  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "chain_used " << chain_used 
      << " remaining " << remaining << endl;
  }


  if( prob.ga_params.use_candidate_relays)
  {
    if( remaining > gene_pool.size() && candidateRelays.size())
    {
      int n_sel_candidates = 0;
      double pCandidate = 1.0 / candidateRelays.size();
      int n_candidates = remaining - (int)gene_pool.size();

      if(prob.ga_params.debug)
      {
	(*g_odebug) <<  "pCandidate set to " << pCandidate << " candidateRelays "
	  << " has size " << candidateRelays.size() 
	  << " gene_pool " << gene_pool.size()
	  << " missing " << n_candidates <<  endl;
      }
      int j=remaining;
      while( gene_pool.size() < remaining && j--)
      {
	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "Selecting candidates: " << n_sel_candidates  << endl;
	}
	if( GAFlipCoin(0.5))
	{
	  /// Iterate forward
	  FOREACH(it, candidateRelays)
	  {
	    if( remaining <= gene_pool.size())
	      break;
	    if(GAFlipCoin(pCandidate))
	    {
	      if( new_genes.find(*it) == new_genes.end())
	      {
		gene_pool.insert(*it);
		n_sel_candidates++;
	      }
	    }
	  }
	} else {
	  /// Iterate backwards
	  RFOREACH(it, candidateRelays)
	  {
	    if( remaining <= gene_pool.size())
	      break;
	    if(GAFlipCoin(pCandidate))
	    {
	      if( new_genes.find(*it) == new_genes.end())
	      {
		gene_pool.insert(*it);
		n_sel_candidates++;
	      }
	    }
	  }
	}
      }
      if(prob.ga_params.debug)
      {
	(*g_odebug) << "Taking " << n_sel_candidates << " candidate relays " << endl;
      }
      g_nUsedCandidates+=n_sel_candidates;
    }
  }

  remaining = n_uniform_genes - (int)new_genes.size();
  /// if still we need more relays, take random
  if( remaining > gene_pool.size())
  {

    if(prob.ga_params.debug)
      (*g_odebug) <<  "still remaining " << remaining 
	<< " > gene_pool: " << gene_pool.size() << endl;
      while( gene_pool.size() < remaining )
      {
	int rid = getRandomRelay();

	if( new_genes.find(rid) == new_genes.end())
	  gene_pool.insert(rid);
      }
      if( remaining > gene_pool.size())
      {
	fprintf(stderr,"WTF???\n");
	exit(-1);
      }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "gene_pool: " << gene_pool.size() 
      << " new_genes: " << new_genes.size() 
      << " my_size " << my_size << endl;
  }
  if(new_genes.size() < my_size)
  {
    /// Now, take remaining genes from gene pool, 
    /// but take those far away from 
    /// genes we have taken
    if(prob.ga_params.debug)
      (*g_odebug) << "Doing " << my_size - new_genes.size() << " by distance " << 
	" gene_pool: " << gene_pool.size() << endl;

    remaining = my_size - (int)new_genes.size();

    if( remaining == gene_pool.size() )
    {
      /// not much to do, add all
      FOREACH(it, gene_pool)
      {

	if( new_genes.find(*it) != new_genes.end())
	{
	  fprintf(stderr, "gene_pool has duplicated %d\n", *it);
	  exit(-1);
	}
	new_genes.insert(*it);
      }
      if( new_genes.size() != my_size)
      {
	fprintf(stderr, "new_genes: %d my_size %d - this should not happen\n",
		new_genes.size(), my_size);
	exit(-1);
      }
    }
    while( new_genes.size() < my_size)
    {
      vector< pair<double, int> > genes_by_distance;
      FOREACH(itb, gene_pool)
      {
	int bg = *itb;
	Node &r1 = ga_network_with_relays->getRelay(bg);
	double sum_dist = 0.0;

	FOREACH(it, new_genes)
	{
	  Node &r2 = ga_network_with_relays->getRelay(*it);
	  double d = ga_network_with_relays->distanceLink(r1,r2);
	  sum_dist += d;
	}
	double avg_relay_dist = sum_dist / (1.0 * new_genes.size());
	genes_by_distance.push_back(make_pair(avg_relay_dist, bg));
      }
      /// Sort bad_genes according to average distance to good_genes
      std::sort( genes_by_distance.begin(), genes_by_distance.end());


      /// Take farthest one - tail of vector
      for(int i=genes_by_distance.size()-1; i>=0; i--)
      {
	if( new_genes.size() >= my_size)
	  break;
	int bg = genes_by_distance[i].second;
	double meandist = genes_by_distance[i].first;
	/// Shouldn't be there, but check just in case
	if(new_genes.find(bg) == new_genes.end())
	{

	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "Selected relay with mean distance " << meandist << endl;
	  }
	  new_genes.insert(bg);
	  gene_pool.erase(bg);
	  break;
	}
      }
    }
  }

  child->gene(0,new_genes.size());
  assert( new_genes.size() <= prob.ga_params.genome_size);

  /// Finally!, set child's genes using new_genes set
  c=1;
  FOREACH(it, new_genes)
  {
    child->gene(c++, *it);
  }

  /// This should work
  //
  if(prob.ga_params.debug)
    (*g_odebug) <<  "Exit crossover " << endl; 
}

void 
MyCustomCrossover(GA1DArrayGenome<int> &mom,
		       GA1DArrayGenome<int> &dad,
		       GA1DArrayGenome<int> *sis,
		       GA1DArrayGenome<int> *bro)
{
  if( sis)
  {
    EnhancedCrossoverV3(mom,dad,sis);
  }

  if( bro)
  {
    EnhancedCrossoverV3(mom,dad,bro);
  }
}




void MyUniformCrossover(GA1DArrayGenome<int> &mom,
			GA1DArrayGenome<int> &dad,
			GA1DArrayGenome<int> *sis,
			GA1DArrayGenome<int> *bro)
{
  int new_size;
  int max_size = max(mom[0],dad[0]);
  int min_size = min(mom[0],dad[0]);

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "Crossing " << mom << " " << dad << endl;
  }

  // Note: try to keep consistency of genome
  set<int> gene_pool;

  for(int i=dad[0];i>0; i--)
  {
    gene_pool.insert(dad.gene(i));
  }

  for(int i=mom[0];i>0;i--)
  {
    gene_pool.insert(mom.gene(i));
  }

  vector<int> genes(gene_pool.size());
  int c=0;
  FOREACH(it, gene_pool)
  {
    genes[c++] = *it;
  }

  if( sis)
  {
    /// Decide size
    new_size = GARandomInt(min_size, max_size);
    sis->gene(0,new_size);
    std::random_shuffle(genes.begin(), genes.end());
    // Lets do it two times, to be fair
    std::random_shuffle(genes.begin(), genes.end());
    c=1;
    for(int i=0; i< new_size; i++)
    {
      int g = genes[i];
      sis->gene(c++,g);
    }
  }

  if( bro)
  {
    /// Decide size
    new_size = GARandomInt(min_size, max_size);
    bro->gene(0,new_size);
    std::random_shuffle(genes.begin(), genes.end());
    std::random_shuffle(genes.begin(), genes.end());
    c=1;
    for(int i=0; i< new_size; i++)
    {
      int g = genes[i];
      bro->gene(c++,g);
    }
  }
}

  int
MyGenome::Cross(const GAGenome& mom, const GAGenome& dad,
		GAGenome* sis, GAGenome* bro)
{
  GA1DArrayGenome<int> &ga_mom = (GA1DArrayGenome<int> &)mom;
  GA1DArrayGenome<int> &ga_dad = (GA1DArrayGenome<int> &)dad;

  GA1DArrayGenome<int> *ga_sis = (GA1DArrayGenome<int> *)sis;
  GA1DArrayGenome<int> *ga_bro = (GA1DArrayGenome<int> *)bro;

  bool done_sis = false;
  bool done_bro = false;


  /// Replace sister before crossover
  if(sis)
  {
    /// First, try expat
    if( expat_set.size() &&  GAFlipCoin(prob.ga_params.pUseExpat))
    {

      MyGenome *mygen_sis = (MyGenome*)sis;
      MyGenome *expat = expat_set.back();
      MyGenome &expat_ref = *expat;

      expat_set.pop_back();

      mygen_sis->copy(*expat);
      done_sis = true;
      ga_sis = NULL;


      if(prob.ga_params.log_population)
      {
	o_population_log << "EXPAT: " << myGenomeToString(*expat) << endl;
      }
      // Now that it is with us, delete expat
      delete expat;

    }

    /// If still not done, try adoption
    if( !done_sis && GAFlipCoin(prob.ga_params.pAdoption))
    {

      if(prob.ga_params.debug)
	(*g_odebug) << "ADOPTED" << endl;
      ((MyGenome*)ga_sis)->initialize();
      done_sis = true;
      ga_sis = NULL;
    }

  }

  /// Replace broter before crossover
  if(bro)
  {
    /// First, try expat
    if( expat_set.size() &&  GAFlipCoin(prob.ga_params.pUseExpat))
    {

      MyGenome *mygen_bro = (MyGenome*)bro;
      MyGenome *expat = expat_set.back();
      MyGenome &expat_ref = *expat;

      expat_set.pop_back();

      mygen_bro->copy(*expat);
      done_bro = true;
      ga_bro = NULL;


      if(prob.ga_params.log_population)
      {
	o_population_log << "EXPAT: " << myGenomeToString(*expat) << endl;
      }
      // Now that it is with us, delete expat
      delete expat;

    }

    /// If still not done, try adoption
    if( !done_bro && GAFlipCoin(prob.ga_params.pAdoption))
    {

      ((MyGenome*)ga_bro)->initialize();

      if(prob.ga_params.debug)
	(*g_odebug) << "ADOPTED" << endl;
      done_bro = true;
      ga_bro =  NULL;
    }
  }

  /// If we are already done, exit
  if(done_sis && done_bro)
    return 2;


  if(mom.compare(dad)==0)
  {
    if(prob.ga_params.debug)
      (*g_odebug) << "doing uniform because equal" << endl;
    MyUniformCrossover(ga_mom, ga_dad, ga_sis, ga_bro);
  }
  else
  {
    if(GAFlipCoin(prob.ga_params.pUniformCross))
    {
      if(prob.ga_params.debug)
	(*g_odebug) << "doing uniform" << endl;
      MyUniformCrossover(ga_mom,ga_dad,ga_sis,ga_bro);
    }else
    {
      MyCustomCrossover(ga_mom,ga_dad,ga_sis,ga_bro);
    }
  }


  int nc=0;
  if(ga_sis || done_sis)
  {
    nc++;
    ((MyGenome*)sis)->repair();
  }
  if(ga_bro || done_bro)
  {
    nc++;
    ((MyGenome*)bro)->repair();
  }
  return nc;
}


#if 0

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DEPRECATED


/// Version 2 of enhanced crossover operator
/// version 1 significantly biases the procedure
  void 
EnhancedCrossoverV2(GA1DArrayGenome<int> &mom, 
		    GA1DArrayGenome<int> &dad, 
		    GA1DArrayGenome<int> *child)
{

  if(prob.ga_params.debug)
    (*g_odebug) <<  "Entering crossover " << endl; 
  set<int> new_genes;
  int max_size = max(mom[0],dad[0]);
  int min_size = min(mom[0],dad[0]);
  double score_dad, score_mom;

  score_dad = dad.score();
  score_mom = mom.score();

  GA1DArrayGenome<int> &best_parent = (score_dad>score_mom?mom:dad);
  GA1DArrayGenome<int> &bad_parent = (score_dad>score_mom?dad:mom);
  /// If parent are the same, create random child
  if( mom.equal(dad))
  {
    (*g_olog) << "Equal parents, creating random child" << endl;
    cout << "Equal parents, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  } 
  else if( score_dad == score_mom)
  {
    (*g_olog) << "Equal scores, creating random child" << endl;
    cout << "Equal scores, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  }

  int my_size = GARandomInt(min_size, max_size);

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "my_size " << my_size << " min_size "
      << min_size << " max_size " << max_size << endl;
  }
  set<int> parent_pool;
  for(int i=mom[0];i>0;i--)
    parent_pool.insert(mom[i]);

  for(int i=dad[0];i>0;i--)
    parent_pool.insert(dad[i]);

  set<int> gene_pool;
#if 0
  /// Now, remove blacklisted
  std::set_difference(parent_pool.begin(), 
		      parent_pool.end(), 
		      g_redList.begin(), 
		      g_redList.end(),
		      std::inserter(gene_pool, gene_pool.end()));
#endif
  gene_pool.insert( parent_pool.begin(), parent_pool.end());


#if 0

  /// some candidate relays -> gene_pool
  /// Now, we fix the number of candidates to get

  int n_candidates = min_size;

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "n_candidates initially set to " << n_candidates 
      << " min_size " << min_size 
      << " my_size " << my_size 
      << " parent_pool had " << parent_pool.size() 
      << " gene_pool currently has " << gene_pool.size() 
      << " diff " << (my_size - (int)gene_pool.size()) << endl;
  }
  if( (my_size - (int)gene_pool.size()) > n_candidates)
  {
    /// We will need more candidates
    n_candidates = my_size - gene_pool.size();

    if(prob.ga_params.debug)
    {
      (*g_odebug) <<  "not enough relays in gene pool. setting n_candidates to " << n_candidates << endl;
    }
  }

  int n_sel_candidates = 0;
  double pCandidate = 1.0 / candidateRelays.size();

  if(prob.ga_params.debug)
  {
    (*g_odebug) <<  "pCandidate set to " << pCandidate << " candidateRelays "
      << " has size " << candidateRelays.size() <<  endl;
  }
  while( n_sel_candidates < n_candidates)
  {
    if(prob.ga_params.debug)
    {
      (*g_odebug) << "Selecting candidates: " << n_sel_candidates  << endl;
    }
    if( GAFlipCoin(0.5))
    {
      /// Iterate forward
      FOREACH(it, candidateRelays)
      {
	if( n_sel_candidates >= n_candidates)
	  break;
	if(GAFlipCoin(pCandidate))
	{
	  gene_pool.insert(*it);
	  n_sel_candidates++;
	}
      }
    } else {
      /// Iterate backwards
      RFOREACH(it, candidateRelays)
      {
	if( n_sel_candidates >= n_candidates)
	  break;
	if(GAFlipCoin(pCandidate))
	{
	  gene_pool.insert(*it);
	  n_sel_candidates++;
	}
      }
    }
  }

  if(prob.ga_params.debug)
  {
    (*g_odebug) << "Taking " << n_sel_candidates << " candidate relays " << endl;
  }
#endif
  /// Put then inside vector
  vector<int> gene_pool_vec(gene_pool.size());
  int k=0;
  FOREACH(it, gene_pool)
  {
    gene_pool_vec[k++] = *it;
  }

  /// ... and shuffle
  std::random_shuffle( gene_pool_vec.begin(), gene_pool_vec.end());
  std::random_shuffle( gene_pool_vec.begin(), gene_pool_vec.end());

  int n_uniform_genes = my_size;

  if(prob.ga_params.debug)
    (*g_odebug) <<  "my_size: " << my_size << endl; 
  if(prob.ga_params.debug)
  {
    (*g_odebug) << "UniformGenes: " << n_uniform_genes  
      << " from genePool of size " << gene_pool.size() << endl;
  }
  int c=0;
  bool check_conflicts = prob.ga_params.use_conflicts;
  while(new_genes.size() < n_uniform_genes && c < gene_pool_vec.size())
  {
    int mgene = gene_pool_vec[c++];
    bool conflict = false;
    if( check_conflicts)
    {
      set<int> &cset = g_conflictMap[mgene];
      FOREACH(rr, new_genes)
      {
	if( cset.find(*rr) != cset.end())
	{
	  /// We have a conflict
	  conflict = true;
	  break;
	}
      }
    }
    if( conflict)
    {
      /// Remove it from gene pool
      gene_pool.erase(mgene);
      continue;
    }

    new_genes.insert(mgene);
    /// Remove it from gene pool
    gene_pool.erase(mgene);
    /// If we have space for another gene, try partner relays
    if(new_genes.size() < n_uniform_genes) 
    {
      if(GAFlipCoin(prob.ga_params.pPartnerRelay))
      {
	/// Choose a partner
	set<int> &partner = partner_relays[mgene];
	if(partner.size())
	{
	  set<int>::const_iterator it(partner.begin());
	  int rp = GARandomInt(0,partner.size()-1);
	  std::advance(it,rp);
	  /// Partner found
	  new_genes.insert(*it);
	  /// Just in case, try to remove it from gene_pool
	  /// if it was there
	  gene_pool.erase(*it);
	  /// log used partner
	  g_nUsedPartners++;
	}
      }
    }
  }

  int remaining = my_size - (int)new_genes.size();
  if( prob.ga_params.use_candidate_relays)
  {
    if( remaining > gene_pool.size() && candidateRelays.size())
    {

      int n_sel_candidates = 0;
      double pCandidate = 1.0 / candidateRelays.size();
      int n_candidates = remaining - (int)gene_pool.size();

      if(prob.ga_params.debug)
      {
	(*g_odebug) <<  "pCandidate set to " << pCandidate << " candidateRelays "
	  << " has size " << candidateRelays.size() 
	  << " gene_pool " << gene_pool.size()
	  << " missing " << n_candidates <<  endl;
      }
      int j=remaining;
      while( gene_pool.size() < remaining && j--)
      {
	if(prob.ga_params.debug)
	{
	  (*g_odebug) << "Selecting candidates: " << n_sel_candidates  << endl;
	}
	if( GAFlipCoin(0.5))
	{
	  /// Iterate forward
	  FOREACH(it, candidateRelays)
	  {
	    if( remaining <= gene_pool.size())
	      break;
	    if(GAFlipCoin(pCandidate))
	    {
	      if( new_genes.find(*it) == new_genes.end())
	      {
		gene_pool.insert(*it);
		n_sel_candidates++;
	      }
	    }
	  }
	} else {
	  /// Iterate backwards
	  RFOREACH(it, candidateRelays)
	  {
	    if( remaining <= gene_pool.size())
	      break;
	    if(GAFlipCoin(pCandidate))
	    {
	      if( new_genes.find(*it) == new_genes.end())
	      {
		gene_pool.insert(*it);
		n_sel_candidates++;
	      }
	    }
	  }
	}
      }
      if(prob.ga_params.debug)
      {
	(*g_odebug) << "Taking " << n_sel_candidates << " candidate relays " << endl;
      }
    }
  }

  remaining = my_size - (int)new_genes.size();
  /// if still we need more relays, take random
  if( remaining > gene_pool.size())
  {

    if(prob.ga_params.debug)
      (*g_odebug) <<  "still remaining " << remaining 
	<< " > gene_pool: " << gene_pool.size() << endl;
      while( gene_pool.size() < remaining )
      {
	int rid = getRandomRelay();

	if( new_genes.find(rid) == new_genes.end())
	  gene_pool.insert(rid);
      }
      if( remaining > gene_pool.size())
      {
	fprintf(stderr,"WTF???\n");
	exit(-1);
      }
  }

  if(prob.ga_params.debug)
    (*g_odebug) <<  "gene_pool: " << gene_pool.size() 
      << "new_genes: " << new_genes.size() << endl;
  if(new_genes.size() < my_size)
  {

    /// Now, take remaining genes from gene pool, 
    /// but take those far away from 
    /// genes we have taken

    if(prob.ga_params.debug)
      (*g_odebug) << "Doing " << my_size - new_genes.size() << " by distance " << 
	" gene_pool: " << gene_pool.size() << endl;

    remaining = my_size - (int)new_genes.size();

    if( remaining == gene_pool.size() )
    {
      /// not much to do, add all
      FOREACH(it, gene_pool)
      {

	if( new_genes.find(*it) != new_genes.end())
	{
	  fprintf(stderr, "gene_pool has duplicated %d\n", *it);
	  exit(-1);
	}
	new_genes.insert(*it);
      }
      if( new_genes.size() != my_size)
      {
	fprintf(stderr, "new_genes: %d my_size %d - this should not happen\n",
		new_genes.size(), my_size);
	exit(-1);
      }
    }
    while( new_genes.size() < my_size)
    {
      vector< pair<double, int> > genes_by_distance;
      FOREACH(itb, gene_pool)
      {
	int bg = *itb;
	Node &r1 = ga_network_with_relays->getRelay(bg);
	double sum_dist = 0.0;

	FOREACH(it, new_genes)
	{
	  Node &r2 = ga_network_with_relays->getRelay(*it);
	  double d = ga_network_with_relays->distanceLink(r1,r2);
	  sum_dist += d;
	}
	double avg_relay_dist = sum_dist / (1.0 * new_genes.size());
	genes_by_distance.push_back(make_pair(avg_relay_dist, bg));
      }
      /// Sort bad_genes according to average distance to good_genes
      std::sort( genes_by_distance.begin(), genes_by_distance.end());


      /// Take farthest one - tail of vector
      for(int i=genes_by_distance.size()-1; i>=0; i--)
      {
	if( new_genes.size() >= my_size)
	  break;
	int bg = genes_by_distance[i].second;
	double meandist = genes_by_distance[i].first;
	/// Shouldn't be there, but check just in case
	if(new_genes.find(bg) == new_genes.end())
	{

	  if(prob.ga_params.debug)
	  {
	    (*g_odebug) << "Selected relay with mean distance " << meandist << endl;
	  }
	  new_genes.insert(bg);
	  gene_pool.erase(bg);
	  break;
	}
      }
    }
  }

  child->gene(0,new_genes.size());
  assert( new_genes.size() <= prob.ga_params.genome_size);

  /// Finally!, set child's genes using new_genes set
  c=1;
  FOREACH(it, new_genes)
  {
    child->gene(c++, *it);
  }

  /// This should work
  //
  if(prob.ga_params.debug)
    (*g_odebug) <<  "Exit crossover " << endl; 
}




  void 
CreateOneChild(GA1DArrayGenome<int> &mom, 
	       GA1DArrayGenome<int> &dad, 
	       GA1DArrayGenome<int> *child)
{
  int max_size = max(mom[0],dad[0]);
  int min_size = min(mom[0],dad[0]);
  double score_dad, score_mom;

  score_dad = dad.score();
  score_mom = mom.score();

  GA1DArrayGenome<int> &best_parent = (score_dad>score_mom?mom:dad);
  GA1DArrayGenome<int> &bad_parent = (score_dad>score_mom?dad:mom);
  /// If parent are the same, create random child
  if( mom.equal(dad))
  {
    (*g_olog) << "Equal parents, creating random child" << endl;
    cout << "Equal parents, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  } else if( score_dad == score_mom)
  {
    (*g_olog) << "Equal scores, creating random child" << endl;
    cout << "Equal scores, creating random child" << endl;
    ((MyGenome*)child)->initialize();
    g_olog->flush();
    return;
  }
  /// Decide size
  int best_size = best_parent[0];
  int bad_size = bad_parent[0];


  /// best parent genes -> good_genes
  set<int> good_genes;
  for(int i=best_size;i>0;i--)
    good_genes.insert(best_parent[i]);

  /// some candidate relays -> good_genes
  FOREACH(it, candidateRelaySet)
  {
    if(GAFlipCoin(prob.ga_params.pCandidateRelay))
    {
      good_genes.insert((*it).second);
    }
  }

  /// bad parent genes -> bad_genes
  set<int> bad_genes;
  for(int i=bad_size;i>0;i--)
    bad_genes.insert(bad_parent[i]);


  set<int> bad_gene_pool;

  /// Compute bad_genes diff good_genes
  /// that will be our bad_gene_pool
  std::set_difference(bad_genes.begin(), 
		      bad_genes.end(), 
		      good_genes.begin(), 
		      good_genes.end(),
		      std::inserter(bad_gene_pool, bad_gene_pool.end()));


#if 0
  /// If there are no bad genes, let the crossover still do the 
  /// cross with the candidate relays
  /// Ok, if bad_gene_pool is empty (i.e. bad_parent is contained in best_parent)
  /// crossover just copies best_parent
  if(!bad_gene_pool.size())
  {
    if(prob.ga_params.debug)
      printf("Parents contained, returning same\n");
    for(int i=best_size;i>0;i--)
    {
      child->gene(i,best_parent[i]);
    }
    child->gene(0, best_size);
    /// and... we are done
    return;
  }
#endif


  /// Select how many genes from best parent
  /// at least half, but taking care of how many bad genes we
  /// will need to take to fill the genome

  /// If best_size=1, take all best_genes
  int lb_good_size = max(best_size/2, 1);

  /// If just few bad genes in bad_gene_pool
  /// take more genes from best parent
  if( (best_size - lb_good_size) > bad_gene_pool.size())
  {
    lb_good_size = best_size - bad_gene_pool.size();
  }

  int n_good_genes = GARandomInt(lb_good_size, best_size);
  if(n_good_genes <=0)
    n_good_genes = 1;

  if(prob.ga_params.debug)
    (*g_odebug) << "Taking " << n_good_genes << " good genes " 
      << " best_size=" << best_size 
      << " bad_size=" << bad_gene_pool.size() << endl;

  /// This will contain our child's genes
  set<int> new_genes;

  /// Lets take good genes first
  vector<int> good_genes_vec(good_genes.size());
  int k=0;
  FOREACH(it, good_genes)
  {
    good_genes_vec[k++] = *it;
  }

  /// Random permutation, take first n_good_genes
  std::random_shuffle( good_genes_vec.begin(), good_genes_vec.end());
  std::random_shuffle( good_genes_vec.begin(), good_genes_vec.end());

  int c=0;
  while(new_genes.size() < n_good_genes)
  {
    int mgene = good_genes_vec[c++];
    new_genes.insert(mgene);
    /// If we have space for another gene, try partner relays
    if(new_genes.size() < n_good_genes) {
      if(GAFlipCoin(prob.ga_params.pPartnerRelay))
      {
	/// Choose a partner
	set<int> &partner = partner_relays[mgene];
	if(partner.size())
	{
	  set<int>::const_iterator it(partner.begin());
	  int rp = GARandomInt(0,partner.size()-1);
	  std::advance(it,rp);
	  /// Partner found
	  new_genes.insert(*it);
	}
      }
    }
  }

  if(new_genes.size() < best_size)
  {

    /// Now, take genes from bad parent, but take those far away from 
    /// genes we have taken

    vector< pair<double, int> > bad_genes_by_distance;
    FOREACH(itb, bad_gene_pool)
    {
      int bg = *itb;
      Node &r1 = ga_network_with_relays->getRelay(bg);
      double sum_dist = 0.0;

      FOREACH(it, new_genes)
      {
	Node &r2 = ga_network_with_relays->getRelay(*it);
	double d = ga_network_with_relays->distanceLink(r1,r2);
	sum_dist += d;
      }
      double avg_relay_dist = sum_dist / (1.0 * new_genes.size());
      bad_genes_by_distance.push_back(make_pair(avg_relay_dist, bg));
    }
    /// Sort bad_genes according to average distance to good_genes
    std::sort( bad_genes_by_distance.begin(), bad_genes_by_distance.end());


    /// Take farthest ones
    for(int i=bad_genes_by_distance.size()-1; i>=0; i--)
    {
      if( new_genes.size() >= best_size)
	break;
      int bg = bad_genes_by_distance[i].second;
      if(new_genes.find(bg) == new_genes.end())
      {
	new_genes.insert(bg);
      }
    }
  }

  child->gene(0,new_genes.size());
  assert( new_genes.size() <= prob.ga_params.genome_size);

  /// Finally!, set child's genes using new_genes set
  c=1;
  FOREACH(it, new_genes)
  {
    child->gene(c++, *it);
  }

  /// GAME OVER
}

void GASol::testRandom()
{
  CpuTime cput;
  cput.start();
  for(int i=10000000; i>=0; i--)
  {
    GAFlipCoin(0.5);
  }
  cout << "Time elapsed " << cput.cpu_time_elapsed() << endl;




}





/** OLD STUFF - DO NOT USE
  void MyOldCrossover(GA1DArrayGenome<int> &mom,
  GA1DArrayGenome<int> &dad,
  GA1DArrayGenome<int> *sis,
  GA1DArrayGenome<int> *bro)
  {
  int new_size;
  int max_size = max(mom[0],dad[0]);
  int min_size = min(mom[0],dad[0]);
  if( sis)
  {
/// Decide size
new_size = GARandomInt(min_size, max_size);
/// Fill 
if( max_size > min_size)
{
if(mom[0] == max_size)
{
for(int i=new_size; i>min_size;i--)
{
sis->gene(i,mom[i]);
}
}else
{
for(int i=new_size; i>min_size;i--)
{
sis->gene(i,dad[i]);
}
}
}

sis->gene(0,new_size);
unsigned int i = min_size;
while(i>0)
{
std::set<int> relaysSelected;

if(GAFlipCoin(0.5))
{
/// choose mom's gene
sis->gene(i,mom[i]);
relaysSelected.insert((*sis)[i]);
i--;
if(i>0)
{
sis->gene(i,findGeneMatch(mom[i], dad,relaysSelected));
relaysSelected.insert((*sis)[i]);
i--;
continue;
}else
{
break;
}
}
else
{
/// choose dad's gene
sis->gene(i,dad[i]);
relaysSelected.insert((*sis)[i]);
i--;
if(i>0)
{
sis->gene(i,findGeneMatch(dad[i], mom,relaysSelected));
relaysSelected.insert((*sis)[i]);
i--;
continue;
}else
{
break;
}
}
}
}

if( bro)
{
  /// Decide size
  new_size = GARandomInt(min_size, max_size);
  /// Fill 
  if( max_size > min_size)
  {
    if(mom[0] == max_size)
    {
      for(int i=new_size; i>min_size;i--)
      {
	bro->gene(i, mom[i]);
      }
    }else
    {
      for(int i=new_size; i>min_size;i--)
      {
	bro->gene(i,dad[i]);
      }
    }
  }

  bro->gene(0,new_size);
  unsigned int i = min_size;
  while(i>0)
  {
    std::set<int> relaysSelected;

    if(GAFlipCoin(0.5))
    {
      /// choose mom's gene
      bro->gene(i,mom[i]);
      relaysSelected.insert((*bro)[i]);
      i--;
      if(i>0)
      {
	bro->gene(i,findGeneMatch(mom[i], dad,relaysSelected));
	relaysSelected.insert((*bro)[i]);
	i--;
	continue;
      }else
      {
	break;
      }
    }
    else
    {
      /// choose dad's gene
      bro->gene(i,dad[i]);
      relaysSelected.insert((*bro)[i]);
      i--;
      if(i>0)
      {
	bro->gene(i,findGeneMatch(dad[i], mom,relaysSelected));
	relaysSelected.insert((*bro)[i]);
	i--;
	continue;
      }else
      {
	break;
      }
    }
  }
}
}
*/

#endif
