/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/15/2010 11:08:37 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include "main.h"
#include "Random.h"
#include "generator.h"
#include "graph.h"
#include "network.h"
#include "lpsolver.h"
#include "ga_solver.h"
#include "GetPot.h"
#include "heuristic.h"
#include "mst.h"
#include "solver.h"

#ifndef MATHEURISTIC_DISABLED
#include "matheuristic.h"
#endif

#ifndef CPLEX_DISABLED
#include "cplex_solver.h"
#endif

/** 
    CPLEX C++ API not supported yet
    #ifdef ILOG_CPLEX_CPP_ENABLED
    #define COMPILE_CPLEX_CPP
    #include "cplex_solver_v2.h"
    #endif
**/


#include "rnpsolution.h"

///NOTES:
/// fixing relays only works with mycplex - 01/2013

extern ofstream cdebug;  // debug output stream

extern param prob;
void generate_grid(Network *net, double dimX, double dimY, double r );
void generate_random_nodes(Network *net, RandomGenerator *, double dimX, double dimY,double,double, int num, nodeType t);
GetPot read_config_file(char *config_file);
void load_defaults(param &);
void load_config_file(GetPot infile, param &);
void checkParameters();
Network * los_band_restriction(Network *, double, double);
void generate_cluster_nodes(Network *net, RandomGenerator *gen, double dimX, double dimY,
			    int num_static, int num_base, clusters_t &clu);
#if WITH_OLD_BOOST_GEOMETRY
// disabled until we fix dependencies error with boost geometry
void generate_grid_h2(Network *net, double dimX, double dimY, double min_res, double max_res);
Network * convex_hull_restriction(Network *);
#endif
void generate_simulation(Graph *G, Network *net, string filename);

int fexist( char *filename );
void init_debug_file(string);
void close_debug_file();
void output_graphs(Graph *G, string filename);



void
testBoostSerialize()
{
  LpSolution *lpsol = new LpSolution();
  lpsol->value["TEST"] = 12345;
  lpsol->objval = 12345;
  int buf_len = 65536;
  int msg_size;
  std::vector<char> buff_vec;

  printf("Attempting to serialize lpsol %f\n", lpsol->objval);
  save_lpsolution(*lpsol, buff_vec, buf_len);
  printf("Serialized BINARY size %d\n", buff_vec.size());

}

int redirect_output(string ofile)
{
  int file = open(ofile.c_str(),O_WRONLY | O_CREAT, 0600);
  if( file < 0) {
    printf("Error opening %s\n", ofile.c_str());
    return 1;
  }

  /* yes, we are about to close standard output */


  /* Now we call dup to create a duplicate file descriptor. Because dup always
     uses the smallest available file descriptor and we just closed file descriptor
     1, the return value on dup will be 1. */

  if(dup2(file, STDOUT_FILENO)!=1)
  {
    printf("Error in dup2\n");
    return 1;
  }

  close(file);
  //if(dup(file) != 1)    return 1;

  //return dup2(file, fileno(stdout));
  /* We don't need our original file descriptor, because now the file descriptor 1
     points to our file. */
  //close(file);
  return 0;

}


void
output_connectivity(string filename, int n)
{
  ofstream connFile(filename.c_str());
  connFile << "CONNECTED_COMPONENTS " << n;
  connFile.close();
}


void
output_plain(Network *net, string filename)
{
  /*  Generate plain file with description */
  string plain_file = prob.output_path + prob.instance_id + "-" + filename + ".net";

  ofstream netFile (plain_file.c_str());

  for(int i=0;i<net->size();i++){
    Node &n = net->getNode(i);
    netFile << "node\t" << n.id << " " << n.x << " " << n.y << endl;
  }
  netFile.close();


}

//TODO: Check if demand must be integer,if we use integer flows
void
setDemandFromFile(string filename, Network *net)
{
  printf("Setting demand from file %s\n", 
	 filename.c_str());
  ifstream dfile( filename.c_str());
  int num;
  string id;
  double demand;
  dfile >> num;
  while(num--)
  {
    dfile >> id >> demand;
    //      printf("Setting d[%s] = %.2f\n",id.c_str(),demand);
    int ix = net->getNodeIndex(id);
    //      printf("Node index for %s is %d\n",id.c_str(), ix);
    net->setDemand(ix, demand);
  }

}
#define USE_CONVEXHULL 0
#define USE_LOS_STRATEGY 0
#define DOTESTS 0

//#define TEST_COMM 1
int main(int argc, char *argv[])
{

  GetPot infile;
  string cmd;
  Network *net_hull,*net_los;
  Graph *G, *sol_graph;
  Heuristic *heuristic = NULL;
  Network *net;

  load_defaults(prob);

  /// Read command line first
  GetPot cl(argc,argv);


  /// Read problem configuration files
  vector<string> configfiles;
  cl.init_multiple_occurrence();
  string configfile = cl.follow("none", 1,"-c");
  if(configfile != "none")
    configfiles.push_back(configfile);
  while(configfile != "none")
  {
    configfile = cl.follow("none", 1,"-c");
    if(configfile != "none")
      configfiles.push_back(configfile);
  }
  cl.enable_loop();


  if( !configfiles.size())
  {
    fprintf(stderr, "No config files. Aborting...\n");
    exit(1);
  }
  FOREACH(fileit, configfiles)
  {
    configfile = *fileit;
    printf("Using config file %s\n",configfile.c_str());
    if( fexist((char *)configfile.c_str()))
    {
      infile = read_config_file((char *)configfile.c_str());
      load_config_file(infile, prob);
    }
    else
    {
      fprintf(stderr, "Invalid config file %s\n", configfile.c_str());
      exit(1);
    }
  }

  /// Check for additional parameters (not in config files)
  /// --sid : suffix_id  : string to append at the end of instance id
  string suffix_id = cl.follow("none", 1,"--sid");
  if(suffix_id != "none")
    prob.instance_id = prob.instance_id + suffix_id;

  prob.instance_path = cl.follow(prob.instance_path.c_str(), 1,"--instance-path");
  prob.output_path = cl.follow(prob.output_path.c_str(), 1,"--output-path");
  prob.generation_path = cl.follow(prob.generation_path.c_str(), 1,"--generation-path");
  /// Finally, check parameters and fix file paths
  checkParameters();

  printf("OUTPUT_PATH %s\n", prob.output_path.c_str());
  /// Initialize debug stream

  DEBUG(1)
  {
    string dbgfile = 
      prob.output_path  + prob.debug_prefix 
      + prob.instance_id + ".txt";
    cout << "Debugging to " << dbgfile << endl;
    init_debug_file(dbgfile);
  }
  if( DOTESTS )
  {
    testBoostSerialize();
    printf("Tests ended\n");
    exit(1);
  }

  if(!prob.generate_instance)
  {
    net = new Network(prob.instance_path + prob.network_file);
    net->print_info();
  }
  else
  {
    /** Display general information */

    VERBOSE(1){
      printf("generating random instance:\n");
      printf("	static nodes: %d\n",prob.st_nodes);
      printf("	  base nodes: %d\n",prob.b_nodes);
      printf(" transmission range: %f \n",prob.tx_range);
    }

    /** Network creation */
    net = new Network(); /* create an empty network structure */
    net->setRange(prob.tx_range); /* set transmission range for the network */
    net->setRegionSize(prob.dimX, prob.dimY); /* set area size */

    UniformRandom *rand = new UniformRandom(prob.seed);


    if(prob.use_clusters){
      printf("Generating clusters...\n");
      generate_cluster_nodes(net,rand,prob.dimX,prob.dimY,
			     prob.st_nodes,prob.b_nodes, prob.clu);
    }else{
      printf("Generating random locations\n");
      //generate random static node locations
      generate_random_nodes(net, rand, prob.dimX, 
			    prob.dimY,0,0,prob.st_nodes, STATIC);   

      //generate random base node locations
      generate_random_nodes(net, rand, prob.dimX, 
			    prob.dimY,0,0,prob.b_nodes, BASE);  
    }
    delete rand;
  }

  if( prob.demand_type == "file")
  {
    /// removed instance_path prefix
    //	    setDemandFromFile(prob.instance_path + prob.demand_file, net);
    setDemandFromFile(prob.demand_file, net);
  }
  if( net->numRelays() > 1000)
  {
    cerr << "Number of relays is restricted to be < 1000" << endl;
    exit(1);
  }


  if(prob.gen_dot)
  {
    G = net->createGraph();
    output_graphs(G,"static");
    delete G;
  }

  if(prob.gen_net)
  {
    output_plain(net,"static");
  }

  //Save a copy of static network
  Network *static_net = new Network(*net);
  //	printf("DEBUG demand %f\n",static_net->getDemand(2));


  /* If asked to generate simulation files
   * try to solve without any relays
   * */
  /*  
      if(prob.output_sim){
      GLPSol *lp = new GLPSol(net);

      lp_solution_t *lp_sol = lp->solve(prob.K,NULL); // solve lp
      Graph *static_graph = lp->create_graph_from_sol(lp_sol);
      output_graphs(static_graph, "static-sol");
      generate_simulation(static_graph, net, "static");
      }
      */
  if(prob.use_heuristic && prob.solve && prob.solution_strategy == "default"){
    if(prob.heuristic == "Low_res")
    {
      /// TODO: Enable use of heuristics
      /*
	printf("Creating heuristic Low_res (factor = %d)\n",prob.Low_res_factor);
	if(prob.Low_res_use_file)
	{
	lp_solution_t *prev_sol = new lp_solution_t(prob.Low_res_solution_file);
	heuristic = new Low_res_heuristic(prev_sol,
	prob.grid_res,prob.Low_res_factor);
	delete prev_sol;
	}
	else{
	heuristic = new Low_res_heuristic(net,prob.grid_res,prob.Low_res_factor);
	}
	*/


    }
  }else
  {
    heuristic = NULL;
  }

  if(prob.generate_instance || 
     (prob.solve && !(prob.lp_params.use_mps && prob.solution_strategy == "default")))
  {
    if(
      (prob.solution_strategy == "default" || prob.solution_strategy == "ga" 
       || prob.solution_strategy == "mth" || prob.solution_strategy == "multimth_cplex" 
       || prob.solution_strategy == "multimth_ga" || prob.solution_strategy=="mycplex") && prob.K != 0)
    {
      if(prob.grid_strategy == "default")
      {
	generate_grid(net, prob.dimX, prob.dimY,prob.grid_res); 
      }
      else if(prob.grid_strategy == "grid_h2")
      {
#if WITH_OLD_BOOST_GEOMETRY
	generate_grid_h2(net,prob.dimX, prob.dimY, prob.gridh2_minr, prob.gridh2_maxr);
#else
	fprintf(stderr,"grid_h2 not compiled\n");
	exit(1);
#endif
      }
      else 
      {
	printf("No relays generated (maybe read from file?)\n");
	//exit(-1);
      }
      printf("number of relays %d\n",net->numRelays());
      if(prob.use_convex_hull && !prob.no_position){

#if WITH_OLD_BOOST_GEOMETRY
	net_hull = convex_hull_restriction(net);
	delete net;
	net = net_hull;
	printf("number of relays (hull) %d\n",net->numRelays());
#else
	fprintf(stderr,"convex_hull not compiled\n");
	exit(1);
#endif
      }

      if(USE_LOS_STRATEGY)
      {
	net_los = los_band_restriction(net_hull, net->getRange(), prob.dimX);
	delete net;
	net = net_los;
	printf("number of relays (LOS) %d\n",net->numRelays());
      }
      if(prob.fixed_relays)
      {
	///NOTE: Theres a problem doing this approach:
	///Metric is affected by the removal of other relays
	///in fact, fixing relays should force the ilp to use
	///those relays, and allowing to use others
	///FIXED: JAN 2013
	/*
	  Network *fixed_network = new Network(*static_net);
	  printf("Using fixed relays set of size %d\n", 
	  prob.relays.size());
	  for(int i=0; i<prob.relays.size();i++)
	  {
	  Node &n = net->getNodeById(prob.relays[i]);
	  if(!fixed_network->addNode(n))
	  {
	  fprintf(stderr,"Error adding node %s\n",
	  n.id.c_str());
	  exit(1);
	  }
	  }
	  delete net;
	  net = fixed_network;
	  */
      }

      if(prob.gen_dot)
      {
	G = net->createGraph();
	output_graphs(G,"static+grid");
	delete G;
      }

      if(prob.gen_net)
      {
	output_plain(net,"static+grid");
      }
      /* plot node locations */
      string plot_filename = prob.output_path + prob.instance_id + ".svg";
      //Plotter *plot = net->plot(plot_filename);
      //plot->write();
    }
  }

  if(prob.generate_instance)
  {
    if(prob.lp_params.save_mps && prob.solution_strategy == "default")
    {
      //TODO Enable mps stand-alone generation

      /** Generate mps and save it to file */
      //GLPSol::generate_and_save_mps(net,prob.generation_path + prob.output_mps_file);
    }
    net->write(prob.generation_path + prob.network_file);
  }


  if(prob.solve)
  {
    if( prob.solution_strategy == "default")
    {
      /** Solve problem */
      MyLpSolver *lp;
      if(prob.lp_params.use_mps){
	//TODO: Enable use of existing mps file
	//lp = new MyLpSolver(net, prob.instance_path + prob.input_mps_file);
      }else
      {
	if(prob.minK)
	  lp = new MyLpSolver(net,prob.K, prob.minK);
	else
	  lp = new MyLpSolver(net,prob.K);

      }

      bool solved = lp->solve();
      RNPSolutionPtr sol = lp->solution();
      //lp_solution_t *lp_sol = lp->solve(prob.K,heuristic); // solve lp
      if(prob.output_sol)
      {
	lp->write_time_to_file(prob.output_path + prob.output_time_file);
	lp->write_solution_to_file(prob.output_path + prob.output_sol_file);
      }

      if(solved)
      {
	sol_graph = sol->create_graph_from_sol(net);

	if(prob.output_sim){
	  generate_simulation(sol_graph, net, "sol");
	}
	/** print output data */

	output_graphs(sol_graph,"solution");
	delete sol_graph;

	if( sol->isOptimal())
	  printf("Optimal value: ");
	else
	  printf("SubOptimal value: ");
	printf("%f\n",sol->getObjVal());

	// Just to check
	/* comment me
	   vector<string> &rtouse = lp->relaysToUse;
	   for(int r=0; r < rtouse.size(); r++)
	   {
	   double rflow = lp->getOutgoingFlow(rtouse[r]);
	   cout << "Relay " << rtouse[r] << ": " << rflow << endl;
	   }
	   */
      }
      else
      {
	printf("Solution not found\n");
      }


      delete lp;

    }
    else if(prob.solution_strategy == "mst")
    {
      MST_heuristic(net);
      G = net->createGraph();
      output_graphs(G,"static+grid");
      output_plain(net,"static+grid");
      delete G;

    }
    else if(prob.solution_strategy == "ga")
    {
      printf("Using GA\n");
      GASol *ga = new GASol(static_net, net);
      double ga_sol = ga->solve();
      //ga->testRandom();
      //double ga_sol = ga->test();

      bool solved  = ga->hasSolution();

      if(solved)
      {
	RNPSolutionPtr sol = ga->solution();
	if(prob.output_sol)
	{
	  sol->writeToFile(prob.output_path + prob.output_sol_file);
	}

	sol_graph = sol->create_graph_from_sol(net);

	output_graphs(sol_graph,"solution");
	delete sol_graph;
      }
      printf("Solution value: %f\n",ga_sol);
      delete ga;
      printf("GA: Terminated\n");
    }
    else if(prob.solution_strategy == "connectivity")
    {
      Graph *G = static_net->createGraph();
      int n = numberOfConnectedComponents(G);
      output_connectivity(prob.output_path + prob.output_conn_file, n);
    }
    else if(prob.solution_strategy == "mycplex")
    {

      RNPSolutionPtr init_sol;
      if( prob.lp_params.use_initial_k0)
      {
	if( prob.K == 0)
	{
	  printf("WARNING: the initial solution would be the optimal DUMB!\n");
	  exit(-1);
	}
	MyCplexSolver *cplex_solver_k0 = new MyCplexSolver(static_net, static_net);
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
      }
      //		  printf("DEBUG demand %f\n",net->getDemand(2));
      MyCplexSolver *cplex_solver = new MyCplexSolver(static_net, net);

      if( prob.lp_params.use_initial_k0)
      {
	printf("Using initial k0\n");
	cplex_solver->useInitial(init_sol);
      }

      ///fixed:jan 2013 see above
      if(prob.fixed_relays)
      {

	printf("Using fixed relays set of size %d\n", 
	       prob.relays.size());
	cplex_solver->disableAllRelays();
	for(int i=0; i<prob.relays.size();i++)
	{
	  cplex_solver->enableRelay(prob.relays[i]);
	}
      }
      bool solved  = cplex_solver->solve();

      if(solved)
      {
	RNPSolutionPtr sol = cplex_solver->solution();

	if( cplex_solver->solution()->isOptimal())
	  printf("Optimal value: ");
	else
	  printf("SubOptimal value: ");
	printf("%f\n",cplex_solver->getObjVal());

	if(prob.output_sol)
	{
	  sol->writeToFile(prob.output_path + prob.output_sol_file);
	}

	sol_graph = sol->create_graph_from_sol(net);

	if(prob.output_sim){
	  generate_simulation(sol_graph, net, "sol");
	}

	output_graphs(sol_graph,"solution");
	delete sol_graph;
      }


    }
    //        else if(prob.solution_strategy == "random


    else if(prob.solution_strategy == "mycplexv2")
    {
#ifdef COMPILE_CPLEXCPP
      MyCplexSolverV2 *cplex_solver = new MyCplexSolverV2(static_net, net);
      bool solved  = cplex_solver->solve();

      if(solved)
      {
	if( cplex_solver->solution()->isOptimal())
	  printf("Optimal value: ");
	else
	  printf("SubOptimal value: ");
	printf("%f\n",cplex_solver->getObjVal());
      }

#else
      printf("CPLEX CPP NOT COMPILED\n");
#endif
    }
    else if(prob.solution_strategy == "mth")
    {
      // Matheuristic (ILP + GA in parallel)
      // Child will be GA
      // Parent will be ILP Solver
      printf("Solving using SIMPLE MATHEURISTIC\n");
      Matheuristic *mth = new Matheuristic();

      // Here's where the stuff begins
      pid_t pID = fork();
      if (pID == 0)  // child
      {
	// Change some parameters 
	prob.instance_id = prob.instance_id + "_MT_GA";
	string stdout_file = prob.output_path + 
	  prob.instance_id + ".stdout";
	prob.ga_params.time_limit = prob.matheuristic_params.time_limit;

	string mth_outfile = prob.output_path + 
	  prob.instance_id + ".mthout";

	printf("FORK OK. Im child\n");
	printf("Child: Redirecting output to %s\n", 
	       stdout_file.c_str());
	if(redirect_output(stdout_file))
	  printf("Error redirecting stdout\n");

	//				prob.lp_params.verbose = 0;
	mth->setServer(false);
	mth->setOutput(mth_outfile);
	GASol *ga = new GASol(static_net, net);
	ga->setMatheuristic(mth);

#ifdef TEST_COMM
	double rcv_val = 12345;
	printf("Testing communication - SENDING %d\n", rcv_val);
	mth->testPush(rcv_val);
	mth->testPull(&rcv_val);
	printf("TEST: Received %f\n", rcv_val);
#else
	double ga_sol = ga->solve();
#endif
      }
      else if (pID < 0)  // failed to fork
      {

      }
      else     // parent
      {
	prob.instance_id = prob.instance_id + "_MT_LP";
	prob.lp_params.strong_timelim = prob.matheuristic_params.time_limit;
	string stdout_file = prob.output_path  + 
	  prob.instance_id + ".stdout";

	string mth_outfile = prob.output_path + 
	  prob.instance_id + ".mthout";
	printf("FORK OK. Im Parent\n");
	printf("Parent: Redirecting output to %s\n", 
	       stdout_file.c_str());

	if(redirect_output(stdout_file))
	  printf("Error redirecting stdout\n");
	// This are the required parameters
	// to make the MTH work
	prob.lp_params.verbose = 1;
	prob.lp_params.noMIP = 0;
	prob.lp_params.force_branch_and_cut = 1;

	MyCplexSolver *cplex_solver = new MyCplexSolver(static_net, net);


	mth->setServer(true);

	mth->setOutput(mth_outfile);
	cplex_solver->setMatheuristic(mth);
	bool solved = false;
#ifdef TEST_COMM
	printf("Testing communication\n");
	double rcv_val;
	mth->testPull(&rcv_val);
	printf("TEST: Received %f\n", rcv_val);
	printf("Testing communication - SENDING %d\n", rcv_val);
	mth->testPush(rcv_val);
#else 
	solved = cplex_solver->solve();
	RNPSolutionPtr sol = cplex_solver->solution();

#endif
	if(solved)
	{
	  RNPSolutionPtr sol = cplex_solver->solution();

	  if( cplex_solver->solution()->isOptimal())
	    printf("Optimal value: ");
	  else
	    printf("SubOptimal value: ");
	  printf("%f\n",cplex_solver->getObjVal());

	  if(prob.output_sol)
	  {
	    sol->writeToFile(prob.output_path + prob.output_sol_file);
	  }

	  sol_graph = sol->create_graph_from_sol(net);

	  if(prob.output_sim){
	    generate_simulation(sol_graph, net, "sol");
	  }
	  output_graphs(sol_graph,"solution");
	  delete sol_graph;
	}
	bool wait_4_all = false;
	if( wait_4_all)
	{
	  printf("PARENT DONE. Waiting for child\n");
	  while (true) {
	    int status;
	    pid_t done = wait(&status);
	    if (done == -1) {
	      //if (errno == ECHLD) break; // no more child processes
	      break;
	    } else {
	      if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		cerr << "pid " << done << " failed" << endl;
		exit(1);
	      }
	    }
	  }
	}
	else
	{
	  printf("PARENT DONE. Killing child process..\n");
	  /// Kill'em all!!
	  kill(pID, SIGTERM);
	}
      }
    }
    else if( prob.solution_strategy == "multimth_cplex" ||
	     prob.solution_strategy == "multimth_ga"  )
    {
      Matheuristic *mth = new Matheuristic();
#define TEST_MULTIMTH 0
#if TEST_MULTIMTH
      while(true)
      {
	printf("Press any key to send test data\n");
	getchar();
	mth->testPush(rand());
	printf("Press any key to receive test data\n");
	getchar();
	double testval;
	int ret;
	while( ret = mth->testPull(&testval) )
	{
	  printf("ret = %d Received test val %d\n", ret, testval);
	}
      }
#endif

      string mth_outfile = prob.output_path + 
	prob.instance_id + ".mthout";
      mth->setOutput(mth_outfile);
      bool solved = false;
      if( prob.solution_strategy == "multimth_cplex" )
      {
	MyCplexSolver *cplex_solver = new MyCplexSolver(static_net, net);

	cplex_solver->setMatheuristic(mth);

	solved = cplex_solver->solve();
	RNPSolutionPtr sol = cplex_solver->solution();
	if(solved)
	{
	  RNPSolutionPtr sol = cplex_solver->solution();

	  if( cplex_solver->solution()->isOptimal())
	    printf("Optimal value: ");
	  else
	    printf("SubOptimal value: ");
	  printf("%f\n",cplex_solver->getObjVal());

	  if(prob.output_sol)
	  {
	    sol->writeToFile(prob.output_path + prob.output_sol_file);
	  }

	  sol_graph = sol->create_graph_from_sol(net);

	  if(prob.output_sim){
	    generate_simulation(sol_graph, net, "sol");
	  }
	  output_graphs(sol_graph,"solution");
	  delete sol_graph;
	}
      }
      else if( prob.solution_strategy == "multimth_ga" )
      {
	GASol *ga = new GASol(static_net, net);
	ga->setMatheuristic(mth);
	double ga_sol = ga->solve();
      }
    }
  }

  delete static_net;
  delete net;


  DEBUG(1)
    close_debug_file();
  return 0;

}
