/*
 * =====================================================================================
 *
 *       Filename:  generator.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/01/2010 01:08:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _GENERATOR
#define _GENERATOR

/** Minimun node distance used when generating random locations 
 *  (set in meters)*/
#define MIN_NODE_DISTANCE 0.6


enum nodeType { BASE, STATIC, RELAY};
static string nodeTypeId[3] = {"static", "base", "relay"};


typedef struct clusters{
  int n_x;
  int n_y;
  double *density;
} clusters_t;

typedef struct grb_params{
  int cuts;
} grb_params_t;

//typedef struct cplex_params{
//char * name;
//int cpx_index;
//char type;
//double db_value;
//int int_value;
//struct cplex_params * next;

//} cplex_params_t;

typedef struct  
{
  double mipGap;
  int trelim;
  int logToFile;
  int workmem;
  int nodefileind;
  int mipemphasis;
  int threads;
  int lpmethod;
  int preind;
  int parallelmode;
  int mipdisplay;
  int probe;
  bool save_cpx_params; /// bool: save cplex params
} mycplex_params_t;


typedef struct prr_metric_params
{
  double ttx;
  double ple;
  double demand;
  double sigma;
}prr_metric_params_t;

typedef struct
{
  int   	seed;
  int           genome_size;     /// if (genome_size <= 0) => equal to prob.K
  int   	populationSize;
  float 	pCrossover;
  float 	pUniformCross;
  float 	pMutation;
  int   	nGenerations;
  float 	pMutSizeChange;
  float 	pMutSizeIncrease;
  float 	pUniformMutate;
  float 	pEnhancedMutationRandomRelay;
  float 	adjacent_txrange;
  int   	verbose;
  int    	debug;
  string 	debug_file;
  string 	algorithm;
  string 	ga_selection_scheme;
  string 	ga_scaling_scheme;
  double 	report_score_interval;
  double 	time_limit;
  bool 		scale_obj;
  double 	obj_scale_factor;
  string 	output_file;
  bool 		log_to_file;
  string 	log_file;
  bool 		log_population;
  string 	log_population_file;
  bool 		do_local_search;
  int 		local_search_time;

  /// whether to do the updates after the evaluation
  /// if updates are not done but there is an attempt to 
  /// use the feature in the crossover (e.g., do_update_conflicmap = 0 AND use_conflicts=1)
  /// there could be errors
  bool 		do_update_candidate;
  bool 		do_update_partner;
  bool 		do_update_conflictmap;
  bool 		do_update_chain;

  bool 		fix_rnpsol;


  int 		candidate_set_size;

  /// version 1
  double 	pCandidateRelay; //< deprecated, -> use_candidate_relays
  bool 		candidate_neighbor_check;
  int 		max_candidate_neighbor;

  /// version 2
  int 		candidate_region_x;
  int 		candidate_region_y;


  /// whether to use conflicts and candidates in the crossover
  /// still, the conflicts and candidates might be updated after each evaluation
  /// to completely disable these features, set do_update_candidates do_update_conflicts
  bool use_conflicts;
  bool use_candidate_relays;
  bool use_partners;
  bool use_chains;

  double pPartnerRelay;
  double pUseExpat;
  double pAdoption;
  double pUseChain;  /// overriden if use_chains == false

  bool log_generation;
  string log_generation_file;
  string which_solver;
  bool always_accept_expat;
  int max_evaluation_time;
  bool evaluate_use_initial_k0;

} ga_solver_params_t;

typedef struct
{
  string comm;
  string socket_address;
  int socket_port;
  int time_limit;
} matheuristic_params_t;

typedef struct
{
  bool flow_relay_constraints;
  bool integer_flow;
  double max_capacity;
  double max_local_flow;
  double max_flow;
  int max_degree;
  double default_demand;
  double flow_neighbor_penalty_w;
  double penalty_relay;
  bool use_node_degree_constraints;
  bool use_flow_neighbor_constraints;
  //! improved version that involves RNs 26.03.2015
  bool use_im_flow_neighbor_constraints;
  int im_flow_neighbor_nh;
  //! use interference range
  bool im_flow_neighbor_use_interference_range;
  double im_flow_neighbor_interference_range;
  //! strict flow neighbor constraints
  bool strict_flow_neighbor_constraints;
  //
  bool use_flow_neighbor_sink;
  //!
  double lost_flow_penalty;
  //!
  bool allow_incomplete_delivery;
  double incomplete_delivery_penalty;
} ilp_model_params;

typedef struct
{
  string output_path; // Directory where to put output files (. )
  string input_mps_file;
  string output_mps_file;
  int use_mps;
  int save_mps;
  string model_file; // LP Model filename (model.mod)
  string data_file; // LP Data filename (model.dat)
  string solver;
  int get_suboptimal;
  int verbose;
  int keep_datafile; // Dont delete model data file
  bool do_postsolve;
  bool allow_parallel;
  bool noMIP;
  bool force_branch_and_cut;
  int weak_timelim;
  int strong_timelim;
  bool log_progress;
  bool use_initial_k0;
} lp_solver_params_t;



/// Problem parameters (default value)

typedef struct
{
  int st_nodes; 		//! Number of static nodes (10)
  int b_nodes; 			//! Number of base station nodes (2)
  int K; 			//! Maximum number of relays to use (+INF)
  int minK; 			//! Minimum number of relays (-1)
  bool fixed_relays; 		//! Use determined set of relays
  vector<string> relays;  	//! predetermined set of relays
  double dimX, dimY; 		//! Width and length of field (10 x 10)
  double tx_range; 		//! Transmission range (1.0)
  bool no_position; 		//! do not consider node positions
  bool fixed_links; 		//! links are fixed and given as input

  string demand_type;
  string demand_file;
  unsigned int seed; //Random seed (time(0))

  double grid_res; //Grid resolution (1.0)

  int verbose; // Verbose level (1)


  bool gen_net; // Generate net files (true)	
  bool gen_dot; // Generate dot files: graphviz (true)	
  bool gen_graphs; // Generate graphs using graphviz (false)	
  int debug; // Debug level (0)
  string instance_id;  // id to identify problem instance (default)
  //     string instance_id_suffix; // string to append to instance_id ("")
  int glpk_out; //Display glpk output (0) 
  int grb_out; //Display gurobi output (0)
  string output_path; // Directory where to put output files (. )

  clusters_t clu;

  lp_solver_params_t lp_params;
  grb_params_t grb;

  //cplex_params *cplex_par;
  mycplex_params_t mycplex_params;

  ga_solver_params_t ga_params;

  ilp_model_params ilp_model;
  matheuristic_params_t matheuristic_params;


  int use_clusters;
  double gridh2_minr; // Grid h2 minimum delta
  double gridh2_maxr; // Grid h2 maximum delta
  string grid_strategy; // Grid strategy to use (default)
  string solution_strategy;
  double iter_res_init;
  int iter_res_niter;
  int use_convex_hull;
  string debug_prefix;
  int use_heuristic;
  string heuristic;

  string instance_path;
  string network_file;

  int output_sol;
  string output_sol_file;
  string output_time_file;

  string output_conn_file;

  int Low_res_factor;
  int Low_res_use_file;
  string Low_res_solution_file;


  int generate_instance;
  string generation_path;
  int solve;

  int output_sim; // Generate sim files (default 1)

  double lambda; // coefficient for metric cost
  string metric; // Which cost metric to use (prr)
  prr_metric_params_t prr_par; // PRR Metric parameters

  double uniform_random_demand_low;
  double uniform_random_demand_high;


} param;
/*
  void generate_grid(Network *net, double dimX, double dimY, double r );
  void generate_random_nodes(Network *net, int num, nodeType t);
  GetPot read_config_file(char *config_file);
  void load_defaults();
  void load_config_file(GetPot infile);
  */





#endif
