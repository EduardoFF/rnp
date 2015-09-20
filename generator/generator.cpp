/*
 * =====================================================================================
 *
 *       Filename:  generator.cpp
 *
 *    Description:  Instance generator
 *
 *        Version:  1.0
 *        Created:  05/31/2010 05:26:46 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eduardo Feo
 *        Company:  
 *
 * =====================================================================================
 */

#include "main.h"
#include "Plotter.h"
#include "Random.h"
#include "generator.h"
#include "graph.h"
#include "network.h"
#include "GetPot.h"

param prob;
//param def; //Default values for parameters
//ofstream cdebug;

#if WITH_OLD_BOOST_GEOMETRY
using namespace boost::geometry;
bool lineIntersection(boost::geometry::linestring_2d &l1,
		      boost::geometry::linestring_2d &l2,
		      double *X,
		      double *Y);
bool pointInBox(point_2d &p,
		point_2d &p1,
		point_2d &p2);
#endif

void init_debug_file(string filename){
		cdebug.open(filename.c_str());
		cdebug << "********* START OF DEBUG FILE *************" << endl;

}

void close_debug_file(){
	cdebug.close();
}



void generate_simulation(Graph *G, Network *net, string filename){

	string sim_config_file 	= prob.output_path + "simconfig-" + prob.instance_id + "-" + filename + ".txt";
	string net_file 	= prob.output_path + "sim-" + prob.instance_id + "-" + filename + ".net";

	
	ofstream netFile (net_file.c_str());

	ofstream simconfig (sim_config_file.c_str());

	simconfig << "network_file\t" << net_file << "\n";
	simconfig << "sim_time\t" << "60" << "\n";
	simconfig << "area_size\t" << net->dimX << "\t" << net->dimY << "\n";
	simconfig << "random_seed\t" << prob.seed << "\n";

	simconfig.close();
	G->toNetFile(netFile);
/*  Not necesary - possible links added in sim script 	
	for(int i=0;i<net->size();i++){
		for(int j=i+1;j<net->size();j++){
			Node &n = net->getNode(i);
			Node &m = net->getNode(j);
			if(net->inRange(n,m)){
				netFile << "link\t" << i << "\t" << j << endl;
				netFile << "link\t" << j << "\t" << i << endl;
				
			}
		}
	}
	*/
	netFile.close();

}


/**
 * @brief Generates grid points with indicated resolution 
 * @param  net network to superpose the relay grid
 * @param  dimX  area length
 * @param  dimY  area width  
 * @param  r resolution to use
 * @return 
 */

  void
generate_grid(Network *net, double dimX, double dimY, double r )
{
  double x, y;
  int cnt =0;
  x = 0;
  DEBUG(2) cdebug << "Generating relays..." << endl;
  int rx = (int)floor(dimX/r) + 1;

  int ry = (int)floor(dimX/r) + 1;

  DEBUG(2) cdebug << "X-Res " << rx << ", Y-Res " << ry << endl;

  for(int i=0;i<rx;i++){
    y = 0;
    for(int j =0;j<ry;j++){
      if(net->addNode(x,y,RELAY,""))
	DEBUG(2) cdebug << "Relay at " << x << ", " << y << endl;
      y += r;
    }
    x += r;
  }

}

  void 
generate_grid(Network *net, double dimX, double offset_x, 
	      double dimY, double offset_y, double r )
{

  DEBUG(3) cdebug << "griding dimX " << dimX << " dimY " << dimY << " offset_x " << offset_x << " offset_y " << offset_y << " r " << r << endl;
  double x, y;

  x = 0;
  int rx = (int)floor(dimX/r) + 1;
  int ry = (int)floor(dimX/r) + 1;
  for(int i=0;i<rx;i++){
    y = 0;
    for(int j=0;j<ry;j++){
      net->addNode(x+offset_x,y+offset_y,RELAY,"");
      y += r;
    }
    x += r;
  }

}


#if WITH_GRID_H2
// depends on old booost geometry
typedef boost::geometry::point_2d point_t;
typedef boost::geometry::linestring_2d segment_t;
//TODO Parametrize these values
#define GRID_H2_X 5
#define GRID_H2_Y 5
#define GRID_H2_RANGE 4
typedef struct gridh2_region
{
  double offset_x;
  double offset_y;
  double density;
  double avg_neighbor_density;
  double res;

} gridh2_region_t;

  void 
generate_grid_h2(Network *net, double dimX, double dimY, double min_res, double max_res)
{
  /**  Divide region in sub-areas and calculate node density
   *  for each region */

  DEBUG(3) cdebug << "grid_h2 min_res " << min_res << " max_res " << max_res << endl;
  vector< point_t> nodes;
  gridh2_region_t regions[GRID_H2_X][GRID_H2_Y];
  double off_x,off_y;
  double size_x, size_y;
  double resolutions[GRID_H2_RANGE];

  for(int i=0;i<GRID_H2_RANGE-1;i++)
  {
    resolutions[i] = min_res + ((max_res-min_res)/(GRID_H2_RANGE-1))*i;
    DEBUG(3) cdebug << "res " << i << " = " << resolutions[i] << endl;
  }
  resolutions[GRID_H2_RANGE-1] = max_res;

  /** Get node locations as points */
  for(int i=0;i<net->size();i++)
  {
    Node &n = net->getNode(i);
    nodes.push_back(boost::geometry::make<point_t>(n.x,n.y));
  }

  size_x = dimX/GRID_H2_X;
  size_y = dimY/GRID_H2_Y;

  double standard_density = 1.0/(GRID_H2_X * GRID_H2_Y);
  DEBUG(2) cdebug << "standard density " << standard_density << endl;


  off_x = 0.0;
  for(int i=0;i<GRID_H2_X;i++){
    off_y = 0.0;
    for(int j=0;j<GRID_H2_Y;j++){
      regions[i][j].offset_x = off_x;
      regions[i][j].offset_y = off_y;

      point_2d p(off_x,off_y);
      point_2d q(off_x+size_x,off_y+size_y);

      DEBUG(3) cdebug << "region " 
	<< i << j << ": " << "(" << off_x 
	<< " " << off_y << ")"	<< endl;

      /** calculate density */
      int count=0;
      for(int k=0;k<nodes.size();k++)
	if( pointInBox(nodes[k],p,q))
	  count++;
      regions[i][j].density = (1.0*count)/(1.0*nodes.size());
      DEBUG(3) cdebug << "region " 
	<< i << j << ": (" << count 
	<< ") " << regions[i][j].density << endl;
      off_y += size_y;
    }
    off_x+= size_x;
  }

  /** calculate min_neighbor_density */
  int dir[8][2] = {{-1,-1},
    {-1, 0},
    {-1,1},
    { 0, -1},
    { 0, 1},
    { 1,-1},
    { 1, 0},
    { 1, 1}};
  for(int i=1;i < GRID_H2_X-1;i++)
  {
    for(int j = 1; j < GRID_H2_Y-1;j++)
    {
      double avg_d = 0;
      for(int k=0;k<8;k++){
	int ii = i+dir[k][0];
	int jj = j+dir[k][1];
	avg_d +=regions[ii][jj].density;
      }
      regions[i][j].avg_neighbor_density = avg_d/8;;
    }
  }

  /** 
   * Assign maximum delta(resolution) to border and 
   * regions surrounded by empty regions 
   * */

  for(int j=0;j<GRID_H2_Y;j++)
  {
    regions[0][j].res = max_res;
    regions[GRID_H2_X-1][j].res = max_res;
  }
  for(int i=0;i<GRID_H2_X;i++)
  {
    regions[i][0].res = max_res;
    regions[i][GRID_H2_Y-1].res = max_res;
  }

  /** For other regions - calculate resolution based on
   *  average resolution of region and surroundings 
   *  a standard density is calculated and region has more or less
   *  resolution based on the distance to the standard
   *  */
  for(int i=1;i<GRID_H2_X-1;i++)
  {
    for(int j=1;j<GRID_H2_Y-1;j++)
    {
      //double avg_density = 
      //(regions[i][j].avg_neighbor_density + regions[i][j].density)/2.0;
      double avg_density = regions[i][j].density;

      //if(avg_density == 1.0) avg_density-=0.001;
      int rang = 
	fabs(standard_density - avg_density - 0.0001)/(standard_density/GRID_H2_RANGE);
      int rang_ix = min(rang,GRID_H2_RANGE-1);
      if( rang_ix < 0 || rang_ix > GRID_H2_RANGE-1)
      {
	fprintf(stderr, "WTF! error in grid_h2\n");
      }
      regions[i][j].res = resolutions[rang_ix];
      DEBUG(3) cdebug << "region " << i << j << " dens " 
	<< regions[i][j].density << " rank " << rang 
	<< " res " << regions[i][j].res << endl;
    }
  }

  /** Generate nodes */
  for(int i=0;i<GRID_H2_X;i++)
  {
    for(int j=0;j<GRID_H2_Y;j++)
    {
      generate_grid(net, size_x, 
		    regions[i][j].offset_x, 
		    size_y, regions[i][j].offset_y, 
		    regions[i][j].res);
    }
  }
}


// depends on old boost geometry located in utils
void 
generate_grid_h1(Network *net, double dimX, double dimY)
{
  vector< point_t> grid;
  vector< segment_t> sgs;
  for(int i=0;i<net->size();i++){
    Node &n = net->getNode(i);
    grid.push_back(boost::geometry::make<point_t>(n.x,n.y));
  }

  /** Create horizontal segments */
  for(int i=0;i<grid.size();i++)
  {
    point_t p1,p2;
    p1 = boost::geometry::make<point_t>(0.0,grid[i].y());
    p2 = boost::geometry::make<point_t>(dimX, grid[i].y());

    segment_t seg;
    seg.push_back(p1);
    seg.push_back(p2);
    sgs.push_back( seg);
  }
  /** Create horizontal segments */
  for(int i=0;i<grid.size();i++)
  {
    point_t p1,p2;
    p1 = boost::geometry::make<point_t>(grid[i].x(),0.0);
    p2 = boost::geometry::make<point_t>(grid[i].x(),dimY);

    segment_t seg;
    seg.push_back(p1);
    seg.push_back(p2);
    sgs.push_back( seg);
  }

  for(int i=0;i<sgs.size();i++)
  {
    for(int j=i+1;j<sgs.size();j++)
    {
      //std::vector<point_t > intersection;
      double x,y;
      if(lineIntersection(sgs[i],sgs[j],&x,&y))
      {
	net->addNode(x,y,RELAY,"");
      }
      //    				 boost::geometry::intersection_inserter<point_t >(sgs[i], sgs[j], 
      //										  std::back_inserter(intersection));
      //	 cout << "Intersection size = " << intersection.size() << endl;

    }
  }

}
#endif

/**
 * @brief  generates random locations and add nodes in these locations to the network
 *
 *         To handle minimum inter-node distance, the function net::addNode returns false
 *         if the location where we want to locate the node violates the minimum distance
 * @param  net network to add the nodes
 * @param  gen random number generator - which gives a random value between [0,1)
 * @param  dimX area length
 * @param  dimY area width
 * @param  offset_x offset in x-axis
 * @param  offset_y offset in y_axis
 * @param  num  number of nodes to add
 * @param  t    type of nodes to add 
 * @return 
 */
void 
generate_random_nodes(Network *net, RandomGenerator *gen, double dimX, double dimY, 
		      double offset_x, double offset_y, int num, nodeType t){
  VERBOSE(3){ 
    printf("generating %d %s nodes...\n",num, nodeTypeId[t].c_str());
  }
  int i = 0;  /* Number of added nodes */
  while(i < num){
    double x, y;
    x = (*gen)() * dimX + offset_x;
    y = (*gen)() * dimY + offset_y;
    VERBOSE(3){
      cout << x << " " << y << endl;
    }
    /** Retry until we have found an appropiate location */
    if(net->addNode(x,y,t,""))
      i++;
  }
}


/**
 * @brief  generates random locations based on clusters settings. Add nodes to the network in these locations
 *
 *         This function presents a particularly situation when considering a minimum inter-node distance.
 *         Fortunately this issue is handled inside the network class, which does not allow a node to be added
 *         in a location which violates the minimum distance
 *
 * @param  net: Network to add the nodes
 * @param  gen: Random number generator
 * @param  dimX: Area width in meters
 * @param  dimY: Area height in meters
 * @param  num_static: Number of static nodes to generate
 * @param  num_base:   Number of sink nodes to generate
 * @param  clu: Cluster settings (see clusters_t structure)
 * @return 
 */

/**
 * @brief 
 * 
 * Detailed description starts here.
 */
  void 
generate_cluster_nodes(Network *net, RandomGenerator *gen, double dimX, double dimY,
		       int num_static, int num_base, clusters_t &clu)
{

  int num_clusters = clu.n_x * clu.n_y;
  double offset_x=0.0, offset_y=0.0;
  double size_x = dimX/clu.n_x;
  double size_y = dimY/clu.n_y;

  int clu_static[clu.n_x][clu.n_y];

  int clu_base[clu.n_x][clu.n_y];

  memset(clu_static, 0,sizeof(clu_static));
  memset(clu_base, 0,sizeof(clu_base));


  /// Generate static nodes 

  /// Calculate number of nodes to be placed in each cluster
  double ds_x=0.0;
  for(int i=0;i<clu.n_x;i++)
    for(int j=0;j<clu.n_y;j++)
      ds_x += clu.density[i*clu.n_y + j];
  /// re-scale densities and set number of nodes
  int c_static = 0;
  int ns=0;
  for(int i=0;i<clu.n_x;i++)
  {
    for(int j=0;j<clu.n_y;j++)
    {
      if(i==clu.n_x-1 && j == clu.n_y-1)
      {
	/** place the remaining nodes in the last cluster */
	clu_static[clu.n_x-1][clu.n_y-1] = num_static - c_static;
	DEBUG(2) cdebug << "placing " << (num_static - c_static) 
	  << " in cluster " << i << " " << j << endl;

      }
      else
      {
	clu.density[i*clu.n_y + j] = clu.density[i*clu.n_y + j]/ds_x;
	if(clu.density[i*clu.n_y + j] == 0)
	{
	  ns = 0;
	}
	else
	{

	  ns = num_static*clu.density[i*clu.n_y + j];
	  c_static+= ns;
	}
	clu_static[i][j] = ns;
	DEBUG(2) cdebug << "placing " << ns 
	  << " in cluster " << i << " " << j << endl;
      }
    }
  }
  offset_x = 0.0;
  for(int i=0;i < clu.n_x;i++){
    offset_y = 0.0;
    for(int j=0;j < clu.n_y;j++){
      generate_random_nodes(net,gen,size_x,size_y,
			    offset_x,offset_y,clu_static[i][j],STATIC);
      offset_y += size_y;
    }
    offset_x += size_x;
  }


  /// Generate base nodes
  for(int b =0; b< num_base;b++){
    /** Select randomly where to place the node 
     *  based on cluster densities      
     *  */

    double r = (*gen)(); // Get random number between 0,1
    double p = 0.0;
    bool found = false;
    offset_x = 0.0;
    // Find where to place it
    for(int i=0;i<clu.n_x && !found;i++){
      offset_y = 0.0;
      for(int j=0;j<clu.n_y && !found;j++){
	if(p <= r && r <= p + clu.density[i*clu.n_y + j]){
	  generate_random_nodes(net,gen,size_x,size_y,
				offset_x,offset_y,1,BASE);
	  found = true;
	}
	offset_y += size_y;
	p+= clu.density[i*clu.n_y + j];
      }
      offset_x += size_x;

    }

  }



}


GetPot read_config_file(char *config_file)
{
  /// Define the syntax for the input configuration file and start reading
  string base_directory  = "./";  
  string comment_start   = "#";  
  string comment_end     = "\n";
  string field_separator = ",";

  string input_file    = config_file;

  GetPot infile((input_file).c_str(), 
		       comment_start.c_str(), comment_end.c_str(), field_separator.c_str()); 
  return (infile);
}


void load_defaults(param &def)
{
  def.st_nodes = 10;
  def.b_nodes = 2;
  def.K = 1000000;
  def.minK = 0;
  def.tx_range = 1.0;
  def.grid_res = 1.0;
  def.dimX = 10.0;
  def.dimY = 10.0;
  def.seed = time(0);
  def.verbose = 1;
  def.debug = 0;
  def.instance_id = "default";
  //	def.instance_id_suffix = "";
  def.gen_graphs = false;
  def.gen_dot = true;
  def.gen_net = true;
  def.no_position = false;
  def.fixed_links = false;
  def.demand_type = "uniform";
  def.demand_file = "default";
  def.uniform_random_demand_low = 1.0;
  def.uniform_random_demand_high =4.0;




  def.glpk_out = 0;
  def.grb_out = 0;
  def.grb.cuts = -1;
  def.use_clusters = 0;
  def.clu.n_x = 0;
  def.clu.n_y = 0;
  def.gridh2_minr = 1.0;
  def.gridh2_maxr = 1.0;
  def.grid_strategy = "default";
  def.solution_strategy = "default";
  def.iter_res_init = 1.0;
  def.iter_res_niter = 1;
  def.use_convex_hull = 1;
  def.debug_prefix = "debug.";
  def.use_heuristic = 0;
  def.heuristic = "no-heuristic";

  def.Low_res_factor = 1;
  def.Low_res_use_file = 0;
  def.Low_res_solution_file = "default";

  def.network_file = "default";
  def.instance_path = "./";


  def.solve = true;
  def.generate_instance = false;
  def.generation_path = "./";
  def.output_sol_file = "default";
  def.output_conn_file = "default";
  def.output_time_file = "default";
  def.output_sol = 0;

  def.metric = "PRR";
  def.output_sim = 1;

  def.lambda = 1.0;
  def.prr_par.ttx = 4.486; // 128 bytes message
  def.prr_par.ple = 3.0;
  def.prr_par.sigma = 3.2;
  def.prr_par.demand = 0.012; // 12 packets per second

  def.lp_params.input_mps_file = "default";
  def.lp_params.output_mps_file = "default";
  def.lp_params.use_mps = 0;
  def.lp_params.model_file = "./model.mod";
  def.lp_params.data_file = "default";	
  def.lp_params.save_mps = false;
  def.lp_params.output_path = "./";
  def.lp_params.solver = "CPLEX";
  def.lp_params.keep_datafile = 0;
  def.lp_params.get_suboptimal = 1;
  def.lp_params.verbose = 0;
  def.lp_params.do_postsolve = 1;
  def.lp_params.allow_parallel = 1;
  def.lp_params.noMIP = 0;
  def.lp_params.weak_timelim = 0; 
  def.lp_params.strong_timelim = 7200;
  def.lp_params.force_branch_and_cut = 0;
  def.lp_params.log_progress = 1;
  def.lp_params.use_initial_k0 = 0;

  def.ga_params.seed = 12345;
  def.ga_params.genome_size = 0;
  def.ga_params.populationSize = 10;
  def.ga_params.pCrossover = 0.9;
  def.ga_params.pUniformCross = 0.2;

  def.ga_params.adjacent_txrange = 2.0;
  def.ga_params.pMutation = 0.1;
  def.ga_params.pMutSizeChange = 0.5;
  def.ga_params.pMutSizeIncrease = 0.5;
  def.ga_params.pUniformMutate = 0.0;
  def.ga_params.pEnhancedMutationRandomRelay = 0.5;
  def.ga_params.candidate_neighbor_check = 0;
  def.ga_params.max_candidate_neighbor = 1;

  def.ga_params.nGenerations = 1;
  def.ga_params.verbose = 0;
  def.ga_params.ga_selection_scheme = "GATournamentSelector";
  def.ga_params.ga_scaling_scheme = "GALinearScaling";
  def.ga_params.algorithm = "GASteadyState";
  def.ga_params.report_score_interval = 30.0; // Every 1/2 minute
  def.ga_params.time_limit = 300; // 5 minutes
  def.ga_params.scale_obj = false;
  def.ga_params.obj_scale_factor = 1.0;
  def.ga_params.output_file = "default";
  def.ga_params.log_to_file = 0;
  def.ga_params.log_file = "default";
  def.ga_params.log_population = 0;
  def.ga_params.log_population_file = "default";

  def.ga_params.log_generation = 0;
  def.ga_params.log_generation_file = "default";

  def.ga_params.fix_rnpsol = 1;
  def.ga_params.do_update_candidate = 1;
  def.ga_params.do_update_partner = 1;
  def.ga_params.do_update_conflictmap = 1;
  def.ga_params.do_update_chain = 1;

  def.ga_params.do_local_search = 1;
  def.ga_params.local_search_time = 60;
  def.ga_params.debug = 0;
  def.ga_params.debug_file = "default";

  def.ga_params.candidate_set_size = 25;
  def.ga_params.pCandidateRelay = 0.25;

  def.ga_params.candidate_region_x = 4;
  def.ga_params.candidate_region_y = 4;

  def.ga_params.pPartnerRelay = 0.25;
  def.ga_params.pAdoption = 0.1;
  def.ga_params.pUseChain = 0.25;

  def.ga_params.which_solver = "cplexsolver";
  def.ga_params.pUseExpat = 1.0;
  def.ga_params.always_accept_expat = 1;
  def.ga_params.max_evaluation_time = 120; // 2 minutes

  def.ga_params.use_conflicts = 1;
  def.ga_params.use_candidate_relays = 1;
  def.ga_params.use_partners = 1;
  def.ga_params.use_chains = 1;

  def.ga_params.evaluate_use_initial_k0 = 0;

  def.ilp_model.flow_relay_constraints = 1;
  def.ilp_model.max_capacity = 70.0;


  def.ilp_model.use_flow_neighbor_constraints = 1;
  def.ilp_model.use_im_flow_neighbor_constraints = 0;
  
  def.ilp_model.use_node_degree_constraints = 1;

  def.ilp_model.max_local_flow = def.ilp_model.max_capacity;
  def.ilp_model.max_flow = 10000.0;
  def.ilp_model.max_degree = 6;
  def.ilp_model.default_demand = 1.0;

  def.ilp_model.integer_flow = false;

  def.ilp_model.flow_neighbor_penalty_w = 0.1;
  def.ilp_model.penalty_relay = 1.0;
  def.ilp_model.im_flow_neighbor_nh = 2;
  def.ilp_model.im_flow_neighbor_use_interference_range = true;
  def.ilp_model.im_flow_neighbor_interference_range = 2*def.tx_range;

  def.ilp_model.allow_incomplete_delivery = false;
  def.ilp_model.incomplete_delivery_penalty = 1.0;

  def.ilp_model.use_flow_neighbor_sink = false;
  def.ilp_model.lost_flow_penalty = 0.1;

  def.mycplex_params.mipGap = 10e-6;
  def.mycplex_params.trelim = 2048;
  def.mycplex_params.workmem = 128;
  def.mycplex_params.nodefileind = 1;
  def.mycplex_params.logToFile = 0;
  def.mycplex_params.mipemphasis = 0;
  def.mycplex_params.threads = 0;
  def.mycplex_params.lpmethod = 0;
  def.mycplex_params.preind = 1;
  def.mycplex_params.parallelmode = -1;
  def.mycplex_params.mipdisplay = 2;
  def.mycplex_params.save_cpx_params = false;

  def.matheuristic_params.comm = "pipes";
  def.matheuristic_params.time_limit = 7200;
  def.matheuristic_params.socket_address = "127.0.0.1";
  def.matheuristic_params.socket_port = 12345;
}

string fixPath(string path)
{
  string ret_path = path;
  if(path[path.size()-1] != '/')
  {
    ret_path += "/";
  }
  return ret_path;
}

void 
load_config_file(GetPot infile, param &def)
{
  prob.st_nodes = infile("ProblemInstance/static",def.st_nodes);
  prob.b_nodes = infile("ProblemInstance/base",def.b_nodes);
  prob.K = infile("ProblemInstance/max_relays",def.K);

  prob.minK = infile("ProblemInstance/min_relays",def.minK);

  prob.tx_range = infile("ProblemInstance/tx_range",def.tx_range);

  prob.grid_res = infile("ProblemInstance/grid_res",def.grid_res);

  prob.no_position = infile("ProblemInstance/no_position", def.no_position);
  prob.fixed_links = infile("ProblemInstance/fixed_links", 
			    def.fixed_links);

  /// the behaviour should be:
  /// (a) if fixed_relays is true in the current config file
  ///     then, load the relays
  /// (b) else, leave everyting as before
  bool fixed_relays = infile("ProblemInstance/fixed_relays",false);
  if( fixed_relays)
  {
    prob.fixed_relays = fixed_relays;
    if(prob.relays.size())
    {
      prob.relays.clear();
    }
    int nrelays = infile.vector_variable_size("ProblemInstance/relays");
    for(int i=0;i<nrelays;i++){
      string rid = infile("ProblemInstance/relays","no-relay",i);
      prob.relays.push_back(rid);
    }
  }

  prob.demand_type = infile("ProblemInstance/demand_type", def.demand_type.c_str());
  prob.demand_file = infile("ProblemInstance/demand_file",def.demand_file.c_str());
  prob.uniform_random_demand_low = infile("ProblemInstace/uniform_random_demand_low", def.uniform_random_demand_low);
  prob.uniform_random_demand_high = infile("ProblemInstace/uniform_random_demand_high", def.uniform_random_demand_high);




  prob.dimX = infile("ProblemInstance/X",def.dimX);
  prob.dimY = infile("ProblemInstance/Y",def.dimY);
  prob.seed = infile("Parameters/seed",(int)def.seed);
  prob.verbose = infile("Parameters/verbose",def.verbose);
  prob.debug = infile("Parameters/debug_level",def.debug);

  prob.gen_dot = infile("Parameters/gen_dot",def.gen_dot);
  prob.gen_net = infile("Parameters/gen_net",def.gen_net);
  prob.gen_graphs = infile("Parameters/gen_graphs",def.gen_graphs);
  prob.instance_id = infile("ProblemInstance/instance_id",def.instance_id.c_str());
  //	prob.instance_id_suffix = infile("ProblemInstance/instance_id_suffix", "");
  //	if( prob.instance_id_suffix !=
  //	prob.instance_id = prob.instance_id + 
  prob.glpk_out = infile("Parameters/glpk_out",def.glpk_out);
  prob.grb_out = infile("Parameters/grb_out",def.grb_out);
  prob.output_path = infile("Parameters/output_path",def.output_path.c_str());

  // Avoiding problem when using several conf files
  // If we have already loaded clusters, use previous values as default
  // and do not allocate new memory
  // Remember, param def is actually a reference to prob variable
  bool def_has_clusters = def.use_clusters;
  int def_cluster_size = 0;
  if(def.use_clusters)
  {
    def_cluster_size = def.clu.n_x * def.clu.n_y;
  }


  prob.use_clusters =  infile("Clusters/use_clusters",def.use_clusters);
  prob.clu.n_x =  infile("Clusters/X",def.clu.n_x);
  prob.clu.n_y =  infile("Clusters/Y",def.clu.n_y);
  if(prob.use_clusters){
    if(!def_has_clusters)
    {
      // If default param had no clusters
      // we can allocate memory
      prob.clu.density = 
	(double*)malloc(prob.clu.n_x*prob.clu.n_y*sizeof(double));
    }else
    {
      // If had cluster, but differente size
      if(def_cluster_size != (prob.clu.n_x * prob.clu.n_y))
      {
	free(prob.clu.density);
	prob.clu.density = 
	  (double*)malloc(prob.clu.n_x*prob.clu.n_y*sizeof(double));
      }
    }

    int ndens = infile.vector_variable_size("Clusters/density");
    for(int i=0;i<ndens;i++){
      prob.clu.density[i] = infile("Clusters/density",1.0,i);
    }
  }
  prob.gridh2_minr =  infile("Grid_H2/min_res",def.gridh2_minr);
  prob.gridh2_maxr =  infile("Grid_H2/max_res",def.gridh2_maxr);
  prob.grid_strategy = infile("Parameters/grid_strategy",def.grid_strategy.c_str());

  prob.solution_strategy = infile("Parameters/solution_strategy",def.solution_strategy.c_str());

  prob.iter_res_init = infile("Parameters/iter_res_init",def.iter_res_init);

  prob.iter_res_niter = infile("Parameters/iter_res_niter",def.iter_res_niter);
  prob.use_convex_hull = infile("Parameters/use_convex_hull",def.use_convex_hull);
  prob.debug_prefix = infile("Parameters/debug_prefix",def.debug_prefix.c_str());

  prob.use_heuristic = infile("Parameters/use_heuristic",def.use_heuristic);
  prob.heuristic = infile("Parameters/heuristic",def.heuristic.c_str());

  /// LP Parameters
  ////TODO This is a bit confusing, which parameters are LP, CPLEX, ILP
  prob.lp_params.output_mps_file = infile("LP/output_mps_file",def.lp_params.output_mps_file.c_str());
  prob.lp_params.input_mps_file = infile("LP/input_mps_file",def.lp_params.input_mps_file.c_str());
  prob.lp_params.model_file = infile("LP/model",def.lp_params.model_file.c_str());
  prob.lp_params.data_file = infile("LP/data",def.lp_params.data_file.c_str());
  prob.lp_params.use_mps = infile("LP/use_mps",def.lp_params.use_mps);
  prob.lp_params.keep_datafile = infile("LP/keep_datafile",def.lp_params.keep_datafile);
  prob.lp_params.save_mps = infile("LP/save_mps",def.lp_params.save_mps);

  prob.lp_params.output_path = infile("LP/output_path",def.lp_params.output_path.c_str());
  prob.lp_params.get_suboptimal = infile("LP/get_suboptimal",def.lp_params.get_suboptimal);
  prob.lp_params.verbose = infile("LP/verbose",def.lp_params.verbose);
  prob.lp_params.solver = infile("LP/solver",def.lp_params.solver.c_str());
  prob.lp_params.do_postsolve = 
    infile("LP/do_postsolve", def.lp_params.do_postsolve);
  prob.lp_params.allow_parallel = 
    infile("LP/allow_parallel", def.lp_params.allow_parallel);
  prob.lp_params.noMIP = 
    infile("LP/noMIP", def.lp_params.noMIP);
  prob.lp_params.strong_timelim = 
    infile("LP/strong_timelim", def.lp_params.strong_timelim);
  prob.lp_params.weak_timelim = 
    infile("LP/weak_timelim", def.lp_params.weak_timelim);
  prob.lp_params.log_progress = 
    infile("LP/log_progress", def.lp_params.log_progress);
  prob.lp_params.force_branch_and_cut = 
    infile("LP/force_branch_and_cut", 
	   def.lp_params.force_branch_and_cut);
  prob.lp_params.use_initial_k0 =  infile("LP/use_initial_k0", 
					  def.lp_params.use_initial_k0);

  prob.output_conn_file = infile("Parameters/output_conn_file", def.output_conn_file.c_str());
  prob.output_sol = infile("Parameters/output_sol", def.output_sol);
  prob.output_sol_file = infile("Parameters/output_sol_file", def.output_sol_file.c_str());

  prob.output_time_file = infile("Parameters/output_time_file", def.output_time_file.c_str());

  prob.output_sim = infile("Parameters/output_sim", def.output_sim);


  prob.instance_path = infile("ProblemInstance/instance_path",def.instance_path.c_str());
  prob.network_file = infile("ProblemInstance/network_file",def.network_file.c_str());


  prob.solve = infile("Main/solve",def.solve);
  prob.generate_instance = infile("Main/generate_instance",def.generate_instance);

  prob.generation_path = infile("Parameters/generation_path",def.generation_path.c_str());

  prob.metric = infile("Parameters/metric",def.metric.c_str());

  prob.Low_res_factor = infile("Low_res_heuristic/factor",def.Low_res_factor);
  prob.Low_res_use_file = infile("Low_res_heuristic/use_file",def.Low_res_use_file);
  prob.Low_res_solution_file = infile("Low_res_heuristic/solution_file",def.Low_res_solution_file.c_str());


  prob.grb.cuts = infile("Gurobi/cuts",def.grb.cuts);

  prob.Low_res_factor = infile("Low_res_heuristic/factor",def.Low_res_factor);

  prob.lambda = infile("Parameters/lambda",def.lambda);
  prob.prr_par.ttx = infile("PRR_METRIC/ttx",def.prr_par.ttx);
  prob.prr_par.ple = infile("PRR_METRIC/ple",def.prr_par.ple);
  prob.prr_par.sigma = infile("PRR_METRIC/sigma",def.prr_par.sigma);
  prob.prr_par.demand = infile("PRR_METRIC/demand",def.prr_par.demand);

  prob.ga_params.genome_size = infile("GA/genome_size",def.ga_params.genome_size);
  prob.ga_params.seed = infile("GA/seed",def.ga_params.seed);
  prob.ga_params.populationSize = infile("GA/populationSize",def.ga_params.populationSize);
  prob.ga_params.pCrossover = infile("GA/pCrossover",def.ga_params.pCrossover);
  prob.ga_params.pUniformCross = infile("GA/pUniformCross",def.ga_params.pUniformCross);

  prob.ga_params.pMutation = infile("GA/pMutation",def.ga_params.pMutation);
  prob.ga_params.adjacent_txrange = infile("GA/adjacent_txrange",def.ga_params.adjacent_txrange);
  prob.ga_params.pMutSizeChange = infile("GA/pMutSizeChange",def.ga_params.pMutSizeChange);
  prob.ga_params.pMutSizeIncrease = infile("GA/pMutSizeIncrease",def.ga_params.pMutSizeIncrease);
  prob.ga_params.pUniformMutate = infile("GA/pUniformMutate",def.ga_params.pUniformMutate);
  prob.ga_params.pEnhancedMutationRandomRelay = infile("GA/pEnhancedMutationRandomRelay",def.ga_params.pEnhancedMutationRandomRelay);

  prob.ga_params.nGenerations = infile("GA/nGenerations",def.ga_params.nGenerations);
  prob.ga_params.verbose = infile("GA/verbose",def.ga_params.verbose);

  prob.ga_params.algorithm = infile("GA/algorithm",def.ga_params.algorithm.c_str());
  prob.ga_params.ga_selection_scheme = infile("GA/GASelectionScheme",def.ga_params.ga_selection_scheme.c_str());
  prob.ga_params.ga_scaling_scheme = infile("GA/GAScalingScheme",def.ga_params.ga_scaling_scheme.c_str());

  prob.ga_params.time_limit = infile("GA/time_limit",def.ga_params.time_limit);
  prob.ga_params.report_score_interval = infile("GA/report_score_interval",def.ga_params.report_score_interval);
  prob.ga_params.scale_obj = infile("GA/scale_obj",def.ga_params.scale_obj);
  prob.ga_params.obj_scale_factor = infile("GA/obj_scale_factor",def.ga_params.obj_scale_factor);
  prob.ga_params.output_file = infile("GA/output_file", def.ga_params.output_file.c_str());
  prob.ga_params.log_to_file = infile("GA/log_to_file",def.ga_params.log_to_file);
  prob.ga_params.log_file = infile("GA/log_file",def.ga_params.log_file.c_str());
  prob.ga_params.debug_file = infile("GA/debug_file",def.ga_params.debug_file.c_str());

  prob.ga_params.log_population = infile("GA/log_population",def.ga_params.log_population);
  prob.ga_params.log_population_file = 
    infile("GA/log_population_file",def.ga_params.log_population_file.c_str());

  prob.ga_params.log_generation = infile("GA/log_generation",def.ga_params.log_generation);
  prob.ga_params.log_generation_file = 
    infile("GA/log_generation_file",def.ga_params.log_generation_file.c_str());


  prob.ga_params.fix_rnpsol = infile("GA/fix_rnpsol",def.ga_params.fix_rnpsol);
  prob.ga_params.do_update_candidate = 
    infile("GA/do_update_candidate",def.ga_params.do_update_candidate);
  prob.ga_params.do_update_partner = 
    infile("GA/do_update_partner",def.ga_params.do_update_partner);
  prob.ga_params.do_update_conflictmap = 
    infile("GA/do_update_conflictmap",def.ga_params.do_update_conflictmap);
  prob.ga_params.do_update_chain = 
    infile("GA/do_update_chain",def.ga_params.do_update_chain);

  prob.ga_params.use_conflicts = infile("GA/use_conflicts",def.ga_params.use_conflicts);
  prob.ga_params.use_candidate_relays = infile("GA/use_candidate_relays",def.ga_params.use_candidate_relays);
  prob.ga_params.use_partners = infile("GA/use_partners",def.ga_params.use_partners);
  prob.ga_params.use_chains = infile("GA/use_chains",def.ga_params.use_chains);

  prob.ga_params.do_local_search = infile("GA/do_local_search",def.ga_params.do_local_search);
  prob.ga_params.local_search_time = infile("GA/local_search_time",def.ga_params.local_search_time);
  prob.ga_params.candidate_set_size = infile("GA/candidate_set_size",def.ga_params.candidate_set_size);
  prob.ga_params.pCandidateRelay = infile("GA/pCandidateRelay",def.ga_params.pCandidateRelay);
  prob.ga_params.candidate_neighbor_check = infile("GA/candidate_neighbor_check",def.ga_params.candidate_neighbor_check);
  prob.ga_params.max_candidate_neighbor = infile("GA/max_candidate_neighbor",def.ga_params.max_candidate_neighbor);

  prob.ga_params.candidate_region_x = infile("GA/candidate_region_x",def.ga_params.candidate_region_x);
  prob.ga_params.candidate_region_y = infile("GA/candidate_region_y",def.ga_params.candidate_region_y);


  prob.ga_params.pPartnerRelay = infile("GA/pPartnerRelay",def.ga_params.pPartnerRelay);
  prob.ga_params.pAdoption = infile("GA/pAdoption",def.ga_params.pAdoption);
  prob.ga_params.pUseChain = infile("GA/pUseChain",def.ga_params.pUseChain);

  prob.ga_params.debug = infile("GA/debug",def.ga_params.debug);

  prob.ga_params.which_solver = infile("GA/which_solver",def.ga_params.which_solver.c_str());

  prob.ga_params.evaluate_use_initial_k0 = infile("GA/evaluate_use_initial_k0",def.ga_params.evaluate_use_initial_k0);

  prob.ga_params.pUseExpat = infile("GA/pUseExpat",def.ga_params.pUseExpat);

  prob.ga_params.always_accept_expat = infile("GA/always_accept_expat",
					      def.ga_params.always_accept_expat);

  /******** ILP ************/ 
  prob.ilp_model.flow_relay_constraints = infile("ILP/flow_relay_constraints",def.ilp_model.flow_relay_constraints);

  prob.ilp_model.use_im_flow_neighbor_constraints =
    infile("ILP/use_im_flow_neighbor_constraints", def.ilp_model.use_im_flow_neighbor_constraints);

  prob.ilp_model.im_flow_neighbor_nh =
    infile("ILP/im_flow_neighbor_nh", def.ilp_model.im_flow_neighbor_nh);
  //!
  prob.ilp_model.im_flow_neighbor_use_interference_range =
    infile("ILP/im_flow_neighbor_use_interference_range",
	   def.ilp_model.im_flow_neighbor_use_interference_range);
  prob.ilp_model.im_flow_neighbor_interference_range =
    infile("ILP/im_flow_neighbor_interference_range",
	   def.ilp_model.im_flow_neighbor_interference_range);

  prob.ilp_model.use_flow_neighbor_sink =
    infile("ILP/use_flow_neighbor_sink",
	   def.ilp_model.use_flow_neighbor_sink);
  
  prob.ilp_model.use_flow_neighbor_constraints =
    infile("ILP/use_flow_neighbor_constraints",
	   def.ilp_model.use_flow_neighbor_constraints);
  //!
  prob.ilp_model.strict_flow_neighbor_constraints =
    infile("ILP/strict_flow_neighbor_constraints",
	   def.ilp_model.strict_flow_neighbor_constraints);
  //!
  prob.ilp_model.lost_flow_penalty =
    infile("ILP/lost_flow_penalty",
	   def.ilp_model.lost_flow_penalty);
  
  prob.ilp_model.use_node_degree_constraints = infile("ILP/use_node_degree_constraints",
						      def.ilp_model.use_node_degree_constraints);

  prob.ilp_model.max_capacity = infile("ILP/max_capacity",def.ilp_model.max_capacity);

  prob.ilp_model.max_flow = infile("ILP/max_flow",def.ilp_model.max_flow);
  prob.ilp_model.max_local_flow = infile("ILP/max_local_flow",def.ilp_model.max_local_flow);

  prob.ilp_model.max_degree = infile("ILP/max_degree",def.ilp_model.max_degree);

  prob.ilp_model.default_demand = infile("ILP/default_demand",def.ilp_model.default_demand);

  prob.ilp_model.flow_neighbor_penalty_w = 
    infile("ILP/flow_neighbor_penalty_w",def.ilp_model.flow_neighbor_penalty_w);

  prob.ilp_model.penalty_relay = 
    infile("ILP/penalty_relay",def.ilp_model.penalty_relay);

  prob.ilp_model.integer_flow =
    infile("ILP/integer_flow", def.ilp_model.integer_flow);
  prob.ilp_model.allow_incomplete_delivery =
    infile("ILP/allow_incomplete_delivery", 
	   def.ilp_model.allow_incomplete_delivery);
  prob.ilp_model.incomplete_delivery_penalty =
    infile("ILP/incomplete_delivery_penalty", 
	   def.ilp_model.incomplete_delivery_penalty);


    /********  CPLEX config *********/
  prob.mycplex_params.mipGap = infile("CPLEX/mipGap", def.mycplex_params.mipGap);
  prob.mycplex_params.trelim = infile("CPLEX/trelim", def.mycplex_params.trelim);
  prob.mycplex_params.logToFile = infile("CPLEX/logToFile", def.mycplex_params.logToFile);
  prob.mycplex_params.workmem = infile("CPLEX/workmem", def.mycplex_params.workmem);
  prob.mycplex_params.nodefileind = 
    infile("CPLEX/nodefileind", def.mycplex_params.nodefileind);
  prob.mycplex_params.mipemphasis = 
    infile("CPLEX/mipemphasis", def.mycplex_params.mipemphasis);
  prob.mycplex_params.threads = infile("CPLEX/threads", def.mycplex_params.threads);
  prob.mycplex_params.lpmethod = infile("CPLEX/lpmethod", def.mycplex_params.lpmethod);
  prob.mycplex_params.preind = infile("CPLEX/preind", def.mycplex_params.preind);
  prob.mycplex_params.parallelmode = 
    infile("CPLEX/parallelmode", def.mycplex_params.parallelmode);
  prob.mycplex_params.mipdisplay = 
    infile("CPLEX/mipdisplay", def.mycplex_params.mipdisplay);
  prob.mycplex_params.probe = infile("CPLEX/probe", def.mycplex_params.probe);
  prob.mycplex_params.save_cpx_params =
    infile("CPLEX/save_cpx_params", def.mycplex_params.save_cpx_params);

  /// Matheuristic
  prob.matheuristic_params.comm = 
    infile("MATHEURISTIC/comm",def.matheuristic_params.comm.c_str());;
  prob.matheuristic_params.time_limit = 
    infile("MATHEURISTIC/time_limit",def.matheuristic_params.time_limit);;
  prob.matheuristic_params.socket_address = 
    infile("MATHEURISTIC/socket_address",def.matheuristic_params.socket_address.c_str());
  prob.matheuristic_params.socket_port = 
    infile("MATHEURISTIC/socket_port",def.matheuristic_params.socket_port);
}

void 
checkParameters()
{
  prob.output_path = fixPath(prob.output_path);
  prob.generation_path = fixPath(prob.generation_path);
  prob.instance_path = fixPath(prob.instance_path);
  prob.lp_params.output_path = fixPath(prob.lp_params.output_path);

  if(prob.ga_params.log_population_file == "default")
  {
    prob.ga_params.log_population_file = prob.instance_id + ".ga_pop_log";
  }
  if(prob.ga_params.log_generation_file == "default")
  {
    prob.ga_params.log_generation_file = prob.instance_id + ".ga_gen_log";
  }

  if(prob.ga_params.log_file == "default")
  {
    prob.ga_params.log_file = prob.instance_id + ".ga_log";
  }

  if(prob.ga_params.debug_file == "default")
  {
    prob.ga_params.debug_file = prob.instance_id + ".ga_debug";
  }
  if(prob.ga_params.output_file == "default")
  {
    prob.ga_params.output_file = prob.instance_id + ".ga_sol";
  }
  if(prob.Low_res_use_file && prob.Low_res_solution_file == "default")
    prob.Low_res_solution_file = prob.instance_id + ".sol";
  if(prob.network_file == "default")
    prob.network_file = prob.instance_id +".net";
  if(prob.output_time_file == "default")
    prob.output_time_file = prob.instance_id + ".time";	
  if(prob.output_sol_file == "default")
    prob.output_sol_file = prob.instance_id + ".sol";	
  if(prob.lp_params.input_mps_file == "default")
    prob.lp_params.input_mps_file = prob.instance_id +".mps.gz";
  if(prob.lp_params.output_mps_file == "default")
    prob.lp_params.output_mps_file = prob.instance_id +".mps";
  if(prob.output_conn_file == "default")
    prob.output_conn_file = prob.instance_id + ".conn";	
}

int 
fexist( char *filename ) 
{
  struct stat buffer ;
  if ( stat( filename, &buffer ) ) return 0 ;
  return 1 ;
}

void 
output_graphs(Graph *G, string filename)
{
  string graph_file,graph_dot_file,cmd;
  /** print graph */
  graph_file = prob.output_path + prob.instance_id + "-" + filename + ".ps";
  graph_dot_file = prob.output_path + prob.instance_id + "-" + filename + ".dot";
  G->toGraphviz(graph_dot_file); 
  if(prob.gen_graphs)
  {
    cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
    DEBUG(2) cdebug << "neato cmd" << cmd << endl;
    if(system(cmd.c_str()))
      printf("Generator: Error executing command - %s\n",cmd.c_str());
  }
}

