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
#include "generator.h"
#include "graph.h"
#include "solver.h"
#include "network.h"


Network net;

param prob;
param def; //Default values for parameters




Graph *G;

/**
 * @brief  generate plot file for the initial random locations
 * @param  none
 * @return none
 */
void plotFile(){

	ofstream staticFile ("nodes.dat");
	ofstream baseFile ("base.dat");
	ofstream relayFile ("relay.dat");


	if (staticFile.is_open() && baseFile.is_open() && relayFile.is_open() )
	{
		for(int i = 0; i < net.numStatic(); i++)
			staticFile << net.getStatic(i).x << " " << net.getStatic(i).y << "\n";
		for(int i = 0; i < net.numBases(); i++)
			baseFile << net.getBase(i).x << " " << net.getBase(i).y << "\n";
		for(int i = 0; i < net.numRelays(); i++)
			relayFile << net.getRelay(i).x << " " << net.getRelay(i).y << "\n";

		staticFile.close();
		baseFile.close();
		relayFile.close();
		///system("gnuplot -persist plot.sh");
	}else printf("generator: (Plot) Unable to open files\n");

}



bool inRange(const Node &m, const Node &n){
	  double dx = m.x - n.x;
	  double dy = m.y - n.y;
	 return (sqrt(dx * dx + dy * dy) < prob.tx_range);
}

void createGraph(){
	G = new Graph(net.size());
	for(int i=0;i<net.size();i++)
		G->set_vertex_id(i,net.getNode(i).id);
	for( int i=0; i < net.size(); i++)
		for(int j = i+1; j < net.size(); j++)
			if( inRange(net.getNode(i),net.getNode(j)))
				G->addEdge(i,j);



}

/**
 * @brief Generates grid points with indicated resolution 
 * @param  r resolution to use
 * @return 
 */
	

void generateGrid(double r){
	double x, y;

	x = 0;
	while(x < prob.dimX){
		y = 0;
		while(y < prob.dimY){
			Node n(x,y);
			n.t = RELAY;
			net.addNode(n);
			y += r;
		}
		x += r;
	}

}


void generateRandom(){
//randgen<boost::uniform_real<> > gen_x(0,1);
  	variate_generator<mt19937, uniform_real<> >
		gen_x(mt19937(prob.seed),uniform_real<>(0,prob.dimX));
  	variate_generator<mt19937, uniform_real<> >
		gen_y(mt19937(prob.seed*prob.seed),uniform_real<>(0,prob.dimY));


//	gen_x.min = 0.0;
//
        VERBOSE(1) printf("generating %d static nodes...\n",prob.st_nodes);
	for(int i = 0; i < prob.st_nodes; i++){

		Node n(gen_x(),gen_y());
		n.t = STATIC;
		net.addNode(n);

		VERBOSE(2)cout << n.x << ", " << n.y << endl;


	}
	VERBOSE(1)printf("generating %d base nodes...\n",prob.b_nodes);
	for(int i = 0; i < prob.b_nodes; i++){

		Node n(gen_x(),gen_y());
		n.t = BASE;
		net.addNode(n);
		VERBOSE(2) printf("( %f, %f)\n",n.x,n.y);


	}

}

GetPot read_config_file(char *config_file){

  // Define the syntax for the input configuration file and start reading
  string base_directory  = "./";  
  string comment_start   = "#";  
  string comment_end     = "\n";
  string field_separator = ",";

  string input_file    = config_file;

  GetPot infile((base_directory + input_file).c_str(), 
		       comment_start.c_str(), comment_end.c_str(), field_separator.c_str()); 
  return (infile);
}


void load_defaults(){
	def.model_file = "./model.mod";
	def.data_file = "./model.dat";	
	def.st_nodes = 10;
	def.b_nodes = 2;
	def.K = 1000000;
	def.tx_range = 1.0;
	def.grid_res = 1.0;
	def.dimX = 10.0;
	def.dimY = 10.0;
	def.seed = time(0);
	def.verbose = 1;
	def.debug = 0;
	def.instance_id = "default";
	def.use_gurobi = false;
	def.glpk_out = 0;
	def.grb_out = 0;

}

void load_config_file(GetPot infile){
	prob.st_nodes = infile("ProblemInstance/static",def.st_nodes);
        prob.b_nodes = infile("ProblemInstance/base",def.b_nodes);
        prob.K = infile("ProblemInstance/max_relays",def.K);

        prob.tx_range = infile("ProblemInstance/tx_range",def.tx_range);

        prob.grid_res = infile("ProblemInstance/grid_res",def.grid_res);

        prob.dimX = infile("ProblemInstance/X",def.dimX);
        prob.dimY = infile("ProblemInstance/Y",def.dimY);
        prob.seed = infile("Parameters/seed",(int)def.seed);
        prob.model_file = infile("ProblemInstance/model",def.model_file.c_str());
        prob.data_file = infile("ProblemInstance/data",def.data_file.c_str());
        prob.verbose = infile("Parameters/verbose_level",def.verbose);
        prob.debug = infile("Parameters/debug_level",def.debug);
	prob.instance_id = infile("ProblemInstance/instance_id",def.instance_id.c_str());
	prob.use_gurobi = infile("Parameters/use_gurobi",def.use_gurobi);
	prob.glpk_out = infile("Parameters/glpk_out",def.glpk_out);
	prob.grb_out = infile("Parameters/grb_out",def.grb_out);




}
int main(int argc, char *argv[]){

	GetPot infile;


	load_defaults();
	if(argc == 1){
		cout << "No config file specified. Using default values" << endl;
		prob = def; // using default values
	}else if(argc == 2){
		printf("Using config file %s\n",argv[1]);
		infile = read_config_file(argv[1]);
		load_config_file(infile);
	}

	VERBOSE(1){
		printf("generating random instance:\n");
		printf("	static nodes: %d\n",prob.st_nodes);
		printf("	  base nodes: %d\n",prob.b_nodes);
		printf(" transmission range: %f \n",prob.tx_range);
	}
	generateRandom();   //generate random locations 
	generateGrid(prob.grid_res); //generate grid points and possible relays
	createGraph();  // create network graph using generated random locations and tx range

	GLPSol lp(prob.tx_range);  // create LP problem
	Graph *sol = lp.solve(prob.K); // solve lp

	 
	/** print output data */
	string graph_file = prob.instance_id + "-graph.ps";
	string cmd;
	string sol_file = prob.instance_id + "-sol.ps";

	G->toGraphviz("graph.dot"); 
	sol->toGraphviz("sol.dot");

	cmd = GRAPHVIZ_CMD(graph.dot) + graph_file;
	DEBUG(1) printf("neato cmd = %s\n",cmd.c_str());
	system(cmd.c_str());
	cmd = GRAPHVIZ_CMD(sol.dot) + sol_file;
	system(cmd.c_str());

	plotFile();




}
