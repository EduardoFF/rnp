/*
 * =====================================================================================
 *
 *       Filename:  test.cpp
 *
 *    Description:  Unit test for project components
 *
 *        Version:  1.0
 *        Created:  06/15/2010 10:56:55 AM
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
#include "Plotter.h"
#include "network.h"
#include "solver.h"
#include "GetPot.h"

extern param prob;
extern param def;

 
void generate_grid(Network *net, double dimX, double dimY, double r );
void generate_cluster_nodes(Network *net, RandomGenerator *gen, double dimX, double dimY,
			    int num_static, int num_base, clusters_t &clu);

void generate_random_nodes(Network *net, RandomGenerator *, double dimX, double dimY,double,double,int num, nodeType t);
GetPot read_config_file(char *config_file);
void load_defaults();
void load_config_file(GetPot infile);

void test_geometry();
Network * convex_hull_restriction(Network *);
Network * los_band_restriction(Network *, double, double);
void generate_grid_h1(Network *net, double dimX, double dimY);
void generate_grid_h2(Network *net, double dimX, double dimY,double,double);






#ifdef _COMPILE_CAIROM
void test_grid_h1(){
	Network *net = new Network();	
	Graph *g;
	string graph_file;
	string cmd;
	int dimX, dimY;
	
	int seed = 123;

	UniformRandom *rand = new UniformRandom(seed);
	double rad = 125.0;
	dimX = 20;
	dimY = 20;	
	clusters_t clu;
	clu.n_x = 5;
	clu.n_y = 5;

	int ds_size = (clu.n_x*clu.n_y*sizeof(double));
	clu.density = (double *)malloc(ds_size);
	double density[] = { 0.1,0.0,0.0,0.1,0.0,
		        0.0,0.1,0.0,0.0,0.1,
			0.04,0.1,0.0,0.1,0.0,
			0.0,0.0,0.1,0.0,0.1,
			0.0,0.1,0.0,0.1,0.06};
	
	memcpy(clu.density,density,ds_size);

	generate_cluster_nodes(net, rand, dimX, dimY,
			    100, 5, clu);


	net->setRange(0.5);
	net->setRegionSize(dimX,dimY);


	generate_grid_h2(net, dimX,dimY,0.15,0.5);
	Network *net_hull = convex_hull_restriction(net);
	net = net_hull;
	
	
	printf("Grid generated - number of relays %d\n",net->numRelays());
	Plotter *plot = net->plot("grid_h1_plot.svg");
//	plot->write();
	g = net->createGraph();
	string graph_dot_file = "testgraph_gridh1.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "testgraph_gridh1.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());
	

}
#endif

void test_clusters(){
	Network *net = new Network();	
	Graph *g;
	string graph_file;
	string cmd;
	int dimX, dimY;
	
	int seed = 1234;

	double rad = 125.0;
	dimX = 5;
	dimY = 5;	


	net->setRange(1.0);
	net->setRegionSize(dimX,dimY);
	UniformRandom *rand = new UniformRandom(seed);
	clusters_t clu;
	clu.n_x = 5;
	clu.n_y = 5;

	int ds_size = (clu.n_x*clu.n_y*sizeof(double));
	clu.density = (double *)malloc(ds_size);
	double density[] = { 0.1,0.0,0.0,0.1,0.0,
		        0.0,0.1,0.0,0.0,0.1,
			0.04,0.1,0.0,0.1,0.0,
			0.0,0.0,0.1,0.0,0.1,
			0.0,0.1,0.0,0.1,0.06};
	
	memcpy(clu.density,density,ds_size);

	generate_cluster_nodes(net, rand, dimX, dimY,
			    100, 5, clu);
	Plotter *plot = net->plot("cluster_plot.svg");
//	plot->write();
	g = net->createGraph();
	string graph_dot_file = "clustergraph.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "clustergraph.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());
	generate_grid(net, dimX,dimY,0.25);

	g = net->createGraph();
	graph_dot_file = "clustergraph-grid.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "clustergraph-grid.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());
	



	


}



/**
 * @brief  Test for graphviz output
 * @param  
 * @return 
 */

#define BAND_FACTOR 30

#ifdef _COMPILE_CAIROM
void test_plot(){
	Network *net = new Network();	
	Graph *g;
	string graph_file;
	string cmd;
	int dimX, dimY;
	
	int seed = 123;

	double rad = 125.0;
	dimX = 500;
	dimY = 500;	


	net->setRange(1.0);
	net->setRegionSize(dimX,dimY);
	UniformRandom *rand = new UniformRandom(seed);
	generate_random_nodes(net, rand, dimX, dimY,0,0, 10, STATIC);
	generate_random_nodes(net,rand, dimX, dimY,0,0, 1 ,BASE);
	generate_grid(net, dimX,dimY,1.0);
	
	printf("Grid generated - number of relays %d\n",net->numRelays());
	
	//net->plot();
	Network *net_hull = convex_hull_restriction(net);
	Network *net_los;
	net_los = los_band_restriction(net_hull, net->getRange()*BAND_FACTOR,rad );
	Plotter *plot_net = net_los->plot("plot.svg");
//	plot_net->write();
	g = net_los->createGraph();
	string graph_dot_file = "testgraph_los.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "testgraph_los.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());
	//
	//Plotter *plot = new Plotter("plottest.svg",200,200);
	//plot->write();
}

#endif
int test_graphviz(){

	// Testing the output for specific network topologies
	Network *net = new Network();	
	Graph *g;
	string graph_file;
	string cmd;
	int dimX, dimY;
	int seed = 123;
	double rad = 5.0;
	dimX = 30;
	dimY = 30;	


	net->setRange(0.3);
	net->setRegionSize(dimX,dimY);
	UniformRandom *rand = new UniformRandom(seed);
	generate_random_nodes(net, rand, dimX, dimY,0,0, 20, STATIC);
	generate_random_nodes(net,rand, dimX, dimY,0,0, 3 ,BASE);
	generate_grid(net, dimX,dimY,0.3);
	
	printf("Grid generated - number of relays %d\n",net->numRelays());
	
	//net->plot();
	Network *net_hull = convex_hull_restriction(net);
	Network *net_los;

	net_los = los_band_restriction(net_hull, net->getRange()*BAND_FACTOR, rad);

	printf("Convex hull - number of relays %d\n",net_hull->numRelays());

	printf("LOS Restriction - number of relays %d\n",net_los->numRelays());

	g = net->createGraph();
	string graph_dot_file = "testgraph.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "testgraph.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());

	g = net_hull->createGraph();
	
	graph_dot_file = "testgraph_hull.dot";
	g->toGraphviz(graph_dot_file);
	graph_file = "testgraph_hull.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());


	//test_geometry();
	
	g = net_los->createGraph();

	graph_dot_file = "testgraph_los.dot";

	g->toGraphviz(graph_dot_file);
	graph_file = "testgraph_los.ps";
	
	cmd = GRAPHVIZ_CMD(graph_dot_file) + graph_file;
	if(system(cmd.c_str()))
		printf("Generator: Error executing command - %s\n",cmd.c_str());



}

void
test_serialization()
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


int main(){
	load_defaults();  //Load default values for problem parameters
	//test_graphviz();
	//test_plot();
//	test_clusters();
//	test_grid_h1();


}
