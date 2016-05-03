/*
 * =====================================================================================
 *
 *       Filename:  mst.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/22/2010 11:32:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */



#include "main.h"
#include "generator.h"
//#include "Plotter.h"
#include "graph.h"
#include "network.h"
//#include "metric.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/config.hpp>

using namespace boost;


void
MST_heuristic(Network *net)
{
  typedef adjacency_list < vecS, vecS, undirectedS,
			   property<vertex_distance_t, int>, property < edge_weight_t, int > > Graph;

  Graph g(net->size());

  for( int i=0; i < net->size(); i++)
    {
      for(int j = i+1; j < net->size(); j++)
	{
	  //		if( net->inRange(net->getNode(i),net->getNode(j)))
	  add_edge(i,j,g);
	}
    }

  graph_traits < Graph >::edge_iterator ei, ei_end;
  property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
  int i = 0;
  for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei, ++i)
    {
      int is = source(*ei, g);
      int js = target(*ei,g);	
      Node &nn = net->getNode(is);
      Node &nm = net->getNode(js);
      double d;
      if(nm.t == BASE && nn.t == BASE)
	{
	  d=0;
	}
      else
	{
	  d = net->distanceLink(nn,nm);
	}
      weightmap[*ei] = d;
      std::cout << "weight for " << is << " " << js << " = " << d << std::endl;
    }
  std::vector < graph_traits < Graph >::vertex_descriptor >
    p(num_vertices(g));
  prim_minimum_spanning_tree(g, &p[0]);
  int r_cnt = 0;
  for (std::size_t i = 0; i != p.size(); ++i)
    {
      Node &n = net->getNode(i);
      std::cout << n.id;
      if (p[i] != i)
	{
	  Node &m = net->getNode(p[i]);
	  if(n.t == BASE && m.t == BASE)
	    continue;
	  double dnm = net->distanceLink(n,m);
	  if(dnm <= net->getRange()){
	    continue;
	  }
	  std::cout << "Adding relays to " << n.id << " " << m.id << " ";
	  int n_r = (int)floor(dnm/net->getRange());
	  std::cout << n_r << endl;

	  double ds = dnm/(n_r+1);
	  std::cout << "S: (" << n.x << "," << n.y << " " << m.x << "," << m.y << ": ";

	  for(int j=1;j<=n_r;j++)
	    {
	      int k1 = j;
	      int k2 = n_r-j+1;
	      double xr = (k1*n.x + k2*m.x)/(k1+k2);
	      double yr = (k1*n.y + k2*m.y)/(k1+k2);
	      std::cout << " (" << xr << "," << yr << ") ";
	      net->addNode(xr,yr,RELAY,"");
	    }
	  std::cout << endl;
	  std::cout << "parent[" << i << "] = " << p[i] << std::endl;
	}
      else
	{
	  std::cout << "parent[" << i << "] = no parent" << std::endl;
	}
    }
}
