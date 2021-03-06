/*
 * =====================================================================================
 *
 *       Filename:  network.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 04:49:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _NETWORK 
#define _NETWORK

class Plotter;


class Node{
 public:
  Node(const Node &n)
    {
      id = n.id;
      t = n.t;
      x = n.x;
      y = n.y;
    }
 Node(double _x, double _y) : x(_x), y(_y){}
  nodeType t;
  string id;
  double x, y;

};



/**
 * @class Network
 * @brief Represents a network topology which consists on static nodes, base stations and relay nodes.
 */

typedef struct metrics{
  int avg_path_len;
} metrics_t;

class Network{
 public:
  vector<Node> nodes;	///< @brief set of nodes of the network
  double tx_range;	///< @brief transmission range
  vector<int> s_ix;	 ///< @brief set of indexes of static nodes
  vector<int> b_ix;	 ///< @brief set of indexes of base stations
  vector<int> r_ix;	 ///< @brief set of indexes of relay nodes

  map<string, int> id_node;
  map<string, int> id_relay;

  std::map< pair<int, int>, double > links; ///< @brief set of links (when no_position is enabled)
  int m_expected_hops; /// (when no_position is enabled)
  int m_maxPathLength;
  enum {
    INVALID_NODE_INDEX=-1,
  };

  map<int, double> m_demand;

  double dimX, dimY;

  int numStatic(){ return s_ix.size();}
  int numBases(){ return b_ix.size();}
  int numRelays(){ return r_ix.size();}
  Node &getNodeById(string id);   
						 
  Network(string);
  Network(){ m_maxPathLength = 0;}
  Network(const Network &n)
    {
      tx_range = n.tx_range;
      nodes = n.nodes;
      s_ix = n.s_ix;
      b_ix = n.b_ix;
      r_ix = n.r_ix;
      id_node = n.id_node;
      dimX = n.dimX;
      dimY = n.dimY;
      m_demand = n.m_demand;
      m_expected_hops = n.m_expected_hops;
      links = n.links;
      m_maxPathLength = n.m_maxPathLength;
    }


  void write(string);

  bool isStatic(int i);
  bool isBase(int i);
  bool isRelay(int i);

  Node &getStatic(int i);
  Node &getBase(int i);
  Node &getRelay(int i);
  void print_info();

  int size();
  Node &getNode(int i);

  void setExpectedHopCount(int h) { 
    printf("Got ExpectedHopCount %d\n", h);
    m_expected_hops = h;}
  int getNodeIndex(string id);
  int getRelayIndex(string id);

  bool addNode(double _x, double _y, nodeType _t, string s); 
  bool addNode(Node &);
  double distanceLink(const Node &n, const Node &m);
		
  int numNeighbors(const Node &n);

  int neighborsLink(const Node &m, const Node &n);

  bool inRange(const Node &n, const Node &m);
  Graph *createGraph();
  void setRegionSize(double, double);
  void setRange(double);
  double getRange();

  void setLink(const Node &n, const Node &m, double w);
  bool areLinked(const Node &n, const Node &m);

  double getLinkWeight(const Node &m, const Node &n);
  double getLinkWeight(string, string );

  void setLink(string s, string t, double w);
  bool areLinked(string s, string t);
  Plotter *plot(string filename);

  double getDemand(int i);
  void setDemand(int i, double d);

  int expectedHopCount();
  int maxPathLength();


  metrics_t *get_metrics();


};
#endif
