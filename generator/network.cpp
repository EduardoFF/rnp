/*
 * =====================================================================================
 *
 *       Filename:  network.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 04:50:13 PM
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

#include "graph.h"

//#include "Plotter.h"
#include "network.h"


extern param prob;

bool Network::isBase(int i){ 
	return (i < b_ix.size());
}

bool Network::isStatic(int i){ 
	return (i < s_ix.size());
}
bool Network::isRelay(int i){ 
	return (i < r_ix.size());
}


Node &Network::getBase(int i){ 
	if(i < b_ix.size())
		return nodes[b_ix[i]];
	else
		printf("Network: Invalid index\n");
}

Node &Network::getStatic(int i){ 
	if(i < s_ix.size())
		return nodes[s_ix[i]];
	else
		printf("Network: Invalid index\n");
}
Node &Network::getRelay(int i){ 
	if(i < r_ix.size())
		return nodes[r_ix[i]];
	else
		printf("Network: Invalid index\n");
}


double Network::getDemand(int i)
{
	assert(0 <= i && i < numStatic());
	return m_demand[i];
}

void Network::setDemand(int i, double d)
{
	assert(0 <= i && i < numStatic());
	m_demand[i] = d;
}
void Network::write(string filename){


	ofstream outfile(filename.c_str());
	
	outfile << "STATIC\t" << numStatic() << endl;
	outfile << "BASES\t" << numBases() << endl;
	outfile << "RELAYS\t" << numRelays() << endl;

	outfile << "X\t" << dimX << endl;
	outfile << "Y\t" << dimY << endl;
	outfile << "TX_RANGE\t" << tx_range << endl;

	/** Calculate number of links
	 *  n_links_all: All links, including grid locations
	 *  n_links_sb:  Only links in initial topology 
	 *  */
	int n_links_all = 0;
	int n_links_sb = 0;
	for(int i=0; i<nodes.size(); i++){
		Node &ni = nodes[i];
		for(int j=0;j<nodes.size();j++){
			if( i == j)
				continue;
			Node &nj = nodes[j];
			if(!inRange(ni,nj))
				continue;
			
			if(ni.t != STATIC){
				n_links_all +=1;
				if(ni.t != RELAY && nj.t != RELAY)
					n_links_sb +=1;
			}
			if(nj.t != STATIC){
				n_links_all +=1;
				if(ni.t != RELAY && nj.t != RELAY)
					n_links_sb +=1;
			}
			
		}
	}
	outfile << "LINKS_ALL\t" << n_links_all << endl;
	outfile << "LINKS_INIT\t" << n_links_sb << endl;

	for(int i=0;i<numStatic();i++){
		Node &n = getStatic(i);
		outfile << n.id << "\t" << n.x << "\t" << n.y << endl; 
	}
	for(int i=0;i<numBases();i++){
		Node &n = getBase(i);
		outfile << n.id << "\t" << n.x << "\t" << n.y << endl; 
	}
	for(int i=0;i<numRelays();i++){
		Node &n = getRelay(i);
		outfile << n.id << "\t" << n.x << "\t" << n.y << endl; 
	}

	outfile.close();

}

/**
 * \brief  loads network toology from file
 * 	   file format is
 * 	   <number static> <number bs> <number rns>
 *
 * \param filename
 */
Network::Network(string filename)
{
  ifstream ifile(filename.c_str());

  m_maxPathLength = 0;
  cout << "Reading network from : " << filename << endl;
  int nums,numb,numr;
  double dx,dy,tx;
  int nlinks, nlinks_relays, ehops;
  string label;
  ifile >> label >> nums;
  ifile >> label >> numb; 
  ifile >> label >> numr;
  ifile >> label >> dx;
  ifile >> label >> dy;
  ifile >> label	>> tx;
  setRange(tx);
  setRegionSize(dx,dy);
  cout << nums << " " << numb << " " << numr << endl;
  cout << dx << " " << dy << " " << tx << endl;

  ifile >> label >> nlinks;
  ifile >> label >> nlinks_relays;

  cout << nlinks << " " << nlinks_relays << endl;
  for(int i=0;i<nums;i++)
  {
    string id;
    double x,y;
    ifile >> id >> x >> y;

    cout << id << " " << x << " " << y <<endl;
    if(!addNode(x,y,STATIC,id))
    {
      cout << "Invalid network (inter-node distance < " << MIN_NODE_DISTANCE << ")\n";
      exit(-1);
    }

  }
  for(int i=0;i<numb;i++){
    string id;
    double x,y;
    ifile >> id >> x >> y;
    if(!addNode(x,y,BASE,id)){
      cout << "Invalid network (inter-node distance < " << MIN_NODE_DISTANCE << ")\n";
      exit(-1);

    }
  }
  for(int i=0;i<numr;i++){
    string id;
    double x,y;
    ifile >> id >> x >> y;
    addNode(x,y,RELAY,id);
  }

  if( prob.no_position || prob.fixed_links )
  {
    /// Read links

    for(int i=0;i<nlinks;i++)
    {
      string ids, idt;
      double w;
      ifile >> ids >> idt >> w;
      setLink(ids,idt,w);
    }
  }

  if( prob.no_position)
  {
    /// Read also expected hop count
    ifile >> label >> ehops;
    setExpectedHopCount(ehops);
  }

  ifile.close();



}

void Network::print_info(){
	cout << "Network info: " << endl;
	cout << numStatic() << " static nodes" << endl;
	cout << numBases() << " sink nodes" << endl;
	cout << numRelays() << " relay nodes" << endl;
	cout << links.size() << " links" << endl;
}

/**
 * @brief  adds a node in the network topology at the specified location
 *
 * 	   Note that the minimum distance condition is assumed for static a nodes and sinks
 * 	   For relays it does not make sense because they are in a grid, which could place them 
 * 	   as close as required
 * @param  x: x-coordinate of location
 * @param  y: y-coordinate of location
 * @param  _t: type of the node (defined in generator.h)
 * @param  _id: identifier of the node. If blank, the function generates one based on its type
 * @return 1 if the node location does not violate the inter-node distance (MIN_NODE_DISTANCE)
 * 	   0 otherwise 
 */
bool Network::addNode(Node &n)
{
	return addNode(n.x,n.y,n.t,n.id);
}

  bool 
Network::addNode(double _x, double _y, nodeType _t, string _id)
{
  Node n(_x,_y);
  n.t = _t;


  if(_id == "")
  {

    // If id is not provided, 
    // define a nice node_id

    stringstream ss;
    if(n.t == STATIC)
    {
      ss << "s" <<  s_ix.size();
    }
    else if(n.t == BASE)
    {
      ss << "b" <<  b_ix.size();
    }
    else if(n.t == RELAY)
    {
      ss << "r" <<  r_ix.size();
    }

    n.id = ss.str();
  }
  else
  {
    n.id = _id;
  }
  // This step is required for GA solver
  // If node already exists in network, quietly ignore it
  // assert(id_node.find(n.id) == id_node.end());
  if( id_node.find(n.id) != id_node.end())
  {
    return true;
  }

  /// If it is a static or sink, check distance between 
  /// other nodes to see if the location is correct
  if(!prob.no_position && !prob.fixed_links)
  {
    for(int i=0;i<numStatic();i++)
    {
      Node &m = getStatic(i);
      if (distanceLink(n,m) < MIN_NODE_DISTANCE)
      {
	fprintf(stderr, "Invalid node distance!\n");
	return false;			
      }
    }
    for(int i=0;i<numBases();i++)
    {
      Node &m = getBase(i);
      if (distanceLink(n,m) < MIN_NODE_DISTANCE)
      {
	fprintf(stderr, "Invalid node distance!\n");
	return false;			
      }
    }
    /// check distance with other RELAYS iff the new node is a RELAY
    if( n.t == RELAY )
    {
      for(int i=0;i<numRelays();i++)
      {
	Node &m = getRelay(i);
	if (distanceLink(n,m) < MIN_NODE_DISTANCE)
	{
	  fprintf(stderr, "Invalid node distance!\n");
	  return false;			
	}
      }
    }
  }
  if (n.t == STATIC)
  {
    s_ix.push_back(nodes.size());
    if(prob.demand_type == "uniform")
    {
      m_demand[s_ix.size()-1] = prob.ilp_model.default_demand;
    } 
    else if( prob.demand_type == "random")
    {
      double delta = prob.uniform_random_demand_high - prob.uniform_random_demand_low; 
      double r_demand = 1.0*rand()/RAND_MAX;
      r_demand *= delta;
      r_demand += prob.uniform_random_demand_low;
      m_demand[s_ix.size()-1] = r_demand;
    }
    else
    {
      /// must be from file, set to 1 by default
      m_demand[s_ix.size()-1] = 1;
    }
    //printf("Setting demand to %f\n", prob.ilp_model.default_demand);
  }
  else if(n.t == RELAY)
  {
    id_relay[n.id] = r_ix.size();
    r_ix.push_back(nodes.size());
  }
  else if(n.t == BASE)
  {
    b_ix.push_back(nodes.size());
  }

  id_node[n.id] = nodes.size();
  nodes.push_back(n);
  return true;
}
/**
 * @brief  Determines if two node locations are in range considering the current transmission range
 * @param  m is a Node
 * @param  n is a Node
 * @return true if location of node m and n are in range
 */
bool Network::inRange(const Node &m, const Node &n){
  if(prob.no_position || prob.fixed_links)
  {
    bool b = areLinked(m,n);
//    printf("%s %s are linked? %d\n",m.id.c_str(), n.id.c_str(), b);
    return b;

  } else 
  {
	double dx = m.x - n.x;
	double dy = m.y - n.y;
	return (sqrt(dx * dx + dy * dy) <= tx_range);
  }
}


double Network::distanceLink(const Node &m, const Node &n){
	double dx = m.x - n.x;
	double dy = m.y - n.y;
	double dist = sqrt(dx * dx + dy * dy);
	return dist;

}


/*  Calculates the number of neighbors a link has 
 *  defined as the sum of the number of neighbors of 
 *  its source and target */


int Network::neighborsLink(const Node &m, const Node &n){
	int cnt = 0;
	for(int i =0;i< size();i++){

		Node &k = getNode(i);
		if(k.t == RELAY)
				continue;
		if(inRange(n,k) || inRange(m,k))
			cnt++;
	}
	return cnt - 2;

}

int Network::numNeighbors(const Node &n){
	int cnt = 0;
/*  	for(int i =0;i< size();i++){

		Node &m = getNode(i);
		if(m.t == RELAY)
				continue;
		if(inRange(n,m))
			cnt++;
	}
*/
	return cnt -1;
}

int 
Network::getRelayIndex(string id)
{
	ITERATOR(id_relay) it = id_relay.find(id);
	if( it != id_relay.end())
		return it->second;
	else
		return INVALID_NODE_INDEX;
}
int 
Network::getNodeIndex(string id)
{
	ITERATOR(id_node) it = id_node.find(id);
	if( it != id_node.end())
		return it->second;
	else
		return INVALID_NODE_INDEX;
}
Node &Network::getNodeById(string id){

	return nodes[id_node[id]];
	if(id[0] == 's'){
		int ix = lexical_cast<int> (id.substr(1,id.size()-1));
		if(ix < s_ix.size())
			return nodes[s_ix[ix]];
		else
			printf("Network: Invalid node id\n");


	}else if(id[0] == 'r'){
		int ix = lexical_cast<int> (id.substr(1,id.size()-1));
		if(ix < r_ix.size())
			return nodes[r_ix[ix]];
		else
			printf("Network: Invalid node id\n");

	}else if(id[0] == 'b'){
		int ix = lexical_cast<int> (id.substr(1,id.size()-1));
		if(ix < b_ix.size())
			return nodes[b_ix[ix]];
		else
			printf("Network: Invalid node id\n");

	}else
		printf("Network: Invalid node id\n");

}
	Node &Network::getNode(int i){
		if( i < nodes.size())
			return nodes[i];
		else{
			printf("Network: Invalid index\n");
			exit(-1);
		}
	}
int
Network::size()
{
  return nodes.size();
}

void
Network::setRange(double r)
{
  tx_range = r;
}

double
Network::getRange(){
  return tx_range;
}

/**
 * @brief  creates an instance of the graph class using the nodes
 * @param  
 * @return 
 */
Graph *Network::createGraph(){
  Graph *G = new Graph(size());
  for(int i=0;i<size();i++){
    G->set_vertex_id(i,getNode(i).id);

    if( !prob.no_position)
    {
      /** node coordinates are normalized  */
      double x = getNode(i).x/dimX;
      double y = getNode(i).y/dimY;
      G->set_vertex_loc(i,x, y);
    }
    G->set_vertex_type(i,getNode(i).t);
  }
  for( int i=0; i < size(); i++)
    for(int j = i+1; j < size(); j++)
      if( inRange(getNode(i),getNode(j)))

        G->addEdge(i,j);
  return G;
}

void Network::setRegionSize(double _x, double _y){
	dimX  =_x;
	dimY = _y;
}
/**
 * @brief  generate plot file for the initial random locations
 * @param  none
 * @return none
 */
Plotter *Network::plot(string filename){
#ifdef _COMPILE_CAIROM

	Plotter *plot = new Plotter(filename,dimX,dimY);
	for(int i = 0; i<numStatic();i++){
		Node &n = getStatic(i);
		plot->draw_point(n.x, n.y);
	}	
	for(int i=0;i<numBases();i++){
		Node &n = getBase(i);
		plot->draw_point(n.x,n.y);
	}
	return plot;
	//plot.write();
#else
	return NULL;
#endif



/*	ofstream staticFile ("nodes.dat");
	ofstream baseFile ("base.dat");
	ofstream relayFile ("relay.dat");


	if (staticFile.is_open() && baseFile.is_open() && relayFile.is_open() )
	{
		for(int i = 0; i < numStatic(); i++)
			staticFile << getStatic(i).x << " " << getStatic(i).y << "\n";
		for(int i = 0; i < numBases(); i++)
			baseFile << getBase(i).x << " " << getBase(i).y << "\n";
		for(int i = 0; i < numRelays(); i++)
			relayFile << getRelay(i).x << " " << getRelay(i).y << "\n";

		staticFile.close();
		baseFile.close();
		relayFile.close();
		///system("gnuplot -persist plot.sh");
	}else printf("Network: (Plot) Unable to open files\n");
	*/

}


metrics_t *Network::get_metrics(){
	metrics_t *m = new metrics_t();
	Graph *g = createGraph();
	int minDist[numStatic()]; 
	

	/** for each base node - call bfs and get distances */
/*  
	for(int b = 0;b < numBases();b++){
		vertex_d b_index = b_ix[b]; // this index also match index in graph
		boost::graph_traits<MyGraphType>::vertices_size_type d[size()];
		std::fill_n(d, size(), 0);

		boost::breadth_first_search(g,b_index,
					    boost::make_bfs_visitor(
						    std::make_pair(boost::record_distances(d,boost::on_tree_edge()),
								   bfs_visitor<null_visitor>())));







		
	}
*/
	return m;
}

double Network::getLinkWeight(const Node &m, const Node &n)
{ 
  int mi = getNodeIndex(m.id);
  int ni = getNodeIndex(n.id);
  assert(mi != INVALID_NODE_INDEX && ni != INVALID_NODE_INDEX);
  ITERATOR(links) it = links.find( make_pair(mi, ni));
  if( it == links.end())
    return 0.0;
  else
    return it->second;
}

double Network::getLinkWeight(string mid, string nid)
{
  Node &m = getNodeById(mid);
  Node &n = getNodeById(nid);
  return getLinkWeight(m,n);
}
void Network::setLink(string mid, string nid, double w)
{
  Node &m = getNodeById(mid);
  Node &n = getNodeById(nid);
  setLink(m,n,w);
}

bool Network::areLinked(string mid, string nid)
{
  Node &m = getNodeById(mid);
  Node &n = getNodeById(nid);
  return areLinked(m,n);
}

void Network::setLink(const Node &m, const Node &n, double w)
{
  int mi = getNodeIndex(m.id);
  int ni = getNodeIndex(n.id);
  //printf("Setting link %s %s %.2f\n", m.id.c_str(), n.id.c_str(),w);
  assert(mi != INVALID_NODE_INDEX && ni != INVALID_NODE_INDEX);
  links[ make_pair(mi, ni)] = w;
}

bool Network::areLinked(const Node &m, const Node &n)
{
  int mi = getNodeIndex(m.id);
  int ni = getNodeIndex(n.id);
  assert(mi != INVALID_NODE_INDEX && ni != INVALID_NODE_INDEX);
  ITERATOR(links) it = links.find( make_pair(mi, ni));
  return it != links.end();
}

int Network::maxPathLength()
{
  if( prob.no_position )
    {
      fprintf(stderr, "maxPathLength not supported for no_position\n");
      exit(-1);
    }
  //! we need to compute it the first time
  if( m_maxPathLength <= 0)
    {
      printf("Network: computing maxPathLength\n");
      for(int i = 0; i< numStatic();i++)
	{
	  Node &ns = getStatic(i);
	  int mpl = -1;
	  for(int j=0;j < numBases();j++)
	    {
	      Node &nb = getBase(j);
	      double d = distanceLink(ns,nb);
	      int hc = (int)ceil(1.0*d/getRange());
	      if( hc > mpl )
		mpl = hc;
	  }
	  if( m_maxPathLength < mpl )
	    m_maxPathLength = mpl;
	}
      printf("Network: maxPathLength %d\n", m_maxPathLength);
    }
  return m_maxPathLength;
}

int Network::expectedHopCount()
{
    if( !prob.no_position)
    {
      int expected_hop_count = 0;
      for(int i = 0; i< numStatic();i++){
        Node &ns = getStatic(i);
        int hop_count = 100000000;
        for(int j=0;j < numBases();j++){
          Node &nb = getBase(j);
          double d = distanceLink(ns,nb);
          int hc = (int)ceil(1.0*d/getRange());
          hop_count = min(hc,hop_count);
        }
        expected_hop_count += hop_count;
      }

      return expected_hop_count;
    } else
    {
      return m_expected_hops;
    }
}
