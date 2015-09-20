/*
 * =====================================================================================
 *
 *       Filename:  graph.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 11:59:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _GRAPH
#define _GRAPH

#define GRAPHVIZ_POINTSCALE 50
#define GRAPHVIZ_NODEWIDTH 1
#define GRAPHVIZ_CMD(x) "neato -Tps -n2 " + x + " > "

//#include "network.h"
//enum nodeType { STATIC, BASE, RELAY};

//struct node_id_t {
//	     typedef vertex_property_tag kind;
//};

/*  struct vertex_location_t {
    typedef vertex_property_tag kind;
};

namespace boost
{
    BOOST_INSTALL_PROPERTY(vertex, location);
}
*/
//enum vertexType {STATIC,BASE,RELAY}; 

struct vertex_info{
	string id;
	nodeType t;
	double x,y;

};

struct edge_info{
	string label;
};

//typedef property<vertex_name_t, string, property<vertex_location_t, pair< double, double> > > VertexProperty;
typedef adjacency_list<vecS, vecS, undirectedS, vertex_info,edge_info> MyGraphType;

typedef boost::graph_traits<MyGraphType>::edge_descriptor edge_d;
typedef graph_traits < MyGraphType >::vertex_descriptor vertex_d;
class Graph : public MyGraphType {	
	public:
		Graph(int n): MyGraphType(n){}
		
		void addEdge(int i, int j);
		void set_edge_label(int i, int j, string l);
		void set_vertex_id(int , string );
		void set_vertex_loc(int, double, double);
		void set_vertex_type(int, nodeType);

		int toGraphviz(string fname);
		int toNetFile(ofstream &f);

	private:
		//MyGraphType G;
};
class edge_writer {
	public:
		map< edge_d, string> idmap;

		edge_writer(map< edge_d, string> vs): idmap(vs) {}
		template <class VertexOrEdge>
			void operator()(std::ostream& out, const VertexOrEdge& v) const;
};


class position_writer {
	public:
		vector<string> idmap;
		vector< pair< double, double> > locmap;
		vector< nodeType> typemap;

		position_writer(vector<string> vs, vector< pair< double, double> > vl, vector<nodeType> vt): idmap(vs), locmap(vl), typemap(vt){}
		template <class VertexOrEdge>
			void operator()(std::ostream& out, const VertexOrEdge& v) const;
};

edge_writer make_edge_writer(map<edge_d, string> m );
position_writer make_position_writer(vector<string> ids, vector< pair<double, double> > locs, vector<nodeType> vt);
int numberOfConnectedComponents(Graph *);





#endif
