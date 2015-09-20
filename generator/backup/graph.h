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

#define GRAPHVIZ_POINTSCALE 100
#define GRAPHVIZ_NODEWIDTH 0.01


//struct node_id_t {
//	     typedef vertex_property_tag kind;
//};
typedef property<vertex_name_t, string > IDProperty;
typedef adjacency_list<vecS, vecS, undirectedS, IDProperty> MyGraphType;


class Graph : public MyGraphType {	
	public:
		Graph(int n): MyGraphType(n){}
		
		void addEdge(int i, int j);
		void set_vertex_id(int , string );

		int toGraphviz(string fname);

	private:
		//MyGraphType G;
};


class position_writer {
	public:
		vector<string> idmap;

		position_writer(vector<string> vs): idmap(vs){}
		template <class VertexOrEdge>
			void operator()(std::ostream& out, const VertexOrEdge& v) const;
};

position_writer make_position_writer(vector<string> ids);

#endif
