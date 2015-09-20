/*
 * =====================================================================================
 *
 *       Filename:  graph.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2010 12:01:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "generator.h"
#include "network.h"
#include "graph.h"

extern Network net;

extern param prob;

void Graph::addEdge(int i, int j){
	add_edge(i,j,*this);
}
void Graph::set_vertex_id(int i, string id){
	property_map<MyGraphType, vertex_name_t>::type
		vertex_name = get(vertex_name_t(), *this);
	//printf("Setting id[%d] = %s\n",i,id.c_str());
	boost::put(vertex_name,i,id ); 
	
}


struct graph_writer {
    void operator()(std::ostream& out) const {
      out << "label = \" " << prob.instance_id <<" \""<< endl;
      out << "fontsize = 30" << endl;
      out << "labelloc = \"t\" " << endl;
      out << "overlap = \"scale\"" << endl;
    }
  };

int Graph::toGraphviz(string filename){
	ofstream dotFile (filename.c_str());
	if( !dotFile.is_open()){
		return -1;
	}

	vector<string> ids;

//	cout << "Num vertices " << num_vertices(*this) << endl;
	property_map<MyGraphType, vertex_name_t>::type idmap = get(vertex_name_t(),*this);
	 graph_traits<MyGraphType>::vertex_iterator vi, vi_end, next;
	   tie(vi, vi_end) = vertices(*this);
	     for (next = vi; next != vi_end; next++) {
	       // printf("Vertex id [%d] = %s\n",*next,idmap[*next].c_str());
		ids.push_back(idmap[*next]);

	     }


	     



write_graphviz(dotFile, *this,make_position_writer(ids),default_writer(),graph_writer());
return 0;

}



position_writer make_position_writer(vector<string> ids){
	return position_writer(ids);
}
template <class VertexOrEdge>
void position_writer::operator()(std::ostream& out, const VertexOrEdge& v) const {


	Node &n = net.getNodeById(idmap[v]);
	out << "[pos=\"" << (n.x)*GRAPHVIZ_POINTSCALE << "," << (n.y)*GRAPHVIZ_POINTSCALE  << "!\"]";

	if(n.t == BASE)
		out << "[shape = triangle, color = blue]";
	else if(n.t == STATIC)
		out << "[shape = circle, color = black]";
	else if(n.t == RELAY)
		out << "[shape = box, color = red]";
	out << "[width = " << GRAPHVIZ_NODEWIDTH << "]";
	out << "[label = " << n.id << "]";
}
