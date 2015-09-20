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

#include "main.h"
#include "generator.h"
//#include "network.h"
#include "graph.h"
#include "boost/graph/connected_components.hpp"

//extern Network net;

extern param prob;

void Graph::addEdge(int i, int j){
	add_edge(i,j,*this);
}
void Graph::set_vertex_id(int i, string _id){
//	property_map<MyGraphType, vertex_name_t>::type
//		vertex_name = get(vertex_name_t(), *this);
	//printf("Setting id[%d] = %s\n",i,id.c_str));
//	boost::put(vertex_name,i,id ); 
	(*this)[i].id = _id;
	

}

void Graph::set_edge_label(int i, int j, string label){
	boost::graph_traits<MyGraphType>::edge_descriptor e1;
	bool found;
	tie(e1, found) = edge(i, j, *this);

	(*this)[e1].label = label;
}
void Graph::set_vertex_type(int i, nodeType _t){
	(*this)[i].t = _t;
}
void Graph::set_vertex_loc(int i, double _x, double _y){
	(*this)[i].x = _x;
	(*this)[i].y = _y;
}


struct graph_writer {
    void operator()(std::ostream& out) const {
      out << "label = \" " << prob.instance_id <<" \""<< endl;
      out << "fontsize = 30" << endl;
      out << "labelloc = \"t\" " << endl;
      out << "overlap = \"scale\"" << endl;
      //out << "ratio = " << prob.dimY / prob.dimY << endl;
     // out << "size = 10" << endl;
    }
  };


int Graph::toNetFile(ofstream &netFile){
	if(!netFile.is_open())
		return -1;
	graph_traits<MyGraphType>::vertex_iterator vi, vi_end, next;
	tie(vi, vi_end) = vertices(*this);
	for (next = vi; next != vi_end; next++) {
		// printf("Vertex id [%d] = %s\n",*next,idmap[*next].c_str());
		netFile << "node\t" << (*this)[*next].id 
			<< "\t" << (*this)[*next].x*prob.dimX
			<< "\t" << (*this)[*next].y*prob.dimY
			<< "\t" << (*this)[*next].t
			<< "\n";

	}



	
	graph_traits<MyGraphType>::edge_iterator ei,e_end,e_next;
	graph_traits<MyGraphType>::edge_descriptor e;
         tie(ei, e_end) = edges( (*this)); 
	 for(e_next = ei; e_next != e_end; e_next++) {
		/* Not necesary - links added in sim script 
		 netFile << "link\t" << source(*e_next, *this) 
			 << "\t" << target(*e_next, *this)
			 << "\n";
	         */
		 netFile << "route\t" << source(*e_next, *this) 
			 << "\t" << target(*e_next, *this)
			 << "\t" << (*this)[*e_next].label
			 << "\n";

	 }
	 return 0;


}

int Graph::toGraphviz(string filename){
	ofstream dotFile (filename.c_str());
	if( !dotFile.is_open()){
		return -1;
	}

	vector<string> ids;
	vector< pair< double, double> > locs;
	vector<nodeType> types;
	map< edge_d,string> edge_labels;

	//	cout << "Num vertices " << num_vertices(*this) << endl;
	//	property_map<MyGraphType, vertex_name_t>::type idmap = get(vertex_name_t(),*this);
	graph_traits<MyGraphType>::vertex_iterator vi, vi_end, next;
	tie(vi, vi_end) = vertices(*this);
	for (next = vi; next != vi_end; next++) {
		// printf("Vertex id [%d] = %s\n",*next,idmap[*next].c_str());
		ids.push_back((*this)[*next].id);
		locs.push_back( make_pair( (*this)[*next].x, (*this)[*next].y));
		types.push_back((*this)[*next].t);

	}
	
	graph_traits<MyGraphType>::edge_iterator ei,e_end,e_next;
	graph_traits<MyGraphType>::edge_descriptor e;
         tie(ei, e_end) = edges( (*this)); 
	 for(e_next = ei; e_next != e_end; e_next++) {
		 edge_labels[*e_next] = (*this)[*e_next].label;
	 }







	write_graphviz(dotFile, *this,make_position_writer(ids,locs,types),make_edge_writer(edge_labels),graph_writer());
	return 0;

}


edge_writer make_edge_writer(map<edge_d, string> _m){
	return edge_writer(_m);
}

position_writer make_position_writer(vector<string> ids, vector< pair< double, double> > locs, vector<nodeType> types){
	return position_writer(ids,locs,types);
}
template <class VertexOrEdge>
void position_writer::operator()(std::ostream& out, const VertexOrEdge& v) const {
    if( !prob.no_position)
    {
      //Node &n = net.getNodeById(idmap[v]);
      double x = locmap[v].first;
      double y = locmap[v].second;
      out << "[pos=\"" << (locmap[v].first)*GRAPHVIZ_POINTSCALE << "," 
        << (locmap[v].second)*GRAPHVIZ_POINTSCALE  << "!\"]";
    }

 	if(typemap[v] == BASE){
		out << "[shape = triangle, color = blue]";
			out << "[label = " << idmap[v] << "]";
			out << "[width = " << GRAPHVIZ_NODEWIDTH << "]";
	}else if(typemap[v] == STATIC){
		out << "[shape = circle, color = black]";
		out << "[label = " << idmap[v] << "]";
		out << "[width = " << GRAPHVIZ_NODEWIDTH << "]";
	}else if(typemap[v] == RELAY){
		out << "[shape = box, color = red]"; 
		out << "[width = " << GRAPHVIZ_NODEWIDTH/2 << "]";
		out << "[style = filled,fillcolor = red]" << endl;
		out << "[label = " << idmap[v] << "]";
	}
}
template <class VertexOrEdge>
void edge_writer::operator()(std::ostream& out, const VertexOrEdge& v) const {
	edge_d e = v;
	map< edge_d, string>::const_iterator it = idmap.find(v);
	//idmap[e] = "a";
	out << "[label = \"" << it->second << "\"]"; 
}

int numberOfConnectedComponents(Graph *G)
{
	typedef std::map<graph_traits<MyGraphType>::vertex_descriptor, 	graph_traits<MyGraphType>::vertices_size_type> component_type;
	component_type component;
	boost::associative_property_map< component_type > component_map(component);

	int num_components = connected_components(*G, component_map);
	return num_components;
}


