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
#include "generator.h"

#include "network.h"

/*std:stream& operator<<(std:stream& s, const Node& n){
    s << n.id;
    
    }
    */

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
void Network::addNode(Node &n){
	Node nn(n);
	 stringstream ss;//create a stringstream
	if(n.t == STATIC){
		ss << "s" <<  s_ix.size();
		s_ix.push_back(nodes.size());
	}else if(n.t == BASE){
		ss << "b" <<  b_ix.size();
		b_ix.push_back(nodes.size());
	}else if(n.t == RELAY){
		ss << "r" <<  r_ix.size();
		r_ix.push_back(nodes.size());
	}
	nn.id = ss.str();
	nodes.push_back(nn);


}

Node &Network::getNodeById(string id){
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
int Network::size(){
	return nodes.size();
}



