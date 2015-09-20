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


enum nodeType { STATIC, BASE, RELAY};


class Node{
public:

	Node(double _x, double _y) : x(_x), y(_y){}
	nodeType t;
	string id;
	double x, y;

};

//std:stream& operator<<(std:stream& s, const Node& n);
class Network{
	public:
		vector<Node> nodes;
		vector<int> s_ix;
		vector<int> b_ix;
		vector<int> r_ix;

		int numStatic(){ return s_ix.size();}
		int numBases(){ return b_ix.size();}
		int numRelays(){ return r_ix.size();}

		Node &getNodeById(string id);

		Node &getStatic(int i);
		Node &getBase(int i);
		Node &getRelay(int i);

		int size();
		Node &getNode(int i);
		void addNode(Node &n);
};
#endif
