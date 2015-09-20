/*
 * =====================================================================================
 *
 *       Filename:  graphsample.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/01/2010 04:42:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
enum family
{ jeanie, debbie, rick, john, amanda, margaret, benjamin, n };
	int
main()
{
	using namespace boost;
	const char *name[] = { "jeanie", "debbie", "rick", "john", "amanda",
		"margaret", "benjamin"
	};

	adjacency_list <> g(n);
	add_edge(jeanie, debbie, g);
	add_edge(jeanie, rick, g);
	add_edge(jeanie, john, g);
	add_edge(debbie, amanda, g);
	add_edge(rick, margaret, g);
	add_edge(john, benjamin, g);

	graph_traits < adjacency_list <> >::vertex_iterator i, end;
	graph_traits < adjacency_list <> >::adjacency_iterator ai, a_end;
	property_map < adjacency_list <>, vertex_index_t >::type
		index_map = get(vertex_index, g);

	for (tie(i, end) = vertices(g); i != end; ++i) {
		std::cout << name[get(index_map, *i)];
		tie(ai, a_end) = adjacent_vertices(*i, g);
		if (ai == a_end)
			std::cout << " has no children";
		else
			std::cout << " is the parent of ";
		for (; ai != a_end; ++ai) {
			std::cout << name[get(index_map, *ai)];
			if (boost::next(ai) != a_end)
				std::cout << ", ";
		}
		std::cout << std::endl;
	}
	return exit_success;
}

}
}
}
