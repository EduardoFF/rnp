/*
 * =====================================================================================
 *
 *       Filename:  heuristic.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/08/2010 06:32:52 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef _HEURISTIC
#define _HEURISTIC

class LpSolution;
class Heuristic{
	public:
		virtual LpSolution *operator()() = 0;

};

//struct lp_solution_t;

class Low_res_heuristic : public Heuristic{
	public:
		LpSolution *hvalue;
		Low_res_heuristic(Network *, double , int, bool silent =0);
		Low_res_heuristic(LpSolution*, double, int, bool silent = 0);
		LpSolution *update_solution(LpSolution* cnt_sol, map<int,int> index_map);
		LpSolution *operator()();
		
};
#endif
