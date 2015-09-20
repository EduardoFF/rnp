/*
 * =====================================================================================
 *
 *       Filename:  wrapper.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/11/2010 05:56:35 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

typedef struct{
	int grb_out;
	int glp_out;
	int grb_cuts;
	int grb_rootmethod;
} wrapper_params;




double solve_glp_grb(glp_prob *mip, wrapper_params *);


