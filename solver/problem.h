///
///      @file  problem.h
///     @brief  
///
/// Abstract class "Problem" to be inherited from for defining particular problems to
/// be solved 
///
///    @author  Eduardo Feo (), eduardo@idsia.ch
///
///  @internal
///    Created  06/20/2011
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  IDSIA (Dalle Molle Institute for Artificial Intelligence)
///  Copyright  Copyright (c) 2011, Eduardo Feo
///
/// This program is free software; you can redistribute it and/or modify it 
/// under the terms of the GNU General Public License as published by the 
/// Free Software Foundation; either version 2 of the License, or (at your option) 
/// any later version. This program is distributed in the hope that it will be useful, 
/// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
/// or FITNESS FOR A PARTICULAR PURPOSE. 
/// See the GNU General Public License for more details at  
/// http://www.gnu.org/copyleft/gpl.html

///=====================================================================================
///




#ifndef PROBLEM_H
#define PROBLEM_H

using namespace std;

class Solution;
class LpSolution;


/**
 * @class Problem
 * @brief Abstract class to define MILP problems
 * 
 */
class Problem{
	public:
	Problem(){}
	virtual bool writeToDataFile(ofstream &)=0;
	virtual Solution* fromLpsol(LpSolution *) = 0;
};
#endif
