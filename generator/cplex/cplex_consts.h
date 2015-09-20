#ifndef __CPX_CONSTS
#define __CPX_CONSTS

#include <stdio.h>
#include <string>
#include <map>

using namespace std;

extern map<std::string, int> cplex_parameters;

void populate_cpx_const();

#endif
