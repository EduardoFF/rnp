#include "problem.h"





Problem::Problem(Grid *g){
	grid = g;
}


int man_distance(locationid_t l1, locationid_t l2){
	return (abs(l1.first-l2.first) + abs(l1.second - l2.second));
}
double Problem::energy_cost(locationid_t lid, agent_t a){
	return 1;
}

/* Initially, same for all agents */
bool Problem::adyacent_loc(locationid_t l1, locationid_t l2, agent_t a){
	if(abs(l1.first - l2.first) <=1 && abs(l1.second - l2.second) <= 1)
		return true;
	else
		return false;

}


double Problem::exploration(locationid_t lid, agent_t ag){
	double v = grid->getVegetation(lid);
	double s = grid->getSlope(lid);
	double a = grid->getAltitude(lid);

	double res;
	if(ag->type == HUMAN){
		/**  */
		res = (1.0/( pow(v,3) * (a < ALTITUDE_TH ? (1.0):(1.0+a/ALTITUDE_TH))*(s<1.0?(1.0):(s*s))));
	}else if(ag->type == DOG){
		res = (1.0/( pow(v,1) * (a < ALTITUDE_TH ? (1.0):sqrt(1.0+a/ALTITUDE_TH))*(s<1.0?(1.0):(s))));

	}

	return res;
}

bool Problem::connected_loc(locationid_t l1, agent_t a1, locationid_t l2, agent_t a2){
	return (man_distance(l1,l2) < 5);
}
/*
string Problem::toPlainText(){
	stringstream s;
	
}

void Problem::generateMILPData(){
	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << toPlainText();
	myfile.close();
	

}

*/
