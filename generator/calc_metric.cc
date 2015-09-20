/*
 * =====================================================================================
 *
 *       Filename:  calc_metric.cc
 *
 *    Description:  Calculates the numerical estimation of PRR based on SINR Physical Model and MAC model
 *
 *        Version:  1.0
 *        Created:  10/22/2010 03:35:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eduardo Feo
 *        Company:  IDSIA
 *
 * =====================================================================================
 */

#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>



using namespace std;

#define MAX_NODES 500

map<int,pair <double,double > > nodes;
//vector< pair < int, double > > distances;
map< pair<int, int>, double > all_pair_dist;

map< int, double> pmac_values;
double **pactive = NULL;
map< pair<int, int>, double> p_active;
double sigma;
double ple;
double demand;
double ttx;

double link_dist;
#define PLE  3.0
#define SIGMA  3.2
#define M 0.230258509 // Log[10]/10
#define P0 -55.0 // Reference power -55 dBm
#define SINRTH 6.0 // SINR threshold
#define MACTH -72 // MAC Threshold
#define NOISE -105.0
#define is_bit_set(X,i) X & (1 << i)
#define LAMBDA 0.05
#define TTX 2.533
#define MAC_R  10.0
#define NET_R  10.0
#define INT_R  15.0


#define USE_PACTIVE 1
#define USE_PMACS 0
#define USE_PMACS_AT_RECEIVE 1
#define USE_MAC_PACTIVE 0
#define UP_LIMIT 16

static double muZ(vector<pair< int, double> > &dist, unsigned int mask, double sZ);
static double sigmaZ(vector<pair <int, double> > &dist,unsigned int mask);
static double probZge(vector< pair<int, double> > &dist, unsigned int mask,double TH);
double prob_ij(int i, int j);
static double muW(double);
static double sigmaW(double);


struct sort_pred {
	    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
		            return left.second < right.second;
			        }
};
#define RCV_RANGE 10.0

double PNOTRCV(vector< pair<int, double> > &dist, int mask, int sender){
	double q = demand*ttx;
	double ret = 0.0;
	bitset<32> bs(mask);
	for(int i=0;i<dist.size();i++){
		if(dist[i].second > RCV_RANGE){
			continue;
		}
		double oneminusq = 1.0;
		for(int j=0;j<dist.size();j++){
			if(dist[j].second > RCV_RANGE)
				continue;
			if(i == j)
				continue;
			//double coef = p_active[make_pair(dist[i].first,dist[j].first)];
			double coef = pactive[dist[i].first][dist[j].first];
			oneminusq *= (1 - coef*q);
		}
		//double sender_pmac = p_active[make_pair(sender,dist[i].first)];
		double sender_pmac = pactive[sender][dist[i].first];
		double pmacs = pmac_values[dist[i].first];
		if(!USE_PMACS_AT_RECEIVE){
			pmacs = 1.0;
		}

		ret += sender_pmac*pmacs*q*oneminusq;
	}
	return 1 - ret;
}

double PRR(vector< pair<int, double> > &dist,int sender){
	sort(dist.begin(),dist.end(),sort_pred());
	if( dist.size() > UP_LIMIT){
		while( dist.size() > UP_LIMIT){
			dist.pop_back();
		}
	}
	double q = demand*ttx;
	double ret = 0.0;
	register unsigned i;
	for(i=(1 << dist.size())-1;i!=0;i--){	
		double pb = 1.0 - probZge(dist,i,SINRTH);
		bitset<32> bs(i);

		int len = bs.count();
	
		/*  Improvement: Probability that contenders transmit at same time
		 *  is equal to probability that any pair of them stops transmitting due to CCA
		 *   Lets call this probability p_activated
		 *   */
		double node_zero_pmac = 1.0;
		double qside = 1.0;
		double oneminusq = 1.0;
		for(int j = 0;j < dist.size();++j){
			if( bs[j]){ 
				// Both j and k are in contending subset
				//double p = p_active[make_pair(sender,dist[j].first)];
				double p = pactive[sender][dist[j].first];
				node_zero_pmac *= p;
			}

		}

		
		if( USE_PACTIVE){			
			qside = 1.0;
			double coef = 1.0;
			double pmacs = 1.0;
			for(int j=dist.size()-1;j>=0;--j){
				if(bs[j]){
					pmacs *= pmac_values[dist[j].first];
					for(int k=dist.size()-1;k>=0;--k){
						if(k == j)
							continue;
						if(bs[k]){
							//double p = p_active[make_pair(dist[j].first,dist[k].first)];
							double p = pactive[dist[j].first][dist[k].first];
							//p *= p;
							coef *= p;
					}
					}
				}
			}
			if( !USE_PMACS){
				pmacs = 1.0;
			}
			qside *= pmacs*pow(q,len)*sqrt(coef);
			oneminusq = 1.0;
			for(int j=dist.size()-1;j>=0;--j){
				if(!bs[j]){
					double coef=1.0;
					for(int k=dist.size()-1;k>=0;--k){
						if(bs[k]){
							//coef = min(coef,p_active[make_pair(dist[j].first,dist[k].first)]);
							coef = min(coef,pactive[dist[j].first][dist[k].first]);
						}
					}
					oneminusq *= (1-coef*q);
				}
			}

		}else{
			qside = pow(q, len);
			oneminusq = pow(1-q,(int)dist.size()-len);
		}

		double a =  qside*oneminusq*node_zero_pmac*pb;
		double ploss =  a;


		ret += ploss;
	}
	return max(1 - ret,0.0);


}


/*  Probability that node i and j do not share same medium */
double prob_ij(int i, int j){

	double d = all_pair_dist[make_pair(i,j)];
	//cout << "r = " << d;
	double ss = P0 - 10.0*ple*log(d)/log(10.0) -MACTH;
	double r = 1.0 - 0.5*(1 + erf(ss/(sqrt(2)*sigma)));
	//cout << "p: " << r << endl;
	return r;
}



static double muW(double d){
	return log(1.0/pow(d,ple));
}

static double sigmaW(double d){
	return M*M*sigma*sigma;
}

static double muZ(vector< pair<int, double> > &dist, unsigned int mask, double sZ){
	double sum1 = 0.0;
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
			double d = dist[i].second;
			double exponent= log(1.0/(pow(d,ple)));
			sum1 += exp(exponent);
		}
	}
	double b = log(sum1);
	b += M*M*sigma*sigma*0.5;
	b -= sZ*0.5;
	//b += log(pow(10.0,P0/10.0));
	//b += pow(10,NOISE/10.0);
	return b;


}

static double sigmaZ(vector<pair <int, double> > &dist,unsigned int mask){
	double sum1 = 0.0;
//	cout << "calc sigmaZ: ";
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
			double d = dist[i].second;
		//	cout << " " << d;
			double exponent= 2.0*log(1.0/(pow(d,ple)));
			sum1 += exp(exponent);
		}

	}
	//cout << endl;
	double sum2 = 0.0;
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
			double d = dist[i].second;
			double exponent= log(1.0/(pow(d,ple)));
			sum2 += exp(exponent);
		}
	}
	sum2 = sum2*sum2;

	double b = sum1/sum2;
	b *= (exp(M*M*sigma*sigma) - 1);
	b += 1;
	return log(b);

}
double probZge(vector<pair <int, double> > &dist, unsigned int mask,double TH){
	double sz = sigmaZ(dist,mask);
	double mz = muZ(dist,mask,sz);
	//cout << "sz " << sz << " mz " << (1/M)*mz << endl;
	double sw = sigmaW(link_dist);
	double mw = muW(link_dist);
	
	mz += (pow(10,NOISE/10.0)/pow(10.0,P0/10.0));
	double sy = sz +sw;
       	double my = mw - mz;
	sy *= (1/M)*(1/M);
	my *= (1/M);
	double zle = 0.5*(1 + erf((TH - my)/(sqrt(2*sy))));
	return (1.0 - zle);

}


/***************** PMAC ***************************/
double muZ_MAC(vector< pair<int, double> > &dist, unsigned int mask, double sZ){
	double sum1 = 0.0;
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
	//		cout << "bitset " << i;
			double d = dist[i].second;
	//		cout << " " << d;
			double exponent= log(1.0/(pow(d,ple)));
			sum1 += exp(exponent);
		}
	}
	//cout << endl;
	double b = log(sum1);
	b += M*M*sigma*sigma*0.5;
	b -= sZ*0.5;
	b += log(pow(10.0,P0/10.0));
	b += pow(10,NOISE/10.0);
	return b;


}

double sigmaZ_MAC(vector<pair <int, double> > &dist,unsigned int mask){
	double sum1 = 0.0;
//	cout << "calc sigmaZ: ";
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
			double d = dist[i].second;
		//	cout << " " << d;
			double exponent= 2.0*log(1.0/(pow(d,ple)));
			sum1 += exp(exponent);
		}

	}
	//cout << endl;
	double sum2 = 0.0;
	for(int i=0;i<dist.size();i++){
		if(is_bit_set(mask,i)){
			double d = dist[i].second;
			double exponent= log(1.0/(pow(d,ple)));
			sum2 += exp(exponent);
		}
	}
	sum2 = sum2*sum2;

	double b = sum1/sum2;
	b *= (exp(M*M*sigma*sigma) - 1);
	b += 1;
	return log(b);

}

double probZge_MAC(vector<pair <int, double> > &dist, unsigned int mask,double TH){
	double sz = sigmaZ_MAC(dist,mask);
	double mz = muZ_MAC(dist,mask,sz);
	//cout << "sz " << sz << " mz " << (1/M)*mz << endl;
	sz *= (1/M)*(1/M);
	mz *= (1/M);
	double zle = 0.5*(1 + erf((TH - mz)/sqrt(2*sz)));
	return (1.0 - zle);

}


/*  Receives vector with CS and return Pbusy */
double PBCA(vector< pair<int, double> > &dist){
	double q = demand*ttx;
	double ret = 0.0;
	for(unsigned long  i=1;i< (1 << dist.size());i++){
		double pb = probZge_MAC(dist,i,MACTH);
		bitset<32> bs(i);

		int len = bs.count();
	
		/*  Improvement: Probability that contenders transmit at same time
		 *  is equal to probability that any pair of them stops transmitting due to CCA
		 *   Lets call this probability p_activated
		 *   */
		double	qside= 1.0;
		double oneminusq = 1.0;
		
		if( USE_MAC_PACTIVE){
			qside = 1.0;
			/* q side part */
			/* relation between active contenders */
			double qcoef = 1.0;
			for(int j=0;j < dist.size();j++){
				if(bs[j] == 1){
					double coef = 1.0;
					for(int k=j+1;k < dist.size();k++){
						if(bs[k]){
							/* Nodes j and k are both active */
							//double p = p_active[make_pair(dist[j].first,dist[k].first)];
							double p = pactive[dist[j].first][dist[k].first];
							qcoef *= (p);

						}
					}
					if(coef > 1.0 || coef < 0.0)
						printf("ERROR ON CALCULATION\n");
			

				}
			}
			qside = pow(q,len)*qcoef;
			/* one minus q side */
			oneminusq = 1.0;
			for(int j=0;j< dist.size();j++){
				if(!bs[j]){
					double coef = 1.0;
					for(int k=0;k<dist.size();k++){
						if(bs[k]){
							//double p = p_active[make_pair(dist[j].first,dist[k].first)];
							double p = pactive[dist[j].first][dist[k].first];
							coef = min(coef,p);
						}
					}
					if(coef > 1.0 || coef < 0.0)
						printf("ERROR ON CALCULATION\n");
					oneminusq *= (1 - (coef)*q);
				}

			}
		}else{
			qside = pow(q,len);
			oneminusq = pow(1-q,(int)dist.size()-len);
		}
		ret += qside*oneminusq*pb;
	}

	return min(ret,1.0);


}


/*  Receives a vector for CS set and returns PMAC */

double PMAC(vector<pair <int, double> > &dist){

	
	/* Remove exceed */
	sort(dist.begin(),dist.end(),sort_pred());
	if( dist.size() > UP_LIMIT){
		while( dist.size() > UP_LIMIT){
			dist.pop_back();
		}
	}

	double pcca = 1 - PBCA(dist);
	double ret = 0.0;
	for(int i=1;i<=5;i++){
		ret += pcca*pow((1-pcca),i-1);

	}
	return ret;
}

/************** END PMAC ***********************/

double distance(double x1, double y1, double x2, double y2){
	double dx = (x1 - x2)*(x1 - x2);
	double dy = (y1 - y2)*(y1 - y2);
	return sqrt(dx + dy);

}

void calc_metric(vector< pair<double, double> > &net, vector<int> &mask_send, vector<int> &mask_rcv,
		 vector<int> &mask_int, map< pair< int, int>, double> &ret_val,
	       double _demand, double _sigma, double _ple, double _ttx	){
	demand = _demand;
	sigma = _sigma;
	ple = _ple;
	ttx = _ttx;
	
	int ni =0;
	cout << "Calc metric v. 1.0" << endl;
	for(int i=0;i<net.size();i++){
		double x = net[i].first;
		double y = net[i].second;
		nodes[ni++] = make_pair(x,y);
	}
	pactive = new double*[ni];
	for(int i=0;i<ni;i++){
		pactive[i] = new double[ni];
	}

	for( int i=0; i < ni; i++){
		for( int j= i+1; j<ni;j++){
			double d = distance(nodes[i].first,nodes[i].second,nodes[j].first,nodes[j].second);
			all_pair_dist[make_pair(i,j)] = max(d,1.0);	
			double pa = prob_ij(i,j);

			pactive[i][j] = pa;
			pactive[j][i] = pa;
			p_active[make_pair(i,j)] = pa;
			p_active[make_pair(j,i)] = p_active[make_pair(i,j)];
		}
	}
 

	cout << "Debug\n";
	/* Calculate PMAC value for each node */

	for(int i=0;i < mask_send.size(); i++){
		int ii = mask_send[i];
		vector< pair < int, double > > mac_distances;
		for(int j=0; j < mask_int.size();j++){
			int jj = mask_int[j];
			if(ii ==jj)
				continue;
			double d = distance(nodes[ii].first,nodes[ii].second,nodes[jj].first, nodes[jj].second);
			if( d < MAC_R){
				mac_distances.push_back(make_pair(jj,d));
			}
		}

		double pmac = PMAC(mac_distances);
		//printf("pmac value for %d = %f n = %d\n",ii,pmac,mac_distances.size());
		pmac_values[ii] = pmac;
	}


	/* Calculate PPR for each link */
	for(int i=0;i<mask_send.size();i++){
		int ii = mask_send[i];
		for(int j=0;j<mask_rcv.size();j++){
			int jj = mask_rcv[j];
			vector< pair < int, double > > int_distances;
			if(ii == jj)
				continue;
			double d = distance(nodes[ii].first,nodes[ii].second,nodes[jj].first, nodes[jj].second);
			if(d > NET_R)
				continue;
			/* Search for interferer nodes around receiver jj */
			for(int k=0;k<mask_int.size();k++){
				int kk = mask_int[k];
				if(kk == ii || kk == jj)
					continue;
				double di = distance(nodes[kk].first,nodes[kk].second,nodes[jj].first, nodes[jj].second);
				if(di > INT_R)
					continue;
				int_distances.push_back(make_pair(kk,di));
			}
			/* We got interference set */
			link_dist = d; /* Global var defining link distance between i -> j */
			double a = PRR(int_distances,ii);
			//cerr << "link " << i << "-> " << j << " ready" << endl;
			double pint = (1 - a*a);
			double pnotrev = PNOTRCV(int_distances,1<<int_distances.size(),i);
			double ploss = (1-pnotrev) + (pnotrev*pint);
			//printf("%d -> %d  ISET = %d - ld %f - a %f - pint %f - pnotrcv %f - ploss %f\n",i,j,int_distances.size(),link_dist,a,pint,pnotrev,ploss);
			double metric = pmac_values[ii]*(1-ploss);
			ret_val[make_pair(ii,jj)] = metric;
			//cout << filename << " " << i << " " << j << " " << metric << endl;
		}
	}
	for(int i=0;i<ni;i++){
		delete [] pactive[i];
	}
	 delete [] pactive;


}
#ifndef _METRIC
int main2(int argc, char *argv[]){
	string filename(argv[1]);

	demand = LAMBDA;
	sigma = SIGMA;
	ple = PLE;
	ttx = TTX;
	for(int i=2;i<argc;i++){

		if(strcmp(argv[i],"-d")==0)
			demand = atof(argv[i+1]);
		if(strcmp(argv[i],"-s")==0)
			sigma = atof(argv[i+1]);
		if(strcmp(argv[i],"-ple")==0)
			ple = atof(argv[i+1]);
		if(strcmp(argv[i],"-ttx")==0)
			ttx = atof(argv[i+1]);



	}
	string line;
	ifstream myfile(filename.c_str());
	
	int ni = 0; // Node counter
	//cout << "Reading network from : " << filename << endl;
	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			getline (myfile,line);
			istringstream iss(line);

			do
			{
				string sub;

				string id,xs,ys;
				iss >> sub;
				if(sub == "node"){
					iss >> id >> xs >> ys;
					int ii = atoi(id.c_str());
					double x = atof(xs.c_str());
					double y = atof(ys.c_str());
					nodes[ni++] = make_pair(x,y);


				}
			} while (iss);

		}
		myfile.close();
	}

	else cout << "Unable to open file"; 

	for( int i=0; i < ni; i++){
		for( int j= i+1; j<ni;j++){
			double d = distance(nodes[i].first,nodes[i].second,nodes[j].first,nodes[j].second);
			all_pair_dist[make_pair(i,j)] = max(d,1.0);	
			double pa = prob_ij(i,j);
			//printf("pactive %d %d d=%f pa %f\n",i,j,d,pa);

			pactive[i][j] = pa;
			pactive[j][i] = pa;
			p_active[make_pair(i,j)] = pa;
			p_active[make_pair(j,i)] = p_active[make_pair(i,j)];


		}
	}
 

	/* Calculate PMAC value for each node */

	for(int i=0;i < ni; i++){
		vector< pair < int, double > > mac_distances;
		for(int j=0;j<ni;j++){
			if(i ==j)
				continue;
			double d = distance(nodes[i].first,nodes[i].second,nodes[j].first, nodes[j].second);
			if( d < MAC_R){
				mac_distances.push_back(make_pair(j,d));
			}
		}

		double pmac = PMAC(mac_distances);
		printf("pmac value for %d = %f n = %d\n",i,pmac,mac_distances.size());
		pmac_values[i] = pmac;
	}


	/* Calculate PPR for each link */
	for(int i=0;i<ni;i++){
		for(int j=0;j<ni;j++){
			vector< pair < int, double > > int_distances;
			if(i == j)
				continue;
			double d = distance(nodes[i].first,nodes[i].second,nodes[j].first, nodes[j].second);
			if(d > NET_R)
				continue;
			/* Search for interferer nodes around receiver j */
			for(int k=0;k<ni;k++){
				if(k == i || k == j)
					continue;
				double di = distance(nodes[k].first,nodes[k].second,nodes[j].first, nodes[j].second);
				if(di > INT_R)
					continue;
				int_distances.push_back(make_pair(k,di));
			}
			/* We got interference set */
			link_dist = d; /* Global var defining link distance between i -> j */
			double a = PRR(int_distances,i);
			cerr << "link " << i << "-> " << j << " ready" << endl;
			double pint = (1 - a*a);
			double pnotrev = PNOTRCV(int_distances,1<<int_distances.size(),i);
			double ploss = (1-pnotrev) + (pnotrev*pint);
			printf("%d -> %d  ISET = %d - ld %f - a %f - pint %f - pnotrcv %f - ploss %f\n",i,j,int_distances.size(),link_dist,a,pint,pnotrev,ploss);
			double metric = pmac_values[i]*(1-ploss);
			cout << filename << " " << i << " " << j << " " << metric << endl;


		}
	}

/*
	for( int i=0; i < ni; i++){
		for( int j= i+1; j<ni;j++){
			double d = distance(nodes[i].first,nodes[i].second,nodes[j].first,nodes[j].second);
			all_pair_dist[make_pair(i,j)] = max(d,1.0);	
			p_active[make_pair(i,j)] = prob_ij(i,j);
			p_active[make_pair(j,i)] = p_active[make_pair(i,j)];
			//double ss = P0 - 10.0*ple*log(d)/log(10.0) -MACTH - NOISE;
			//cout << "r: " << d << " p_a: " << prob_ij(i,j) << endl; //" ss: " << ss << endl;
		}
	}
	
	  
	//double sz = sigmaZ(distances);
	//double mz = muZ(distances,sz);
	//cout << "sigmaz " << sz << endl;
	//cout << "muz " << mz*1/M << endl;
*/
//	cout << "ProbZGE " << probZge(distances,(1 << distances.size()),SINRTH) << endl;

	
/*
	double q = demand*ttx;
	double a = PRR(distances);
	double pnotrev = PNOTRCV(distances,1<<distances.size());
*/
	//double pnotrev = (1 - distances.size()*q*pow(1-q,distances.size()-1));
//	double pnotrev = 1.0;

	//   Previous 
//	double pint = ((1 -a) + 2*a*(1-a) + (1-a)*(1-a))/2.0;

//	double pint = (1 - a*a);
	//double pint = 1-a;
//	double ploss = (1-pnotrev) + (pnotrev*pint);
	//cout  << filename << "\t" << max((1-ploss),0.0) << endl;
//	cout << "prev: " << (1 - pnotrev) << " pint: " << pint << " 1-prr: " << 1-a << endl;


	return 0;

}
#endif
