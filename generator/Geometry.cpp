	/*
	 * =====================================================================================
	 *
	 *       Filename:  Geometry.cpp
	 *
	 *    Description:  
	 *
	 *        Version:  1.0
	 *        Created:  06/16/2010 04:31:03 PM
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
#include  "Plotter.h"
#include "graph.h"
#include "network.h"

	using namespace boost::geometry;


	bool equal_points(point_2d& p1, point_2d &p2)
	{

	  return (fabs(p1.y() -p2.y()) < EPSILON) && (fabs(p1.x() - p2.x()) < EPSILON);
	}


/** Returns true if point p is inside box
 *  p1 ----- x
 *  |        |
 *  x ----- p2 */

bool pointInBox(point_2d &p,
		point_2d &p1,
		point_2d &p2){


	return ((p1.x() <= p.x() && p.x() <= p2.x()) && (p1.y() <= p.y() && p.y() <= p2.y()));
}







//  Determines the intersection point of the line defined by points A and B with the
//  line defined by points C and D.
//
//  Returns 1 if the intersection point was found, and stores that point in X,Y.
//  Returns 0 if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

bool lineIntersection(linestring_2d &l1,
		      linestring_2d &l2,
		      double *X,
		      double *Y){
double Ax = l1[0].x();
       
double Ay = l1[0].y();

double Bx = l1[1].x();
double By = l1[1].y();

double Cx = l2[0].x();
double Cy = l2[0].y();

double Dx = l2[1].x();
double Dy = l2[1].y();

  double  distAB, theCos, theSin, newX, ABpos ;

  //  Fail if either line is undefined.
  if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy) return 0;

  //  (1) Translate the system so that point A is on the origin.
  Bx-=Ax; By-=Ay;
  Cx-=Ax; Cy-=Ay;
  Dx-=Ax; Dy-=Ay;

  //  Discover the length of segment A-B.
  distAB=sqrt(Bx*Bx+By*By);

  //  (2) Rotate the system so that point B is on the positive X axis.
  theCos=Bx/distAB;
  theSin=By/distAB;
  newX=Cx*theCos+Cy*theSin;
  Cy  =Cy*theCos-Cx*theSin; Cx=newX;
  newX=Dx*theCos+Dy*theSin;
  Dy  =Dy*theCos-Dx*theSin; Dx=newX;

  //  Fail if the lines are parallel.
  if (Cy==Dy) return 0;

  //  (3) Discover the position of the intersection point along line A-B.
  ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

  //  (4) Apply the discovered position to line A-B in the original coordinate system.
  *X=Ax+ABpos*theCos;
  *Y=Ay+ABpos*theSin;

  //  Success.
  return 1; 

}

	bool less_eq_points(const point_2d &p1, const point_2d &p2){
		if(p1.x() < p2.x())
			return true;
		if(p1.x() > p2.x())
			return false;
		return (p1.y() < p2.y());
	}


	

	double dist(const point_2d &p1, const point_2d &p2){
		return sqrt(pow(p1.x() - p2.x(),2) + pow(p1.y() - p2.y(),2));

	}
Network *convex_hull_restriction(Network *net){

  Network *net_hull = new Network();
  linestring_2d points;
  //	point_2d *points = new point_2d[net->numStatic()];

  for(int i=0; i < net->numStatic();i++){
    Node &n = net->getStatic(i);
    net_hull->addNode(n.x, n.y, STATIC,"");
    net_hull->setDemand(i,net->getDemand(i));
    append(points,make<point_2d>(n.x, n.y));
  }
  for(int i=0; i< net->numBases(); i++){

    Node &n = net->getBase(i);
    net_hull->addNode(n.x,n.y,BASE,"");
    append(points,make<point_2d>(n.x,n.y));
  }

  //	std::cout << dsv(points) << std::endl;

  polygon_2d hull;
  convex_hull(points,hull);

  //	cout << "Hull" << endl;
  //	cout << dsv(hull) << endl;


  /*  Create new network - removing relays outside hull of static nodes */
  net_hull->setRange(net->getRange());
  net_hull->setRegionSize(net->dimX, net->dimY);

  for(int i=0;i<net->numRelays();i++){
    Node &n = net->getRelay(i);
    point_2d p;
    assign(p,n.x,n.y);
    if(within(p,hull)){
      net_hull->addNode(n.x,n.y,RELAY,n.id);
    }
  }
  return net_hull;



}

	//
	//
	//


	polygon_2d *los_band(point_2d a, point_2d b, double band){
		point_2d p;
		point_2d q;
		if(a.x() <= b.x()){
			assign(p,a.x(),a.y());
			assign(q,b.x(),b.y());
		}else{
			assign(p,b.x(),b.y());
			assign(q,a.x(),a.y());
		}

		point_2d p1,p2,p3,p4;
		if(q.x() == p.x()){ // vertical straight line case
			if(q.y() > p.y()){
				assign(p1, p.x() - band/2, p.y());
				assign(p2, p.x() - band/2, q.y());
				assign(p3, q.x() + band/2, q.y());
				assign(p4, q.x() + band/2, p.y());
		}else{
			assign(p1, p.x() - band/2, q.y());
			assign(p2, p.x() - band/2, p.y());
			assign(p3, p.x() + band/2, p.y());
			assign(p4, p.x() + band/2, q.y());
		}
	}else if(q.y() == p.y()){ //horizontal line case
			assign(p1, p.x(), p.y() - band/2);
			assign(p2, p.x(), p.y() + band/2);
			assign(p3, q.x(), q.y() + band/2);
			assign(p4, q.x(), q.y() - band/2);

	}else{ // general cases
		if(q.y() > p.y()){ //positive slope
			double alpha = atan((q.y() - p.y())/(q.x() - p.x()));
			double dy = cos(alpha) * (band/2.0);
			double dx = sin(alpha) * (band/2.0);
			assign(p1, p.x() + dx, p.y() - dy);
			assign(p2, p.x() - dx, p.y() + dy);
			assign(p3, q.x() - dx, q.y() + dy);
			assign(p4, q.x() + dx, q.y() - dy);
		}else{
			double alpha = atan((q.y() - p.y())/(q.x() - p.x()));
			if(alpha < 0)
				alpha *= -1.0;
			else{
				printf("upss: angle positive when should be negative!\n");
			}
			double dy = cos(alpha) * band/2.0;
			double dx = sin(alpha) * band/2.0;
			assign(p1, q.x() - dx, q.y() - dy);
			assign(p2, p.x() - dx, p.y() - dy);
			assign(p3, p.x() + dx, p.y() + dy);
			assign(p4, q.x() + dx, q.y() + dy);
		}


	}
	polygon_2d *ret = new polygon_2d();
	append(*ret, p1);
	append(*ret, p2);
	append(*ret, p3);

	append(*ret, p4);
	append(*ret, p1);
	correct(*ret);
	//cout << dsv(*ret) << endl;
	return ret;

}

Network *los_band_restriction(Network *net, double band, double rad){
	Network *net_los = new Network();

	vector<point_2d> points;
//	point_2d *points = new point_2d[net->numStatic()];

	for(int i=0; i < net->numStatic();i++){
		Node &n = net->getStatic(i);
		net_los->addNode(n.x, n.y, STATIC,"");
		points.push_back(make<point_2d>(n.x, n.y));
	}
	for(int i=0; i< net->numBases(); i++){

		Node &n = net->getBase(i);
		net_los->addNode(n.x,n.y,BASE,"");
		points.push_back(make<point_2d>(n.x,n.y));
	}
	/*  Create new network - removing relays outside hull of static nodes */
	net_los->setRange(net->getRange());
	net_los->setRegionSize(net->dimX, net->dimY);

	/* Plot resulting areas */
#ifdef _COMPILE_CAIROM
	Plotter *plot = net_los->plot("plot_los.svg");

	linestring_2d relaysToUse;
	for(int i=0;i<points.size();i++){
		for(int j=i+1;j<points.size();j++){
			if(dist(points[i],points[j]) < rad){
				polygon_2d * poly = 
					los_band(points[i],points[j], band);

				/*  plot resulting polygons */
				linear_ring<point_2d> &outer = poly->outer();
				for(int i=0;i<4;i++){
					plot->draw_line(outer[i].x(),outer[i].y(),
							outer[i+1].x(),outer[i+1].y());
				}


				

				for(int k = 0; k < net->numRelays();k++){
					Node &n = net->getRelay(k);
					point_2d p;
					assign(p,n.x,n.y);
					if(within(p,*poly)){
						append(relaysToUse,p);
					}
				}
			}
		}
	}
	plot->write();
	sort(relaysToUse.begin(),relaysToUse.end(),less_eq_points);
	//cout << dsv(relaysToUse) << endl;
	linestring_2d::iterator p_end = unique(relaysToUse.begin(),relaysToUse.end(), equal_points);

	//relaysToUse.resize(it - relaysToUse.begin());
	for(linestring_2d::iterator it = relaysToUse.begin(); it != p_end ;it++){
		net_los->addNode(it->x(),it->y(),RELAY,"");
	}
	

#endif



	return net_los;
	

}
void test_geometry(){

	los_band(make<point_2d>(2.0, 3.0), make<point_2d>(10.0,9.0), 0.5);
	los_band(make<point_2d>(2.0, 9.0), make<point_2d>(10.0,3.0), 0.5);
	los_band(make<point_2d>(2.0, 9.0), make<point_2d>(10.0,9.0), 0.5);
	los_band(make<point_2d>(2.0, 2.0), make<point_2d>(2.0,10.0), 0.5);
}

