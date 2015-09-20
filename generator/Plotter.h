/*
 * =====================================================================================
 *
 *       Filename:  Plotter.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/22/2010 05:08:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifdef _COMPILE_CAIROM
#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#endif

class Plotter{
	public:
#ifdef _COMPILE_CAIROM	
	Plotter(string _fn, double dx, double dy);
	Cairo::RefPtr<Cairo::SvgSurface> surface;
	Cairo::RefPtr<Cairo::Context> cr; 
	double height, width,dimX, dimY;
	string svg_filename;

	void draw_point(double ,double);
	void draw_point(double x, double y, double alpha);

	void draw_line(double, double, double,double);
	double scale_x(double);
	double scale_y(double);
	void write();
#else
	Plotter(string _fn, double dx, double dy){}
#endif
};
