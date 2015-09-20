/*
 * =====================================================================================
 *
 *       Filename:  Plotter.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/22/2010 05:08:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#define POINT_RES 5.0
#define DIM_RES 500

#include "main.h"
#include "Plotter.h"
#ifdef _COMPILE_CAIROM
double Plotter::scale_x(double x){
	return ((x/dimX) * width);
}

double Plotter::scale_y(double y){
	return height - ((y/dimY) * height);
}

Plotter::Plotter(string fname, double dx, double dy){
	svg_filename = fname;
	dimX = dx;
	dimY = dy;
	/* Scale dimensions */
	if(dx > dy){
		width = DIM_RES;
		height = DIM_RES * (dy/dx);
	}else{
		height = DIM_RES;
		width = DIM_RES * (dx/dy);
	}

	Cairo::RefPtr<Cairo::SvgSurface> sf = 
		Cairo::SvgSurface::create(svg_filename, width, height);
	surface = Cairo::RefPtr<Cairo::SvgSurface>::cast_static(sf);
	
	//Cairo::RefPtr<Cairo::Context> 
	cr = Cairo::Context::create(surface);
		
	cr->save(); // save the state of the context
	cr->set_source_rgb(1.00, 1.00, 1.00);
	cr->paint();    // fill image with the color
	cr->restore();  // color is back to black now

	cr->save();
	// draw a border around the image
	//cr->set_line_width(20.0);    // make the line wider
	//cr->rectangle(0.0, 0.0, cairo_image_surface_get_width(surface->cobj()), height);
	cr->stroke();

	cr->set_source_rgba(0.0, 0.0, 0.0, 0.7);
//	cr->show_page();
	surface->flush();
	//surface->finish();
}

void Plotter::draw_point(double x, double y, double alpha){

  //Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    cr->arc(scale_x(x), scale_y(y), 
            alpha*POINT_RES, 0.0, 2.0 * M_PI);
    cr->stroke(); 
}


void Plotter::draw_point(double x, double y){

  //Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    cr->arc(scale_x(x), scale_y(y), 
            POINT_RES, 0.0, 2.0 * M_PI);
    cr->stroke(); 
}

void Plotter::draw_line(double x1, double y1, double x2, double y2){
	//Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);

  cr->move_to (scale_x(x1), scale_y(y1));
    cr->line_to (scale_x(x2), scale_y(y2));
    
    cr->stroke();
    cr->save();

}

void Plotter::write(){
//	Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
  //  cr->show_page();
   // cr->flush();
   // cr->finish();
  surface->flush();
surface->finish();
    std::cout << "Wrote SVG file \"" << svg_filename << "\"" << std::endl;

}
#endif
