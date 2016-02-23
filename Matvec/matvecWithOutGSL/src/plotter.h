/*
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Library General Public
 License as published by the Free Software Foundation; either
 version 2 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Library General Public License for more details.
 
 You should have received a copy of the GNU Library General Public
 License along with this library; if not, write to the Free
 Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 MA 02111-1307, USA 
 */

#ifndef MATVEC_PLOTTER_H
#define MATVEC_PLOTTER_H

#include <cstdio>
#include <iostream>
#include <string>
#include "doublematrix.h"
namespace matvec {
/*!
   \class Plotter  plotter.h
   \brief an object to plot curves
*/
class Plotter
{
   public:
      Plotter(const std::string &pp = "/usr/local/bin/gnuplot");
      Plotter(const Plotter &p);
      ~Plotter(void){close();}

      enum plot_cmd {plt, replt};

      int       plot(const std::string &cmd,Plotter::plot_cmd pc = Plotter::plt);
      int       plot3D(const std::string &cmd);
      int       plot(const doubleMatrix& dat,Plotter::plot_cmd pc = Plotter::plt,
                     const std::string &options = "title 'y' with lines");
      int       plot3D(const doubleMatrix& dat,
                     const std::string &options = "title 'z' with lines");
      int       replot(void);
      int       save(const std::string &fname, const std::string &device = "postscript");
      int       set(const std::string &para);
      int       display(void) const {std::cout<<"\tan object of Plotter\n";return 1;}
      int       open(void);
      void      close(void);

   protected:
      FILE*     plot_pipe;
      std::string    plot_program,data_fname;
};
}
#endif
