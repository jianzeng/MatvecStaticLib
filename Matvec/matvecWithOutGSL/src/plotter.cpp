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

#include "session.h"
#include "util.h"
#include "plotter.h"
#include "matrix.h"

namespace matvec {

Plotter::Plotter(const std::string &pp)
{
    plot_program = pp;
    plot_pipe = 0;
    data_fname = SESSION.mktemp();
}

Plotter::Plotter(const Plotter &p)
{
   plot_program = p.plot_program;
   plot_pipe    = p.plot_pipe;
   data_fname   = p.data_fname;
}

int Plotter::plot(const std::string &cmd, Plotter::plot_cmd pc)
{
   this->open();
   if (pc == Plotter::plt) {
      fprintf(plot_pipe,"plot %s\n",cmd.c_str());
      fflush(plot_pipe);
   }
   else if (pc == Plotter::replt) {
      fprintf(plot_pipe,"replot %s\n",cmd.c_str());
      fflush(plot_pipe);
   }
   return 1;
}

int Plotter::plot3D(const std::string &cmd)
{
   this->open();
   fprintf(plot_pipe,"splot %s\n",cmd.c_str());
   fflush(plot_pipe);
   return 1;
}

int Plotter::plot(const doubleMatrix& dat,Plotter::plot_cmd pc,const std::string &options)
{
   if (pc == Plotter::replt) { data_fname = SESSION.mktemp();}
   dat.save(data_fname,std::ios::out);
   this->open();
   if (dat.num_cols() >= 1) {
      if (pc == Plotter::plt) {
         fprintf(plot_pipe,"plot \"%s\" %s\n",data_fname.c_str(),options.c_str());
      }
      else if (pc == Plotter::replt) {
         fprintf(plot_pipe,"replot \"%s\" %s\n",data_fname.c_str(),options.c_str());
      }
   }
   else {
      warning("Plotter::plot(data): data is empty");
   }
   fflush(plot_pipe);
   return 1;
}

int Plotter::plot3D(const doubleMatrix& dat,const std::string &options)
{
   dat.save(data_fname.c_str(),std::ios::out);
   this->open();
   if (dat.num_cols() >= 3) {
      fprintf(plot_pipe,"splot \"%s\" %s\n",data_fname.c_str(),options.c_str());
   }
   else {
      warning("Plotter::plot(data): data needs at least three columns");
   }
   fflush(plot_pipe);
   return 1;
}

int Plotter::replot( void)
{
   if (plot_pipe) {
      fprintf(plot_pipe,"set terminal x11\n");
      fprintf(plot_pipe,"replot\n");  fflush(plot_pipe);
   }
   else {
      warning("Plotter::replot(): no  previous plotting has been found");
   }
   return 1;
}

int Plotter::set(const std::string &settings)
{
   if (!plot_pipe) this->open();
   fprintf(plot_pipe,"set %s\n",settings.c_str());
   fflush(plot_pipe);
   return 1;
}

int Plotter::save(const std::string &fname,const std::string &device)
{
   fprintf(plot_pipe,"set terminal %s\n",device.c_str());
   fprintf(plot_pipe,"set output \"%s\" \n",fname.c_str());
   fprintf(plot_pipe,"replot\n");
   fflush(plot_pipe);
   return 1;
}

int Plotter::open(void)
{
   if (! plot_pipe) plot_pipe = popen(plot_program.c_str(),"w");
   return 1;
}

void Plotter::close(void)
{
   if (plot_pipe) {
      pclose(plot_pipe);
      plot_pipe = 0;
   }
}
}

