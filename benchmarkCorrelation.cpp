/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: benchmarkCorrelation.cpp

Author: Moritz Fink <fink.moritz@gmail.com>

This program uses the following library:



Grid physics library, www.github.com/paboyle/Grid

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>


using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int nData;
int nLoops;
std::vector<int> latt_size;
int nThreads;

bool overlapComms = false;


double average(double* array, int len)
{
  double av = 0.0;
  for(int i=0; i<len; i++) {
    av += array[i];
  }
  return av/len;
}

double standardDeviation(double* array, int len)
{
  double squares = 0.0;
  for(int i=0; i<len; i++) {
    squares += array[i]*array[i];
  }
  squares /= len;

  double av = average(array,len);
  double var = squares - av*av;

  return Grid::sqrt(var);
}

bool processCmdLineArgs(int argc, char ** argv)
{
  nData  = 2;
  nLoops = 10;
  latt_size = {4,4,4,8};
  nThreads = omp_get_max_threads();

  for(int i=1; i<argc; i++) {
    std::string option = std::string(argv[i]);
    if(option == "--lattice") {
      if(i+5 == argc) { //--lattice must be last option
        for(int j=0; j<4; j++) {
          latt_size[j] = atoi(argv[i+j+1]);
        }
        i+=4;
      } else {
        std::cerr << "--lattice x y z t must be the last option." << std::endl;
        return false;
      }
    }
    else if(option == "--nData") {
      if(i+1 < argc) {
        nData = atoi(argv[++i]);
      } else {
        std::cerr << "--nData option requires one argument." << std::endl;
        return false;
      }
    }
    else if(option == "--nLoops") {
      if(i+1 < argc) {
        nLoops = atoi(argv[++i]);
      } else {
        std::cerr << "--nLoops option requires one argument." << std::endl;
        return false;
      }
    }
    else if(option == "--nThreads") {
      if(i+1 < argc) {
        nThreads = atoi(argv[++i]);
        omp_set_num_threads(nThreads);
      } else {
        std::cerr << "--nThreads option requires one argument." << std::endl;
        return false;
      }
    }
  }
  std::cout << "Lattice = " << latt_size[0] << " " << latt_size[1] << " " << latt_size[2] << " " << latt_size[3] << std::endl
            << "Measurements = " << nData << std::endl
            << "Loops per measurement = " << nLoops << std::endl
            << "Threads = " << omp_get_max_threads() << std::endl << std::endl;
  return true;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  if( GridCmdOptionExists(argv,argv+argc,"--asynch") ){
    overlapComms = true;
  }

  if(!processCmdLineArgs(argc,argv)) return 1;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers({1,2,3,4});
  //pRNG.SeedRandomDevice();

  Gamma gamma5(Gamma::Gamma5);
  LatticePropagator quark_propagator(&Grid);
  random(pRNG,quark_propagator);

  for(int i0=0; i0<quark_propagator._odata.size(); i0++) {
  for(int i1=0; i1<Ns; i1++) {
  for(int i2=0; i2<Ns; i2++) {
  for(int i3=0; i3<Nc; i3++) {
  for(int i4=0; i4<Nc; i4++) {
    quark_propagator._odata[i0]._internal._internal[i1][i2]._internal[i3][i4] = 3.14;
  }}}}}

  LatticePropagator anti_quark = gamma5 * quark_propagator * gamma5;
  anti_quark = adj(anti_quark);

  double timeData[nData];
  double flopData[nData];

  int latticeSize = latt_size[0] * latt_size[1] * latt_size[2] * latt_size[3];
  unsigned long flopsPerLoop = 2 * 30262 * latticeSize;

  for(int j=0; j<nData; j++) {
    double start = usecond();

    for(int i=0; i<nLoops; i++) {
      LatticeComplex corr_fn = trace(anti_quark * gamma5 * quark_propagator * gamma5);
      //std::cout << corr_fn[0] << std::endl; break; //for test purposes
    }

    double stop = usecond();

    timeData[j] = (stop-start)/1000000.0;
    flopData[j] = flopsPerLoop/timeData[j]/1000000000.0*nLoops;
  }

  std::cout << endl << "time for " << nLoops << " loops = "
            << average(timeData,nData) << " +/- " << standardDeviation(timeData,nData) << " secs" << endl
            << "GFlops/second = " << average(flopData,nData) << " +/- " << standardDeviation(flopData,nData) << endl
            << "(Average over " << nData << " measurements)" << endl;

  Grid_finalize();
}
