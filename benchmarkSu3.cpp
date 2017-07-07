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

#define WARM_UP 10

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int nLoops;
std::vector<int> latt_size;
std::vector<int> mpi_layout(4);
int nThreads;
std::string outFileName;

bool overlapComms = false;


void error(double* array, int len, double& average, double& error)
{
  average = 0.0;
  double square = 0.0;

  for(int i=0; i<len; i++) {
    average += array[i];
    square += array[i]*array[i];
  }

  average = average/len;
  square = square/len;

  error = std::sqrt(square - average*average);
  error /= std::sqrt(len);
}

bool processCmdLineArgs(int argc, char ** argv)
{
  nLoops = 1000;
  latt_size = {4,4,4,4};
  nThreads = omp_get_max_threads();
  outFileName = "output.txt";
  mpi_layout = {1,1,1,1};

  for(int i=1; i<argc; i++) {
    std::string option = std::string(argv[i]);
    if(option == "--lattice") {
      if(i+5 == argc) { //--lattice must be last argument
        for(int j=0; j<4; j++) {
          latt_size[j] = atoi(argv[i+j+1]);
        }
        i+=4;
      } else {
        std::cerr << "--lattice x y z t must be the last option." << std::endl;
        return false;
      }
    } else if(option == "--nLoops") {
        if(i+1 < argc) {
          nLoops = atoi(argv[++i]);
        } else {
          std::cerr << "--nLoops option requires one argument." << std::endl;
          return false;
        }
      } else if(option == "--nThreads") {
        if(i+1 < argc) {
          nThreads = atoi(argv[++i]);
          omp_set_num_threads(nThreads);
        } else {
          std::cerr << "--nThreads option requires one argument." << std::endl;
          return false;
        }
      } else if(option == "--outFile") {
        if(i+1 < argc) {
          outFileName = argv[++i];
        } else {
          std::cerr << "--outFile option requires one argument." << std::endl;
          return false;
        }
      } else if(option == "--mpiLayout") {
        if(i+4 < argc) {
          for(int j=0; j<4; j++) {
            mpi_layout[j] = atoi(argv[i+j+1]);
          }
        i+=4;
        } else {
          std::cerr << "--mpiLayout option requires four arguments." << std::endl;
          return false;
        }
      }
    }
    std::cout << "Lattice = " << latt_size[0] << " " << latt_size[1] << " " << latt_size[2] << " " << latt_size[3] << std::endl
              << "Loops = " << nLoops << std::endl
              << "Threads = " << omp_get_max_threads() << std::endl
              << "Mpi Layout = " << mpi_layout[0] << " " << mpi_layout[1] << " " << mpi_layout[2] << " " << mpi_layout[3] << std::endl
              << "Output file = " << outFileName << std::endl << std::endl;
    return true;
  }


  int main (int argc, char ** argv)
  {
    Grid_init(&argc,&argv);

    if( GridCmdOptionExists(argv,argv+argc,"--asynch") ){
      overlapComms = true;
    }

    if(!processCmdLineArgs(argc,argv)) {
       Grid_finalize();
       return 1;
    }

    /*//////////////////
    // Initialization //
    //////////////////*/

    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);

    LatticePropagator prop1(&Grid);
    LatticePropagator prop2(&Grid);

    for(int i0=0; i0<quark_propagator._odata.size(); i0++) {
    for(int i1=0; i1<Ns; i1++) {
    for(int i2=0; i2<Ns; i2++) {
    for(int i3=0; i3<Nc; i3++) {
    for(int i4=0; i4<Nc; i4++) {
      prop1._odata[i0]._internal._internal[i1][i2]._internal[i3][i4] = 3.14;
      prop2._odata[i0]._internal._internal[i1][i2]._internal[i3][i4] = 2.5;
    }}}}}

    LatticeComplex corr(&Grid);

    /*///////////////
    // Calculation //
    //   Warm up   //
    ///////////////*/

    for(int i=0; i<WARM_UP; i++) {
      corr = trace(prop1 * prop2);
    }

    /*///////////////
    // Calculation //
    // Measurement //
    ///////////////*/

    double timeData[nLoops];

    for(int i=0; i<nLoops; i++) {
      double start = usecond();

      corr = trace(prop1 * prop2);

      double stop = usecond();
      timeData[i] = stop-start;
    }

    /*//////////////
    // Evaluation //
    //////////////*/

    double time, timeError;
    error(timeData,nLoops,time,timeError);

    time /= 1000000.0;
    timeError /= 1000000.0;

    int vol = latt_size[0] * latt_size[1] * latt_size[2] * latt_size[3];
    unsigned long flopsPerLoop = 2 * (10080 + 22);
    double flops = flopsPerLoop/1000000000.0*vol;

    double flopsPerSec = flops/time;
    double flopsPerSecError = timeError/time * flopsPerSec;

    /*/////////////////
    // Print results //
    /////////////////*/

    if(Grid.IsBoss()) {
      ofstream file;
      file.open(outFileName,ios::app);
      if(file.is_open()) {
        file << nThreads << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
             << vol << "\t" << time << "\t" << timeError << "\t" << flopsPerSec << "\t" << flopsPerSecError << std::endl;
        file.close();
      } else {
          std::cerr << "Unable to open file!" << std::endl;
        }
      }

      Grid_finalize();
    }
