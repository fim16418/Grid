/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: benchmarkDerivative.cpp

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
std::vector<int> latt_size(4);
std::vector<int> mpi_layout(4);
int nThreads;
int mu;
int length;
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
  mu = 0;
  length = 1;
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
    } else if(option == "--mu") {
      if(i+1 < argc) {
        mu = atoi(argv[++i]);
      } else {
        std::cerr << "--mu option requires one argument." << std::endl;
        return false;
      }
    } else if(option == "--length") {
      if(i+1 < argc) {
        length = atoi(argv[++i]);
      } else {
        std::cerr << "--length option requires one argument." << std::endl;
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
            << "Derivative in direction " << mu << " with length " << length << std::endl
            << "Mpi Layout = " << mpi_layout[0] << " " << mpi_layout[1] << " " << mpi_layout[2] << " " << mpi_layout[3] << std::endl
            << "Output file = " << outFileName << std::endl << std::endl;
  return true;
}


int main (int argc, char ** argv)
{
  if(!processCmdLineArgs(argc,argv)) {
    return 1;
  }
  
  Grid_init(&argc,&argv);

  if( GridCmdOptionExists(argv,argv+argc,"--asynch") ){
    overlapComms = true;
  }

  /*//////////////////
  // Initialization //
  //////////////////*/

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(latt_size,simd_layout,mpi_layout);
  std::vector<LatticeColourMatrix> U(Nd,UGrid);

  LatticeGaugeField Umu(&Grid);
  for(int i=0; i<Umu._odata.size(); i++) {
  for(int dir=0; dir<Nd; dir++) {
  for(int c1=0; c1<Nc; c1++) {
  for(int c2=0; c2<Nc; c2++) {
    Umu._odata[i]._internal[dir]._internal._internal[c1][c2] = dir+c1+c2;
  }}}}

  for(int dir=0; dir<Nd; dir++) {
    U[dir] = PeekIndex<LorentzIndex>(Umu,dir);
  }

  Gamma gamma5(Gamma::Algebra::Gamma5);
  LatticePropagator quark_propagator(&Grid);

  for(int index=0; index<quark_propagator._odata.size(); index++) {
  for(int s1=0; s1<Ns; s1++) {
  for(int s2=0; s2<Ns; s2++) {
  for(int c1=0; c1<Nc; c1++) {
  for(int c2=0; c2<Nc; c2++) {
    quark_propagator._odata[index]._internal._internal[s1][s2]._internal[c1][c2] = index+s1+s2+c1+c2;
  }}}}}

  LatticePropagator anti_quark = gamma5 * quark_propagator * gamma5;
  anti_quark = adj(anti_quark);

  LatticeColourMatrix gField = U[mu];
  LatticePropagator tmp(&Grid);

  LatticeComplex corr(&Grid);

  /*///////////////
  // Calculation //
  //   Warm up   //
  ///////////////*/

  for(int i=0; i<WARM_UP; i++) {
    tmp = adj(gField) * quark_propagator;
    corr = trace(anti_quark * gamma5 * (gField*Cshift(quark_propagator,mu,length) - Cshift(tmp,mu,-length))*gamma5);
  }

  /*///////////////
  // Calculation //
  // Measurement //
  ///////////////*/

  double timeData[nLoops];

  for(int i=0; i<nLoops; i++) {
    double start = usecond();

    tmp = adj(gField) * quark_propagator;
    corr = trace(anti_quark * gamma5 * (gField*Cshift(quark_propagator,mu,length) - Cshift(tmp,mu,-length))*gamma5);
    //std::cout << corr[0] << std::endl; break; //for test purposes

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
  unsigned long flopsPerLoop = 2 * (5*10080 + 22);
  double flops = flopsPerLoop/1000000000.0*vol;

  double flopsPerSec = flops/time;
  double flopsPerSec_error = timeError/time * flopsPerSec;

  /*/////////////////
  // Print results //
  /////////////////*/

  if(Grid.IsBoss()) {
    ofstream file;
    file.open(outFileName,ios::app);
    if(file.is_open()) {
      file << omp_get_max_threads() << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
           << vol << "\t" << time << "\t" << timeError << "\t" << flopsPerSec << "\t" << flopsPerSec_error << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file!" << std::endl;
    }
  }

  Grid_finalize();
}
