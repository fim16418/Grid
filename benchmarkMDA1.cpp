/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: benchmarkMDA1.cpp

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
#include <iostream>
#include <fstream>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int nLoops;
std::vector<int> latt_size(4);
std::vector<int> mpi_layout(4);
int nThreads;
std::string outFileName;

double average(double* array, int len)
{
  double av = 0.0;
  for(int i=0; i<len; i++) {
    av += array[i];
  }
  return av/len;
}

bool processCmdLineArgs(int argc,char** argv)
{
  nLoops = 1000;
  nThreads = omp_get_max_threads();
  outFileName = "output.txt";
  mpi_layout = {1,1,1,1};

  for(int i=1; i<argc; i++) {
    std::string option = std::string(argv[i]);
    if(option == "--nLoops") {
      if(i+1 < argc) {
        nLoops = atoi(argv[++i]);
      } else {
        std::cerr << "--nLoops option requires one argument." << std::endl;
        return false;
      }
    } else if(option == "--nThreads") {
      if(i+1 < argc) {
        nThreads = atoi(argv[++i]);
        GridThread::SetThreads(nThreads);
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
  std::cout << "Loops per measurement = " << nLoops << std::endl
            << "Threads = " << omp_get_max_threads() << std::endl
            << "Mpi Layout = " << mpi_layout[0] << " " << mpi_layout[1] << " " << mpi_layout[2] << " " << mpi_layout[3] << std::endl
            << "Output file = " << outFileName << std::endl << std::endl;
  return true;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  if(!processCmdLineArgs(argc,argv)) {
    Grid_finalize();
    return 1;
  }

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG rng(&Grid);
  rng.SeedRandomDevice();

  int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  LatticePropagator p1(&Grid); random(rng,p1);
  LatticePropagator p2(&Grid); random(rng,p2);

  LatticeColourMatrix colMat1[Ns*Ns](&Grid);
  LatticeColourMatrix colMat2[Ns*Ns](&Grid);

  for(int s1=0; s1<Ns; s1++) {
  for(int s2=0; s2<Ns; s2++) {
    colMat1[s1*Ns+s2] = peekSpin(p1,s1,s2);
    colMat2[s1*Ns+s2] = peekSpin(p2,s1,s2);
  }}

  LatticeColourMatrix tmp(&Grid);
  LatticeComplex mda[Ns*Ns*Ns*Ns](&Grid);

  for(int s1=0; s1<Ns; s1++) {
  for(int s2=0; s2<Ns; s2++) {
  for(int s3=0; s3<Ns; s3++) {
  for(int s4=0; s4<Ns; s4++) {
    //mult(tmp,colMat1[s1*Ns+s2],colMat2[s3*Ns+s4]);
    tmp = colMat1[s1*Ns+s2] * colMat2[s3*Ns+s4];
    mda[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = traceColour(tmp); //TensorIndexRecursion<ColourIndex>::traceIndex(tmp);
  }}}}

  Grid_finalize();
}
