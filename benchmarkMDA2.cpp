/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: benchmarkMDA2.cpp

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

#define WARM_UP 10

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int nLoops;
std::vector<int> latt_size(4);
std::vector<int> mpi_layout(4);
int nThreads;
std::string outFileName;


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

bool processCmdLineArgs(int argc,char** argv)
{
  nLoops = 1000;
  nThreads = omp_get_max_threads();
  outFileName = "output.txt";
  mpi_layout = {1,1,1,1};
  latt_size = {8,8,8,8};

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
    } else if(option == "--lattice") {
      if(i+4 < argc) {
        for(int j=0; j<4; j++) {
          latt_size[j] = atoi(argv[i+j+1]);
        }
        i+=4;
      } else {
        std::cerr << "--lattice option requires four arguments." << std::endl;
        return false;
      }
    }
  }
  std::cout << "Loops = " << nLoops << std::endl
            << "Threads = " << omp_get_max_threads() << std::endl
            << "Lattice = " << latt_size[0] << " " << latt_size[1] << " " << latt_size[2] << " " << latt_size[3] << std::endl
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

  /*//////////////////
  // Initialization //
  //////////////////*/

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG rng(&Grid);
  rng.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  LatticePropagator p1(&Grid); //random(rng,p1);
  LatticePropagator p2(&Grid); //random(rng,p2);

  for(int x=0; x<p1._odata.size(); x++) {
  for(int s1=0; s1<Ns; s1++) {
  for(int s2=0; s2<Ns; s2++) {
  for(int c1=0; c1<Nc; c1++) {
  for(int c2=0; c2<Nc; c2++) {
    p1._odata[x]._internal._internal[s1][s2]._internal[c1][c2] = x*10000 + s1*1000 + s2*100 + c1*10 + c2;
    p2._odata[x]._internal._internal[s1][s2]._internal[c1][c2] = x*10000 + s1*1000 + s2*100 + c1*10 + c2;
  }}}}}

  LatticeSpinMatrix sMat1(&Grid);
  LatticeSpinMatrix sMat2(&Grid);

  // Work-around for LatticeComplex mda[Ns*Ns*Ns*Ns](&Grid)
  void* raw_memory = operator new[](Ns*Ns*Ns*Ns * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* mda = static_cast<LatticeComplex*>( raw_memory );
  for(int i=0; i<Ns*Ns*Ns*Ns; i++) new( &mda[i] )LatticeComplex(&Grid);
 
  // Work-around for LatticeComplex a[Ns*Ns*Nc*Nc](&Grid)
  void* raw_memory2 = operator new[](Ns*Ns*Nc*Nc * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* a = static_cast<LatticeComplex*>( raw_memory2 );
  for(int i=0; i<Ns*Ns*Nc*Nc; i++) new( &a[i] )LatticeComplex(&Grid);

  // Work-around for LatticeComplex b[Ns*Ns*Nc*Nc](&Grid)
  void* raw_memory3 = operator new[](Ns*Ns*Nc*Nc * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* b = static_cast<LatticeComplex*>( raw_memory3 );
  for(int i=0; i<Ns*Ns*Nc*Nc; i++) new( &b[i] )LatticeComplex(&Grid);

  /*///////////////
  // Preparation //
  //   Warm up   //
  ///////////////*/

  for(int i=0; i<WARM_UP; i++) {
    for(int c1=0; c1<Nc; c1++) {
    for(int c2=0; c2<Nc; c2++) {
      sMat1 = peekColour(p1,c1,c2);
      sMat2 = peekColour(p2,c1,c2);

      for(int s1=0; s1<Ns; s1++) {
      for(int s2=0; s2<Ns; s2++) {
        a[c1*Nc*Ns*Ns+c2*Ns*Ns+s1*Ns+s2] = peekSpin(sMat1,s1,s2);
        b[c1*Nc*Ns*Ns+c2*Ns*Ns+s1*Ns+s2] = peekSpin(sMat2,s1,s2);
      }}
    }}
  }

  /*///////////////
  // Preparation //
  // Measurement //
  ///////////////*/

  double timePrep[nLoops];
  double start, stop;

  for(int i=0; i<nLoops; i++) {
    start = usecond();

    for(int c1=0; c1<Nc; c1++) {
    for(int c2=0; c2<Nc; c2++) {
      sMat1 = peekColour(p1,c1,c2);
      sMat2 = peekColour(p2,c1,c2);

      for(int s1=0; s1<Ns; s1++) {
      for(int s2=0; s2<Ns; s2++) {
        a[c1*Nc*Ns*Ns+c2*Ns*Ns+s1*Ns+s2] = peekSpin(sMat1,s1,s2);
        b[c1*Nc*Ns*Ns+c2*Ns*Ns+s1*Ns+s2] = peekSpin(sMat2,s1,s2);
      }}
    }}

    stop = usecond();
    timePrep[i] = stop-start;
  }

  /*///////////////
  // Calculation //
  //   Warm up   //
  ///////////////*/

  for(int i=0; i<WARM_UP; i++) {
    
  #pragma omp parallel for collapse(4)
    for(int s1=0; s1<Ns; s1++) {
    for(int s2=0; s2<Ns; s2++) {
    for(int s3=0; s3<Ns; s3++) {
    for(int s4=0; s4<Ns; s4++) {

      mda[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = a[0*Ns*Ns+s1*Ns+s2] * b[0*Ns*Ns+s3*Ns+s4] +
                                           a[3*Ns*Ns+s1*Ns+s2] * b[1*Ns*Ns+s3*Ns+s4] +
                                           a[6*Ns*Ns+s1*Ns+s2] * b[2*Ns*Ns+s3*Ns+s4] +
                                           a[1*Ns*Ns+s1*Ns+s2] * b[3*Ns*Ns+s3*Ns+s4] +
                                           a[4*Ns*Ns+s1*Ns+s2] * b[4*Ns*Ns+s3*Ns+s4] +
                                           a[7*Ns*Ns+s1*Ns+s2] * b[5*Ns*Ns+s3*Ns+s4] +
                                           a[2*Ns*Ns+s1*Ns+s2] * b[6*Ns*Ns+s3*Ns+s4] +
                                           a[5*Ns*Ns+s1*Ns+s2] * b[7*Ns*Ns+s3*Ns+s4] +
                                           a[8*Ns*Ns+s1*Ns+s2] * b[8*Ns*Ns+s3*Ns+s4];
    }}}}
  }

  /*///////////////
  // Calculation //
  // Measurement //
  ///////////////*/

  double timeComp[nLoops];

  for(int i=0; i<nLoops; i++) {
    start = usecond();

#pragma omp parallel for collapse(4)
    for(int s1=0; s1<Ns; s1++) {
    for(int s2=0; s2<Ns; s2++) {
    for(int s3=0; s3<Ns; s3++) {
    for(int s4=0; s4<Ns; s4++) {

      mda[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = a[0*Ns*Ns+s1*Ns+s2] * b[0*Ns*Ns+s3*Ns+s4] +
                                           a[3*Ns*Ns+s1*Ns+s2] * b[1*Ns*Ns+s3*Ns+s4] +
                                           a[6*Ns*Ns+s1*Ns+s2] * b[2*Ns*Ns+s3*Ns+s4] +
                                           a[1*Ns*Ns+s1*Ns+s2] * b[3*Ns*Ns+s3*Ns+s4] +
                                           a[4*Ns*Ns+s1*Ns+s2] * b[4*Ns*Ns+s3*Ns+s4] +
                                           a[7*Ns*Ns+s1*Ns+s2] * b[5*Ns*Ns+s3*Ns+s4] +
                                           a[2*Ns*Ns+s1*Ns+s2] * b[6*Ns*Ns+s3*Ns+s4] +
                                           a[5*Ns*Ns+s1*Ns+s2] * b[7*Ns*Ns+s3*Ns+s4] +
                                           a[8*Ns*Ns+s1*Ns+s2] * b[8*Ns*Ns+s3*Ns+s4];
    }}}}

    stop = usecond();
    timeComp[i] = stop-start;
  }

  /*//////////////
  // Evaluation //
  //////////////*/

  unsigned long flopsPerLoop = (9*4+8*2)*Ns*Ns*Ns*Ns; // vol placed below
  double flops = flopsPerLoop/1000000000.0*vol;

  double tPrep, tPrepError, tComp, tCompError;
  error(timePrep,nLoops,tPrep,tPrepError);
  error(timeComp,nLoops,tComp,tCompError);

  tPrep /= 1000000.0;
  tPrepError /= 1000000.0;
  tComp /= 1000000.0;
  tCompError /= 1000000.0;

  double flopsPerSec = flops/tComp;
  double flopsPerSecError = tCompError/tComp * flopsPerSec;

  /*/////////////////
  // Print results //
  /////////////////*/

  if(Grid.IsBoss()) {
    std::cout << "mda[0] = " << mda[0]._odata[0] << std::endl;

    ofstream file;
    file.open(outFileName,ios::app);
    if(file.is_open()) {
      file << omp_get_max_threads() << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
           << vol << "\t" << tPrep << "\t" << tPrepError << "\t" << tComp << "\t" << tCompError << "\t"
           << flopsPerSec << "\t" << flopsPerSecError << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file!" << std::endl;
    }
  }
  
  /*///////////////
  // Destructors //
  ///////////////*/

  for(int i=Ns*Ns*Ns*Ns-1; i>=0; i--) {
    mda[i].~LatticeComplex();
  }
  for(int i=Ns*Ns*Nc*Nc-1; i>=0; i--) {
    a[i].~LatticeComplex();
    b[i].~LatticeComplex();
  }
  operator delete[]( raw_memory );
  operator delete[]( raw_memory2 );
  operator delete[]( raw_memory3 );

  Grid_finalize();
}
