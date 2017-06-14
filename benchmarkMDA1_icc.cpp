/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: benchmarkMDA1_icc.cpp

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
#include <vector>

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
  std::cout << "Loops per measurement = " << nLoops << std::endl
            << "Threads = " << omp_get_max_threads() << std::endl
            << "Lattice = " << latt_size[0] << " " << latt_size[1] << " " << latt_size[2] << " " << latt_size[3] << std::endl
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
  //rng.SeedRandomDevice();
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

  //LatticeColourMatrix colMat1[Ns*Ns](&Grid);
  //LatticeColourMatrix colMat2[Ns*Ns](&Grid);

  // Work-around for LatticeColourMatrix colMat1[Ns*Ns](&Grid)
  void* raw_memory = operator new[](Ns*Ns * sizeof(LatticeColourMatrix(&Grid)));
  LatticeColourMatrix* colMat1 = static_cast<LatticeColourMatrix*>( raw_memory );
  for(int i=0; i<Ns*Ns; i++) new( &colMat1[i] )LatticeColourMatrix(&Grid);

  // Work-around for LatticeColourMatrix colMat2[Ns*Ns](&Grid)
  void* raw_memory2 = operator new[](Ns*Ns * sizeof(LatticeColourMatrix(&Grid)));
  LatticeColourMatrix* colMat2 = static_cast<LatticeColourMatrix*>( raw_memory2 );
  for(int i=0; i<Ns*Ns; i++) new( &colMat2[i] )LatticeColourMatrix(&Grid);

  LatticeColourMatrix tmp(&Grid);
  //LatticeComplex mda[Ns*Ns*Ns*Ns](&Grid);

  // Work-around for LatticeComplex mda[Ns*Ns*Ns*Ns](&Grid)
  void* raw_memory3 = operator new[](Ns*Ns*Ns*Ns * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* mda = static_cast<LatticeComplex*>( raw_memory3 );
  for(int i=0; i<Ns*Ns*Ns*Ns; i++) new( &mda[i] )LatticeComplex(&Grid);

  for(int s1=0; s1<Ns; s1++) {
  for(int s2=0; s2<Ns; s2++) {
    colMat1[s1*Ns+s2] = peekSpin(p1,s1,s2);
    colMat2[s1*Ns+s2] = peekSpin(p2,s1,s2);
  }}

  double start = usecond();

  for(int i=0; i<nLoops; i++) {

    for(int s1=0; s1<Ns; s1++) {
    for(int s2=0; s2<Ns; s2++) {
    for(int s3=0; s3<Ns; s3++) {
    for(int s4=0; s4<Ns; s4++) {
      tmp = colMat1[s1*Ns+s2] * colMat2[s3*Ns+s4];
      mda[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = trace(tmp);
    }}}}
  }

  double stop = usecond();
  double time = (stop-start)/1000000.0;

  Grid.Barrier();

  double sumTime;
  MPI_Reduce(&time,&sumTime,1,MPI_DOUBLE,MPI_SUM,Grid.BossRank(),Grid.communicator);
  int nProc = Grid.ProcessorCount();
  time = sumTime/nProc;

  unsigned long flopsPerLoop = (Nc*Nc*16+4)*Ns*Ns*Ns*Ns; // vol placed below
  double flops = flopsPerLoop/1000000000.0*vol*nLoops;

  if(Grid.IsBoss()) {
    std::cout << "mda[0] = " << mda[0]._odata[0] << std::endl; // check the result

    ofstream file;
    file.open(outFileName,ios::app);
    if(file.is_open()) {
      file << nThreads << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
           << vol << "\t" << time << "\t" << flops/time << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file!" << std::endl;
    }
  }

  for(int i=Ns*Ns-1; i>=0; i--) {
    colMat1[i].~LatticeColourMatrix();
    colMat2[i].~LatticeColourMatrix();
  }
  for(int i=Ns*Ns*Ns*Ns-1; i>=0; i--) {
    mda[i].~LatticeComplex();
  }
  operator delete[]( raw_memory);
  operator delete[]( raw_memory2);
  operator delete[]( raw_memory3);

  Grid_finalize();
}
