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
  rng.SeedRandomDevice();

  int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  LatticePropagator p1(&Grid); random(rng,p1);
  LatticePropagator p2(&Grid); random(rng,p2);

  LatticeSpinMatrix sMat1[Nc*Nc](&Grid);
  LatticeSpinMatrix sMat2[Nc*Nc](&Grid);

  LatticeComplex mda[Ns*Ns*Ns*Ns](&Grid);

  for(int c1=0; c1<Nc; c1++) {
  for(int c2=0; c2<Nc; c2++) {
    sMat1[c1*Nc+c2] = peekColour(p1,c1,c2);
    sMat2[c1*Nc+c2] = peekColour(p2,c1,c2);
  }}

  double start = usecond();

  for(int i=0; i<nLoops; i++) {

    for(int s1=0; s1<Ns; s1++) {
    for(int s2=0; s2<Ns; s2++) {
    for(int s3=0; s3<Ns; s3++) {
    for(int s4=0; s4<Ns; s4++) {

      LatticeComplex a[Nc*Nc](&Grid);
      LatticeComplex b[Nc*Nc](&Grid);

      for(int color=0; color<Nc*Nc; color++) {
        a[color] = peekSpin(sMat1[color],s1,s2);
        b[color] = peekSpin(sMat2[color],s3,s4);
      }

      mda[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2] + a[1]*b[3] + a[4]*b[4] +
                                           a[7]*b[5] + a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
    }}}}
  }

  double stop = usecond();
  double time = (stop-start)/1000000.0;

  Grid.Barrier();

  double sumTime;
  MPI_Reduce(&time,&sumTime,1,MPI_DOUBLE,MPI_SUM,Grid.BossRank(),Grid.communicator);
  int nProc = Grid.ProcessorCount();
  time = sumTime/nProc;

  unsigned long flopsPerLoop = (Nc*Nc*16+4)*Ns*Ns*Ns*Ns;
  double flops = flopsPerLoop/1000000000.0*vol*nLoops; //FIXME check me

  if(Grid.IsBoss()) {
    ofstream file;
    file.open(outFileName,ios::app);
    if(file.is_open()) {
      file << nThreads << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
           << time << "\t" << flops/time << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file!" << std::endl;
    }
  }

  Grid_finalize();
}
