/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: mdaFull_icc.cc

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

#define MYNAMESPACE TIMING // or MDA

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

//int nLoops;
std::vector<int> latt_size(4);
std::vector<int> mpi_layout(4);
int nThreads;
std::string outFileName;
int nProps;

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
  //nLoops = 1000;
  nThreads = omp_get_max_threads();
  outFileName = "output.txt";
  mpi_layout = {1,1,1,1};
  latt_size = {8,8,8,8};
  nProps = 4;

  for(int i=1; i<argc; i++) {
    std::string option = std::string(argv[i]);
    /*if(option == "--nLoops") {
      if(i+1 < argc) {
        nLoops = atoi(argv[++i]);
      } else {
        std::cerr << "--nLoops option requires one argument." << std::endl;
        return false;
      }
    } else */if(option == "--nThreads") {
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
    } else if(option == "--nProps") {
      if(i+1 < argc) {
        nProps = atoi(argv[++i]);
      } else {
        std::cerr << "--nProps option requires one argument." << std::endl;
        return false;
      }
    }
  }
  std::cout /*<< "Loops per measurement = " << nLoops << std::endl*/
            << "Threads = " << omp_get_max_threads() << std::endl
            << "Lattice = " << latt_size[0] << " " << latt_size[1] << " " << latt_size[2] << " " << latt_size[3] << std::endl
            << "Mpi Layout = " << mpi_layout[0] << " " << mpi_layout[1] << " " << mpi_layout[2] << " " << mpi_layout[3] << std::endl
            << "Output file = " << outFileName << std::endl
            << "Number of propagators = " << nProps << std::endl << std::endl;
  return true;
}

namespace MDA {

  inline void derivative(const LatticePropagator& prop, const LatticeGaugeField& u, int dir, int len, LatticePropagator& ret)
  {
    //ret = 0.5 * (Cshift(prop,dir,len) - Cshift(prop,dir,-len));
    LatticeColourMatrix u_mu = PeekIndex<LorentzIndex>(u,dir);
    LatticePropagator tmp = adj(u_mu)*prop;
    ret = u_mu*Cshift(prop,dir,len) - Cshift(tmp,dir,-len);
  }

  inline void arrangeData(const LatticePropagator& prop, LatticeComplex* data)
  {
    LatticeSpinMatrix sMat(prop._grid);

    for(int c1=0; c1<Nc; c1++) {
    for(int c2=0; c2<Nc; c2++) {
      sMat = peekColour(prop,c1,c2);

      for(int s1=0; s1<Ns; s1++) {
      for(int s2=0; s2<Ns; s2++) {
        data[c1*Nc*Ns*Ns+c2*Ns*Ns+s1*Ns+s2] = peekSpin(sMat,s1,s2);
      }}
    }}
  }

  inline void mda(LatticeComplex* a, LatticeComplex* b, LatticeComplex* ret)
  {
    for(int s1=0; s1<Ns; s1++) {
    for(int s2=0; s2<Ns; s2++) {
    for(int s3=0; s3<Ns; s3++) {
    for(int s4=0; s4<Ns; s4++) {

      ret[s1*Ns*Ns*Ns+s2*Ns*Ns+s3*Ns+s4] = a[0*Ns*Ns+s1*Ns+s2] * b[0*Ns*Ns+s3*Ns+s4] +
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

}

namespace TIMING {

  inline void derivative(const LatticePropagator& prop, const LatticeGaugeField& u, int dir, int len, LatticePropagator& ret)
  {
    double start = usecond();

    MDA::derivative(prop,u,dir,len,ret);

    double stop = usecond();
    std::cout << std::endl << "derivative time = "
              << (stop - start)/1000000.0 << " secs" << std::endl;
  }

  inline void arrangeData(const LatticePropagator& prop, LatticeComplex* data)
  {
    double start = usecond();

    MDA::arrangeData(prop,data);

    double stop = usecond();
    std::cout << std::endl << "arrangement time = "
              << (stop - start)/1000000.0 << " secs" << std::endl;
  }

  inline void mda(LatticeComplex* a, LatticeComplex* b, LatticeComplex* ret)
  {
    double start = usecond();

    MDA::mda(a,b,ret);

    double stop = usecond();
    std::cout << std::endl << "computation time = "
              << (stop - start)/1000000.0 << " secs" << std::endl;
  }

}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  if(!processCmdLineArgs(argc,argv)) {
    Grid_finalize();
    return 1;
  }

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  GridCartesian    Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG rng(&Grid);
  //rng.SeedRandomDevice();
  rng.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  int derivativeLen = 1;

  //LatticePropagator props[nProps](&Grid);
  // Work-around for LatticeColourMatrix colMat1[Ns*Ns](&Grid)
  void* raw_memory = operator new[](nProps * sizeof(LatticePropagator(&Grid)));
  LatticePropagator* props = static_cast<LatticePropagator*>( raw_memory );
  for(int i=0; i<nProps; i++) new( &props[i] )LatticePropagator(&Grid);

  for(int i=0; i<nProps; i++) {
    random(rng,props[i]);
  }

  LatticeGaugeField U(&Grid);
  random(rng,U);

  LatticePropagator dProp1(props[0]._grid);
  LatticePropagator dProp2(props[0]._grid);

  // Work-around for LatticeComplex data1[Nc*Nc*Ns*Ns](props[0]._grid);
  void* raw_memory2 = operator new[](Nc*Nc*Ns*Ns * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* data1 = static_cast<LatticeComplex*>( raw_memory2 );
  for(int i=0; i<Nc*Nc*Ns*Ns; i++) new( &data1[i] )LatticeComplex(&Grid);

  // Work-around for LatticeComplex data2[Nc*Nc*Ns*Ns](props[0]._grid);
  void* raw_memory3 = operator new[](Nc*Nc*Ns*Ns * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* data2 = static_cast<LatticeComplex*>( raw_memory3 );
  for(int i=0; i<Nc*Nc*Ns*Ns; i++) new( &data2[i] )LatticeComplex(&Grid);

  // Work-around for LatticeComplex resultBuffer[Ns*Ns*Ns*Ns](props[0]._grid);
  void* raw_memory4 = operator new[](Ns*Ns*Ns*Ns * sizeof(LatticeComplex(&Grid)));
  LatticeComplex* resultBuffer = static_cast<LatticeComplex*>( raw_memory4 );
  for(int i=0; i<Ns*Ns*Ns*Ns; i++) new( &resultBuffer[i] )LatticeComplex(&Grid);

  double timerStart = usecond();

  for(int p1=0; p1<nProps; p1++) {

    const LatticePropagator& prop1 = props[p1];

    for(int mu=0; mu<Nd; mu++) {

      MYNAMESPACE::derivative(prop1,U,mu,derivativeLen,dProp1);

      MYNAMESPACE::arrangeData(dProp1,data1);

      for(int p2=0; p2<nProps; p2++) {

        const LatticePropagator& prop2 = props[p2];

        MYNAMESPACE::arrangeData(prop2,data2);

        MYNAMESPACE::mda(data1,data2,resultBuffer);

      }

      for(int nu=0; nu<Nd; nu++) {

        if(mu == nu) continue;

        MYNAMESPACE::derivative(dProp1,U,nu,derivativeLen,dProp2);

        MYNAMESPACE::arrangeData(dProp2,data1);

        for(int p2=0; p2<nProps; p2++) {

          const LatticePropagator& prop2 = props[p2];

          MYNAMESPACE::arrangeData(prop2,data2);

          MYNAMESPACE::mda(data1,data2,resultBuffer);

        }
      }
    }

    for(int p2=0; p2<nProps; p2++) {

      const LatticePropagator& prop2 = props[p2];

      for(int mu=0; mu<Nd; mu++) {

        MYNAMESPACE::derivative(prop1,U,mu,derivativeLen,dProp1);

        MYNAMESPACE::arrangeData(dProp1,data1);

        for(int nu=0; nu<Nd; nu++) {

          if(mu >= nu) continue;

          MYNAMESPACE::derivative(prop2,U,nu,derivativeLen,dProp2);

          MYNAMESPACE::arrangeData(dProp2,data2);

          MYNAMESPACE::mda(data1,data2,resultBuffer);

        }
      }
    }
  }

  double timerStop = usecond();
  double timerTime = (timerStop - timerStart) / 1000000.0;

  if(Grid.IsBoss()) {
    ofstream file;
    file.open(outFileName,ios::app);
    if(file.is_open()) {
      file << nThreads << "\t" << latt_size[0] << latt_size[1] << latt_size[2] << latt_size[3] << "\t"
           << vol << "\t" << timerTime << std::endl;
      file.close();
    } else {
      std::cerr << "Unable to open file!" << std::endl;
    }
  }

  for(int i=nProps-1; i>=0; i--) {
    props[i].~LatticePropagator();
  }
  for(int i=Ns*Ns*Nc*Nc-1; i>=0; i--) {
    data1[i].~LatticeComplex();
    data2[i].~LatticeComplex();
  }
  for(int i=Ns*Ns*Ns*Ns-1; i>=0; i--) {
    resultBuffer[i].~LatticeComplex();
  }
  operator delete[]( raw_memory );
  operator delete[]( raw_memory2 );
  operator delete[]( raw_memory3 );
  operator delete[]( raw_memory4 );

  Grid_finalize();
}
