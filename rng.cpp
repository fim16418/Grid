/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: rng.cpp

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

int main(int argc, char **argv)
/* Demo for the use of the
 * Random Number Generator (RNG)
 * on a LatticePropagator */
{
    Grid_init(&argc,&argv);

    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);

    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers({1,2,3,4}); //fixed seed
    //pRNG.SeedRandomDevice(); //random seed

    LatticePropagator prop(&Grid);
    random(pRNG,prop);

    std::cout << prop[0] << std::endl;

    Grid_finalize();
}
