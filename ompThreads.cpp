/*************************************************************************************

Grid examples, www.github.com/fim16418/Grid

Copyright (C) 2017

Source code: ompThreads.cpp

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
/* Test of the thread
 * manipulation via
 * omp functions */
{
    Grid_init(&argc,&argv);

    std::cout << "#threads=" << omp_get_num_threads() << std::endl;
    //#threads=1 because this is a sequential section

    std::cout << "#maxThreads=" << omp_get_max_threads() << std::endl;
    omp_set_num_threads(1);
    std::cout << "#maxThreads=" << omp_get_max_threads() << std::endl;

    omp_set_num_threads(2);
    std::cout << "#maxThreads=" << omp_get_max_threads() << std::endl;

    Grid_finalize();
}
