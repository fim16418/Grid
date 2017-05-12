# Grid
### C++ code examples for the Grid library (www.github.com/paboyle/Grid)

## Abstract
This repository contains some example programs for the Grid library. This repo addresses people that like to do Lattice Quantum Chromodynamic (LQCD) and want to try out this new library. The examples show how to use Grid's features. They also provide some benchmarks whose corresponding Chroma files can be found in my Chroma repo (www.github.com/fim16418/Chroma) and used for comparison of the two libraries. For the compilation I used the release/v0.6.0 branch of the Grid repository but other branches should work too.

## The files
### rng.cpp
This example shows how to set up an instance of GridCartesian and how to fill it with random numbers from the GridParallelRNG (Random Number Generator).

### gridThreads.cpp & ompThreads.cpp
When it comes to using multiple threads, one can either use the OpenMP (omp.h) or Grid functions. As can be seen in gridThreads.cpp, the number of threads cannot be increased once set (which should not be done anyways). I would suggest users to use the functions from gridThreads.cpp.

### peekPoke.cpp & peekPoke2.cpp
Accessing tensors for read (peek) and write (poke) operations can be done via two different ways. Either the peek & poke functions are used in the way presented in peekPoke.cpp or the elements are accessed directly as done in peekPoke2.cpp

### benchmarkCorrelation.cpp & benchmarkDerivative.cpp
Those programs serve as benchmarks for Grid's performance while computing the correlation (two-point) function (with/without derivative).

### Benchmark_su3.cpp
This file is a modified version of the Benchmark_su3.cc file from the Grid repository. Command line arguments where added to make it more comfortable to use (e.g. with bash scripts).

### loop.sh & slurm.sh
In order to run the benchmarks for various lattice sizes and other parameters, this bash script can be modified and executed.

### benchmarkMDA[x].cpp
Those benchmarks aim to test different ways of calculating the Meson Distribution Amplitude (MDA) with Grid.
