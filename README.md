# F2PY augmented Machine Learning (FML)

FML is a Python library to efficiently carry out machine learning calculation. Everything is accessed via a simple Python interface, while efficient, parallel code (written in F90 with OpenMP) is executed behind the scenes.


## 1) Installation

By default FML compiles with GCC's gfortran and your system's standard BLAS+LAPACK libraries (-lblas -llapack). I recommend you switch to MKL for the math libraries, but the difference between ifort and gfortran shouldn't be to significant. BLAS and LAPACK implementations are the only libraries required for FML. Additionally you need a functional Python 2.x interpreter and NumPy (with F2PY) installed. Most Linux systems can install BLAS and LAPACK like this:

    sudo apt-get install libblas-dev liblapack-dev

Ok, on to the installation instructions:

1.1) First you clone this repository: 

    git clone https://github.com/andersx/fml.git

1.2) Then you simply compile by typing make in the fml folder:

    make

1.3) To make everything accessible to your Python export the fml root-folder to your PYTHONPATH. E.g. this is what I export:

    export PYTHONPATH=/home/andersx/dev/fml:$PYTHONPATH

## 2) How to use:

Will be updated soon! :D
