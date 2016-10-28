OBJECTS = fml/math/fcho_solve.so fml/kernels/fkernels.so

# Flags for Ifort and MKL
IFORT_COMPILER_FLAGS = --opt='-xHost -O3' --fcompiler=intelem --f90flags='-qopenmp -I${MKLROOT}/include'  
IFORT_LINKER_FLAGS = -liomp5 -lpthread -lm -ldl
IFORT_MKL_LINKER_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_rt

# Flags for GCC compilers and MKL
GNU_COMPILER_FLAGS = --f90flags='-fopenmp -O3 -m64 -I${MKLROOT}/include'
GNU_LINKER_FLAGS = -lgomp -lpthread -lm -ldl
GNU_MKL_LINKER_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_rt

# Uncomment to use GCC with MKL
# COMPILER_FLAGS = $(GNU_COMPILER_FLAGS)
# LINKER_FLAGS = $(GNU_LINKER_FLAGS)
# MKL_LINKER_FLAGS = $(GNU_MKL_LINKER_FLAGS)

# Uncomment to use Intel comilers with MKL
COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
LINKER_FLAGS = $(IFORT_LINKER_FLAGS)
MKL_LINKER_FLAGS = $(IFORT_MKL_LINKER_FLAGS)

all: $(OBJECTS)

fml/math/fcho_solve.so: fml/math/fcho_solve.f90
	f2py -c -m fcho_solve fml/math/fcho_solve.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS) $(MKL_LINKER_FLAGS)
	mv fcho_solve.so fml/math/

fml/kernels/fkernels.so: fml/kernels/fkernels.f90
	f2py -c -m fkernels fml/kernels/fkernels.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)
	mv fkernels.so fml/kernels/

clean:
	rm -f fml/*.pyc
	rm -f fml/math/*.so
	rm -f fml/math/*.pyc
	rm -f fml/kernels/*.so
	rm -f fml/kernels/*.pyc
