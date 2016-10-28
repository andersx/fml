OBJECTS = fml/math/fcho_solve.so

IFORT_COMPILER_FLAGS = --opt='-xHost -O3' --fcompiler=intelem --f90flags='-qopenmp -I${MKLROOT}/include'  
IFORT_MKL_LINKER_FLAGS = -liomp5  -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

GNU_COMPILER_FLAGS = --f90flags='-fopenmp -O3 -m64 -I${MKLROOT}/include'
GNU_MKL_LINKER_FLAGS = -lgomp -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl


# Uncomment to use GCC with MKL
# COMPILER_FLAGS = $(GNU_COMPILER_FLAGS)
# LINKER_FLAGS = $(GNU_MKL_LINKER_FLAGS)

# Uncomment to use Intel comilers with MKL
COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
LINKER_FLAGS = $(IFORT_MKL_LINKER_FLAGS)

all: $(OBJECTS)

fml/math/fcho_solve.so: fml/math/fcho_solve.f90
	f2py -c -m fcho_solve fml/math/fcho_solve.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)
	mv fcho_solve.so fml/math/

clean:
	rm -f fml/*.pyc
	rm -f fml/math/fcho_solve.so
	rm -f fml/math/*.pyc
