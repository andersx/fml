OBJECTS = \
			fml/math/fcho_solve.so \
			fml/math/fdistance.so \
			fml/math/farad.so \
			fml/kernels/farad_kernels.so \
			fml/kernels/fkernels.so \
			fml/kernels/fforce_kernels.so \
			fml/representations/frepresentations.so

# Flags for Ifort and MKL
IFORT_COMPILER_FLAGS = --fcompiler=intelem \
					   --opt='-xHost -O3 -axAVX -qopenmp -funroll-loops -qopt-prefetch' \
					   --f77flags='-warn unused -I${MKLROOT}/include'
					   --f90flags='-warn unused -I${MKLROOT}/include'
IFORT_LINKER_FLAGS = -liomp5 -lpthread -lm -ldl
IFORT_MKL_LINKER_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_rt

# Flags for GCC compilers and MKL
GNU_COMPILER_FLAGS = --opt='-O3 -fopenmp -O3 -m64 -march=native' --f90flags=''
# GNU_COMPILER_FLAGS = --opt='-O3 -fopenmp -O3 -m64 -march=native' --f90flags='-I${MKLROOT}/include'
GNU_LINKER_FLAGS = -lgomp -lpthread -lm -ldl
# GNU_MKL_LINKER_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_rt
GNU_MKL_LINKER_FLAGS = -lblas -llapack

# Uncomment to use GCC with MKL
COMPILER_FLAGS = $(GNU_COMPILER_FLAGS)
LINKER_FLAGS = $(GNU_LINKER_FLAGS)
MKL_LINKER_FLAGS = $(GNU_MKL_LINKER_FLAGS)

# Uncomment to use Intel comilers with MKL
# COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
# LINKER_FLAGS = $(IFORT_LINKER_FLAGS)
# MKL_LINKER_FLAGS = $(IFORT_MKL_LINKER_FLAGS)

F2PY_EXEC = f2py2.7

all: $(OBJECTS)

fml/math/fcho_solve.so: fml/math/fcho_solve.f90
	$(F2PY_EXEC) -c -m fcho_solve fml/math/fcho_solve.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS) $(MKL_LINKER_FLAGS)
	mv fcho_solve*.so fml/math/fcho_solve.so

fml/math/fdistance.so: fml/math/fdistance.f90
	$(F2PY_EXEC) -c -m fdistance fml/math/fdistance.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)
	mv fdistance*.so fml/math/fdistance.so

fml/math/farad.so: fml/math/farad.f90
	$(F2PY_EXEC) -c -m farad fml/math/farad.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)
	mv farad*.so fml/math/farad.so

fml/kernels/farad_kernels.so: fml/kernels/farad_kernels.f90
	$(F2PY_EXEC) -c -m farad_kernels fml/kernels/farad_kernels.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)  $(MKL_LINKER_FLAGS)
	mv farad_kernels*.so fml/kernels/farad_kernels.so

fml/kernels/fkernels.so: fml/kernels/fkernels.f90
	$(F2PY_EXEC) -c -m fkernels fml/kernels/fkernels.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS) $(MKL_LINKER_FLAGS)
	mv fkernels*.so fml/kernels/fkernels.so

fml/kernels/fforce_kernels.so: fml/kernels/fforce_kernels.f90
	$(F2PY_EXEC) -c -m fforce_kernels fml/kernels/fforce_kernels.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS) $(MKL_LINKER_FLAGS)
	mv fforce_kernels*.so fml/kernels/fforce_kernels.so

fml/representations/frepresentations.so: fml/representations/frepresentations.f90
	$(F2PY_EXEC) -c -m frepresentations fml/representations/frepresentations.f90 $(COMPILER_FLAGS) $(LINKER_FLAGS)
	mv frepresentations*.so fml/representations/frepresentations.so
clean:
	rm -f fml/*.pyc
	rm -f fml/math/*.so
	rm -f fml/math/*.pyc
	rm -f fml/kernels/*.so
	rm -f fml/kernels/*.pyc
	rm -f fml/representations/*.so
	rm -f fml/representations/*.pyc
	rm -f pgopti.dpi pgopti.dpi.lock *.dyn

test:
	tests/test_arad.py
