all:
	/usr/local/bin/f2py -c -m faras faras.f90 --fcompiler=intelem --opt='-xHost -O3 -axAVX -qopenmp -funroll-loops -qopt-prefetch' -liomp5 -lpthread -lm -ldl
