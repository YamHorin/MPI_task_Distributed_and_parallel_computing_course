finally:
	mpicc dsqrt.c -o dsqrt

	mpicc ssqrt.c -o ssqrt

	mpicc sqrt.c -o sqrt
clean:
	rm -f dsqrt ssqrt sqrt