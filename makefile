#LINUX
#
#MKL := /opt/intel.190611/mkl/lib/intel64
#1DSolver.x:	1DSolver.o
#	ifort 1DSolver.o  -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_sequential.a $(MKL)/libmkl_core.a $(MKL)/libmkl_blas95_ilp64.a $(MKL)/libmkl_lapack95_ilp64.a -Wl,--end-group -lpthread -lm -parallel -o Solver.x

#1DSolver.x:	1DSolver.f90
#	ifort -O4 -o Solver.x 1DSolver.f90 -mkl -check arg_temp_created -g -traceback -warn all -check all -heap-arrays -debug extended
#1DSolver.o:	1DSolver.f90
#	ifort 1DSolver.f90 -mkl
#


#WINDOWS
#
#1DSolver.exe:	1DSolver.f90 
#	ifort /O3 /fp:precise /fpscomp:ioformat /fpscomp:logicals  /debug:full /traceback /Qsave /Qzero /gen-interfaces /warn:interfaces /check /fpe0 1DSolver.f90 -Qmkl -link -libpath:c:/FEAST/lib/x64 libfeast_sparse.a

1DSolver.exe:	1DSolver.f90 
	ifort 1DSolver.f90 -Qmkl
clean:
	del *.x *.exe *.out *.obj *.un~ *.o *__genmod.mod *__genmod.f90 *.pdb
