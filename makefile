#LINUX
#
#1DSolver.x:	1DSolver.o
#	ifort 1DSolver.o  -L/opt/ARPACK/ -larpack_Intel -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a -Wl,--end-group -lpthread -lm -parallel -o Solver.x

#1DSolver.o:	1DSolver.f90
#	ifort -O4 -extend_source -c -parallel 1DSolver.f90
#
#WINDOWS
#
1DSolver.exe:	1DSolver.f90 
	ifort /O3 /fp:precise /fpscomp:ioformat /fpscomp:logicals  /debug:full /traceback /Qsave /Qzero /gen-interfaces /warn:interfaces /check /fpe0 1DSolver.f90 -Qmkl -link -libpath:c:/FEAST/lib/x64 libfeast_sparse.a
