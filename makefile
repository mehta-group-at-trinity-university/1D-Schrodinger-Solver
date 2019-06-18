#version https://git-lfs.github.com/spec/v1
#oid sha256:7ce4eecacbc6182f09620163fedfd03c0c8ae0207b73d4585cffdac82b2c71b5
#size 152
1DSolver.exe:	1DSolver.f90 matrix_stuff.f90
	ifort 1DSolver.f90 matrix_stuff.f90 -Qmkl
