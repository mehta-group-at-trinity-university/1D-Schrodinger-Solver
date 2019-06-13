      subroutine CompSqrMatInv(A, N)
        implicit none
        integer N,info,lwk
        integer, allocatable :: ipiv(:)
        double precision, allocatable :: work(:)
        double complex A(N,N)
        allocate(ipiv(N))
        call zgetrf(N, N, A, N, ipiv, info)
        allocate(work(1))
        lwk = -1
        call zgetri(N, A, N, ipiv, work, lwk, info)
        lwk = work(1)
        deallocate(work)
        allocate(work(lwk))
        call zgetri(N, A, N, ipiv, work, lwk, info)
        deallocate(ipiv,work)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SqrMatInv(A, N)
      implicit none
      integer N,info,lwk
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: work(:)
      double precision A(N,N)
      allocate(ipiv(N))
      call dgetrf(N, N, A, N, ipiv, info)
      allocate(work(1))
      lwk = -1
      call dgetri(N, A, N, ipiv, work, lwk, info)
      lwk = work(1)
      deallocate(work)
      allocate(work(lwk))
      call dgetri(N, A, N, ipiv, work, lwk, info)
      deallocate(ipiv,work)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Computes the eigenvalues and eigenvectors of a generalized eigenvalue problem (not banded form)
!     this wrapper preserves the contents of G and L so that they are not overwritten
      subroutine Mydggev(N,G,LDG,L,LDL,eval,evec)
      implicit none
      integer LDG,N,LDL,info
      double precision G(LDG,N),L(LDL,N),eval(N),evec(N,N),GL(LDG,N),LL(LDL,N)
      double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:)
      integer lwork,i,im,in,j

      allocate(alphar(N),alphai(N),beta(N))
      GL=G
      LL=L
      info = 0
      lwork = -1
      allocate(work(1))
      call dggev('N','V',N,GL,LDG,LL,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)
      !write(6,*) "info = ", info, int(work(1))
      alphar=0d0
      alphai=0d0
      beta=0d0
      evec=0d0
      eval=0d0

      lwork=2*int(work(1))
      deallocate(work)
      allocate(work(lwork))
      GL=G
      LL=L
      call dggev('N','V',N,GL,LDG,LL,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)

      do i = 1, N
      if (abs(alphai(i)).ge.1d-15) then
         print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i)
      endif
      if(abs(beta(i)).ge.1d-15) then
         eval(i) = alphar(i)/beta(i)
      endif
      enddo
      deallocate(alphar,alphai,beta)
      deallocate(work)

      end subroutine Mydggev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes the eigenvalues and eigenvectors of a square matrix A of dimension N.
      subroutine MyDSYEV(A,N,EVAL,EVEC)
      implicit none
      double precision A(N,N), EVAL(N), EVEC(N,N), Atemp(N,N)
      character*1 JOBZ, UPLO
      integer  N, LDA, LWORK, INFO
      double precision, allocatable :: WORK(:)
      Atemp=A ! store the matrix A for safe keeping
      JOBZ = 'V' ! compute both eigenvals and eigenvectors
      UPLO = 'U' ! indicates that the upper triangle of of the matrix A is stored.
      allocate(WORK(2))
      LDA = N ! for a symmetric matrix
      LWORK = -1 ! perform a querry to determine the optimal size of WORK
      call DSYEV(JOBZ, UPLO, N, A, LDA, EVAL, WORK, LWORK, INFO)
      LWORK = WORK(1)
      deallocate(WORK)
      allocate(WORK(LWORK))
      call DSYEV(JOBZ, UPLO, N, A, LDA, EVAL, WORK, LWORK, INFO)
      EVEC = A ! store the eigenvectors in EVEC
      A=Atemp ! restore the original matrix in A

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcEigenErrors(info,iparam,MatrixDim,H,LeadDim,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

      integer info,iparam(*),MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
      double precision H(LeadDim,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

      integer j
      double precision dnrm2
      double precision, allocatable :: HPsi(:),SPsi(:)

      if ( info .eq. 0) then

       if (iparam(5) .lt. NumStates) write(6,*) 'Not all states found'

c Compute the residual norm: ||  A*x - lambda*x ||

       allocate(HPsi(MatrixDim))
       allocate(SPsi(MatrixDim))
       do j = 1,NumStates
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,LeadDim,Psi(1,j),1,0.0d0,HPsi,1)
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
        call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
        Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
        Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
       enddo
       deallocate(HPsi)
       deallocate(SPsi)
      else
       write(6,*) ' '
       write(6,*) ' Error with _sband, info= ', info
       write(6,*) ' Check the documentation of _sband '
       write(6,*) ' '
      end if

      return
      end
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style band form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost
c  of n*nev if both are requested (in this case a separate array Z must be
c  supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are called Ritz values and Ritz vectors respectively.  They are referred
c  to as such in the comments that follow.  The computed orthonormal basis
c  for the invariant subspace corresponding to these Ritz values is referred
c  to as a Lanczos basis.
c
c \Usage
c   call dsband
c      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA,
c        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V,
c        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
c
c \Arguments
c
c  RVEC    Logical (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the associated Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute all Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' , SELECT is not referenced.
c
c  D       Double precision array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the
c          Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  AB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC.
c
c  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c          It is not referenced when IPARAM(7) = 1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LA' -> want the NEV eigenvalues of largest REAL part.
c            'SA' -> want the NEV eigenvalues of smallest REAL part.
c            'BE' -> Compute NEV eigenvalues, half from each end of the
c                    spectrum.  When NEV is odd, compute one more from
c                    the high end than from the low end.
c
c          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of OP to be computed.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Lanczos basis constructed
c          by dsaupd for OP.
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Upon INPUT: the NCV columns of V contain the Lanczos basis
c                      vectors as constructed by dsaupd for OP.
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       represent the Ritz vectors that span the desired
c                       invariant subspace.
c          NOTE: The array Z may be set equal to first NEV columns of the
c          Lanczos basis vector array V computed by dsaupd. In this case
c          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c          the desired Ritz vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT:
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c          ------------------------------------------------------------
c          ISHIFT = 1: exact shifts with respect to the reduced
c                      tridiagonal matrix T.  This is equivalent to
c                      restarting the iteration with a starting vector
c                      that is a linear combination of Ritz vectors
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsband for the
c          five modes available.
c
c          IPARAM(8) = NP
c          Not referenced.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the
c          factorization of M or (A-SIGMA*M).
c
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 in DSAUPD.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Author
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
!subroutine setUpDSYVE(A,sz,EigenVals,EigenVecs)
!      subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,
!     >                  tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)
!
!      character        which*2, bmat, howmny
c      integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info
c      Double precision tol, sigma
c      logical          rvec
c
c      integer          iparam(*), iwork(*)
c      logical          select(*)
c      Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*)
c
c      integer          ipntr(14)
c
c      integer          ido, i, j, Row, Col, type, ierr
c
c      Double precision one, zero
c      parameter        (one = 1.0, zero = 0.0)
c
c      Double precision ddot, dnrm2, dlapy2
c      external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy
c
cc iparam(3) : Max number of Arnoldi iterations
c      iparam(3) = 300
c      iparam(7) = 3
c      rvec = .TRUE.
c      howmny = 'A'
c      which = 'LM'
c      bmat = 'G'
c      type = 4
c      ido = 0
c      iparam(1) = 1
c
c      rfac = 0.0d0
c      do i = 1,n
c       do j = i,min(i+k,n)
c        Row = k+1+i-j
c        Col = j
c        rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
c       enddo
c       do j = max(1,i-k),i-1
c        Row = 2*k+1
c        Col = j
c        rfac(Row+i-j,j) = rfac(Row+j-i,i)
c       enddo
c      enddo
c
c      call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
c      if ( ierr .ne. 0 )  then
c       print*, ' '
c       print*, '_SBAND: Error with _gbtrf:',ierr
c       print*, ' '
c       go to 9000
c      end if
c
c  90  continue
c
c      call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
c
c      if (ido .eq. -1) then
c       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
c       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
c       if (ierr .ne. 0) then
c        print*, ' '
c        print*, '_SBAND: Error with _gbtrs.'
c        print*, ' '
c        go to 9000
c       end if
c      else if (ido .eq. 1) then
c       call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
c       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
c       if (ierr .ne. 0) then
c        print*, ' '
c        print*, '_SBAND: Error with _gbtrs.'
c        print*, ' '
c        go to 9000
c       end if
c      else if (ido .eq. 2) then
c       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
c      else
c       if ( info .lt. 0) then
c        print *, ' '
c        print *, ' Error with _saupd info = ',info
c        print *, ' Check the documentation of _saupd '
c        print *, ' '
c        go to 9000
c       else
c        if ( info .eq. 1) then
c         print *, ' '
c         print *, ' Maximum number of iterations reached.'
c         print *, ' '
c        else if ( info .eq. 3) then
c         print *, ' '
c         print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
c         print *, ' '
c        end if
c        if (iparam(5) .gt. 0) then
c         call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
c         if ( info .ne. 0) then
c          print *, ' '
c          print *, ' Error with _neupd = ', info
c          print *, ' Check the documentation of _neupd '
c          print *, ' '
c          go to 9000
c         endif
c        endif
c       endif
c       go to 9000
c      endif
c
c      go to 90
c
c 9000 continue
c
c      end
cc \BeginDoc
cc
cc \Name: dsband
cc
cc \Description:
cc
cc  This subroutine returns the converged approximations to eigenvalues
cc  of A*z = lambda*B*z and (optionally):
cc
cc      (1) The corresponding approximate eigenvectors;
cc
cc      (2) An orthonormal (Lanczos) basis for the associated approximate
cc          invariant subspace;
cc
cc      (3) Both.
cc
cc  Matrices A and B are stored in LAPACK-style band form.
cc
cc  There is negligible additional cost to obtain eigenvectors.  An orthonormal
cc  (Lanczos) basis is always computed.  There is an additional storage cost
cc  of n*nev if both are requested (in this case a separate array Z must be
cc  supplied).
cc
cc  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
cc  are called Ritz values and Ritz vectors respectively.  They are referred
cc  to as such in the comments that follow.  The computed orthonormal basis
cc  for the invariant subspace corresponding to these Ritz values is referred
cc  to as a Lanczos basis.
cc
cc  dsband can be called with one of the following modes:
cc
cc  Mode 1:  A*x = lambda*x, A symmetric
cc           ===> OP = A  and  B = I.
cc
cc  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
cc           ===> OP = inv[M]*A  and  B = M.
cc           ===> (If M can be factored see remark 3 in DSAUPD)
cc
cc  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
cc           ===> OP = (inv[K - sigma*M])*M  and  B = M.
cc           ===> Shift-and-Invert mode
cc
cc  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
cc           KG symmetric indefinite
cc           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
cc           ===> Buckling mode
cc
cc  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
cc           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
cc           ===> Cayley transformed mode
cc
cc  The choice of mode must be specified in IPARAM(7) defined below.
cc
cc \Usage
cc   call dsband
cc      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA,
cc        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V,
cc        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
cc
cc \Arguments
cc
cc  RVEC    Logical (INPUT)
cc          Specifies whether Ritz vectors corresponding to the Ritz value
cc          approximations to the eigenproblem A*z = lambda*B*z are computed.
cc
cc             RVEC = .FALSE.     Compute Ritz values only.
cc
cc             RVEC = .TRUE.      Compute the associated Ritz vectors.
cc
cc  HOWMNY  Character*1  (INPUT)
cc          Specifies how many Ritz vectors are wanted and the form of Z
cc          the matrix of Ritz vectors. See remark 1 below.
cc          = 'A': compute all Ritz vectors;
cc          = 'S': compute some of the Ritz vectors, specified
cc                 by the logical array SELECT.
cc
cc  SELECT  Logical array of dimension NCV.  (INPUT)
cc          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
cc          computed. To select the Ritz vector corresponding to a
cc          Ritz value D(j), SELECT(j) must be set to .TRUE..
cc          If HOWMNY = 'A' , SELECT is not referenced.
cc
cc  D       Double precision array of dimension NEV.  (OUTPUT)
cc          On exit, D contains the Ritz value approximations to the
cc          eigenvalues of A*z = lambda*B*z. The values are returned
cc          in ascending order. If IPARAM(7) = 3,4,5 then D represents
cc          the Ritz values of OP computed by dsaupd transformed to
cc          those of the original eigensystem A*z = lambda*B*z. If
cc          IPARAM(7) = 1,2 then the Ritz values of OP are the same
cc          as the those of A*z = lambda*B*z.
cc
cc  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
cc          On exit, Z contains the B-orthonormal Ritz vectors of the
cc          eigensystem A*z = lambda*B*z corresponding to the Ritz
cc          value approximations.
cc
cc          If  RVEC = .FALSE. then Z is not referenced.
cc          NOTE: The array Z may be set equal to first NEV columns of the
cc          Lanczos basis array V computed by DSAUPD.
cc
cc  LDZ     Integer.  (INPUT)
cc          The leading dimension of the array Z.  If Ritz vectors are
cc          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
cc
cc  SIGMA   Double precision  (INPUT)
cc          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
cc          IPARAM(7) = 1 or 2.
cc
cc  N       Integer.  (INPUT)
cc          Dimension of the eigenproblem.
cc
cc  AB      Double precision array of dimension LDA by N. (INPUT)
cc          The matrix A in band storage, in rows KL+1 to
cc          2*KL+KU+1; rows 1 to KL of the array need not be set.
cc          The j-th column of A is stored in the j-th column of the
cc          array AB as follows:
cc          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
cc
cc  MB      Double precision array of dimension LDA by N. (INPUT)
cc          The matrix M in band storage, in rows KL+1 to
cc          2*KL+KU+1; rows 1 to KL of the array need not be set.
cc          The j-th column of M is stored in the j-th column of the
cc          array AB as follows:
cc          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
cc          Not referenced if IPARAM(7) = 1
cc
cc  LDA     Integer. (INPUT)
cc          Leading dimension of AB, MB, RFAC.
cc
cc  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
cc          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
cc          is invoked.  It is used to store the LU factors of
cc          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
cc          It is not referenced when IPARAM(7) = 1.
cc
cc  KL      Integer. (INPUT)
cc          Max(number of subdiagonals of A, number of subdiagonals of M)
cc
cc  KU      Integer. (OUTPUT)
cc          Max(number of superdiagonals of A, number of superdiagonals of M)
cc
cc  WHICH   Character*2.  (INPUT)
cc          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
cc          the following.
cc
cc            'LM' -> want the NEV eigenvalues of largest magnitude.
cc            'SM' -> want the NEV eigenvalues of smallest magnitude.
cc            'LA' -> want the NEV eigenvalues of largest REAL part.
cc            'SA' -> want the NEV eigenvalues of smallest REAL part.
cc            'BE' -> Compute NEV eigenvalues, half from each end of the
cc                    spectrum.  When NEV is odd, compute one more from
cc                    the high end than from the low end.
cc
cc          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.
cc
cc  BMAT    Character*1.  (INPUT)
cc          BMAT specifies the type of the matrix B that defines the
cc          semi-inner product for the operator OP.
cc          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
cc          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
cc  NEV     Integer. (INPUT)
cc          Number of eigenvalues of OP to be computed.
cc
cc  TOL     Double precision scalar.  (INPUT)
cc          Stopping criterion: the relative accuracy of the Ritz value
cc          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
cc          If TOL .LE. 0. is passed a default is set:
cc          DEFAULT = DLAMCH('EPS')  (machine precision as computed
cc                    by the LAPACK auxiliary subroutine DLAMCH).
cc
cc  RESID   Double precision array of length N.  (INPUT/OUTPUT)
cc          On INPUT:
cc          If INFO .EQ. 0, a random initial residual vector is used.
cc          If INFO .NE. 0, RESID contains the initial residual vector,
cc                          possibly from a previous run.
cc          On OUTPUT:
cc          RESID contains the final residual vector.
cc
cc  NCV     Integer.  (INPUT)
cc          Number of columns of the matrix V (less than or equal to N).
cc          Represents the dimension of the Lanczos basis constructed
cc          by dsaupd for OP.
cc
cc  V       Double precision array N by NCV.  (OUTPUT)
cc          Upon INPUT: the NCV columns of V contain the Lanczos basis
cc                      vectors as constructed by dsaupd for OP.
cc          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
cc                       represent the Ritz vectors that span the desired
cc                       invariant subspace.
cc          NOTE: The array Z may be set equal to first NEV columns of the
cc          Lanczos basis vector array V computed by dsaupd. In this case
cc          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
cc          the desired Ritz vectors.
cc
cc  LDV     Integer.  (INPUT)
cc          Leading dimension of V exactly as declared in the calling
cc          program.
cc
cc  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
cc          IPARAM(1) = ISHIFT:
cc          The shifts selected at each iteration are used to restart
cc          the Arnoldi iteration in an implicit fashion.
cc          It is set to 1 in this subroutine.  The user do not need
cc          to set this parameter.
cc          ------------------------------------------------------------
cc          ISHIFT = 1: exact shifts with respect to the reduced
cc                      tridiagonal matrix T.  This is equivalent to
cc                      restarting the iteration with a starting vector
cc                      that is a linear combination of Ritz vectors
cc                      associated with the "wanted" Ritz values.
cc          -------------------------------------------------------------
cc
cc          IPARAM(2) = No longer referenced.
cc
cc          IPARAM(3) = MXITER
cc          On INPUT:  max number of Arnoldi update iterations allowed.
cc          On OUTPUT: actual number of Arnoldi update iterations taken.
cc
cc          IPARAM(4) = NB: blocksize to be used in the recurrence.
cc          The code currently works only for NB = 1.
cc
cc          IPARAM(5) = NCONV: number of "converged" eigenvalues.
cc          This represents the number of Ritz values that satisfy
cc          the convergence criterion.
cc
cc          IPARAM(6) = IUPD
cc          No longer referenced. Implicit restarting is ALWAYS used.
cc
cc          IPARAM(7) = MODE
cc          On INPUT determines what type of eigenproblem is being solved.
cc          Must be 1,2,3,4,5; See under \Description of dsband for the
cc          five modes available.
cc
cc          IPARAM(8) = NP
cc          Not referenced.
cc
cc          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
cc          OUTPUT: NUMOP  = total number of OP*x operations,
cc                  NUMOPB = total number of B*x operations if BMAT='G',
cc                  NUMREO = total number of steps of re-orthogonalization.
cc
cc WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
cc
cc WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
cc
cc LWORKL   Integer.  (INPUT)
cc          LWORKL must be at least NCV**2 + 8*NCV.
cc
cc IWORK    Integer array of dimension at least N. (WORKSPACE)
cc          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the
cc          factorization of M or (A-SIGMA*M).
cc
cc INFO     Integer.  (INPUT/OUTPUT)
cc          Error flag on output.
cc          =  0: Normal exit.
cc          =  1: Maximum number of iterations taken.
cc                All possible eigenvalues of OP has been found. IPARAM(5)
cc                returns the number of wanted converged Ritz values.
cc          =  3: No shifts could be applied during a cycle of the
cc                Implicitly restarted Arnoldi iteration. One possibility
cc                is to increase the size of NCV relative to NEV.
cc                See remark 4 in DSAUPD.
cc
cc          = -1: N must be positive.
cc          = -2: NEV must be positive.
cc          = -3: NCV-NEV >= 2 and less than or equal to N.
cc          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
cc          = -6: BMAT must be one of 'I' or 'G'.
cc          = -7: Length of private work WORKL array is not sufficient.
cc          = -8: Error return from trid. eigenvalue calculation;
cc                Informational error from LAPACK routine dsteqr.
cc          = -9: Starting vector is zero.
cc          = -10: IPARAM(7) must be 1,2,3,4,5.
cc          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
cc          = -12: NEV and WHICH = 'BE' are incompatible.
cc          = -13: HOWMNY must be one of 'A' or 'P'
cc          = -14: DSAUPD did not find any eigenvalues to sufficient
cc                 accuracy.
cc          = -9999: Could not build an Arnoldi factorization.
cc                   IPARAM(5) returns the size of the current
cc                   Arnoldi factorization.
cc
cc \EndDoc
cc
cc------------------------------------------------------------------------
cc
cc\BeginLib
cc
cc\References:
cc  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
cc     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
cc     pp 357-385.
cc
cc  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
cc     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
cc     May 1995.
cc
cc\Routines called:
cc     dsaupd  ARPACK reverse communication interface routine.
cc     dseupd  ARPACK routine that returns Ritz values and (optionally)
cc             Ritz vectors.
cc     dgbtrf  LAPACK band matrix factorization routine.
cc     dgbtrs  LAPACK band linear system solve routine.
cc     dlacpy  LAPACK matrix copy routine.
cc     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
cc     dcopy   Level 1 BLAS that copies one vector to another.
cc     ddot    Level 1 BLAS that computes the dot product of two vectors.
cc     dnrm2   Level 1 BLAS that computes the norm of a vector.
cc     dgbmv   Level 2 BLAS that computes the band matrix vector product.
cc
cc\Remarks
cc  1. The converged Ritz values are always returned in increasing
cc     (algebraic) order.
cc
cc  2. Currently only HOWMNY = 'A' is implemented. It is included at this
cc     stage for the user who wants to incorporate it.
cc
cc\Author
cc     Danny Sorensen
cc     Richard Lehoucq
cc     Chao Yang
cc     Dept. of Computational &
cc     Applied Mathematics
cc     Rice University
cc     Houston, Texas
cc
cc\SCCS Information: @(#)
cc FILE: sband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
cc
cc\EndLib
cc
cc---------------------------------------------------------------------
cc
c      subroutine dsband( rvec, howmny, select, d, z, ldz, sigma,
c     &           n, ab, mb, lda, rfac, kl, ku, which, bmat, nev,
c     &           tol, resid, ncv, v, ldv, iparam, workd, workl,
c     &           lworkl, iwork, info)
cc
cc     %------------------%
cc     | Scalar Arguments |
cc     %------------------%
cc
c      character        which*2, bmat, howmny
c      integer          n, lda, kl, ku, nev, ncv, ldv,
c     &                 ldz, lworkl, info
c      Double precision
c     &                 tol, sigma
c      logical          rvec
cc
cc     %-----------------%
cc     | Array Arguments |
cc     %-----------------%
cc
c      integer          iparam(*), iwork(*)
c      logical          select(*)
c      Double precision
c     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
c     &                 ab(lda,*), mb(lda,*), rfac(lda,*),
c     &                 workd(*), workl(*)
cc
cc     %--------------%
cc     | Local Arrays |
cc     %--------------%
cc
c      integer          ipntr(14)
cc
cc     %---------------%
cc     | Local Scalars |
cc     %---------------%
cc
c      integer          ido, i, j, type, imid, itop, ibot, ierr
cc
cc     %------------%
cc     | Parameters |
cc     %------------%
cc
c      Double precision
c     &                  one, zero
c      parameter        (one = 1.0, zero = 0.0)
cc
cc
cc     %-----------------------------%
cc     | LAPACK & BLAS routines used |
cc     %-----------------------------%
cc
c      Double precision
c     &                 ddot, dnrm2, dlapy2
c      external         ddot, dcopy, dgbmv, dgbtrf,
c     &                 dgbtrs, dnrm2, dlapy2, dlacpy
cc
cc     %-----------------------%
cc     | Executable Statements |
cc     %-----------------------%
cc
cc     %----------------------------------------------------------------%
cc     | Set type of the problem to be solved. Check consistency        |
cc     | between BMAT and IPARAM(7).                                    |
cc     | type = 1 --> Solving standard problem in regular mode.         |
cc     | type = 2 --> Solving standard problem in shift-invert mode.    |
cc     | type = 3 --> Solving generalized problem in regular mode.      |
cc     | type = 4 --> Solving generalized problem in shift-invert mode. |
cc     | type = 5 --> Solving generalized problem in Buckling mode.     |
cc     | type = 6 --> Solving generalized problem in Cayley mode.       |
cc     %----------------------------------------------------------------%
cc
c      if ( iparam(7) .eq. 1 ) then
c         type = 1
c      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'I') then
c         type = 2
c      else if ( iparam(7) .eq. 2 ) then
c         type = 3
c      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
c         type = 4
c      else if ( iparam(7) .eq. 4 ) then
c         type = 5
c      else if ( iparam(7) .eq. 5 ) then
c         type = 6
c      else
c         print*, ' '
c         print*, 'BMAT is inconsistent with IPARAM(7).'
c         print*, ' '
c         go to 9000
c      end if
cc
cc     %------------------------%
cc     | Initialize the reverse |
cc     | communication flag.    |
cc     %------------------------%
cc
c      ido   = 0
cc
cc     %----------------%
cc     | Exact shift is |
cc     | used.          |
cc     %----------------%
cc
c      iparam(1) = 1
cc
cc     %-----------------------------------%
cc     | Both matrices A and M are stored  |
cc     | between rows itop and ibot.  Imid |
cc     | is the index of the row that      |
cc     | stores the diagonal elements.     |
cc     %-----------------------------------%
cc
c      itop = kl + 1
c      imid = kl + ku + 1
c      ibot = 2*kl + ku + 1
cc
c      if ( type .eq. 2 .or. type .eq. 6 .and. bmat .eq. 'I' ) then
cc
cc         %----------------------------------%
cc         | Solving a standard eigenvalue    |
cc         | problem in shift-invert or       |
cc         | Cayley mode. Factor (A-sigma*I). |
cc         %----------------------------------%
cc
c          call dlacpy ('A', ibot, n, ab, lda, rfac, lda )
c          do 10 j = 1, n
c             rfac(imid,j) =  ab(imid,j) - sigma
c  10      continue
c          call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr )
c          if (ierr .ne. 0) then
c             print*, ' '
c             print*, ' _SBAND: Error with _gbtrf. '
c             print*, ' '
c             go to  9000
c          end if
cc
c      else if ( type .eq. 3 ) then
cc
cc        %----------------------------------------------%
cc        | Solving generalized eigenvalue problem in    |
cc        | regular mode. Copy M to rfac and Call LAPACK |
cc        | routine dgbtrf to factor M.                  |
cc        %----------------------------------------------%
cc
c         call dlacpy ('A', ibot, n, mb, lda, rfac, lda )
c         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
c         if (ierr .ne. 0) then
c             print*, ' '
c             print*,'_SBAND:  Error with _gbtrf.'
c             print*, ' '
c             go to 9000
c         end if
cc
c      else if ( type .eq. 4 .or. type .eq. 5 .or. type .eq. 6
c     &         .and. bmat .eq. 'G' ) then
cc
cc        %-------------------------------------------%
cc        | Solving generalized eigenvalue problem in |
cc        | shift-invert, Buckling, or Cayley mode.   |
cc        %-------------------------------------------%
cc
cc        %-------------------------------------%
cc        | Construct and factor (A - sigma*M). |
cc        %-------------------------------------%
cc
c         do 60 j = 1,n
c            do 50 i = itop, ibot
c               rfac(i,j) = ab(i,j) - sigma*mb(i,j)
c  50        continue
c  60     continue
cc
c         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
c         if ( ierr .ne. 0 )  then
c             print*, ' '
c             print*, '_SBAND: Error with _gbtrf.'
c             print*, ' '
c             go to 9000
c         end if
cc
c      end if
cc
cc     %--------------------------------------------%
cc     |  M A I N   L O O P (reverse communication) |
cc     %--------------------------------------------%
cc
c  90  continue
cc
c      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
c     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
c     &              info )
cc
c      if (ido .eq. -1) then
cc
c         if ( type .eq. 1) then
cc
cc           %----------------------------%
cc           | Perform  y <--- OP*x = A*x |
cc           %----------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                 lda, workd(ipntr(1)), 1, zero,
c     &                 workd(ipntr(2)), 1)
cc
c         else if ( type .eq. 2 ) then
cc
cc           %----------------------------------%
cc           |             Perform              |
cc           | y <--- OP*x = inv[A-sigma*I]*x   |
cc           | to force the starting vector     |
cc           | into the range of OP.            |
cc           %----------------------------------%
cc
c            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                    iwork, workd(ipntr(2)), n, ierr)
c            if ( ierr .ne. 0 ) then
c               print*, ' '
c               print*, ' _SBAND: Error with _bgtrs. '
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 3 ) then
cc
cc           %-----------------------------------%
cc           | Perform  y <--- OP*x = inv[M]*A*x |
cc           | to force the starting vector into |
cc           | the range of OP.                  |
cc           %-----------------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                  lda, workd(ipntr(1)), 1, zero,
c     &                  workd(ipntr(2)), 1)
c            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                    iwork, workd(ipntr(2)), n, ierr)
c            if (ierr .ne. 0) then
c               print*, ' '
c               print*, '_SBAND: Error with sbgtrs.'
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 4 ) then
cc
cc           %-----------------------------------------%
cc           | Perform y <-- OP*x                      |
cc           |           = inv[A-SIGMA*M]*M            |
cc           | to force the starting vector into the   |
cc           | range of OP.                            |
cc           %-----------------------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1),
c     &                 lda, workd(ipntr(1)), 1, zero,
c     &                 workd(ipntr(2)), 1)
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                   iwork, workd(ipntr(2)), n, ierr)
c            if (ierr .ne. 0) then
c               print*, ' '
c               print*, '_SBAND: Error with _gbtrs.'
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 5) then
cc
cc           %---------------------------------------%
cc           | Perform y <-- OP*x                    |
cc           |    = inv[A-SIGMA*M]*A                 |
cc           | to force the starting vector into the |
cc           | range of OP.                          |
cc           %---------------------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                 lda, workd(ipntr(1)), 1, zero,
c     &                 workd(ipntr(2)), 1)
c            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                   iwork, workd(ipntr(2)), n, ierr)
cc
c            if ( ierr .ne. 0 ) then
c               print*, ' '
c               print*, ' _SBAND: Error with _gbtrs. '
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 6 ) then
cc
cc           %---------------------------------------%
cc           | Perform y <-- OP*x                    |
cc           | = (inv[A-SIGMA*M])*(A+SIGMA*M)*x      |
cc           | to force the starting vector into the |
cc           | range of OP.                          |
cc           %---------------------------------------%
cc
c            if ( bmat .eq. 'G' ) then
c               call dgbmv('Notranspose', n, n, kl, ku, one,
c     &                    ab(itop,1), lda, workd(ipntr(1)), 1,
c     &                    zero, workd(ipntr(2)), 1)
c               call dgbmv('Notranspose', n, n, kl, ku, sigma,
c     &                    mb(itop,1), lda, workd(ipntr(1)), 1,
c     &                    one, workd(ipntr(2)), 1)
c            else
c               call dcopy(n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                    lda, workd(ipntr(1)), 1, sigma,
c     &                    workd(ipntr(2)), 1)
c            end if
cc
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                   iwork, workd(ipntr(2)), n, ierr)
cc
c            if (ierr .ne. 0) then
c               print*, ' '
c               print*, '_SBAND: Error with _gbtrs.'
c               print*, ' '
c               go to 9000
c            end if
cc
c         end if
cc
c      else if (ido .eq. 1) then
cc
c         if ( type .eq. 1) then
cc
cc           %----------------------------%
cc           | Perform  y <--- OP*x = A*x |
cc           %----------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                 lda, workd(ipntr(1)), 1, zero,
c     &                 workd(ipntr(2)), 1)
cc
c         else if ( type .eq. 2) then
cc
cc              %----------------------------------%
cc              |             Perform              |
cc              | y <--- OP*x = inv[A-sigma*I]*x.  |
cc              %----------------------------------%
cc
c               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                       iwork, workd(ipntr(2)), n, ierr)
c               if ( ierr .ne. 0 ) then
c                  print*, ' '
c                  print*, '_SBAND: Error with _gbtrs.'
c                  print*, ' '
c                  go to 9000
c               end if
cc
c         else if ( type .eq. 3 ) then
cc
cc           %-----------------------------------%
cc           | Perform  y <--- OP*x = inv[M]*A*x |
cc           %-----------------------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                  lda, workd(ipntr(1)), 1, zero,
c     &                  workd(ipntr(2)), 1)
c            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                    iwork, workd(ipntr(2)), n, ierr)
c            if (ierr .ne. 0) then
c               print*, ' '
c               print*, '_SBAND: error with _bgtrs.'
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 4 ) then
cc
cc           %-------------------------------------%
cc           | Perform y <-- inv(A-sigma*M)*(M*x). |
cc           | (M*x) has been computed and stored  |
cc           | in workd(ipntr(3)).                 |
cc           %-------------------------------------%
cc
c            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
c            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                    iwork, workd(ipntr(2)), n, ierr)
c            if (ierr .ne. 0) then
c               print*, ' '
c               print*, '_SBAND: Error with _gbtrs.'
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 5 ) then
cc
cc           %-------------------------------%
cc           | Perform y <-- OP*x            |
cc           |    = inv[A-SIGMA*M]*A*x       |
cc           | B*x = A*x has been computed   |
cc           | and saved in workd(ipntr(3)). |
cc           %-------------------------------%
cc
c            call dcopy (n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
c            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                   iwork, workd(ipntr(2)), n, ierr)
c            if ( ierr .ne. 0 ) then
c               print*, ' '
c               print*, ' _SBAND: Error with _gbtrs. '
c               print*, ' '
c               go to 9000
c            end if
cc
c         else if ( type .eq. 6) then
cc
cc           %---------------------------------%
cc           | Perform y <-- OP*x              |
cc           | = inv[A-SIGMA*M]*(A+SIGMA*M)*x. |
cc           | (M*x) has been saved in         |
cc           | workd(ipntr(3)).                |
cc           %---------------------------------%
cc
c            if ( bmat .eq. 'G' ) then
c               call dgbmv('Notranspose', n, n, kl, ku, one,
c     &                    ab(itop,1), lda, workd(ipntr(1)), 1,
c     &                    zero, workd(ipntr(2)), 1)
c               call daxpy( n, sigma, workd(ipntr(3)), 1,
c     &                    workd(ipntr(2)), 1 )
c            else
c               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1),
c     &                    lda, workd(ipntr(1)), 1, sigma,
c     &                    workd(ipntr(2)), 1)
c            end if
c            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda,
c     &                   iwork, workd(ipntr(2)), n, ierr)
cc
c         end if
cc
c      else if (ido .eq. 2) then
cc
cc        %----------------------------------%
cc        |        Perform y <-- B*x         |
cc        | Note when Buckling mode is used, |
cc        | B = A, otherwise B=M.            |
cc        %----------------------------------%
cc
c         if (type .eq. 5) then
cc
cc           %---------------------%
cc           | Buckling Mode, B=A. |
cc           %---------------------%
cc
c            call dgbmv('Notranspose', n, n, kl, ku, one,
c     &                ab(itop,1), lda, workd(ipntr(1)), 1,
c     &                zero, workd(ipntr(2)), 1)
c         else
c
c            call dgbmv('Notranspose', n, n, kl, ku, one,
c     &                mb(itop,1), lda, workd(ipntr(1)), 1,
c     &                zero, workd(ipntr(2)), 1)
c         end if
cc
c      else
cc
cc        %-----------------------------------------%
cc        | Either we have convergence, or there is |
cc        | error.                                  |
cc        %-----------------------------------------%
cc
c         if ( info .lt. 0) then
cc
cc           %--------------------------%
cc           | Error message, check the |
cc           | documentation in DSAUPD  |
cc           %--------------------------%
cc
c            print *, ' '
c            print *, ' Error with _saupd info = ',info
c            print *, ' Check the documentation of _saupd '
c            print *, ' '
c            go to 9000
cc
c         else
cc
c            if ( info .eq. 1) then
c               print *, ' '
c               print *, ' Maximum number of iterations reached.'
c               print *, ' '
c            else if ( info .eq. 3) then
c               print *, ' '
c               print *, ' No shifts could be applied during implicit',
c     &                  ' Arnoldi update, try increasing NCV.'
c               print *, ' '
c            end if
cc
c            if (iparam(5) .gt. 0) then
cc
c               call dseupd ( rvec, 'A', select, d, z, ldz, sigma,
c     &                  bmat, n, which, nev, tol, resid, ncv, v, ldv,
c     &                  iparam, ipntr, workd, workl, lworkl, info )
cc
c               if ( info .ne. 0) then
cc
cc                 %------------------------------------%
cc                 | Check the documentation of dneupd. |
cc                 %------------------------------------%
cc
c                  print *, ' '
c                  print *, ' Error with _neupd = ', info
c                  print *, ' Check the documentation of _neupd '
c                  print *, ' '
c                  go to 9000
cc
c               end if
cc
c            end if
cc
c         end if
cc
c         go to 9000
cc
c      end if
cc
cc     %----------------------------------------%
cc     | L O O P  B A C K to call DSAUPD again. |
cc     %----------------------------------------%
cc
c      go to 90
cc
c 9000 continue
cc
c      end
