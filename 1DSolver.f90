!this program uses Guassian-Lobatto quadrature and the DVR method to solve the schrodinger equation in 3D. 
!Uses sparce matricies for T and V
!Outputs can be eigenvalues  
program Solver
use mkl_spblas
use IFLPORT
use, intrinsic :: ISO_C_BINDING, ONLY : C_PTR, C_INT, C_DOUBLE
implicit none
        real*8, allocatable :: nodesX(:), weightsX(:), nodesY(:), weightsY(:), nodesZ(:), &
                weightsZ(:),holder(:),test(:),EigenVals(:),TmatS(:),TmatCol(:),TmatInd(:), &
                Vsparse(:),jobs(:), res(:), & 
                Hsparse(:)!, Hsparsee(:)
        real*8, allocatable :: Xx(:,:),dXx(:,:),Xy(:,:),dXy(:,:),Xz(:,:),dXz(:,:), &
                weightsDMX(:,:),weightsDMY(:,:),weightsDMZ(:,:),Hmat(:,:),EigenVecs(:,:), &
                Tmat(:,:), TmatX(:,:),TmatY(:,:),TmatZ(:,:),dXxtrim(:,:), dXytrim(:,:), dXztrim(:,:)
        real*8, allocatable :: indexOf(:,:,:)
        real*8 :: n,w,valX,valXp,a,b,Tval,epsout,temp
        integer :: i,j,k,ip,jp,kp,m,sX,sY,sZ,sMax,numN0,ijk,ijkp,sXc,sYc,sZc,info,loop,infoh,rCount,row,m0
        integer, allocatable :: feastparam(:), Hrow(:), Hcol(:), Vrow(:), Vcol(:)!, Hroww(:), Hcoll(:)
        integer :: HnumRow, HnumCol, indexing
        !integer(C_INT) :: op
        !real(C_DOUBLE) :: alpha
        !type(C_PTR) :: Hb, He, Hcol
        !type(C_PTR) :: Hsparse
        type(sparse_matrix_t) :: HH, VH, TH
 
        open(unit=1,file="nodesAndWeights10.dat")
        open(unit=2,file="GLnodesOut.dat")
        open(unit=3,file="basisFuncsX.dat")
        open(unit=4,file="eigenValsOut.dat")

        read(1,*) sX
        sY=sX
        sZ=sX
        sXc=sX-2
        sYc=sY-2
        sZc=sZ-2
        sMax=sXc*sYc*sZc

        allocate(nodesX(sX),weightsX(sX),nodesY(sY),weightsY(sY),nodesZ(sZ),weightsZ(sZ),holder(sX),test(200), &
                Xx(sX,sX),dXx(sX,sX),Xy(sY,sY),dXy(sY,sY),Xz(sZ,sZ),dXz(sZ,sZ), &
                weightsDMX(sX,sX),weightsDMY(sY,sY),weightsDMZ(sZ,sZ), &
                indexOf(sX-1,sY-1,sZ-1),Vsparse(sMax),Vrow(sMax),Vcol(sMax), &
                Hmat(sMax,sMax), &
                Tmat(sMax,sMax), TmatX(sXc,sXc), TmatY(sYc,sYc),TmatZ(sZc,sZc), &
                dXxtrim(sXc,sX),dXytrim(sYc,sY),dXztrim(sZc,sZ))

        !read in data for each basis
        do i=1,sX
                read(1,*) n,w
                nodesX(i)=n
                weightsX(i)=w
        end do
      !  do i=1,sY
      !          read(1,*) n,w
      !          nodesY(i)=n
      !          weightsY(i)=w
      !  end do
      !  do i=1,sZ
      !          read(1,*) n,w
      !          nodesZ(i)=n
      !          weightsZ(i)=w
      !  end do
        nodesY=nodesX
        weightsY=weightsX
        nodesZ=nodesX
        weightsZ=weightsX
        
        !rescale the basis functions to the desired range
        write(*,*) "enter rescale range"
        read(*,*) a,b 
        call rescale(sX,a,b,nodesX)
        call rescale(sY,a,b,nodesY)
        call rescale(sZ,a,b,nodesZ)
        write(*,*) "rescaled"
        
        !fill test array with 200 between a and b
        !do i=1,200
        !        test(i)=a+((b-a)/199)*(i-1)
        !end do

        !get matrix of nodes and derivatives
        do i=1,sX
                do j=1,sX
                        call DVRFunctions(holder,nodesX,nodesX(j),i,sX,valX,weightsX)
                        Xx(i,j)=valX

                        call DifferentiateChi(holder,nodesX,nodesX(j),i,sX,valXp,weightsX)
                        dXx(i,j)=valXp
                end do
        end do
        do i=1,sY
                do j=1,sY
                        call DVRFunctions(holder,nodesY,nodesY(j),i,sY,valX,weightsY)
                        Xy(i,j)=valX

                        call DifferentiateChi(holder,nodesY,nodesY(j),i,sY,valXp,weightsY)
                        dXy(i,j)=valXp
                end do
        end do
        do i=1,sZ
                do j=1,sZ
                        call DVRFunctions(holder,nodesZ,nodesZ(j),i,sZ,valX,weightsZ)
                        Xz(i,j)=valX

                        call DifferentiateChi(holder,nodesZ,nodesZ(j),i,sZ,valXp,weightsZ)
                        dXz(i,j)=valXp
                end do
        end do

        print *, "got Chi's"

        !Xxtrim = Xx(2:sX-1,1:sX)
        !Xytrim = Xy(2:sY-1,1:sY)
        !Xztrim = Xz(2:sZ-1,1:sZ)
        dXxtrim = dXx(2:sX-1,1:sX)
        dXytrim = dXy(2:sY-1,1:sY)
        dXztrim = dXz(2:sZ-1,1:sZ)



        !get matris of derivitaves for test array
        !do i=1,sz
        !        do j=1,200
        !                call DifferentiateChi(holder,nodes,test(j),i,sz,valXp,weights)
        !                call LobattoProduct(holder,nodes,test(j),i,sz,valX,weights)
!
!                        dXGL(i,j)=valXp
!                        write(2,*) test(j), valXp
!                        write(3,*) test(j), valX
!                end do
!                write(2,*)
!                write(3,*)
!        end do
        
        !make index array and potential diag matrix
        row = 0
        do i=1,sXc
                do j=1,sYc
                        do k=1,sZc
                                row=row+1
                                indexOf(i,j,k)=row
                                Call V3d(nodesX(i),nodesY(j),nodesZ(k),valX)
                                Vsparse(row)=valX
                                Vrow(row)=row
                                Vcol(row)=row
                        end do
                end do
        end do

        print *, "Got V"

        !make T kenetic evergy matrix for each dimenstion
        !first make diagonal weight matrices
        do i=1,sX
                do ip=1,sX
                        if (i == ip) then
                                weightsDMX(i,ip)=weightsX(i)
                        else
                                weightsDMX(i,ip) = 0
                        end if
                end do
        end do
        do i=1,sY
                do ip=1,sY
                        if (i == ip) then
                                weightsDMY(i,ip)=weightsY(i)
                        else
                                weightsDMY(i,ip) = 0
                        end if
                end do
        end do
        do i=1,sZ
                do ip=1,sZ
                        if (i == ip) then
                                weightsDMZ(i,ip)=weightsZ(i)
                        else
                                weightsDMZ(i,ip) = 0
                        end if
                end do
        end do
        !now make the Tmat for each dimention
        TmatX = (0.5d0)*MATMUL(dXytrim,MATMUL(weightsDMX,TRANSPOSE(dXxtrim)))
        TmatY = (0.5d0)*MATMUL(dXytrim,MATMUL(weightsDMY,TRANSPOSE(dXytrim)))
        TmatZ = (0.5d0)*MATMUL(dXztrim,MATMUL(weightsDMZ,TRANSPOSE(dXztrim)))
        !now combine them using delta properties
        numN0=0
        row=0
        rCount=1
        do m=1,2
                do i=1,sXc
                        do j=1,sYc
                                do k=1,sZc
                                        row=row+1
                                        if(m==2) then
                                                !Hrow(row)=rCount
                                        end if
                                        do ip=1,sXc
                                                do jp=1,sYc
                                                        do kp=1,sZc
                                                                ijk=indexOf(i,j,k)
                                                                ijkp=indexOf(ip,jp,kp)
                                                                temp = 0d0
                                                                !Tmat(ijk,ijkp)=0d0
                                                                if((j==jp).and.(k==kp)) then
                                                                        !Tmat(ijk,ijkp)=Tmat(ijk,ijkp)+TmatX(i,ip)
                                                                        temp=temp+TmatX(i,ip)
                                                                end if
                                                                if((i==ip).and.(k==kp)) then
                                                                        !Tmat(ijk,ijkp)=Tmat(ijk,ijkp)+TmatY(j,jp)
                                                                        temp=temp+TmatY(j,jp)
                                                                end if
                                                                if((j==jp).and.(i==ip)) then
                                                                        !Tmat(ijk,ijkp)=Tmat(ijk,ijkp)+TmatZ(k,kp)
                                                                        temp=temp+TmatZ(k,kp)
                                                                end if
                                                                if((i==ip).and.(j==jp).and.(k==kp)) then
                                                                        temp=temp+Vsparse(ijk)
                                                                end if
                                                                if(temp /= 0) then
                                                                        numN0=numN0+1
                                                                        rCount=rCount+1
                                                                        if(m==2) then
                                                                                Hsparse(numN0)=temp
                                                                                Hcol(numN0)=ijkp
                                                                                Hrow(ijk+1)=Hrow(ijk+1)+1
                                                                        end if
                                                                end if
                                                        end do
                                                end do
                                        end do
                                end do
                        end do
                end do
                if(m==1) then
                        allocate(Hsparse(1:numN0),Hcol(1:numN0),Hrow(1:sMax+2))
                        Hrow=0
                        Hrow(1)=1
                        numN0=0
                        row=0
                        rCount=1
                end if
        end do
        do i=2,sMax+1
                Hrow(i)=Hrow(i)+Hrow(i-1)
        end do

        print *, "Got H"
        m0=100
        allocate(EigenVals(m0),EigenVecs(sMAx,m0),res(m0))
        
        print *, Hsparse(1)
        
        !make handles for T and V
        !info = mkl_sparse_d_create_csr(TH,1,sMax,sMax,Trow,Trow+1,Tcol,Tsparse)
        !print *, "TH: ",info
        !info = mkl_sparse_d_create_csr(VH,1,sMax,sMax,Vrow,Vrow+1,Vcol,Vsparse)
        !print *, "VH: ",info

        !allocate(Hroww(numN0+sMax),Hcoll(numN0+sMax),Hsparsee(numN0+sMax))
        !info = mkl_sparse_d_create_csr(HH,1,sMax,sMax,Hroww,Hroww+1,Hcoll,Hsparsee)
        !print *, "HH: ",info
        !allocate(Hb(numN0),He(numN0),Hcol(numN0),Hsparse(numN0))
        !call mkl_sparse_d_export_csr(TH,indexing,HnumRow,HnumCol,Hb,He,Hcol,Hsparse)
        !print *, Hsparse(1), HnumRow, Hcol(1)

        !add T and V to get H
        !call mkl_dcsradd('n',2,3,sMax,sMax,Tsparse,Tcol,Trow,1d0,Vsparse,Vcol,Vrow,Hsparse,Hcol,Hind,10000000000000,infoh)
        !op=10
        !alpha = 1d0
        !print *, TH, VH, HH
        !info = mkl_sparse_d_add(op,TH,alpha,VH,HH)
        !print *, "added: ", info

        !allocate(Hb(numN0),He(numN0),Hcol(numN0),Hsparse(numN0))
        !get Hsparse from the handle HH
        !info = mkl_sparse_d_export_csr(HH,indexing,HnumRow,HnumCol,Hb,He,Hcol,Hsparse)
        !print *, "Hsparse: ",info
       ! if(associated(Hsparse))then
       !         print *, "true"
       ! else 
       !         print*, "fals"
       ! end if

        allocate(feastparam(64))
        !setup solver
        call feastinit(feastparam)
        !print *, fpm
        
        !call solver
        print *, size(Hsparse)
        print *, Hrow(sMax+1)-1
        print *, size(Hcol)

        feastparam(1)=1
        feastparam(2)=20
        feastparam(4)=3
        feastparam(17)=0
        call dfeast_scsrev('F',sMax,Hsparse,Hrow,Hcol,feastparam,epsout,loop,0d0,10d0,m0,EigenVals,EigenVecs,m,res,info)
        
        print *, "info: ",info

!       call MyDSYEV(Hmat,SZ,EigenVals,EigenVecs)
        

 !       do i=1,sz
         print *, EigenVals(1)
         print *, EigenVals(2)
         !print *, Tmat(1,1)
   !     end do

        close(1)
        close(2)
        close(3)
        close(4)
end program Solver


!rescales the nodes of size n from a to b
subroutine rescale(n,a,b,nodes)
implicit none
        integer :: n,i
        real*8 :: nodes(*)
        real*8 :: x1,xN,a,b
        x1= nodes(1)
        xN= nodes(n)
        do i=1,n
                if((xN-x1)*(b-a)==0) then
                        nodes(i)=0d0
                else 
                        nodes(i)=a+(nodes(i)-x1)/(xN-x1)*(b-a)
                end if
        end do
end subroutine rescale


!Does the lobatto product into array
subroutine DVRFunctions(holder,nodes,node,i,sz,val,weights)
implicit none
        integer :: i,j,sz
        real*8 :: nodes(sz),holder(sz),weights(sz)
        real*8 :: val,node

        do j=1,sz
                if((nodes(i)-nodes(j))==0d0) then
                    holder(j)=0d0
                else 
                    holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
                end if
        end do
        holder(i)=1d0
        val=PRODUCT(holder)*(weights(i)**(-.5d0))
end subroutine DVRFunctions
                

!find derivitave of the nodes
subroutine DifferentiateChi(holder,nodes,node,i,sz,Xsum,weights)
implicit none
        integer :: i,sz,j,k
        real*8 nodes(sz),holder(sz),weights(sz)
        real*8 val,node,Xsum

        Xsum=0
        do k=1,sz
                if(k /= i) then
                        do j=1,sz
                                if((nodes(i)-nodes(j))==0d0) then
                                        holder(j)=0d0
                                else 
                                        holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
                                end if
                        end do
                        holder(i)=1d0
                        holder(k)=1d0
                        if((nodes(i)-nodes(k))==0) then
                                val=0d0
                        else 
                                val=PRODUCT(holder)/((nodes(i)-nodes(k)))
                        end if
                        Xsum = Xsum + val
                end if
        end do
        Xsum=Xsum*(weights(i)**(-.5d0))
end subroutine DifferentiateChi


subroutine VPotHO(x,res)
implicit none
        real*8 x,res
        
        res = (0.5d0)*x**2
end subroutine VPotHO

subroutine Vsech2(x,res)
implicit none
        real*8 x,res

        res = -((1d0)/cosh(x))**2
end subroutine Vsech2

subroutine V3D(x,y,z,res)
implicit none
        real*8 x,y,z,res
        
        res = (0.5d0)*(x**2+y**2+z**2)
end subroutine V3D
