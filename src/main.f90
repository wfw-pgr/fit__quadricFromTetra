program main
  implicit none
  integer                       :: itet, iv
  integer, parameter            :: cLen = 300
  integer, parameter            :: lun  = 50
  integer                       :: nDim, nVert, nTetra
  character(cLen)               :: cmt
  character(cLen)               :: inpFile="dat/tetra.dat"
  character(cLen)               :: outFile="dat/coeff.dat"
  double precision, allocatable :: tetra(:,:,:), coeff(:,:)
  
  ! ------------------------------------------------------ !
  ! --- [1] main                                       --- !
  ! ------------------------------------------------------ !
  open (lun,file=trim(inpFile),status="old")
  read(lun,*) cmt
  read(lun,*) cmt
  read(lun,*) cmt, nTetra, nVert, nDim
  allocate( tetra(nDim,nVert,nTetra), coeff(6,nTetra) )
  tetra(:,:,:) = 0.d0
  coeff(:,:)   = 0.d0
  write(6,*) nDim, nVert, nTetra
  do itet=1, nTetra
     do iv=1, nVert
        read (lun,*) tetra(1:3,iv,itet)
     enddo
  enddo
  close(lun)
  
  ! ------------------------------------------------------ !
  ! --- [2] calculate coefficient                      --- !
  ! ------------------------------------------------------ !
  call fit__quadricFromTetra( tetra, coeff, nTetra )

  ! ------------------------------------------------------ !
  ! --- [3] return                                     --- !
  ! ------------------------------------------------------ !
  open (lun,file=trim(outFile),status="replace")
  do itet=1, nTetra
     write(lun,"(6(e15.8,1x))") coeff(1:6,itet)
  enddo
  close(lun)

contains

  
  ! ====================================================== !
  ! === calculate coefficient                          === !
  ! ====================================================== !
  subroutine fit__quadricFromTetra( tetra, coeff, nTetra )
    implicit none
    integer         , intent(in)  :: nTetra
    double precision, intent(in)  :: tetra(3,4,nTetra)
    double precision, intent(out) :: coeff(6,  nTetra)
    integer                       :: i, j, itet
    double precision, allocatable :: Amat(:,:,:), bvec(:,:)
    integer         , parameter   :: xx_=1, yy_=2, xy_=3, xm_=4, ym_=5, cm_=6
    integer         , parameter   :: x_=1, y_=2, z_=3, iG_=1
    integer         , parameter   :: nSize=6

    
    ! ------------------------------------------------------ !
    ! --- [1] setting of coefficient matrix              --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] preparation                              --  !
    allocate( Amat(6,6,nTetra), bvec(6,nTetra) )
    Amat(:,:,:) = 0.d0
    bvec(:,:)   = 0.d0
    !  -- [1-2]  Fix points                              --  !
    do itet=1, nTetra
       do i=1, 4
          Amat(i,xx_,itet) = tetra(x_,i,itet)**2
          Amat(i,yy_,itet) = tetra(y_,i,itet)**2
          Amat(i,xy_,itet) = tetra(x_,i,itet)*tetra(y_,i,itet)
          Amat(i,xm_,itet) = tetra(x_,i,itet)
          Amat(i,ym_,itet) = tetra(y_,i,itet)
          Amat(i,cm_,itet) = 1.d0
       enddo
    enddo
    !  -- [1-3]  gradiant = 0  @ top                     --  !
    do itet=1, nTetra
       Amat(5,xx_,itet) = 2.d0 * tetra(x_,iG_,itet)
       Amat(5,xy_,itet) =        tetra(y_,iG_,itet)
       Amat(5,xm_,itet) = 1.d0
       Amat(6,yy_,itet) = 2.d0 * tetra(y_,iG_,itet)
       Amat(6,xy_,itet) =        tetra(x_,iG_,itet)
       Amat(6,ym_,itet) = 1.d0
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] setting of the right hand side             --- !
    ! ------------------------------------------------------ !
    do itet=1, nTetra
       do i=1, 4
          bvec(i,itet) = tetra(z_,i,itet)
       enddo
       bvec(5:6,itet)  = 0.d0
    enddo
    ! ------------------------------------------------------ !
    ! --- [3] solve equation                             --- !
    ! ------------------------------------------------------ !
    coeff(:,:) = 0.d0
    do itet=1, nTetra
       call solve__gaussElimin( Amat(:,:,itet), coeff(:,itet), bvec(:,itet), nSize )
    enddo
    return
  end subroutine fit__quadricFromTetra


  ! ========================================================== !
  ! === Gauss Elimination Solver                           === !
  ! ========================================================== !
  subroutine solve__gaussElimin( Amat, xvec, bvec, nSize )
    implicit none
    integer         , intent(in)  :: nSize
    double precision, intent(in)  :: Amat(nSize,nSize)
    double precision, intent(in)  :: bvec(nSize)
    double precision, intent(out) :: xvec(nSize)
    integer                       :: i, j, k, ipivot
    double precision              :: Dinv, buff, vpivot
    double precision, parameter   :: eps = 1.d-10
    double precision, allocatable :: Umat(:,:), vvec(:,:)

    ! ----------------------------------------- !
    ! --- [1] Preparation                   --- !
    ! ----------------------------------------- !
    !  -- [1-1] allocate Umat               --  !
    allocate( Umat(nSize,nSize) )
    !  -- [1-2] Copy AMat & bvec            --  !
    Umat(:,:) = Amat(:,:)
    xvec(:)   = bvec(:)
    
    ! ----------------------------------------- !
    ! --- [2] Forward Ellimination          --- !
    ! ----------------------------------------- !
    do k=1, nSize

       !  -- [2-1] Pivoting                 --  !
       vpivot = abs( Umat(k,k) )
       ipivot = k
       do j=k+1, nSize
          if ( abs( Umat(j,k) ).gt.vpivot ) then
             vpivot = abs( Umat(j,k) )
             ipivot = j
          endif
       end do
       if ( ipivot.ne.k ) then
          do j=k, nSize
             buff           = Umat(ipivot,j)
             Umat(ipivot,j) = Umat(k     ,j)
             Umat(k     ,j) = buff
          enddo
          buff         = xvec(ipivot)
          xvec(ipivot) = xvec(k)
          xvec(k)      = buff
       end if
       if ( abs( Umat(k,k) ).lt.eps ) then
          write(6,*) '[gaussElimin] Amat :: Singular Matrix :: No Solution End :: @ k= ', k
          stop
       endif
       !  -- [2-2] Diagonal Component       --  !
       Dinv      = 1.d0 / Umat(k,k)
       Umat(k,k) = 1.d0

       !  -- [2-3] Non-Diagonal Component   --  !
       if ( k.eq.nSize ) then
          ! -- [    Last Row :: k == nSize ] -- !
          xvec(k) = Dinv * xvec(k)
       else
          ! -- [Not Last Row :: k != nSize ] -- !
          !  - Division    -  !
          Umat(k,k+1:nSize) = Dinv * Umat(k,k+1:nSize)
          xvec(k)           = Dinv * xvec(k)
          !  - subtraction -  !
          do j=k+1,nSize
             Umat(j,k+1:nSize) = Umat(j,k+1:nSize) - Umat(j,k) * Umat(k,k+1:nSize)
             xvec(j)           = xvec(j)           - Umat(j,k) * xvec(k)
             Umat(j,k)         = 0.d0
          enddo
       endif
       
    end do

    ! ----------------------------------------- !
    ! --- [3] Backward Substituition        --- !
    ! ----------------------------------------- !
    do k=nSize-1, 1, -1
       do i=nSize, k+1, -1
          xvec(k) = xvec(k) - Umat(k,i)*xvec(i)
       enddo
    enddo

    return
  end subroutine solve__gaussElimin

    
end program main
