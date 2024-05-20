! --------------------------------------------------------------------
!     Q:    CARTESIAN IN ANGSTROM
!     VOUT: ENERGY IN (eV)
!     DVDQ: DERIVATIVES IN eV/ANGSTROM
! --------------------------------------------------------------------

! --------------------------------------------------------------------
subroutine pes_ag404(Q,VOUT,DVDQ)
  use nnmod_pes1, only: init_pes1,evaluate_pes1,deallocate_pes1
  use nnmod_pes2, only: init_pes2,evaluate_pes2,deallocate_pes2
  use nnmod_pes3, only: init_pes3,evaluate_pes3,deallocate_pes3
  use nnmod_pes12, only: init_pes12,evaluate_pes12,deallocate_pes12
  use nnmod_pes13, only: init_pes13,evaluate_pes13,deallocate_pes13
  use nnmod_pes22, only: init_pes22,evaluate_pes22,deallocate_pes22
  use nnmod_pes33, only: init_pes33,evaluate_pes33,deallocate_pes33
  use nnmod_pes23, only: init_pes23,evaluate_pes23,deallocate_pes23

  IMPLICIT NONE
  integer,parameter :: print_energy=0

  integer,parameter :: NATOM=404,natom1=44,natom2=12,natom3=12
  integer,parameter :: natom12=56,natom13=56
  integer,parameter :: natom22=24,natom33=24,natom23=24
  integer,parameter :: lig1=30,lig2=0  ! 4ftp30 jyh 9/7
  REAL*8,INTENT(IN) :: Q(3,NATOM)
  REAL*8,INTENT(OUT) :: VOUT,DVDQ(3,NATOM)

  REAL*8,PARAMETER :: EV=27.21138505D0
  REAL*8,PARAMETER :: TOEE=27.21138505D0*23.0605D0*0.04184D0
  REAL*8,PARAMETER :: TOANG=0.5291772083D0

  integer,parameter :: ifrac=0,iforce=1
  integer :: i,j,k,l,m,n
  integer,save :: ifirst=0
  character(kind=1,len=2),save :: symbols(NATOM),symbols1(natom1)       
  character(kind=1,len=2),save :: symbols2(natom2),symbols3(natom3)
  character(kind=1,len=2),save :: symbols12(natom12),symbols13(natom13)
  character(kind=1,len=2),save :: symbols22(natom22),symbols33(natom33)
  character(kind=1,len=2),save :: symbols23(natom23)
  real*8 :: q1(3,natom1),q2(3,natom2),q3(3,natom3)
  real*8 :: q12(3,natom12),q13(3,natom13)
  real*8 :: q22(3,natom22),q33(3,natom33),q23(3,natom23)
  real*8 :: f1(3,natom1),f2(3,natom2),f3(3,natom3) 
  real*8 :: f12(3,natom12),f13(3,natom13) 
  real*8 :: f22(3,natom22),f33(3,natom33),f23(3,natom23)    
  real*8 :: e1,e2,e3,e12,e13,e22,e33,e23
  real*8 :: e2tot,e3tot,e12tot,e13tot,e22tot,e33tot,e23tot
  real*8 :: f1every(3,natom1),f2every(3,natom2,lig1),f3every(3,natom3,lig2)
  real*8 :: f12every(3,natom12,lig1),f13every(3,natom13,lig2)
  real*8 :: f22every(3,natom22,lig1,lig1)
  real*8 :: f33every(3,natom33,lig2,lig2)
  real*8 :: f23every(3,natom23,lig1,lig2)
  real*8 :: q2all(3,natom2,lig1),q3all(3,natom3,lig2)
  real*8 :: q12all(3,natom12,lig1),q13all(3,natom13,lig2)
  real*8 :: q22all(3,natom22,lig1,lig1) !jyh 22-4-6 /22-4-20 lig1-1
  real*8 :: q33all(3,natom33,lig2,lig2) !jyh 22-4-6 /22-4-20 lig2-1
  real*8 :: q23all(3,natom23,lig1,lig2)
  real*8 :: f12tot_ag44(3,natom1),f13tot_ag44(3,natom1) !jyh 22-4-21
  real*8 :: ftot_ag44_0(3,natom1),ftot_ag44(3,natom1) !jyh 22-4-20
  real*8 :: f22tot_4ftp(3,natom2,lig1),f23tot_4ftp(3,natom2,lig1)
  real*8 :: ftot_4ftp_1(3,natom2,lig1),ftot_4ftp_2(3,natom2,lig1)
  real*8 :: ftot_4ftp(3,natom2,lig1)
  real*8 :: f33tot_2ftp(3,natom3,lig2),f23tot_2ftp(3,natom3,lig2)
  real*8 :: ftot_2ftp_1(3,natom3,lig2),ftot_2ftp_2(3,natom3,lig2)
  real*8 :: ftot_2ftp(3,natom3,lig2)
  real*8 :: f(3,NATOM)
  real*8 :: e2_every(lig1),e3_every(lig2),e12_every(lig1),e13_every(lig2)
  real*8 :: e22_every(lig1,lig1),e33_every(lig2,lig2),e23_every(lig1,lig2) !jyh
  real*8 :: m1(natom1),m1_sum,q1_com(3)  ! jyh 22/6/13
  real*8 :: m2(natom2),m2_sum,q2_com(3)
  real*8 :: m3(natom2),m3_sum,q3_com(3)
  real*8 :: m21(natom2),m21_sum,q21_com(3)
  real*8 :: m22(natom2),m22_sum,q22_com(3)
  real*8 :: m31(natom2),m31_sum,q31_com(3)
  real*8 :: m32(natom2),m32_sum,q32_com(3)
  real*8 :: rx,w0,w1     
  real*8,parameter :: pi=dacos(-1.d0)
  real*8,parameter :: em_lig1=-75.629226d0,em_cluster=-2460.71677066935d0 !4ftp30
  real*8,parameter :: em_lig2=-75.540533d0
  real*8 :: de1,de2,emax_lig1_lig1,emax_lig1_lig2,emax_lig2_lig2 ! jyh 9/6
  real*8 :: rc1,rc2  ! jyh 9/6

! 8 symbol fuzhi 22/6/15 
  if(ifirst.eq.0) then
    symbols1='Ag'
    symbols2=['H','H','H','H','C','C','C','C','C','C','F','S']
    symbols3=['H','H','H','H','C','C','C','C','C','C','F','S']
    symbols12(1:12)=symbols2; symbols12(13:56)=symbols1
    symbols13(1:12)=symbols3; symbols13(13:56)=symbols1
    symbols22(1:23:2)=symbols2; symbols22(2:24:2)=symbols2
    symbols33(1:23:2)=symbols3; symbols33(2:24:2)=symbols3
    symbols23(1:23:2)=symbols2; symbols23(2:24:2)=symbols3
    ifirst=1
  endif

! 8m fuzhi 
  m1(1:natom1)=107.8682
  m2(1:natom2)=[1.0079,1.0079,1.0079,1.0079,12.011,12.011,12.011,12.011,12.011,12.011,18.998,32.06] 
  m3(1:natom3)=[1.0079,1.0079,1.0079,1.0079,12.011,12.011,12.011,12.011,12.011,12.011,18.998,32.06]
  m21(1:natom2)=m2; m22(1:natom2)=m2
  m31(1:natom3)=m3; m32(1:natom3)=m3
  m1_sum=sum(m1(1:natom1)); m2_sum=sum(m2(1:natom2)); m3_sum=sum(m3(1:natom3))
  m21_sum=m2_sum; m22_sum=m2_sum
  m31_sum=m3_sum; m32_sum=m3_sum

! ag44 q->pes1 e1,fevery(f1)
  q1=Q(1:3,1:natom1)
  call evaluate_pes1(ifrac,iforce,natom1,symbols1,q1,e1,f1)
  f1every(1:3,1:natom1)=f1(1:3,1:natom1)

  if (print_energy.eq.1) write(320,*) e1

! 4ftp 28ge q->pes2 e2tot,f2every
  e2tot=0.d0
  q2all=reshape(Q(1:3,(natom1+1):(natom1+natom2*lig1)),[3,natom2,lig1])
  do i=1,lig1
    call evaluate_pes2(ifrac,iforce,natom2,symbols2,q2all(1:3,1:natom2,i),e2,f2)
    e2_every(i)=e2
    e2tot=e2tot+e2
    f2every(1:3,1:natom2,i)=f2

    if (print_energy.eq.1) then
       write(3300+i,*) e2
      !if(i.eq.18) then
      !  write(33000+i*10,*) natom2
      !  write(33000+i*10,*) e2
      !  do j=1,natom2
      !     write(33000+i*10,*) symbols2(j),q2all(1:3,j,i)
      !  enddo
      !endif
    endif
  end do

! 2ftp 2ge q->pes3 e3tot,f3every
  e3tot=0.d0
  q3all=reshape(Q(1:3,(natom1+natom2*lig1+1):NATOM),[3,natom3,lig2])
  do i=1,lig2    
    call evaluate_pes3(ifrac,iforce,natom3,symbols3,q3all(1:3,1:natom3,i),e3,f3)
    e3_every(i)=e3
    e3tot=e3tot+e3
    f3every(1:3,1:natom3,i)=f3

    if (print_energy.eq.1) write(3400+i,*) e3
  end do

! ag44_4ftp 1*28ge q->pes12 e12tot,f12every(3,56,28)
  e12tot=0.d0
  do i=1,lig1
     q12all(1:3,1:natom2,i)=q2all(1:3,1:natom2,i)
     q12all(1:3,(natom2+1):natom12,i)=q1(1:3,1:natom1)

     !! cut off ! jyh 22/6/13
     q1_com(1)=dot_product(q1(1,1:natom1),m1(1:natom1))/m1_sum
     q1_com(2)=dot_product(q1(2,1:natom1),m1(1:natom1))/m1_sum
     q1_com(3)=dot_product(q1(3,1:natom1),m1(1:natom1))/m1_sum
     q2_com(1)=dot_product(q2all(1,1:natom2,i),m2(1:natom2))/m2_sum
     q2_com(2)=dot_product(q2all(2,1:natom2,i),m2(1:natom2))/m2_sum
     q2_com(3)=dot_product(q2all(3,1:natom2,i),m2(1:natom2))/m2_sum

     rc1=8.d0; rc2=11.d0
     rx=dsqrt(dot_product(q1_com-q2_com,q1_com-q2_com))
     if(rx<rc1) then
       w0=1.d0
     else if(rx>rc2) then
       w0=0.d0
     else 
       w0=cos((rx-rc1)*pi/(rc2-rc1))*0.5d0+0.5d0
     end if

     de1=5.d0; de2=6.d0
     if(e2_every(i)<de1) then
       w1=1.d0
     else if(e2_every(i)>de2) then
       w1=0.d0
     else
       w1=cos((e2_every(i)-de1)*pi/(de2-de1))*0.5d0+0.5d0
     end if

     e12=0.d0; f12=0.d0
     if(w0*w1.gt.1.d-8) call evaluate_pes12(ifrac,iforce,natom12,symbols12,q12all(1:3,1:natom12,i),e12,f12)

     e12_every(i)=e12*w0*w1   !!  energy w1 jyh 22/6/16
     e12tot=e12tot+e12*w0*w1
     f12every(1:3,1:natom12,i)=f12*w0*w1

     if (print_energy.eq.1) write(3500+i,*) e12*w0*w1
  end do

! ag44_2ftp 1*2ge q->pes13 e13tot,f13every(3,56,2)
  e13tot=0.d0
  do i=1,lig2
     q13all(1:3,1:natom3,i)=q3all(1:3,1:natom3,i)
     q13all(1:3,(natom3+1):natom13,i)=q1(1:3,1:natom1)

     !! cut off ! jyh 22/6/13
     q1_com(1)=dot_product(q1(1,1:natom1),m1(1:natom1))/m1_sum
     q1_com(2)=dot_product(q1(2,1:natom1),m1(1:natom1))/m1_sum
     q1_com(3)=dot_product(q1(3,1:natom1),m1(1:natom1))/m1_sum
     q3_com(1)=dot_product(q3all(1,1:natom3,i),m3(1:natom3))/m3_sum
     q3_com(2)=dot_product(q3all(2,1:natom3,i),m3(1:natom3))/m3_sum
     q3_com(3)=dot_product(q3all(3,1:natom3,i),m3(1:natom3))/m3_sum

     rc1=8.d0; rc2=11.d0
     rx=dsqrt(dot_product(q1_com-q3_com,q1_com-q3_com))
     if(rx<rc1) then
       w0=1.d0
     else if(rx>rc2) then
       w0=0.d0
     else
       w0=cos((rx-rc1)*pi/(rc2-rc1))*0.5d0+0.5d0
     end if

     de1=5.d0; de2=6.d0
     if(e3_every(i)<de1) then
       w1=1.d0
     else if(e3_every(i)>de2) then
       w1=0.d0
     else
       w1=cos((e3_every(i)-de1)*pi/(de2-de1))*0.5d0+0.5d0
     end if

     e13=0.d0; f13=0.d0
     if(w0*w1.gt.1.d-8) call evaluate_pes13(ifrac,iforce,natom13,symbols13,q13all(1:3,1:natom13,i),e13,f13)

     e13_every(i)=e13*w0*w1  !!  energy w1 jyh 22/9/6
     e13tot=e13tot+e13*w0*w1
     f13every(1:3,1:natom13,i)=f13*w0*w1

     if (print_energy.eq.1) write(3600+i,*) e13*w0*w1
  end do

! 4ftp4ftp q->pes22 e:28*27/2 e22tot f:heng27*shu28 f22every(3,24,28) 
  do i=1,lig1-1 ! 2022-4-6 jyh
  do m=i+1,lig1 ! 2022-4-6 jyh
     q22all(1:3,1:natom2,m,i)=q2all(1:3,1:natom2,i)
     q22all(1:3,(natom2+1):natom22,m,i)=q2all(1:3,1:natom2,m)
  end do
  end do
  
  e22tot=0.d0
  do i=1,lig1-1 ! 2022-4-6 jyh
  do m=i+1,lig1 ! 2022-4-6 jyh
     q22=q22all(1:3,[1:4,13:16,5:10,17:22,11,23,12,24],m,i) 

     !! cut off ! jyh 22/6/13
     q21_com(1)=dot_product(q2all(1,1:natom2,i),m21(1:natom2))/m21_sum
     q21_com(2)=dot_product(q2all(2,1:natom2,i),m21(1:natom2))/m21_sum
     q21_com(3)=dot_product(q2all(3,1:natom2,i),m21(1:natom2))/m21_sum
     q22_com(1)=dot_product(q2all(1,1:natom2,m),m22(1:natom2))/m22_sum
     q22_com(2)=dot_product(q2all(2,1:natom2,m),m22(1:natom2))/m22_sum
     q22_com(3)=dot_product(q2all(3,1:natom2,m),m22(1:natom2))/m22_sum

     rc1=5.d0; rc2=7.d0
     rx=dsqrt(dot_product(q21_com-q22_com,q21_com-q22_com))
     if(rx<rc1) then
       w0=1.d0
     else if(rx>rc2) then
       w0=0.d0
     else
       w0=cos((rx-rc1)*pi/(rc2-rc1))*0.5d0+0.5d0
     end if

     emax_lig1_lig1=max(e2_every(i),e2_every(m))
     de1=4.d0; de2=5.d0
     if(emax_lig1_lig1<de1) then
       w1=1.d0
     else if(emax_lig1_lig1>de2) then
       w1=0.d0
     else
       w1=cos((emax_lig1_lig1-de1)*pi/(de2-de1))*0.5d0+0.5d0
     end if

     e22=0.d0; f22=0.d0
     if(w0*w1.gt.1.d-8) call evaluate_pes22(ifrac,iforce,natom22,symbols22,q22,e22,f22)

     e22_every(m,i)=e22*w0*w1
     e22tot=e22tot+e22*w0*w1  ! jide /2 ! 22-4-6 buyong/
     f22every(1:3,1:natom22,m,i)=f22(1:3,[1:4,9:14,21,23,5:8,15:20,22,24])*w0*w1

     if (print_energy.eq.1) write(440000+i*100+m,*) e22*w0*w1
  end do  
  end do

! 2ftp2ftp q->pes33 e:2*1/2 e33tot f:heng1*shu2 f33every(3,24,2) 
  do i=1,lig2-1 ! 22-4-6 jyh
  do m=i+1,lig2 ! 22-4-6 jyh
     q33all(1:3,1:natom3,m,i)=q3all(1:3,1:natom3,i)
     q33all(1:3,(natom3+1):natom33,m,i)=q3all(1:3,1:natom3,m)
  end do
  end do

  e33tot=0.d0
  do i=1,lig2-1 ! 22-4-6 jyh
  do m=i+1,lig2 ! 22-4-6 jyh
     q33=q33all(1:3,[1:4,13:16,5:10,17:22,11,23,12,24],m,i)       

     !! cut off ! jyh 22/6/13
     q31_com(1)=dot_product(q3all(1,1:natom3,i),m31(1:natom3))/m31_sum
     q31_com(2)=dot_product(q3all(2,1:natom3,i),m31(1:natom3))/m31_sum
     q31_com(3)=dot_product(q3all(3,1:natom3,i),m31(1:natom3))/m31_sum
     q32_com(1)=dot_product(q3all(1,1:natom3,m),m32(1:natom3))/m32_sum
     q32_com(2)=dot_product(q3all(2,1:natom3,m),m32(1:natom3))/m32_sum
     q32_com(3)=dot_product(q3all(3,1:natom3,m),m32(1:natom3))/m32_sum

     rc1=5.d0; rc2=7.d0
     rx=dsqrt(dot_product(q31_com-q32_com,q31_com-q32_com))
     if(rx<rc1) then
       w0=1.d0
     else if(rx>rc2) then
       w0=0.d0
     else
       w0=cos((rx-rc1)*pi/(rc2-rc1))*0.5d0+0.5d0
     end if

     emax_lig2_lig2=max(e3_every(i),e3_every(m))
     de1=4.d0; de2=5.d0
     if(emax_lig2_lig2<de1) then
       w1=1.d0
     else if(emax_lig2_lig2>de2) then
       w1=0.d0
     else
       w1=cos((emax_lig2_lig2-de1)*pi/(de2-de1))*0.5d0+0.5d0
     end if

     e33=0.d0; f33=0.d0
     if(w0*w1.gt.1.d-8) call evaluate_pes33(ifrac,iforce,natom33,symbols33,q33,e33,f33)

     e33_every(m,i)=e33*w0*w1
     e33tot=e33tot+e33*w0*w1  ! jide /2 ! 22-4-6 buyong/
     f33every(1:3,1:natom33,m,i)=f33(1:3,[1:4,9:14,21,23,5:8,15:20,22,24])*w0*w1

     if (print_energy.eq.1) write(450000+i*100+m,*) e33*w0*w1
  end do
  end do

! 2ftp4ftp q->pes23 e:2*28 e23tot f:heng28*shu2 f23every(3,24,28,2) 
  do i=1,lig2
  do m=1,lig1
     q23all(1:3,1:natom3,m,i)=q3all(1:3,1:natom3,i)
     q23all(1:3,(natom3+1):natom23,m,i)=q2all(1:3,1:natom2,m)
  end do
  end do

  e23tot=0.d0
  do i=1,lig2
  do m=1,lig1
     q23=q23all(1:3,[1:4,13:16,5:10,17:22,11,23,12,24],m,i)       

     !! cut off ! jyh 22/6/13
     q31_com(1)=dot_product(q3all(1,1:natom3,i),m31(1:natom3))/m31_sum
     q31_com(2)=dot_product(q3all(2,1:natom3,i),m31(1:natom3))/m31_sum
     q31_com(3)=dot_product(q3all(3,1:natom3,i),m31(1:natom3))/m31_sum
     q21_com(1)=dot_product(q2all(1,1:natom2,m),m21(1:natom2))/m21_sum
     q21_com(2)=dot_product(q2all(2,1:natom2,m),m21(1:natom2))/m21_sum
     q21_com(3)=dot_product(q2all(3,1:natom2,m),m21(1:natom2))/m21_sum

     rc1=5.d0; rc2=7.d0
     rx=dsqrt(dot_product(q31_com-q21_com,q31_com-q21_com))
     if(rx<rc1) then
       w0=1.d0
     else if(rx>rc2) then
       w0=0.d0
     else
       w0=cos((rx-rc1)*pi/(rc2-rc1))*0.5d0+0.5d0
     end if

     emax_lig1_lig2=max(e2_every(m),e3_every(i))
     de1=4.d0; de2=5.d0
     if(emax_lig1_lig2<de1) then
       w1=1.d0
     else if(emax_lig1_lig2>de2) then
       w1=0.d0 
     else
       w1=cos((emax_lig1_lig2-de1)*pi/(de2-de1))*0.5d0+0.5d0
     end if

     e23=0.d0; f23=0.d0
     if(w0*w1.gt.1.d-8) call evaluate_pes23(ifrac,iforce,natom23,symbols23,q23,e23,f23)

     e23_every(m,i)=e23*w0*w1
     e23tot=e23tot+e23*w0*w1
     f23every(1:3,1:natom23,m,i)=f23(1:3,[1:4,9:14,21,23,5:8,15:20,22,24])*w0*w1 

     if (print_energy.eq.1) write(460000+i*100+m,*) e23*w0*w1
  end do
  end do

! e:eight part
  VOUT=e1+e2tot+e3tot+e12tot+e13tot+e22tot+e33tot+e23tot ! jyh 22/6/16
  VOUT=VOUT-(em_cluster-em_lig1*lig1-em_lig2*lig2)    ! 73step energy jyh 22/6/15

  if (print_energy.eq.1) write(310,*) vout

 !VOUT=VOUT/EV*TOEE

! ag44 f huizong
 ! 28ge f12every_ag44 huizong 
  f12tot_ag44=0.d0
  do i=1,lig1 
    f12tot_ag44(1:3,1:natom1)=f12tot_ag44(1:3,1:natom1)+f12every(1:3,(natom2+1):natom12,i)
  end do
     
 ! 2ge f13every_ag44 huizong 
  f13tot_ag44=0.d0
  do i=1,lig2
    f13tot_ag44(1:3,1:natom1)=f13tot_ag44(1:3,1:natom1)+f13every(1:3,(natom3+1):natom13,i)
  end do

  ftot_ag44_0(1:3,1:natom1)=f12tot_ag44(1:3,1:natom1)+f13tot_ag44(1:3,1:natom1)
 !ftot_ag44_0(1:3,1:natom1)=f12tot_ag44(1:3,(natom2+1):natom12)+f13tot_ag44(1:3,(natom3+1):natom13)
  ftot_ag44(1:3,1:natom1)=f1every(1:3,1:natom1)+ftot_ag44_0(1:3,1:natom1)!cuo jyh

! 28ge 4ftp f(f12_4ftp f22_4ftp f23_4ftp) huizong
  do i=2,lig1
    do m=1,i-1
      f22every(1:3,1:natom2,m,i)=f22every(1:3,(natom2+1):natom22,i,m)
      f22every(1:3,(natom2+1):natom22,m,i)=f22every(1:3,1:natom2,i,m)
    end do
  end do
  
 !f22tot_4ftp=0.d0
  do i=1,lig1
   f22tot_4ftp(1:3,1:natom2,i)=0.d0 ! jyh 22-4-21 
   do m=1,lig1 
   if(m==i) cycle ! jyh 22-4-9
   f22tot_4ftp(1:3,1:natom2,i)=f22tot_4ftp(1:3,1:natom2,i)+f22every(1:3,1:natom2,m,i)
   end do
  end do

 !f23tot_4ftp=0.d0
  do i=1,lig1
   f23tot_4ftp(1:3,1:natom2,i)=0.d0 ! jyh 22-4-21
   do m=1,lig2
     f23tot_4ftp(1:3,1:natom2,i)=f23tot_4ftp(1:3,1:natom2,i)+f23every(1:3,(natom3+1):natom23,i,m)
   end do
  end do

  do i=1,lig1
   ftot_4ftp_1(1:3,1:natom2,i)=f2every(1:3,1:natom2,i)+f12every(1:3,1:natom2,i)
  ftot_4ftp_2(1:3,1:natom2,i)=f22tot_4ftp(1:3,1:natom2,i)+f23tot_4ftp(1:3,1:natom2,i)
   ftot_4ftp(1:3,1:natom2,i)=ftot_4ftp_1(1:3,1:natom2,i)+ftot_4ftp_2(1:3,1:natom2,i)
  end do

! 2ge 2ftp f(f13_2ftp f33_2ftp f23_2ftp) huizong
  do i=2,lig2
    do m=1,i-1
      f33every(1:3,1:natom3,m,i)=f33every(1:3,(natom3+1):natom33,i,m)
      f33every(1:3,(natom3+1):natom33,m,i)=f33every(1:3,1:natom3,i,m)
    end do
  end do
 
 !f33tot_2ftp=0.d0
  do i=1,lig2
   f33tot_2ftp(1:3,1:natom3,i)=0.d0 ! jyh 22-4-21
   do m=1,lig2
   if(m==i) cycle ! jyh 22-4-9
   f33tot_2ftp(1:3,1:natom3,i)=f33tot_2ftp(1:3,1:natom3,i)+f33every(1:3,1:natom3,m,i)
   end do
  end do

 !f23tot_2ftp=0.d0
  do i=1,lig2
   f23tot_2ftp(1:3,1:natom3,i)=0.d0 ! jyh 22-4-21
   do m=1,lig1
   f23tot_2ftp(1:3,1:natom3,i)=f23tot_2ftp(1:3,1:natom3,i)+f23every(1:3,1:natom3,m,i)
   end do
  end do

  do i=1,lig2
   ftot_2ftp_1(1:3,1:natom3,i)=f3every(1:3,1:natom3,i)+f13every(1:3,1:natom3,i)
  ftot_2ftp_2(1:3,1:natom3,i)=f33tot_2ftp(1:3,1:natom3,i)+f23tot_2ftp(1:3,1:natom3,i)
   ftot_2ftp(1:3,1:natom3,i)=ftot_2ftp_1(1:3,1:natom3,i)+ftot_2ftp_2(1:3,1:natom3,i)
  end do

! 404atom f
  f(1:3,1:natom1)=ftot_ag44(1:3,1:natom1)

  f(1:3,(natom1+1):(natom1+natom2*lig1))=reshape(ftot_4ftp(1:3,1:natom2,1:lig1),[3,natom2*lig1]) ! jyh 22-4-21 

  f(1:3,(natom1+natom2*lig1+1):NATOM)=reshape(ftot_2ftp(1:3,1:natom3,1:lig2),[3,natom3*lig2])  ! jyh 22-4-21 

  DVDQ=f
 !DVDQ=-DVDQ/(EV/TOANG)*toee/TOANG

  return

  END SUBROUTINE pes_ag404  

! call init_pes
  subroutine INIT_PES()
    use nnmod_pes1, only: init_pes1
    use nnmod_pes2, only: init_pes2
    use nnmod_pes3, only: init_pes3
    use nnmod_pes12, only: init_pes12
    use nnmod_pes13, only: init_pes13
    use nnmod_pes22, only: init_pes22
    use nnmod_pes33, only: init_pes33
    use nnmod_pes23, only: init_pes23

    implicit none
    call init_pes1
    call init_pes2
    call init_pes3
    call init_pes12
    call init_pes13
    call init_pes22
    call init_pes33
    call init_pes23 
    return
  end subroutine INIT_PES

! call deallocate_pes
  subroutine DEALLOCATE_PES
    use nnmod_pes1, only: deallocate_pes1
    use nnmod_pes2, only: deallocate_pes2
    use nnmod_pes3, only: deallocate_pes3
    use nnmod_pes12, only: deallocate_pes12
    use nnmod_pes13, only: deallocate_pes13
    use nnmod_pes22, only: deallocate_pes22
    use nnmod_pes33, only: deallocate_pes33
    use nnmod_pes23, only: deallocate_pes23

    implicit none
    call deallocate_pes1
    call deallocate_pes2
    call deallocate_pes3
    call deallocate_pes12
    call deallocate_pes13
    call deallocate_pes22
    call deallocate_pes33
    call deallocate_pes23  
    return
  end subroutine DEALLOCATE_PES
