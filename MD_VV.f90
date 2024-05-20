  Subroutine MD_VV(i_step)

  Use ctrl

  Implicit None
  integer,intent(in) :: i_step
  Integer :: iii,jjj,kkk
  Double precision :: rrr
  Double precision , dimension(nfr,1) :: cvec,fvec
  Double precision , dimension(3,natom) :: dedx
!  Double precision , dimension(1) :: egy
  Double precision :: egy
  Logical :: flag_reverse

! First stage
! Update velocity
  if (entype == 1) then
    Call NHC_VV(1)
  else if (entype == 0) then
    Call MD_VV_half()
  end if

! Update coordinates
  do iii = 1, natom
    do jjj = 1, 3
      Cmat(jjj,iii) = Cmat(jjj,iii) + Vmat(jjj,iii) * timestep
    end do
  end do

! Update force !! updated 2022/2/17 @@@@
  if(natom.eq.404) then
  Call pes_ag404(Cmat,egy,dedx)
  else
  stop "natom = 404"
  endif

  Ev = egy * eeVtoint

  do iii = 1, natom
    do jjj = 1, 3
      Fmat(jjj,iii) = dedx(jjj,iii) * eeVtoint
    end do
  end do

! Second stage
! Update velocity
  if (entype == 1) then
    Call NHC_VV(2)
  else if (entype == 0) then
    Call MD_VV_half()
  end if

  Call CalcTemp()

! Constrain !! updated 2022/2/17 @@@@
! !  flag_reverse = .false.
!   do iii = 1, natom
!     if ((Cmat(3,iii) > zupper).and.(Vmat(3,iii) > 0.0D0)) then
! !      flag_reverse = .true.
!       Vmat(3,iii) = -Vmat(3,iii)
!     end if
!   end do
  if ( Cmat(3,natom-1) .gt. zupper .and. Vmat(3,natom-1).gt.0.d0) then
  ! flag_reverse = .true.
    Vmat(3,natom-1) = -Vmat(3,natom-1) ! C
    Vmat(3,natom) = -Vmat(3,natom)     ! O
    write(uout,*) "Desorption (ps) ",i_step*timestep
  endif


  End Subroutine MD_VV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine MD_VV_half()

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision :: hstep

  hstep = 0.5D0 * timestep

  do iii = 1, natom
    do jjj = 1, 3
      Vmat(jjj,iii) = Vmat(jjj,iii) + Fmat(jjj,iii) / mass(iii) * hstep
    end do
  end do

  End Subroutine MD_VV_half

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine Calctemp()

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision :: rrr

  Ek = 0.0D0
  do iii = 1, natom
    do jjj = 1, 3
      Ek = Ek + Mass(iii) * Vmat(jjj,iii) ** 2
    end do
  end do

  tempnow = Ek / ektoint / Dble(nfr)
  Ek = Ek * 0.5D0

  End Subroutine Calctemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine RescaleV(cscale)

  Use ctrl

  Implicit None

  Double precision , intent(in) :: cscale

  Integer :: iii,jjj,kkk

  do iii = 1, natom
    do jjj = 1, 3
      Vmat(jjj,iii) = Vmat(jjj,iii) * cscale
    end do
  end do

  End Subroutine RescaleV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


