!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : df_dump_grnf
!!!           df_dump_sigf
!!!           df_dump_hybf <<<---
!!!           df_dump_grnd
!!!           df_dump_sigd
!!!           df_dump_wssd <<<---
!!!           df_dump_grnk
!!!           df_dump_sigk <<<---
!!!           df_dump_v4_d
!!!           df_dump_v4_m
!!!           df_dump_v4_f <<<---
!!!           df_dump_schi
!!!           df_dump_cchi <<<---
!!! source  : df_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/17/2019 by li huang (last modified)
!!! purpose : dump key observables produced by the diagrammatic framework
!!!           for dynamical mean field theory to external files.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> dump data of local variables                                     <<<
!!========================================================================

!!
!! @sub df_dump_grnf
!!
!! write out impurity green's function in matsubara frequency space
!!
  subroutine df_dump_grnf(rmesh, grnf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity green's function
     complex(dp), intent(in) :: grnf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: df.dmft_g.dat
     open(mytmp, file='df.dmft_g.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), grnf(j,i), czero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_grnf

!!
!! @sub df_dump_sigf
!!
!! write out impurity self-energy function in matsubara frequency space
!!
  subroutine df_dump_sigf(rmesh, sigf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity self-energy function
     complex(dp), intent(in) :: sigf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: df.dmft_s.dat
     open(mytmp, file='df.dmft_s.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), sigf(j,i), czero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_sigf

!!
!! @sub df_dump_hybf
!!
!! write out impurity hybridization function in matsubara frequency space
!!
  subroutine df_dump_hybf(rmesh, hybf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity hybridization function
     complex(dp), intent(in) :: hybf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: dt.dmft_h.dat
     open(mytmp, file='dt.dmft_h.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), hybf(j,i), czero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_hybf

!!========================================================================
!!>>> dump data of dual variables                                      <<<
!!========================================================================

!!
!! @sub df_dump_grnd
!!
!! write out dual green's function in matsubara frequency space
!!
  subroutine df_dump_grnd(rmesh, grnd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! dual green's function
     complex(dp), intent(in) :: grnd(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! open data file: dt.dual_g.dat
     open(mytmp, file='dt.dual_g.dat', form='formatted', status='unknown')

! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nffrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), grnd(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_grnd

!!
!! @sub df_dump_sigd
!!
!! write out dual self-energy function in matsubara frequency space
!!
  subroutine df_dump_sigd(rmesh, sigd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! dual self-energy function
     complex(dp), intent(in) :: sigd(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! open data file: dt.dual_s.dat
     open(mytmp, file='dt.dual_s.dat', form='formatted', status='unknown')

! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nffrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), sigd(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_sigd

!!
!! @sub df_dump_wssd
!!
!! write out dual bare green's function in matsubara frequency space
!!
  subroutine df_dump_wssd(rmesh, wssd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! dual bare green's function
     complex(dp), intent(in) :: wssd(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! open data file: dt.dual_b.dat
     open(mytmp, file='dt.dual_b.dat', form='formatted', status='unknown')

! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nffrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), wssd(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_wssd

!!========================================================================
!!>>> dump data of lattice variables                                   <<<
!!========================================================================

!!
!! @sub df_dump_grnk
!!
!! write out lattice green's function in matsubara frequency space
!!
  subroutine df_dump_grnk(rmesh, grnk)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! lattice green's function
     complex(dp), intent(in) :: grnk(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! open data file: dt.latt_g.dat
     open(mytmp, file='dt.latt_g.dat', form='formatted', status='unknown')

! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nffrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), grnk(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_grnk

!!
!! @sub df_dump_sigk
!!
!! write out lattice self-energy function in matsubara frequency space
!!
  subroutine df_dump_sigk(rmesh, sigk)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! lattice self-energy function
     complex(dp), intent(in) :: sigk(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! open data file: dt.latt_s.dat
     open(mytmp, file='dt.latt_s.dat', form='formatted', status='unknown')

! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nffrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), sigk(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_sigk

!!========================================================================
!!>>> dump data of vertex functions                                    <<<
!!========================================================================

!!
!! @sub df_dump_v4_d
!!
!! write out vertex function (density channel) 
!!
  subroutine df_dump_v4_d()
     implicit none

     return
  end subroutine df_dump_v4_d

!!
!! @sub df_dump_v4_m
!!
!! write out vertex function (magnetic channel)
!!
  subroutine df_dump_v4_m()
     implicit none

     return
  end subroutine df_dump_v4_m

!!
!! @sub df_dump_v4_f
!!
!! write out vertex function (full vertex)
!!
  subroutine df_dump_v4_f()
     implicit none

     return
  end subroutine df_dump_v4_f

!!========================================================================
!!>>> dump data of dynamical quantities                                <<<
!!========================================================================

!!
!! @sub df_dump_schi
!!
!! write out spin susceptibility
!!
  subroutine df_dump_schi()
     implicit none

     return
  end subroutine df_dump_schi

!!
!! @sub df_dump_cchi
!!
!! write out charge susceptibility
!!
  subroutine df_dump_cchi()
     implicit none

     return
  end subroutine df_dump_cchi
