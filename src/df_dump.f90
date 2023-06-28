!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : df_dump_bz_1d
!!!           df_dump_bz_2d
!!!           df_dump_bz_3d
!!!           df_dump_fmesh
!!!           df_dump_bmesh
!!!           df_dump_dmft_g
!!!           df_dump_dmft_s
!!!           df_dump_dmft_h
!!!           df_dump_dual_g
!!!           df_dump_dual_s
!!!           df_dump_dual_b
!!!           df_dump_latt_g
!!!           df_dump_latt_s
!!!           df_dump_susc_c
!!!           df_dump_susc_s
!!! source  : df_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           06/29/2023 by li huang (last modified)
!!! purpose : dump key observables produced by the diagrammatic framework
!!!           for dynamical mean field theory to external files.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> dump data of k-mesh in brillouin zone                            <<<
!!========================================================================

!!
!! @sub df_dump_bz_1d
!!
!! write out the k-mesh for 1d lattice model.
!!
  subroutine df_dump_bz_1d(kx)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpts
     use control, only : nkp_x

     implicit none

!! external arguments
     ! k-mesh along x-axis
     real(dp), intent(in) :: kx(nkp_x)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! open data file: df.bz_1d.dat
     open(mytmp, file='df.bz_1d.dat', form='formatted', status='unknown')

     ! write it
     do i=1,nkp_x
         write(mytmp,'(i6,f16.8)') i, kx(i)
     enddo ! over i={1,nkp_x} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine df_dump_bz_1d

!!
!! @sub df_dump_bz_2d
!!
!! write out the k-mesh for 2d lattice model.
!!
  subroutine df_dump_bz_2d(kx, ky)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpts
     use control, only : nkp_x, nkp_y

     implicit none

!! external arguments
     ! k-mesh along x-axis
     real(dp), intent(in) :: kx(nkp_x)

     ! k-mesh along y-axis
     real(dp), intent(in) :: ky(nkp_y)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

     ! open data file: df.bz_2d.dat
     open(mytmp, file='df.bz_2d.dat', form='formatted', status='unknown')

     ! write it
     do i=1,nkp_x
         do j=1,nkp_y
             write(mytmp,'(2i6,f16.8)') i, j, kx(i), ky(j)
         enddo ! over j={1,nkp_y} loop
     enddo ! over i={1,nkp_x} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine df_dump_bz_2d

!!
!! @sub df_dump_bz_3d
!!
!! write out the k-mesh for 3d lattice model.
!!
  subroutine df_dump_bz_3d(kx, ky, kz)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpts
     use control, only : nkp_x, nkp_y, nkp_z

     implicit none

!! external arguments
     ! k-mesh along x-axis
     real(dp), intent(in) :: kx(nkp_x)

     ! k-mesh along y-axis
     real(dp), intent(in) :: ky(nkp_y)

     ! k-mesh along z-axis
     real(dp), intent(in) :: kz(nkp_z)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.bz_2d.dat
     open(mytmp, file='df.bz_2d.dat', form='formatted', status='unknown')

     ! write it
     do i=1,nkp_x
         do j=1,nkp_y
             write(mytmp,'(2i6,f16.8)') i, j, kx(i), ky(j)
         enddo ! over j={1,nkp_y} loop
     enddo ! over i={1,nkp_x} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine df_dump_bz_3d

!!========================================================================
!!>>> dump data of matsubara grid                                      <<<
!!========================================================================

!!
!! @sub df_dump_fmesh
!!
!! write out the fermionic matsubara grid.
!!
  subroutine df_dump_fmesh()
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nffrq

     implicit none

!! [body
!! body]

     return
  end subroutine df_dump_fmesh

!!
!! @sub df_dump_bmesh
!!
!! write out the bosonic matsubara grid.
!!
  subroutine df_dump_bmesh()
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nbfrq

     implicit none

!! [body
!! body]

     return
  end subroutine df_dump_bmesh

!!========================================================================
!!>>> dump data of local variables                                     <<<
!!========================================================================

!!
!! @sub df_dump_dmft_g
!!
!! write out impurity green's function in matsubara frequency space
!!
  subroutine df_dump_dmft_g(rmesh, grnf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! impurity green's function
     complex(dp), intent(in) :: grnf(nffrq,norbs)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

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

!! body]

     return
  end subroutine df_dump_dmft_g

!!
!! @sub df_dump_dmft_s
!!
!! write out impurity self-energy function in matsubara frequency space
!!
  subroutine df_dump_dmft_s(rmesh, sigf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! impurity self-energy function
     complex(dp), intent(in) :: sigf(nffrq,norbs)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

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

!! body]

     return
  end subroutine df_dump_dmft_s

!!
!! @sub df_dump_dmft_h
!!
!! write out impurity hybridization function in matsubara frequency space
!!
  subroutine df_dump_dmft_h(rmesh, hybf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! impurity hybridization function
     complex(dp), intent(in) :: hybf(nffrq,norbs)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

     ! open data file: df.dmft_h.dat
     open(mytmp, file='df.dmft_h.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_dmft_h

!!========================================================================
!!>>> dump data of dual variables                                      <<<
!!========================================================================

!!
!! @sub df_dump_dual_g
!!
!! write out dual green's function in matsubara frequency space
!!
  subroutine df_dump_dual_g(rmesh, grnd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! dual green's function
     complex(dp), intent(in) :: grnd(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.dual_g.dat
     open(mytmp, file='df.dual_g.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_dual_g

!!
!! @sub df_dump_dual_s
!!
!! write out dual self-energy function in matsubara frequency space
!!
  subroutine df_dump_dual_s(rmesh, sigd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! dual self-energy function
     complex(dp), intent(in) :: sigd(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.dual_s.dat
     open(mytmp, file='df.dual_s.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_dual_s

!!
!! @sub df_dump_dual_b
!!
!! write out dual bare green's function in matsubara frequency space
!!
  subroutine df_dump_dual_b(rmesh, wssd)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! dual bare green's function
     complex(dp), intent(in) :: wssd(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.dual_b.dat
     open(mytmp, file='df.dual_b.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_dual_b

!!========================================================================
!!>>> dump data of lattice variables                                   <<<
!!========================================================================

!!
!! @sub df_dump_latt_g
!!
!! write out lattice green's function in matsubara frequency space
!!
  subroutine df_dump_latt_g(rmesh, grnk)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! lattice green's function
     complex(dp), intent(in) :: grnk(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.latt_g.dat
     open(mytmp, file='df.latt_g.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_latt_g

!!
!! @sub df_dump_latt_s
!!
!! write out lattice self-energy function in matsubara frequency space
!!
  subroutine df_dump_latt_s(rmesh, sigk)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

     ! lattice self-energy function
     complex(dp), intent(in) :: sigk(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.latt_s.dat
     open(mytmp, file='df.latt_s.dat', form='formatted', status='unknown')

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

!! body]

     return
  end subroutine df_dump_latt_s

!!========================================================================
!!>>> dump data of dynamical quantities                                <<<
!!========================================================================

!!
!! @sub df_dump_susc_c
!!
!! write out charge susceptibility
!!
  subroutine df_dump_susc_c(rmesh, susc)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nbfrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh (bosonic)
     real(dp), intent(in)    :: rmesh(nbfrq)

     ! charge susceptibility
     complex(dp), intent(in) :: susc(nbfrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.susc_c.dat
     open(mytmp, file='df.susc_c.dat', form='formatted', status='unknown')

     ! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nbfrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), susc(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine df_dump_susc_c

!!
!! @sub df_dump_susc_s
!!
!! write out spin susceptibility
!!
  subroutine df_dump_susc_s(rmesh, susc)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nbfrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! matsubara frequency mesh (bosonic)
     real(dp), intent(in)    :: rmesh(nbfrq)

     ! spin susceptibility
     complex(dp), intent(in) :: susc(nbfrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! open data file: df.susc_s.dat
     open(mytmp, file='df.susc_s.dat', form='formatted', status='unknown')

     ! write it
     do k=1,nkpts
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# kpt:', k, '  orb:', j
             do i=1,nbfrq
                 write(mytmp,'(i6,5f16.8)') i, rmesh(i), susc(i,j,k), czero
             enddo ! over i={1,nffrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine df_dump_susc_s
