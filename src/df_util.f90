!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : cat_fill_l
!!!           cat_fill_k <<<---
!!!           cat_fft_1d
!!!           cat_fft_2d
!!!           cat_fft_3d <<<---
!!!           cat_dia_1d
!!!           cat_dia_2d
!!!           cat_dia_3d <<<---
!!!           cat_bse_solver
!!!           cat_bse_iterator
!!! source  : dt_util.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           01/02/2018 by li huang (last modified)
!!! purpose : provide some utility subroutines, such as FFT, convolution,
!!!           and Bethe-Salpter equation solver, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> shift matsubara frequency                                        <<<
!!========================================================================

!!
!! @sub cat_fill_l
!!
!! try to fill G(\nu + \omega) by G(\nu), momentum-independent version
!!
  subroutine cat_fill_l(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : beta

     use context, only : fmesh

     implicit none

! external arguments
! shifted frequency, \omega
     real(dp), intent(in) :: shift

! input array, G(\nu)
     complex(dp), intent(in)  :: gin(nffrq,norbs)

! filled array, G(\nu + \omega)
     complex(dp), intent(out) :: gout(nffrq,norbs)

! local variables
! loop index
     integer  :: i

! resultant index for \nu + \omega
     integer  :: k

! resultant frequency, \nu + \omega
     real(dp) :: w

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:) = gin(k,:)
         else
             gout(i,:) = czero
         endif ! back if ( k >= 1 .and. k <= nffrq ) block
     enddo ! over i={1,nffrq} loop

     return
  end subroutine cat_fill_l

!!
!! @sub cat_fill_k
!!
!! try to fill G(\nu + \omega, K) by G(\nu, K), momentum-dependent version
!!
  subroutine cat_fill_k(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use control, only : beta

     use context, only : fmesh

     implicit none

! external arguments
! shifted frequency, \omega
     real(dp), intent(in) :: shift

! input array, G(\nu, K)
     complex(dp), intent(in)  :: gin(nffrq,norbs,nkpts)

! filled array, G(\nu + \omega, K)
     complex(dp), intent(out) :: gout(nffrq,norbs,nkpts)

! local variables
! loop index
     integer  :: i

! resultant index for \nu + \omega
     integer  :: k

! resultant frequency, \nu + \omega
     real(dp) :: w

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:,:) = gin(k,:,:)
         else
             gout(i,:,:) = czero
         endif ! back if ( k >= 1 .and. k <= nffrq ) block
     enddo ! over i={1,nffrq} loop

     return
  end subroutine cat_fill_k
