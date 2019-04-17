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
!!>>> solve bethe-salpeter equation                                    <<<
!!========================================================================

!!
!! note:
!!
!! the bethe-salpeter equation reads
!!
!!     \Gamma = \gamma + \gamma \chi \Gamma
!!
!! or equivalently
!!
!!     1 / \gamma - 1 / \Gamma = \chi
!!
!! here \Gamma is called the fully dressed vertex function, \gamma is the
!! impurity vertex function, and \chi is the two-particle bubble. we can
!! solve it directly by matrix inversion, or iterately.
!!

!!
!! @sub cat_bse_solver
!!
!! try to solve the bethe-salpeter equation directly by matrix inversion
!!
!!     \Gamma = \gamma / ( I - \gamma \chi )
!!
  subroutine cat_bse_solver(chiM, vrtM, GamM)
     use constants, only : dp

     use control, only : nffrq

     implicit none

! external arguments
! two-particle bubble, \chi
     complex(dp), intent(in)  :: chiM(nffrq,nffrq)

! impurity vertex function, \gamma
     complex(dp), intent(in)  :: vrtM(nffrq,nffrq)

! fully dressed vertex function, \Gamma
     complex(dp), intent(out) :: GamM(nffrq,nffrq)

! local variables
! unit matrix
     complex(dp) :: Imat(nffrq,nffrq)

! build unit matrix
     call s_identity_z(nffrq, Imat)

! eval I - \gamma \chi
     GamM = Imat - matmul(vrtM,chiM)

! eval ( I - \gamma \chi )^-1 
     call s_inv_z(nffrq, GamM)

! eval ( I - \gamma \chi )^-1 \gamma
     GamM = matmul(GamM, vrtM)

     return
  end subroutine cat_bse_solver

!!
!! @sub cat_bse_iterator
!!
!! try to solve the bethe-salpeter equation iterately
!!
  subroutine cat_bse_iterator(niter, mix, chiM, vrtM, GamM)
     use constants, only : dp

     use control, only : nffrq

     implicit none

! external arguments
! number of iteration
     integer, intent(in)  :: niter

! mixing factor
     real(dp), intent(in) :: mix

! two-particle bubble, \chi
     complex(dp), intent(in)  :: chiM(nffrq,nffrq)

! impurity vertex function, \gamma
     complex(dp), intent(in)  :: vrtM(nffrq,nffrq)

! fully dressed vertex function, \Gamma
     complex(dp), intent(out) :: GamM(nffrq,nffrq)

! local variables
! loop index
     integer :: it

! difference between two successive iterations
     complex(dp) :: diff

! used to save old GamM
     complex(dp) :: Vold(nffrq,nffrq)

! used to save \gamma * \chi
     complex(dp) :: Vchi(nffrq,nffrq)

! init Vold and Vchi
     Vold = vrtM
     Vchi = matmul(vrtM, chiM)

     BSE_ITERATOR: do it=1,niter

! calculate new \Gamma, and then mix it with Vold
         GamM = ( vrtM + matmul(Vchi,Vold) ) * mix + (1.0_dp - mix) * Vold

! calculate the difference
         diff = abs(sum(GamM - Vold)) / real(nffrq * nffrq)

! update Vold with new GamM
         Vold = GamM

     enddo BSE_ITERATOR ! over it={1,niter} loop 

     return
  end subroutine cat_bse_iterator
