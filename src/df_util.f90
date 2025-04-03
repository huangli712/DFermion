!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : cat_fill_gl
!!!           cat_fill_gk
!!! source  : df_util.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 10/01/2008 by li huang (created)
!!!           04/03/2025 by li huang (last modified)
!!! purpose : provide some utility subroutines to deal with the GFs.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub cat_fill_gl
!!
!! try to fill G(\nu + \omega) by G(\nu), for local green's function only.
!!
  subroutine cat_fill_gl(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : beta

     use context, only : fmesh

     implicit none

!! external arguments
     ! shifted frequency, \omega
     real(dp), intent(in) :: shift

     ! input array, G(\nu)
     complex(dp), intent(in)  :: gin(nffrq,norbs)

     ! filled array, G(\nu + \omega)
     complex(dp), intent(out) :: gout(nffrq,norbs)

!! local variables
     ! loop index
     integer  :: i

     ! resultant index for \nu + \omega
     integer  :: k

     ! resultant frequency, \nu + \omega
     real(dp) :: w

!! [body

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:) = gin(k,:)
         else
             gout(i,:) = czero
         endif ! back if ( k >= 1 .and. k <= nffrq ) block
     enddo ! over i={1,nffrq} loop

!! body]

     return
  end subroutine cat_fill_gl

!!
!! @sub cat_fill_gk
!!
!! try to fill G(\nu + \omega, K) by G(\nu, K), for lattice green's
!! function only.
!!
  subroutine cat_fill_gk(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use control, only : beta

     use context, only : fmesh

     implicit none

!! external arguments
     ! shifted frequency, \omega
     real(dp), intent(in) :: shift

     ! input array, G(\nu, K)
     complex(dp), intent(in)  :: gin(nffrq,norbs,nkpts)

     ! filled array, G(\nu + \omega, K)
     complex(dp), intent(out) :: gout(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer  :: i

     ! resultant index for \nu + \omega
     integer  :: k

     ! resultant frequency, \nu + \omega
     real(dp) :: w

!! [body

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:,:) = gin(k,:,:)
         else
             gout(i,:,:) = czero
         endif ! back if ( k >= 1 .and. k < nffrq ) block
     enddo ! over i={1,nffrq} loop

!! body]

     return
  end subroutine cat_fill_gk
