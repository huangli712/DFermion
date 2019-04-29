!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : df_eval_latt_g
!!!           df_eval_latt_s
!!!           df_eval_susc_c
!!!           df_eval_susc_s
!!! source  : df_eval.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 04/29/2009 by li huang (created)
!!!           04/29/2019 by li huang (last modified)
!!! purpose : main subroutines for the dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  subroutine df_eval_latt_g()
     use constants, only : one

     use control, only : norbs

     use context, only : fmesh

     implicit none

     integer :: i, j, k, w

     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 latt_g(i,j,k) =  ( one / ( dmft_h(i,j) - ek(k) ) ) + &
                                  ( one / ( dmft_h(i,j) - ek(k) ) ) / &
                                  dmft_g(i,j) * dual_g(i,j,k) / dmft_g(i,j) * &
                                  ( one / ( dmft_h(i,j) - ek(k) ) ) 
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     do j=1,norbs
         do i=1,nffrq
             dmft_h(i,j) = dmft_h(i,j) + one / dmft_g(i,j) * sum(dual_g(i,j,:)) / sum(latt_g(i,j,:))
         enddo
     enddo

     do w=1,nffrq
         print *, w, fmesh(w), dmft_h(w,1)
     enddo

     STOP 'in df_eval_latt_g'

     return
  end subroutine df_eval_latt_g

!!
!! @sub df_eval_susc_c
!!
!! calculate the charge susceptibility within the dual fermion framework
!!
  subroutine df_eval_susc_c()
     implicit none

     return
  end subroutine df_eval_susc_c

!!
!! @sub df_eval_susc_s
!!
!! calculate the spin susceptibility within the dual fermion framework
!!
  subroutine df_eval_susc_s()
     implicit none

     return
  end subroutine df_eval_susc_s
