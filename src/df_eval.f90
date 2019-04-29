
!!
!! @sub df_schi
!!
!! calculate the spin susceptibility within the dual fermion framework
!!
  subroutine df_schi()
     implicit none

     return
  end subroutine df_schi

!!
!! @sub df_cchi
!!
!! calculate the charge susceptibility within the dual fermion framework
!!
  subroutine df_cchi()
     implicit none

     return
  end subroutine df_cchi

  subroutine df_eval_latt()
     use constants
     use control
     use context

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

     !!do w=1,nffrq
     !!    print *, w, fmesh(w)
     !!    print *, latt_g(w,1,:)
     !!enddo

     do j=1,norbs
         do i=1,nffrq
             dmft_h(i,j) = dmft_h(i,j) + one / dmft_g(i,j) * sum(dual_g(i,j,:)) / sum(latt_g(i,j,:))
         enddo
     enddo

     do w=1,nffrq
         print *, w, fmesh(w), dmft_h(w,1)
     enddo

     STOP 'in df_eval_latt'
     return
  end subroutine df_eval_latt
