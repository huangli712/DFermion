!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : df_eval_dmft_h
!!!           df_eval_latt_g
!!!           df_eval_latt_s
!!!           df_eval_susc_c
!!!           df_eval_susc_s
!!!           cat_susc_lwq
!!!           cat_susc_conv
!!!           cat_susc_value
!!! source  : df_eval.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 04/29/2009 by li huang (created)
!!!           05/31/2019 by li huang (last modified)
!!! purpose : try to evaluate some key observables.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> evaluate local observables                                       <<<
!!========================================================================

!!
!! @sub df_eval_dmft_h
!!
!! calculate the local hybridization function within the dual fermion framework
!!
  subroutine df_eval_dmft_h()
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq

     use context, only : dmft_g, dmft_h, dmft_y
     use context, only : dual_g
     use context, only : latt_g

     implicit none

! local variables
! loop index for fermionic frequency \omega
     integer :: i

! loop index for orbitals
     integer :: j

!!
!! note:
!!
!! dual_g and latt_g must be updated ahead of time.
!!
!! actually, dmft_h will not be changed, the new hybridization function
!! is stored in dmft_y.
!!
     do j=1,norbs
         do i=1,nffrq
             associate ( zh => ( sum( dual_g(i,j,:) ) / sum( latt_g(i,j,:) ) ) )
                 dmft_y(i,j) = dmft_h(i,j) + one / dmft_g(i,j) * zh 
             end associate
         enddo ! over i={1,nffrq} loop
     enddo ! over j={1,norbs} loop

     return
  end subroutine df_eval_dmft_h

!!========================================================================
!!>>> evaluate lattice observables                                     <<<
!!========================================================================

!!
!! @sub df_eval_latt_g
!!
!! calculate the lattice green's function within the dual fermion framework
!!
  subroutine df_eval_latt_g()
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : ek
     use context, only : dmft_g, dmft_h
     use context, only : dual_g
     use context, only : latt_g

     implicit none

! local variables
! loop index for fermionic frequency \omega
     integer :: i

! loop index for orbitals
     integer :: j

! loop index for k-points
     integer :: k

!!
!! note:
!!
!! dual_g must be updated ahead of time. however, dmft_h is old.
!!
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 associate ( zh => ( one / ( dmft_h(i,j) - ek(k) ) ) )
                     latt_g(i,j,k) =  zh + zh**2 / dmft_g(i,j)**2 * dual_g(i,j,k)
                 end associate
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine df_eval_latt_g

!!
!! @sub df_eval_latt_s
!!
!! calculate the lattice self-energy function within the dual fermion framework
!!
  subroutine df_eval_latt_s()
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : dmft_g, dmft_s
     use context, only : dual_s
     use context, only : latt_s

     implicit none

! local variables
! loop index for fermionic frequency \omega
     integer :: i

! loop index for orbitals
     integer :: j

! loop index for k-points
     integer :: k

!!
!! note:
!!
!! dual_s must be updated ahead of time. however, dmft_g and dmft_s are old.
!!
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 associate ( val => ( dual_s(i,j,k) * dmft_g(i,j) + one ) )
                     latt_s(i,j,k) = dmft_s(i,j) + dual_s(i,j,k) / val
                 end associate
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine df_eval_latt_s

!!========================================================================
!!>>> evaluate q-dependent susceptibilities                            <<<
!!========================================================================

!!
!! @sub df_eval_susc_c
!!
!! calculate the charge susceptibility within the dual fermion framework
!!
  subroutine df_eval_susc_c()
     use constants, only : dp, one

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts

     use context, only : bmesh
     use context, only : susc_c
     use context, only : vert_d

     implicit none

! local variables
     integer :: i

     complex(dp), allocatable :: Lwq(:,:,:)
     complex(dp), allocatable :: gd2(:,:,:)
     complex(dp), allocatable :: gt2(:,:,:)
     complex(dp), allocatable :: gl2(:,:,:)

     allocate(Lwq(nffrq,norbs,nkpts))
     allocate(gd2(nffrq,norbs,nkpts))
     allocate(gt2(nffrq,norbs,nkpts))
     allocate(gl2(nffrq,norbs,nkpts))

     call cat_susc_lwq(Lwq)

     do i=1,nbfrq 
         call cat_susc_conv( bmesh(i), Lwq, gd2, gt2, gl2 )
         call cat_susc_value( susc_c(i,:,:), vert_d(:,:,i), gd2, gt2, gl2 )
         susc_c(i,:,:) = susc_c(i,:,:) * one
     enddo ! over i={1,nbfrq} loop

     do i=1,nkpts
         print *, i, susc_c(1,1,i)
     enddo

     deallocate(Lwq)
     deallocate(gd2)
     deallocate(gt2)
     deallocate(gl2)

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

!!
!! @sub cat_susc_lwq
!!
!! try to calculate L_{\Omega,\omega}(q)
!!
  subroutine cat_susc_lwq(Lwq)
     use constants, only : dp, one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : ek
     use context, only : dmft_g, dmft_h
     use context, only : dual_g, dual_b

     implicit none

! external arguments
     complex(dp), intent(out) :: Lwq(nffrq,norbs,nkpts)

! local variables
     integer :: i
     integer :: j
     integer :: k

!!
!! L_{\Omega,\omega}(q) = G_{dual} G_{latt} / G^{0}_{dual}
!!

! try to calculate G_{latt}
! since latt_g was already updated at df_eval_latt_g(), here we have to
! recompute it. note that dmft_g and dmft_h were old.
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 Lwq(i,j,k) = one / ( one / dmft_g(i,j) + dmft_h(i,j) - ek(k) ) 
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     Lwq = Lwq / dual_b
     Lwq = Lwq * dual_g

     return
  end subroutine cat_susc_lwq

!!
!! @sub cat_susc_conv
!!
!! try to calculate convolution between some lattice quantities
!!
  subroutine cat_susc_conv(omega, Lwq, gd2, gt2, gl2)
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : dual_g
     use context, only : latt_g

     implicit none

! external arguments
     real(dp), intent(in) :: omega
     complex(dp), intent(in)  :: Lwq(nffrq,norbs,nkpts)
     complex(dp), intent(out) :: gd2(nffrq,norbs,nkpts)
     complex(dp), intent(out) :: gt2(nffrq,norbs,nkpts)
     complex(dp), intent(out) :: gl2(nffrq,norbs,nkpts)

! local variables
     complex(dp), allocatable :: gstp(:,:,:)

! allocate memory
     allocate(gstp(nffrq,norbs,nkpts))

! gd2 means the convolution of two dual green's functions
     if ( omega == zero ) then
         gstp = dual_g
     else
         call cat_fill_k(dual_g, gstp, omega)
     endif
     call cat_dia_2d(dual_g, gstp, gd2)

! gt2 means the convolution of two L_{\Omega,\omega}{q}
     if ( omega == zero ) then
         gstp = Lwq
     else
         call cat_fill_k(Lwq, gstp, omega)
     endif
     call cat_dia_2d(Lwq, gstp, gt2)

! gl2 means the convolution of two lattice green's functions
     if ( omega == zero ) then
         gstp = latt_g
     else
         call cat_fill_k(latt_g, gstp, omega)
     endif
     call cat_dia_2d(latt_g, gstp, gl2)

! deallocate memory
     deallocate(gstp)

     return
  end subroutine cat_susc_conv

!!
!! @sub cat_susc_value
!!
!! try to calculate the susceptibility for a given bosonic frequency
!!
  subroutine cat_susc_value(susc, vert, gd2, gt2, gl2)
     use constants, only : dp, cone, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
     complex(dp), intent(in)  :: vert(nffrq,nffrq)
     complex(dp), intent(out) :: susc(norbs,nkpts)

     complex(dp), intent(in) :: gd2(nffrq,norbs,nkpts)
     complex(dp), intent(in) :: gt2(nffrq,norbs,nkpts)
     complex(dp), intent(in) :: gl2(nffrq,norbs,nkpts)

! local variables
     integer :: i
     integer :: k

     complex(dp), allocatable :: yvec(:)
     complex(dp), allocatable :: imat(:,:)
     complex(dp), allocatable :: Gmat(:,:)

! allocate memory
     allocate(yvec(nffrq))
     allocate(imat(nffrq,nffrq))
     allocate(Gmat(nffrq,nffrq))

     do i=1,norbs
         do k=1,nkpts
             call s_diag_z(nffrq, gd2(:,i,k), imat)
             call cat_bse_solver(imat, vert, Gmat)

             yvec = czero
             call zgemv('N', nffrq, nffrq, cone, Gmat, nffrq, gt2(:,i,k), 1, czero, yvec, 1)
             susc(i,k) = dot_product(gt2(:,i,k), yvec)
         enddo ! over k={1,nkpts} loop
     enddo ! over i={1,norbs} loop

! deallocate memory
     deallocate(yvec)
     deallocate(imat)
     deallocate(Gmat)

     return
  end subroutine cat_susc_value
