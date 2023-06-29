!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
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
!!!           06/29/2026 by li huang (last modified)
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
!! note:
!!
!! dual_g and latt_g must be updated ahead of time. actually, dmft_h will
!! not be changed, the new hybridization function is stored in dmft_d.
!!
  subroutine df_eval_dmft_h()
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq

     use context, only : dmft_g, dmft_h, dmft_d
     use context, only : dual_g
     use context, only : latt_g

     implicit none

!! local variables
     ! loop index for fermionic frequency \omega
     integer :: i

     ! loop index for orbitals
     integer :: j

!! [body

     do j=1,norbs
         do i=1,nffrq
             associate ( zh => ( sum( dual_g(i,j,:) ) / sum( latt_g(i,j,:) ) ) )
                 dmft_d(i,j) = dmft_h(i,j) + one / dmft_g(i,j) * zh
             end associate
         enddo ! over i={1,nffrq} loop
     enddo ! over j={1,norbs} loop

!! body]

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
!! note:
!!
!! dual_g must be updated ahead of time. however, dmft_h is old.
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

!! local variables
     ! loop index for fermionic frequency \omega
     integer :: i

     ! loop index for orbitals
     integer :: j

     ! loop index for k-points
     integer :: k

!! [body

     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 associate ( zh => ( one / ( dmft_h(i,j) - ek(k) ) ) )
                     latt_g(i,j,k) =  zh + zh**2 / dmft_g(i,j)**2 * dual_g(i,j,k)
                 end associate
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

!! body]

     return
  end subroutine df_eval_latt_g

!!
!! @sub df_eval_latt_s
!!
!! calculate the lattice self-energy function within the dual fermion framework
!!
!! note:
!!
!! dual_s must be updated ahead of time. however, dmft_g and dmft_s are old.
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

!! local variables
     ! loop index for fermionic frequency \omega
     integer :: i

     ! loop index for orbitals
     integer :: j

     ! loop index for k-points
     integer :: k

!! [body

     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 associate ( val => ( dual_s(i,j,k) * dmft_g(i,j) + one ) )
                     latt_s(i,j,k) = dmft_s(i,j) + dual_s(i,j,k) / val
                 end associate
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

!! body]

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
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts

     use context, only : bmesh
     use context, only : susc_c
     use context, only : vert_d

     implicit none

!! local variables
     ! loop index for bosonic frequencies
     integer :: i

     ! status flag
     integer :: istat

     ! L(\omega, k)
     complex(dp), allocatable :: Lwq(:,:,:)

     ! convolution of dual green's function:
     ! --> \sum_{k} G_{d}(\omgea, k) G_{d}(\omega + \Omega, k + q)
     complex(dp), allocatable :: gd2(:,:,:)

     ! convolution of Lwq:
     ! --> \sum_{k} L(\omega, k) L(\omega + \Omega, k + q)
     complex(dp), allocatable :: gt2(:,:,:)

     ! convolution of lattice green's function:
     ! --> \sum_{k} G_{d}(\omgea, k) G_{d}(\omega + \Omega, k + q)
     complex(dp), allocatable :: gl2(:,:,:)

!! [body

     ! allocate memory
     allocate(Lwq(nffrq,norbs,nkpts), stat=istat)
     allocate(gd2(nffrq,norbs,nkpts), stat=istat)
     allocate(gt2(nffrq,norbs,nkpts), stat=istat)
     allocate(gl2(nffrq,norbs,nkpts), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_eval_susc_c','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! try to calculate L(\omega,k) at first
     call cat_susc_lwq(Lwq)

     V_LOOP: do i=1,nbfrq

         ! for a given bosonic frequency, try to calculate all of the
         ! necessary convolutions. in order to calculate gt2, Lwq is
         ! needed. on the other hand, gd2 and gl2 can be calculated with
         ! dual_g and latt_g, respectively.
         call cat_susc_conv( bmesh(i), Lwq, gd2, gt2, gl2 )

         ! try to calculate the orbital-resolved and k-resolved
         ! susceptibilities.
         call cat_susc_value( susc_c(i,:,:), vert_d(:,:,i), gd2, gt2, gl2 )

         ! save the susceptibilies, here one is the normalization factor
         ! note for charge susceptibility, the factor is one.
         susc_c(i,:,:) = susc_c(i,:,:) * one

     enddo V_LOOP ! over i={1,nbfrq} loop

     ! deallocate memory
     deallocate(Lwq)
     deallocate(gd2)
     deallocate(gt2)
     deallocate(gl2)

!! body]

     return
  end subroutine df_eval_susc_c

!!
!! @sub df_eval_susc_s
!!
!! calculate the spin susceptibility within the dual fermion framework
!!
  subroutine df_eval_susc_s()
     use constants, only : dp
     use constants, only : half

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts

     use context, only : bmesh
     use context, only : susc_s
     use context, only : vert_m

     implicit none

!! local variables
     ! loop index for bosonic frequencies
     integer :: i

     ! status flag
     integer :: istat

     ! L(\omega, k)
     complex(dp), allocatable :: Lwq(:,:,:)

     ! convolution of dual green's function:
     ! --> \sum_{k} G_{d}(\omgea, k) G_{d}(\omega + \Omega, k + q)
     complex(dp), allocatable :: gd2(:,:,:)

     ! convolution of Lwq:
     ! --> \sum_{k} L(\omega, k) L(\omega + \Omega, k + q)
     complex(dp), allocatable :: gt2(:,:,:)

     ! convolution of lattice green's function:
     ! --> \sum_{k} G_{d}(\omgea, k) G_{d}(\omega + \Omega, k + q)
     complex(dp), allocatable :: gl2(:,:,:)

!! [body

     ! allocate memory
     allocate(Lwq(nffrq,norbs,nkpts), stat=istat)
     allocate(gd2(nffrq,norbs,nkpts), stat=istat)
     allocate(gt2(nffrq,norbs,nkpts), stat=istat)
     allocate(gl2(nffrq,norbs,nkpts), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_eval_susc_s','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! try to calculate L(\omega,k) at first
     call cat_susc_lwq(Lwq)

     V_LOOP: do i=1,nbfrq

         ! for a given bosonic frequency, try to calculate all of the
         ! necessary convolutions. in order to calculate gt2, Lwq is
         ! needed. on the other hand, gd2 and gl2 can be calculated with
         ! dual_g and latt_g, respectively.
         call cat_susc_conv( bmesh(i), Lwq, gd2, gt2, gl2 )

         ! try to calculate the orbital-resolved and k-resolved
         ! susceptibilities.
         call cat_susc_value( susc_s(i,:,:), vert_m(:,:,i), gd2, gt2, gl2 )

         ! save the susceptibilies, here one is the normalization factor
         ! note for spin susceptibility, the factor is half.
         susc_s(i,:,:) = susc_s(i,:,:) * half

     enddo V_LOOP ! over i={1,nbfrq} loop

     ! deallocate memory
     deallocate(Lwq)
     deallocate(gd2)
     deallocate(gt2)
     deallocate(gl2)

!! body]

     return
  end subroutine df_eval_susc_s

!!========================================================================
!!>>> service layer                                                    <<<
!!========================================================================

!!
!! @sub cat_susc_lwq
!!
!! try to calculate L(\omega, k), a key quantity for the calculation of
!! spin and charge susceptibilities:
!!
!! L(\omega,k) = G_{dual} G_{latt} / G^{0}_{dual}
!!
  subroutine cat_susc_lwq(Lwq)
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : ek
     use context, only : dmft_g, dmft_h
     use context, only : dual_g, dual_b

     implicit none

!! external arguments
     ! the quantity what we want to calculate, L(\omega,k)
     complex(dp), intent(out) :: Lwq(nffrq,norbs,nkpts)

!! local variables
     ! loop indices
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! try to calculate OLD G_{latt}
     ! since latt_g was already updated at df_eval_latt_g(), here we have
     ! to recompute it. note that dmft_g and dmft_h were old.
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 Lwq(i,j,k) = one / ( one / dmft_g(i,j) + dmft_h(i,j) - ek(k) )
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     Lwq = Lwq / dual_b
     Lwq = Lwq * dual_g

!! body]

     return
  end subroutine cat_susc_lwq

!!
!! @sub cat_susc_conv
!!
!! try to calculate convolution between some lattice quantities
!!
  subroutine cat_susc_conv(omega, Lwq, gd2, gt2, gl2)
     use constants, only : dp

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : dual_g
     use context, only : latt_g

     implicit none

!! external arguments
     ! given bosonic frequency
     real(dp), intent(in) :: omega

     ! L(\omega, k)
     complex(dp), intent(in)  :: Lwq(nffrq,norbs,nkpts)

     ! convolution of dual green's function
     complex(dp), intent(out) :: gd2(nffrq,norbs,nkpts)

     ! convolution of Lwq
     complex(dp), intent(out) :: gt2(nffrq,norbs,nkpts)

     ! convolution of lattice green's function
     complex(dp), intent(out) :: gl2(nffrq,norbs,nkpts)

!! local variables
     ! status flag
     integer :: istat

     ! shifted dual green's function
     complex(dp), allocatable :: gstp(:,:,:)

!! [body

     ! allocate memory
     allocate(gstp(nffrq,norbs,nkpts), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('cat_susc_conv','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! gd2 means the convolution of two dual green's functions
     call cat_fill_gk(dual_g, gstp, omega)
     call cat_dia_2d(dual_g, gstp, gd2)

     ! gt2 means the convolution of two L_{\Omega,\omega}{q}
     call cat_fill_gk(Lwq, gstp, omega)
     call cat_dia_2d(Lwq, gstp, gt2)

     ! gl2 means the convolution of two lattice green's functions
     call cat_fill_gk(latt_g, gstp, omega)
     call cat_dia_2d(latt_g, gstp, gl2)

     ! deallocate memory
     deallocate(gstp)

!! body]

     return
  end subroutine cat_susc_conv

!!
!! @sub cat_susc_value
!!
!! try to calculate the susceptibility for a given bosonic frequency
!!
  subroutine cat_susc_value(susc, vert, gd2, gt2, gl2)
     use constants, only : dp
     use constants, only : cone, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! orbital- and momentum-resolved spin or charge susceptibility
     complex(dp), intent(out) :: susc(norbs,nkpts)

     ! matrix form for vertex function
     complex(dp), intent(in)  :: vert(nffrq,nffrq)

     ! convolution of dual green's function
     complex(dp), intent(in)  :: gd2(nffrq,norbs,nkpts)

     ! convolution of Lwq
     complex(dp), intent(in)  :: gt2(nffrq,norbs,nkpts)

     ! convolution of lattice green's function
     complex(dp), intent(in)  :: gl2(nffrq,norbs,nkpts)

!! local parameters
     ! a flag, it denotes whether the contribution from the lattice
     ! bubble will be included in the calculation of susceptibility.
     complex(dp), parameter :: add_lattice_bubble = czero

!! local variables
     ! loop indices
     integer :: i
     integer :: k

     ! status flag
     integer :: istat

     ! dummy complex(dp) vector
     complex(dp), allocatable :: yvec(:)

     ! matrix form for bubble function (convolution of dual green's function)
     complex(dp), allocatable :: imat(:,:)

     ! fully dressed vertex function, \Gamma
     complex(dp), allocatable :: Gmat(:,:)

!! [body

     ! allocate memory
     allocate(yvec(nffrq)      , stat=istat)
     allocate(imat(nffrq,nffrq), stat=istat)
     allocate(Gmat(nffrq,nffrq), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('cat_susc_value','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     O_LOOP: do i=1,norbs
         K_LOOP: do k=1,nkpts

             ! step 1: try to solve the bethe-salpeter equation to
             ! obtain \Gamma_{m/d}.
             call s_diag_z(nffrq, gd2(:,i,k), imat)
             call cat_bse_solver(imat, vert, Gmat)

             ! step 2: get susceptibility from \Gamma_{m/d} and
             ! convolution of dual green's function.
             yvec = czero
             call zgemv('N', nffrq, nffrq, cone, Gmat, nffrq, gt2(:,i,k), 1, czero, yvec, 1)
             susc(i,k) = dot_product(gt2(:,i,k), yvec)

             ! step 3: add up the contribution from lattice bubble
             ! to susc.
             susc(i,k) = susc(i,k) + sum(gl2(:,i,k)) * add_lattice_bubble

         enddo K_LOOP ! over k={1,nkpts} loop
     enddo O_LOOP ! over i={1,norbs} loop

     ! deallocate memory
     deallocate(yvec)
     deallocate(imat)
     deallocate(Gmat)

!! body]

     return
  end subroutine cat_susc_value
