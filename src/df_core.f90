!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : df_run
!!!           df_std
!!!           df_ladder
!!!           df_dyson
!!!           df_schi
!!!           df_cchi
!!! source  : df_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/17/2019 by li huang (last modified)
!!! purpose : main subroutines for the dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub df_run
!!
!! core computational engine, it is used to dispatch the jobs
!!
  subroutine df_run()
     use control, only : isdia

     implicit none

     DF_CORE: &
     select case ( isdia )

         case (1)
             call df_std()

         case (2)
             call df_ladder()

         case default
             call s_print_error('df_run','this feature is not implemented')

     end select DF_CORE

     return
  end subroutine df_run

!!
!! @sub df_std
!!
!! implement the standard dual fermion approximation framework
!!
  subroutine df_std()
     implicit none

     return
  end subroutine df_std

!!
!! @sub df_ladder
!!
!! implement the ladder dual fermion approximation framework
!!
  subroutine df_ladder()
     use constants, only : dp
     use constants, only : one, half, czero
     use constants, only : mystd

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : ndfit, dfmix
     use control, only : beta
     use control, only : myid, master

     use context, only : fmesh, bmesh
     use context, only : dual_g, dual_s, dual_b
     use context, only : vert_d, vert_m

     implicit none

! local variables
! loop index for dual fermion iterations
     integer  :: it

! loop index for k-points
     integer  :: k

! loop index for orbitals
     integer  :: o

! loop index for bosonic frequency \nu
     integer  :: v

! loop index for fermionic frequency \omega
     integer  :: w

! status flag
     integer  :: istat

! current bosonic frequency
     real(dp) :: om

! dummy complex(dp) arrays, used to do fourier transformation
     complex(dp) :: vr(nkpts)
     complex(dp) :: gr(nkpts)

! two-particle bubble function
     complex(dp), allocatable :: g2  (:,:,:)

! shifted dual green's function
     complex(dp), allocatable :: gstp(:,:,:)

! new dual green's function
     complex(dp), allocatable :: gnew(:,:,:)

! ladder green's function, used to calculate dual self-energy function
     complex(dp), allocatable :: gvrt(:,:,:)

! matrix form for bubble function, \chi
     complex(dp), allocatable :: imat(:,:)

! matrix form for vertex function (magnetic channel, \gamma^m)
     complex(dp), allocatable :: mmat(:,:)

! matrix form for vertex function (density channel, \gamma^d)
     complex(dp), allocatable :: dmat(:,:)

! fully dressed vertex function, \Gamma 
     complex(dp), allocatable :: Gmat(:,:)

! allocate memory
     allocate(g2  (nffrq,norbs,nkpts), stat=istat)
     allocate(gstp(nffrq,norbs,nkpts), stat=istat)
     allocate(gnew(nffrq,norbs,nkpts), stat=istat)
     allocate(gvrt(nffrq,norbs,nkpts), stat=istat)

     allocate(imat(nffrq,nffrq),       stat=istat)
     allocate(mmat(nffrq,nffrq),       stat=istat)
     allocate(dmat(nffrq,nffrq),       stat=istat)
     allocate(Gmat(nffrq,nffrq),       stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_core','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     !!print *, nkpts
     !!print *, fmesh
     !!print *, bmesh
     !!do it=1,nffrq
     !!    print *, it, fmesh(it)
     !!    print *, dual_s(it,1,:)
     !!enddo
     !!
     !!vr(1) = czero
     !!do it=1,nffrq
     !!    print *, it, fmesh(it), abs(sum(dual_g(it,1,:))) / nkpts
     !!    vr(1) = vr(1) + abs(sum(dual_g(it,1,:))) / nkpts
     !!enddo
     !!print *, vr(1) / nffrq 
     !!
     !!STOP

!!========================================================================
!!>>> starting ladder dual fermion iteration                           <<<
!!========================================================================

     DF_LOOP: do it=1,ndfit

         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', it
         endif ! back if ( myid == master ) block

         V_LOOP: do v=1,nbfrq

             om = bmesh(v)
             if ( myid == master ) then ! only master node can do it
                 write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', om
             endif ! back if ( myid == master ) block

             if ( om == 0.0_dp ) then
                 gstp = dual_g
             else
                 call cat_fill_k(dual_g, gstp, om)
             endif
             call cat_dia_2d(dual_g, gstp, g2)
             gvrt = czero

             !!do w=1,nffrq
             !!    print *, w, fmesh(w)
             !!    print *, g2(w,1,:)
             !!enddo
             !!print *
             !!print *
             !!print *
             !!print *
             !!
             !!om = bmesh(4)
             !!print *, om
             !!call cat_fill_k(dual_g, gstp, om)
             !!call cat_dia_2d(dual_g, gstp, g2) 
             !!call cat_dia_2d(dual_g, dual_g, g2) 
             !!do w=1,nffrq
             !!    print *, w, fmesh(w)
             !!    print *, g2(w,1,:)
             !!enddo
             !!STOP

         O_LOOP: do o=1,norbs

             mmat = vert_m(:,:,v)
             dmat = vert_d(:,:,v)

             K_LOOP: do k=1,nkpts

                 call s_diag_z(nffrq, g2(:,o,k), imat)

                 call cat_bse_solver(imat, mmat, Gmat)
                 call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, half * 3.0_dp)
                 call cat_bse_iterator(1, one, imat, mmat, Gmat)
                 call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, -half * half * 3.0_dp)

                 call cat_bse_solver(imat, dmat, Gmat)
                 call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, half * 1.0_dp)
                 call cat_bse_iterator(1, one, imat, dmat, Gmat)
                 call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, -half * half * 1.0_dp)

             enddo K_LOOP

             do w=1,nffrq
                 call cat_fft_2d(+1, nkp_x, nkp_y, gvrt(w,o,:), vr)
                 call cat_fft_2d(-1, nkp_x, nkp_y, gstp(w,o,:), gr)
                 gr = vr * gr / real(nkpts * nkpts)
                 call cat_fft_2d(+1, nkp_x, nkp_y, gr, vr)
                 dual_s(w,o,:) = dual_s(w,o,:) + vr / beta
             enddo

         enddo O_LOOP

         if ( v == 7 ) then
         do w=1,nffrq
             print *, w, fmesh(w)
             print *, dual_s(w,1,:)
         enddo

         STOP
         endif

         enddo V_LOOP

         STOP

         call df_dyson(+1, gnew, dual_s, dual_b)
         call s_mix_z( size(gnew), dual_g, gnew, dfmix)

         dual_g = gnew
         dual_s = czero

         write(mystd,*)
     enddo DF_LOOP

!!========================================================================
!!>>> finishing ladder dual fermion iteration                          <<<
!!========================================================================

     call df_dyson(-1, dual_g, dual_s, dual_b)

     deallocate(gstp)
     deallocate(gnew)
     deallocate(gvrt)
     deallocate(g2)
     deallocate(imat)
     deallocate(mmat)
     deallocate(dmat)
     deallocate(Gmat)

     return
  end subroutine df_ladder

!!
!! @sub df_dyson
!!
!! try to calculate the dual green's function or self-energy function by
!! using the dyson equation
!!
  subroutine df_dyson(op, dual_g, dual_s, dual_b)
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
     integer, intent(in) :: op

     complex(dp), intent(inout) :: dual_g(nffrq,norbs,nkpts)
     complex(dp), intent(inout) :: dual_s(nffrq,norbs,nkpts)
     complex(dp), intent(inout) :: dual_b(nffrq,norbs,nkpts)

     if ( op == 1 ) then
         dual_g = one / ( one / dual_b - dual_s )
     else
         dual_s = one / dual_b - one / dual_g
     endif ! back if ( op == 1 ) block

     return
  end subroutine df_dyson

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
