!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : df_run
!!!           df_std
!!!           df_ladder
!!!           df_dyson
!!!           cat_bse_solver
!!!           cat_bse_iterator
!!!           cat_dia_1d
!!!           cat_dia_2d
!!!           cat_dia_3d
!!!           cat_fft_1d
!!!           cat_fft_2d
!!!           cat_fft_3d
!!! source  : df_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           06/30/2023 by li huang (last modified)
!!! purpose : main subroutines for the dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> core subroutines                                                 <<<
!!========================================================================

!!
!! @sub df_run
!!
!! core computational engine, it is used to dispatch the jobs.
!!
  subroutine df_run()
     use control, only : isdia
     use control, only : myid, master

     use context, only : kx, ky
     use context, only : fmesh, bmesh
     use context, only : dmft_d
     use context, only : dual_g, dual_s, dual_b
     use context, only : latt_g, latt_s
     use context, only : susc_c, susc_s

     implicit none

!! [body

     ! dispatch the jobs, decide which dual fermion engine should be used
     DF_CORE: &
     select case ( isdia )

         case (1) ! only standard 2nd diagrams are considered
             call df_std()

         case (2) ! only ladder diagrams are considered
             call df_ladder()

         case default
             call s_print_error('df_run','this feature is not implemented')

     end select DF_CORE

     ! now dual_s (dual self-energy function)
     ! and dual_g (dual green's function)
     ! are already updated, we can try to evaluate the other quantities.

     ! try to update lattice quantities
     call df_eval_latt_g()
     call df_eval_latt_s()

     ! try to update local hybridization function
     ! it can be fed back to the quantum impurity solver
     call df_eval_dmft_h()

     ! try to calculate charge susceptibility and spin susceptibility
     call df_eval_susc_c()
     call df_eval_susc_s()

     ! save the relevant data to external files. they are the local
     ! hybridization function, dual green's function, dual self-energy
     ! function, dual bath green's function, charge susceptibility, and
     ! spin susceptibility. only the master node can do this
     !
     if ( myid == master ) then
         call df_dump_bz_2d(kx, ky)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_fmesh(fmesh)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_bmesh(bmesh)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_dmft_h(fmesh, dmft_d)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_dual_g(fmesh, dual_g)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_dual_s(fmesh, dual_s)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_dual_b(fmesh, dual_b)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_latt_g(fmesh, latt_g)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_latt_s(fmesh, latt_s)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_susc_c(bmesh, susc_c)
     endif ! back if ( myid == master ) block
     !
     if ( myid == master ) then
         call df_dump_susc_s(bmesh, susc_s)
     endif ! back if ( myid == master ) block

!! body]

     return
  end subroutine df_run

!!
!! @sub df_std
!!
!! implement the standard dual fermion approximation framework. here, only
!! the standard second-order diagrams are taken into considerations.
!!
!!
  subroutine df_std()
     implicit none

!! [body

     CONTINUE

!! body]

     return
  end subroutine df_std

!!
!! @sub df_ladder
!!
!! implement the ladder dual fermion approximation framework. here, only
!! the ladder-type diagrams are taken into considerations.
!!
  subroutine df_ladder()
     use constants, only : dp
     use constants, only : zero, one, half, czero
     use constants, only : mystd

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : ndfit, dfmix
     use control, only : beta
     use control, only : myid, master

     use context, only : bmesh
     use context, only : dual_g, dual_s, dual_b
     use context, only : vert_d, vert_m

     implicit none

!! local variables
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

     ! difference between two dual green's functions
     real(dp) :: gdiff

     ! complex(dp) dummy arrays, for fast fourier transformation
     complex(dp) :: vr(nkpts)
     complex(dp) :: gr(nkpts)

     ! matrix form for bubble function, \chi
     complex(dp), allocatable :: imat(:,:)

     ! matrix form for vertex function (magnetic channel, \gamma^m)
     complex(dp), allocatable :: mmat(:,:)

     ! matrix form for vertex function (density channel, \gamma^d)
     complex(dp), allocatable :: dmat(:,:)

     ! fully dressed vertex function, \Gamma
     complex(dp), allocatable :: Gmat(:,:)

     ! two-particle bubble function
     complex(dp), allocatable :: g2  (:,:,:)

     ! shifted dual green's function
     complex(dp), allocatable :: gstp(:,:,:)

     ! new dual green's function
     complex(dp), allocatable :: gnew(:,:,:)

     ! ladder green's function, used to calculate dual self-energy function
     complex(dp), allocatable :: gvrt(:,:,:)

!! [body

     ! allocate memory
     allocate(imat(nffrq,nffrq),       stat=istat)
     allocate(mmat(nffrq,nffrq),       stat=istat)
     allocate(dmat(nffrq,nffrq),       stat=istat)
     allocate(Gmat(nffrq,nffrq),       stat=istat)

     allocate(g2  (nffrq,norbs,nkpts), stat=istat)
     allocate(gstp(nffrq,norbs,nkpts), stat=istat)
     allocate(gnew(nffrq,norbs,nkpts), stat=istat)
     allocate(gvrt(nffrq,norbs,nkpts), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_ladder','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

!!========================================================================
!!>>> starting ladder dual fermion iteration                           <<<
!!========================================================================

     DF_LOOP: do it=1,ndfit

         ! write out iteration information
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,2(a,i3))') &
                 'dual fermion iteration (ladder):', it, ' /', ndfit
             write(mystd,'(2X,a)') &
                 '-----------------------------------------------------'
         endif ! back if ( myid == master ) block

         V_LOOP: do v=1,nbfrq

             ! get bosonic frequency
             om = bmesh(v)
             !
             if ( myid == master ) then ! only master node can do it
                 write(mystd,'(2X,2(a,i3),a,f12.8,a)') &
                     '> bosonic frequency => ', v, ' /', nbfrq, ' (', om, ')'
             endif ! back if ( myid == master ) block

             ! calculate two-particle bubbles, g2
             call cat_fill_gk(dual_g, gstp, om)
             call cat_dia_2d(dual_g, gstp, g2)
             !
             if ( myid == master ) then ! only master node can do it
                 write(mystd,'(4X,a)') 'calculate two-particle bubbles'
             endif ! back if ( myid == master ) block

             ! extract impurity vertex functions for density and magnetic
             ! channels. in principles, they should be orbital-dependent.
             ! we will fix it later.
             mmat = vert_m(:,:,v)
             dmat = vert_d(:,:,v)

             ! actually, gvrt denotes
             !
             ! 3 ( \Gamma^{m} - 0.5\Gamma^{(2),m} ) +
             ! 1 ( \Gamma^{d} - 0.5\Gamma^{(2),d} )
             !
             ! it is used to calcuate dual self-energy function. next we
             ! try to calculate the fully dressed vertex function \Gamma
             ! and use it to build gvrt.
             gvrt = czero
             K_LOOP: do k=1,nkpts
                 O_LOOP1: do o=1,norbs

                     ! build diagonal matrix for the bubbles
                     call s_diag_z(nffrq, g2(:,o,k), imat)

                     ! evaluate the first term, magnetic channel
                     call cat_bse_solver(imat, mmat, Gmat)
                     call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, half * 3.0_dp)
                     call cat_bse_iterator(1, one, imat, mmat, Gmat)
                     call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, -half * half * 3.0_dp)

                     ! evaluate the second term, density channel
                     call cat_bse_solver(imat, dmat, Gmat)
                     call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, half * 1.0_dp)
                     call cat_bse_iterator(1, one, imat, dmat, Gmat)
                     call s_vecadd_z(nffrq, gvrt(:,o,k), Gmat, -half * half * 1.0_dp)

                 enddo O_LOOP1 ! over o={1,norbs} loop
             enddo K_LOOP ! over k={1,nkpts} loop
             !
             if ( myid == master ) then ! only master node can do it
                 write(mystd,'(4X,a)') 'solve bethe-salpeter equations'
             endif ! back if ( myid == master ) block

             ! now gvrt and gstp are used to calculate dual self-energy
             ! function via fast fourier transformation
             O_LOOP2: do o=1,norbs
                 W_LOOP: do w=1,nffrq

                     call cat_fft_2d(+1, nkp_x, nkp_y, gvrt(w,o,:), vr)
                     call cat_fft_2d(-1, nkp_x, nkp_y, gstp(w,o,:), gr)
                     !
                     gr = vr * gr / real(nkpts * nkpts)
                     !
                     call cat_fft_2d(+1, nkp_x, nkp_y, gr, vr)
                     !
                     dual_s(w,o,:) = dual_s(w,o,:) + vr / beta

                 enddo W_LOOP ! over w={1,nffrq} loop
             enddo O_LOOP2 ! over o={1,norbs} loop
             !
             if ( myid == master ) then ! only master node can do it
                 write(mystd,'(4X,a)') 'solve schwinger-dyson equations'
             endif ! back if ( myid == master ) block

         enddo V_LOOP ! over v={1,nbfrq} loop

         ! write out iteration information
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a)') &
                 '-----------------------------------------------------'
         endif ! back if ( myid == master ) block

         ! determine new dual green's function
         call df_dyson(+1, gnew, dual_s, dual_b)
         !
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a)') 'solve dyson equations'
         endif ! back if ( myid == master ) block

         ! try to mix old and new dual green's function
         call s_mix_z( size(gnew), dual_g, gnew, dfmix)
         gdiff = sum(abs(gnew - dual_g)) / size(gnew)
         !
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a)') "mix dual green's functions"
             write(mystd,'(2X,a,e24.16)') 'gdiff =', gdiff
         endif ! back if ( myid == master ) block

         ! reset dual green's function and dual self-energy function
         dual_g = gnew
         dual_s = czero

         ! write out iteration information
         if ( myid == master ) then
             write(mystd,*)
         endif ! back if ( myid == master ) block

     enddo DF_LOOP ! over it={1,ndfit} loop

!!========================================================================
!!>>> finishing ladder dual fermion iteration                          <<<
!!========================================================================

     ! solve dyson equation again to get final dual self-energy function
     call df_dyson(-1, dual_g, dual_s, dual_b)

     ! deallocate memory
     deallocate(imat)
     deallocate(mmat)
     deallocate(dmat)
     deallocate(Gmat)

     deallocate(g2)
     deallocate(gstp)
     deallocate(gnew)
     deallocate(gvrt)

!! body]

     return
  end subroutine df_ladder

!!
!! @sub df_dyson
!!
!! try to calculate the dual green's function or self-energy function by
!! solving the dyson equation
!!
  subroutine df_dyson(op, dual_g, dual_s, dual_b)
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

!! external arguments
     ! a flag, it indicates the one that should be calculated.
     ! if op = +1, calculate dual green's function
     ! if op = -1, calculate dual self-energy function
     integer, intent(in) :: op

     ! dual green's function
     complex(dp), intent(inout) :: dual_g(nffrq,norbs,nkpts)

     ! dual self-energy function
     complex(dp), intent(inout) :: dual_s(nffrq,norbs,nkpts)

     ! dual bath green's function
     complex(dp), intent(in)    :: dual_b(nffrq,norbs,nkpts)

!! [body

     if ( op == 1 ) then
         dual_g = one / ( one / dual_b - dual_s )
     else
         dual_s = one / dual_b - one / dual_g
     endif ! back if ( op == 1 ) block

!! body]

     return
  end subroutine df_dyson

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

!! external arguments
     ! two-particle bubble, \chi
     complex(dp), intent(in)  :: chiM(nffrq,nffrq)

     ! impurity vertex function, \gamma
     complex(dp), intent(in)  :: vrtM(nffrq,nffrq)

     ! fully dressed vertex function, \Gamma
     complex(dp), intent(out) :: GamM(nffrq,nffrq)

!! local variables
     ! unitary matrix
     complex(dp) :: Imat(nffrq,nffrq)

!! [body

     ! build unitary matrix
     call s_identity_z(nffrq, Imat)

     ! eval I - \gamma \chi
     GamM = Imat - matmul(vrtM,chiM)

     ! eval ( I - \gamma \chi )^-1
     call s_inv_z(nffrq, GamM)

     ! eval ( I - \gamma \chi )^-1 \gamma
     GamM = matmul(GamM, vrtM)

!! body]

     return
  end subroutine cat_bse_solver

!!
!! @sub cat_bse_iterator
!!
!! try to solve the bethe-salpeter equation iterately
!!
  subroutine cat_bse_iterator(niter, mix, chiM, vrtM, GamM)
     use constants, only : dp, one

     use control, only : nffrq

     implicit none

!! external arguments
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

!! local variables
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
         GamM = ( vrtM + matmul(Vchi,Vold) ) * mix + (one - mix) * Vold

         ! calculate the difference
         diff = abs(sum(GamM - Vold)) / real(nffrq * nffrq)

         ! update Vold with new GamM
         Vold = GamM

     enddo BSE_ITERATOR ! over it={1,niter} loop

!! body]

     return
  end subroutine cat_bse_iterator

!!========================================================================
!!>>> calculate bubble diagrams                                        <<<
!!========================================================================

!!
!! @sub cat_dia_1d
!!
!! calculate the two-particle bubble diagram, 1d version.
!!
  subroutine cat_dia_1d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x
     use control, only : beta

     implicit none

!! external arguments
     ! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

     ! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

     ! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

!! [body

     ! we have to make sure nkpts == nkp_x
     call s_assert2(nkpts == nkp_x, 'nkpts != nkp_x')

     do i=1,norbs
         do j=1,nffrq
             gk = gin(j,i,:)
             g1 = czero
             call cat_fft_1d(+1, nkp_x, gk, g1) ! gk -> gr

             gk = ginp(j,i,:)
             g2 = czero
             call cat_fft_1d(+1, nkp_x, gk, g2) ! gk -> gr

             gr = g1 * g2
             gk = czero
             call cat_fft_1d(-1, nkp_x, gr, gk) ! gr -> gk
             chiq(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     chiq = chiq / real(nkpts * beta)

!! body]

     return
  end subroutine cat_dia_1d

!!
!! @sub cat_dia_2d
!!
!! calculate the two-particle bubble diagram, 2d version.
!!
  subroutine cat_dia_2d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : beta

     implicit none

!! external arguments
     ! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

     ! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

     ! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

!! [body

     ! we have to make sure nkpts == nkp_x * nkp_y
     call s_assert2(nkpts == (nkp_x * nkp_y), 'nkpts != (nkp_x * nkp_y)')

     do i=1,norbs
         do j=1,nffrq
             gk = gin(j,i,:)
             g1 = czero
             call cat_fft_2d(+1, nkp_x, nkp_y, gk, g1) ! gk -> gr

             gk = ginp(j,i,:)
             g2 = czero
             call cat_fft_2d(+1, nkp_x, nkp_y, gk, g2) ! gk -> gr

             gr = g1 * g2
             gk = czero
             call cat_fft_2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             chiq(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     chiq = chiq / real(nkpts * nkpts * beta)

!! body]

     return
  end subroutine cat_dia_2d

!!
!! @sub cat_dia_3d
!!
!! calculate the two-particle bubble diagram, 3d version.
!!
  subroutine cat_dia_3d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y, nkp_z
     use control, only : beta

     implicit none

!! external arguments
     ! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

     ! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

     ! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

!! [body

     ! we have to make sure nkpts == nkp_x * nkp_y * nkp_z
     call s_assert2(nkpts == (nkp_x * nkp_y * nkp_z), 'nkpts != (nkp_x * nkp_y * nkp_z)')

     do i=1,norbs
         do j=1,nffrq
             gk = gin(j,i,:)
             g1 = czero
             call cat_fft_3d(+1, nkp_x, nkp_y, nkp_z, gk, g1) ! gk -> gr

             gk = ginp(j,i,:)
             g2 = czero
             call cat_fft_3d(+1, nkp_x, nkp_y, nkp_z, gk, g2) ! gk -> gr

             gr = g1 * g2
             gk = czero
             call cat_fft_3d(-1, nkp_x, nkp_y, nkp_z, gr, gk) ! gr -> gk
             chiq(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     chiq = chiq / real(nkpts * nkpts * nkpts * beta)

!! body]

     return
  end subroutine cat_dia_3d

!!========================================================================
!!>>> fast fourier transformations                                     <<<
!!========================================================================

!!
!! note:
!!
!! require fftw3 software package
!!

!!
!! @sub cat_fft_1d
!!
!! conduct fast fourier transformation in 1d
!!
  subroutine cat_fft_1d(op, nx, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     ! import fftw header file
     include 'fftw3.f03'

!! external arguments
     ! fft direction, forward or backward
     integer, intent(in) :: op

     ! size of operand
     integer, intent(in) :: nx

     ! operand
     complex(dp), intent(inout) :: fin(nx)
     complex(dp), intent(inout) :: fout(nx)

!! local variables
     ! fftw descriptor handler
     type(c_ptr) :: plan

!! [body

     select case (op)

         case (+1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_1d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

!! body]

     return
  end subroutine cat_fft_1d

!!
!! @sub cat_fft_2d
!!
!! conduct fast fourier transformation in 2d
!!
  subroutine cat_fft_2d(op, nx, ny, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     ! import fftw header file
     include 'fftw3.f03'

!! external arguments
     ! fft direction, forward or backward
     integer, intent(in) :: op

     ! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny

     ! operand
     complex(dp), intent(inout) :: fin(nx,ny)
     complex(dp), intent(inout) :: fout(nx,ny)

!! local variables
     ! fftw descriptor handler
     type(c_ptr) :: plan

!! [body

     select case (op)

         case (+1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_2d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

!! body]

     return
  end subroutine cat_fft_2d

!!
!! @sub cat_fft_3d
!!
!! conduct fast fourier transformation in 3d
!!
  subroutine cat_fft_3d(op, nx, ny, nz, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     ! import fftw header file
     include 'fftw3.f03'

!! external arguments
     ! fft direction, forward or backward
     integer, intent(in) :: op

     ! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny
     integer, intent(in) :: nz

     ! operand
     complex(dp), intent(inout) :: fin(nx,ny,nz)
     complex(dp), intent(inout) :: fout(nx,ny,nz)

!! local variables
     ! fftw descriptor handler
     type(c_ptr) :: plan

!! [body

     select case (op)

         case (+1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_3d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

!! body]

     return
  end subroutine cat_fft_3d
