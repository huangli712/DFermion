!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : cat_dia_1d
!!!           cat_dia_2d
!!!           cat_dia_3d
!!! source  : df_diagram.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           04/17/2019 by li huang (last modified)
!!! purpose : try to calculate the bubble diagram (convolution).
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub cat_dia_1d
!!
!! calculate the two-particle bubble diagram, 1d version
!!
  subroutine cat_dia_1d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x
     use control, only : beta

     implicit none

! external arguments
! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

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

     return
  end subroutine cat_dia_1d

!!
!! @sub cat_dia_2d
!!
!! calculate the two-particle bubble diagram, 2d version
!!
  subroutine cat_dia_2d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : beta

     implicit none

! external arguments
! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

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

     return
  end subroutine cat_dia_2d

!!
!! @sub cat_dia_3d
!!
!! calculate the two-particle bubble diagram, 3d version
!!
  subroutine cat_dia_3d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y, nkp_z
     use control, only : beta

     implicit none

! external arguments
! G(\nu, K)
     complex(dp), intent(in)  :: gin (nffrq,norbs,nkpts)

! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

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

     return
  end subroutine cat_dia_3d

!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : cat_fft_1d
!!!           cat_fft_2d
!!!           cat_fft_3d
!!! source  : df_fft.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           06/24/2023 by li huang (last modified)
!!! purpose : try to implement the fast fourier transformation.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

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
