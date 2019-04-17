!!========================================================================
!!>>> fast fourier transformation                                      <<<
!!========================================================================

!!
!! note:
!!
!! need fftw3 software package
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

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx

! operand
     complex(dp), intent(inout) :: fin(nx)
     complex(dp), intent(inout) :: fout(nx)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

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

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny

! operand
     complex(dp), intent(inout) :: fin(nx,ny)
     complex(dp), intent(inout) :: fout(nx,ny)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

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

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny
     integer, intent(in) :: nz

! operand
     complex(dp), intent(inout) :: fin(nx,ny,nz)
     complex(dp), intent(inout) :: fout(nx,ny,nz)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

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

     return
  end subroutine cat_fft_3d
