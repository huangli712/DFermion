!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : df_print_header
!!!           df_print_footer
!!!           df_print_summary
!!!           df_print_control
!!!           df_print_runtime
!!!           df_print_it_info
!!! source  : df_print.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/15/2009 by li huang (created)
!!!           04/03/2025 by li huang (last modified)
!!! purpose : provide printing infrastructure for dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub df_print_header
!!
!! print the startup information for dual fermion framework
!!
  subroutine df_print_header()
     use constants, only : mystd

     use version, only : V_FULL
     use version, only : V_AUTH
     use version, only : V_INST
     use version, only : V_MAIL
     use version, only : V_GPL3

     use control, only : cname
     use control, only : nprocs

     implicit none

!! local variables
     ! string for current date and time
     character (len = 20) :: date_time_string

!! [body

     ! obtain current date and time
     call s_time_builder(date_time_string)

# if defined (MPI)

     write(mystd,'(2X,a)') cname//' (Parallelized Edition)'

# else   /* MPI */

     write(mystd,'(2X,a)') cname//' (Sequential Edition)'

# endif  /* MPI */

     write(mystd,'(2X,a)') 'A Modern Dual Fermion Framework For Quantum Lattice Models'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH//' ('//V_INST//')'
     write(mystd,'(2X,a)') 'Support: '//V_MAIL
     write(mystd,'(2X,a)') 'License: '//V_GPL3
     write(mystd,*)

     write(mystd,'(2X,a)') 'start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', 1

# endif  /* MPI */

!! body]

     return
  end subroutine df_print_header

!!
!! @sub df_print_footer
!!
!! print the ending information for dual fermion framework
!!
  subroutine df_print_footer()
     use constants, only : dp
     use constants, only : mystd

     use control, only : cname

     implicit none

!! local variables
     ! string for current date and time
     character (len = 20) :: date_time_string

     ! used to record the time usage information
     real(dp) :: tot_time

!! [body

     ! obtain time usage information
     call cpu_time(tot_time)

     ! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') cname//' >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') cname//' >>> happy ending at '//date_time_string

!! body]

     return
  end subroutine df_print_footer

!!
!! @sub df_print_summary
!!
!! print the running parameters, only for reference
!!
  subroutine df_print_summary()
     use constants, only : mystd

     use control ! ALL

     implicit none

!! [body

     write(mystd,'(2X,a)') '[configuration parameters] -> core control'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'isdia  / value :', isdia , 'type : i'

     write(mystd,'(2X,a)') '[configuration parameters] -> lattice model'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nband  / value :', nband , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nspin  / value :', nspin , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'norbs  / value :', norbs , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nkpts  / value :', nkpts , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nkp_x  / value :', nkp_x , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nkp_y  / value :', nkp_y , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nkp_z  / value :', nkp_z , 'type : i'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'mune   / value :', mune  , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'beta   / value :', beta  , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'part   / value :', part  , 'type : d'

     write(mystd,'(2X,a)') '[configuration parameters] -> dual fermion framework'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nffrq  / value :', nffrq , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nbfrq  / value :', nbfrq , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ndfit  / value :', ndfit , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nbsit  / value :', nbsit , 'type : i'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'dfmix  / value :', dfmix , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'bsmix  / value :', bsmix , 'type : d'

     write(mystd,*)

!! body]

     return
  end subroutine df_print_summary

!!
!! @sub df_print_control
!!
!! print the control parameters, only for reference
!!
  subroutine df_print_control()
     implicit none

!! [body

     CONTINUE

!! body]

     return
  end subroutine df_print_control

!!
!! @sub df_print_runtime
!!
!! print the runtime information, including some physical observables and
!! statistic data, only for reference
!!
  subroutine df_print_runtime()
     implicit none

!! [body

     CONTINUE

!! body]

     return
  end subroutine df_print_runtime

!!
!! @sub df_print_it_info
!!
!! print the iteration information to the screen
!!
  subroutine df_print_it_info()
     implicit none

!! [body

     CONTINUE

!! body]

     return
  end subroutine df_print_it_info
