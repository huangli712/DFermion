!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : control module
!!!           version module
!!! source  : df_control.f90
!!! type    : module
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/15/2009 by li huang (created)
!!!           04/03/2025 by li huang (last modified)
!!! purpose : define global control parameters for dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module control                                                   <<<
!!========================================================================

!!
!! @mod control
!!
!! define the control parameters and dimensional parameters.
!!
  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! code name of the dual fermion framework
!!
     character(len = 06), public, save :: cname = 'AZALEA'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

!!
!! @var isdia
!!
!! control flag, define the running scheme of the code
!!
!! if isdia == 1:
!!     standard dual fermion approximation scheme, only the second order
!!     diagrams are included.
!!
!! if isdia == 2:
!!     ladder dual fermion approximation scheme, only the ladder diagrams
!!     are included.
!!
!! if isdia == 3:
!!     reserved
!!
     integer, public, save :: isdia  = 2

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nband
!!
!! number of correlated bands
!!
     integer, public, save :: nband  = 1

!!
!! @var nspin
!!
!! number of spin projections
!!
     integer, public, save :: nspin  = 2

!!
!! @var norbs
!!
!! number of correlated orbitals (= nband * nspin)
!!
     integer, public, save :: norbs  = 2

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nffrq
!!
!! number of fermionic frequencies for the two-particle green's function
!!
     integer, public, save :: nffrq  = 16

!!
!! @var nbfrq
!!
!! number of bosonic frequncies for the two-particle green's function
!!
     integer, public, save :: nbfrq  = 7

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nkpts
!!
!! number of k-points (totally)
!!
     integer, public, save :: nkpts  = 64

!!
!! @var nkp_x
!!
!! number of k-points (along x-axis)
!!
     integer, public, save :: nkp_x  = 8

!!
!! @var nkp_y
!!
!! number of k-points (along y-axis)
!!
     integer, public, save :: nkp_y  = 8

!!
!! @var nkp_z
!!
!! number of k-points (along z-axis)
!!
     integer, public, save :: nkp_z  = 8

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var ndfit
!!
!! number of dual fermion iterations
!!
     integer, public, save :: ndfit  = 10

!!
!! @var nbsit
!!
!! number of iterations for solving the Bethe-Salpter equation
!!
     integer, public, save :: nbsit  = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

!!
!! @var mune
!!
!! chemical potential or fermi level
!!
     real(dp), public, save :: mune  = 0.00_dp

!!
!! @var beta
!!
!! inversion of temperature
!!
     real(dp), public, save :: beta  = 1.00_dp

!!
!! @var part
!!
!! hopping parameter t for lattice model
!!
     real(dp), public, save :: part  = 1.00_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var dfmix
!!
!! mixing parameter for dual fermion iteration. it is used to mix the
!! dual green's function or dual self-energy function. 
!!
     real(dp), public, save :: dfmix = 1.00_dp

!!
!! @var bsmix
!!
!! mixing parameter for solving the Bethe-Salpter equation via iteration
!! approach. that is to say it is only useful in cat_bse_iterator().
!!
     real(dp), public, save :: bsmix = 0.70_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

!!
!! @var nprocs
!!
!! number of processors: default value 1
!!
     integer, public, save :: nprocs = 1

!!
!! @var myid
!!
!! the id of current process: default value 0
!!
     integer, public, save :: myid   = 0

!!
!! @var master
!!
!! denote as the controller process: default value 0
!!
     integer, public, save :: master = 0

!!
!! @var cid
!!
!! the id of current process in cartesian topology (cid == myid)
!!
     integer, public, save :: cid    = 0

!!
!! @var cx
!!
!! the x coordinates of current process in cartesian topology
!!
     integer, public, save :: cx     = 0

!!
!! @var cy
!!
!! the y coordinates of current process in cartesian topology
!!
     integer, public, save :: cy     = 0

  end module control

!!========================================================================
!!>>> module version                                                   <<<
!!========================================================================

!!
!! @mod version
!!
!! define the semantic version string.
!!
  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL = 'v0.5.0 @ 2025.04.03D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.5.0'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2025.04.03'

!!
!! @var V_STAT
!!
!! version string, only status info., D means devel, T testing, R released.
!!
     character(len=01), public, parameter :: V_STAT = 'D'

!!
!! @var V_AUTH
!!
!! version string, author info.
!!
     character(len=11), public, parameter :: V_AUTH = 'by li huang'

!!
!! @var V_INST
!!
!! version string, affiliation info.
!!
     character(len=36), public, parameter :: V_INST = 'China Academy of Engineering Physics'

!!
!! @var V_MAIL
!!
!! version string, email info.
!!
     character(len=22), public, parameter :: V_MAIL = 'huangli@caep.cn'

!!
!! @var V_GPL3
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: V_GPL3 = 'GNU General Public License version 3'

  end module version
