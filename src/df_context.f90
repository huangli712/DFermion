!!!-----------------------------------------------------------------------
!!! project : dfermion @ azalea
!!! program : df_mesh module
!!!           df_dmft module
!!!           df_dual module
!!!           df_latt module
!!!           df_susc module
!!!           df_vert module
!!!           context module
!!! source  : df_context.f90
!!! type    : modules
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/16/2009 by li huang (created)
!!!           04/03/2025 by li huang (last modified)
!!! purpose : define the key data structure and global arrays/variables
!!!           for dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module df_mesh                                                   <<<
!!========================================================================

!!
!! @mod df_mesh
!!
!! define some essential k-meshes and frequency meshes
!!
  module df_mesh
     use constants, only : dp

!!
!! @var kx
!!
!! k-mesh along x-axis
!!
     real(dp), public, save, allocatable :: kx(:)

!!
!! @var ky
!!
!! k-mesh along y-axis
!!
     real(dp), public, save, allocatable :: ky(:)

!!
!! @var kz
!!
!! k-mesh along z-axis
!!
     real(dp), public, save, allocatable :: kz(:)

!!
!! @var ek
!!
!! band dispersion
!!
     real(dp), public, save, allocatable :: ek(:)

!!
!! @var fmesh
!!
!! fermionic frequencies
!!
     real(dp), public, save, allocatable :: fmesh(:)

!!
!! @var bmesh
!!
!! bosonic frequencies
!!
     real(dp), public, save, allocatable :: bmesh(:)

  end module df_mesh

!!========================================================================
!!>>> module df_dmft                                                   <<<
!!========================================================================

!!
!! @mod df_dmft
!!
!! define some local variables from the output of quantum impurity solver
!!
  module df_dmft
     use constants, only : dp

     implicit none

!!
!! @var dmft_g
!!
!! local impurity green's function
!!
     complex(dp), public, save, allocatable :: dmft_g(:,:)

!!
!! @var dmft_s
!!
!! local self-energy function
!!
     complex(dp), public, save, allocatable :: dmft_s(:,:)

!!
!! @var dmft_h
!!
!! local hybridization function (old, input)
!!
     complex(dp), public, save, allocatable :: dmft_h(:,:)

!!
!! @var dmft_d
!!
!! local hybridization function (new, output)
!!
     complex(dp), public, save, allocatable :: dmft_d(:,:)

  end module df_dmft

!!========================================================================
!!>>> module df_dual                                                   <<<
!!========================================================================

!!
!! @mod df_dual
!!
!! define dual fermion variables
!!
  module df_dual
     use constants, only : dp

     implicit none

!!
!! @var dual_g
!!
!! dual green's function
!!
     complex(dp), public, save, allocatable :: dual_g(:,:,:)

!!
!! @var dual_s
!!
!! dual self-energy function
!!
     complex(dp), public, save, allocatable :: dual_s(:,:,:)

!!
!! @var dual_b
!!
!! dual bare green's function
!!
     complex(dp), public, save, allocatable :: dual_b(:,:,:)

  end module df_dual

!!========================================================================
!!>>> module df_latt                                                   <<<
!!========================================================================

!!
!! @mod latt
!!
!! define some momentum-dependent (lattice) variables
!!
  module df_latt
     use constants, only : dp

     implicit none

!!
!! @var latt_g
!!
!! lattice green's function
!!
     complex(dp), public, save, allocatable :: latt_g(:,:,:)

!!
!! @var latt_s
!!
!! lattice self-energy function
!!
     complex(dp), public, save, allocatable :: latt_s(:,:,:)

  end module df_latt

!!========================================================================
!!>>> module df_susc                                                   <<<
!!========================================================================

!!
!! @mod df_susc
!!
!! define charge and spin susceptibilities
!!
  module df_susc
     use constants, only : dp

     implicit none

!!
!! @var susc_c
!!
!! charge susceptibility
!!
     complex(dp), public, save, allocatable :: susc_c(:,:,:)

!!
!! @var susc_s
!!
!! spin susceptibility
!!
     complex(dp), public, save, allocatable :: susc_s(:,:,:)

  end module df_susc

!!========================================================================
!!>>> module df_vert                                                   <<<
!!========================================================================

!!
!! @mod df_vert
!!
!! define some vertex functions from the output of quantum impurity solver
!!
  module df_vert
     use constants, only : dp

     implicit none

!!
!! @var vert_d
!!
!! density vertex
!!
     complex(dp), public, save, allocatable :: vert_d(:,:,:)

!!
!! @var vert_m
!!
!! magnetic vertex
!!
     complex(dp), public, save, allocatable :: vert_m(:,:,:)

  end module df_vert

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

!!
!! @mod context
!!
!! containing memory management subroutines, which initialize all of the
!! global variables and arrays
!!
  module context
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : nkp_x, nkp_y, nkp_z, nkpts
     use control, only : nffrq, nbfrq
     use control, only : norbs

     use df_mesh
     use df_dmft
     use df_dual
     use df_latt
     use df_susc
     use df_vert

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

     ! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! declaration of module procedures: allocate memory
     public :: cat_alloc_mesh
     public :: cat_alloc_dmft
     public :: cat_alloc_dual
     public :: cat_alloc_latt
     public :: cat_alloc_susc
     public :: cat_alloc_vert

     ! declaration of module procedures: deallocate memory
     public :: cat_free_mesh
     public :: cat_free_dmft
     public :: cat_free_dual
     public :: cat_free_latt
     public :: cat_free_susc
     public :: cat_free_vert

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_mesh
!!
!! allocate memory for mesh-related variables
!!
  subroutine cat_alloc_mesh()
     implicit none

!! [body

     ! allocate memory
     allocate(kx(nkp_x), stat=istat)
     allocate(ky(nkp_y), stat=istat)
     allocate(kz(nkp_z), stat=istat)
     allocate(ek(nkpts), stat=istat)

     allocate(fmesh(nffrq), stat=istat)
     allocate(bmesh(nbfrq), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     kx = zero
     ky = zero
     kz = zero
     ek = zero

     fmesh = zero
     bmesh = zero

!! body]

     return
  end subroutine cat_alloc_mesh

!!
!! @sub cat_alloc_dmft
!!
!! allocate memory for dmft-related variables
!!
  subroutine cat_alloc_dmft()
     implicit none

!! [body

     ! allocate memory
     allocate(dmft_g(nffrq,norbs), stat=istat)
     allocate(dmft_s(nffrq,norbs), stat=istat)
     allocate(dmft_h(nffrq,norbs), stat=istat)
     allocate(dmft_d(nffrq,norbs), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_dmft','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     dmft_g = czero
     dmft_s = czero
     dmft_h = czero
     dmft_d = czero

!! body]

     return
  end subroutine cat_alloc_dmft

!!
!! @sub cat_alloc_dual
!!
!! allocate memory for dual-related variables
!!
  subroutine cat_alloc_dual()
     implicit none

!! [body

     ! allocate memory
     allocate(dual_g(nffrq,norbs,nkpts), stat=istat)
     allocate(dual_s(nffrq,norbs,nkpts), stat=istat)
     allocate(dual_b(nffrq,norbs,nkpts), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_dual','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     dual_g = czero
     dual_s = czero
     dual_b = czero

!! body]

     return
  end subroutine cat_alloc_dual

!!
!! @sub cat_alloc_latt
!!
!! allocate memory for latt-related variables
!!
  subroutine cat_alloc_latt()
     implicit none

!! [body

     ! allocate memory
     allocate(latt_g(nffrq,norbs,nkpts), stat=istat)
     allocate(latt_s(nffrq,norbs,nkpts), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_latt','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     latt_g = czero
     latt_s = czero

!! body]

     return
  end subroutine cat_alloc_latt

!!
!! @sub cat_alloc_susc
!!
!! allocate memory for susc-related variables
!!
  subroutine cat_alloc_susc()
     implicit none

!! [body

     ! allocate memory
     allocate(susc_c(nbfrq,norbs,nkpts), stat=istat)
     allocate(susc_s(nbfrq,norbs,nkpts), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_susc','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     susc_c = czero
     susc_s = czero

!! body]

     return
  end subroutine cat_alloc_susc

!!
!! @sub cat_alloc_vert
!!
!! allocate memory for vert-related variables
!!
  subroutine cat_alloc_vert()
     implicit none

!! [body

     ! allocate memory
     allocate(vert_d(nffrq,nffrq,nbfrq), stat=istat)
     allocate(vert_m(nffrq,nffrq,nbfrq), stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_vert','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     vert_d = czero
     vert_m = czero

!! body]

     return
  end subroutine cat_alloc_vert

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_mesh
!!
!! deallocate memory for mesh-related variables
!!
  subroutine cat_free_mesh()
     implicit none

!! [body

     if ( allocated(kx)    )  deallocate(kx   )
     if ( allocated(ky)    )  deallocate(ky   )
     if ( allocated(kz)    )  deallocate(kz   )
     if ( allocated(ek)    )  deallocate(ek   )

     if ( allocated(fmesh) )  deallocate(fmesh)
     if ( allocated(bmesh) )  deallocate(bmesh)

!! body]

     return
  end subroutine cat_free_mesh

!!
!! @sub cat_free_dmft
!!
!! deallocate memory for dmft-related variables
!!
  subroutine cat_free_dmft()
     implicit none

!! [body

     if ( allocated(dmft_g) ) deallocate(dmft_g)
     if ( allocated(dmft_s) ) deallocate(dmft_s)
     if ( allocated(dmft_h) ) deallocate(dmft_h)
     if ( allocated(dmft_d) ) deallocate(dmft_d)

!! body]

     return
  end subroutine cat_free_dmft

!!
!! @sub cat_free_dual
!!
!! deallocate memory for dual-related variables
!!
  subroutine cat_free_dual()
     implicit none

!! [body

     if ( allocated(dual_g) ) deallocate(dual_g)
     if ( allocated(dual_s) ) deallocate(dual_s)
     if ( allocated(dual_b) ) deallocate(dual_b)

!! body]

     return
  end subroutine cat_free_dual

!!
!! @sub cat_free_latt
!!
!! deallocate memory for latt-related variables
!!
  subroutine cat_free_latt()
     implicit none

!! [body

     if ( allocated(latt_g) ) deallocate(latt_g)
     if ( allocated(latt_s) ) deallocate(latt_s)

!! body]

     return
  end subroutine cat_free_latt

!!
!! @sub cat_free_susc
!!
!! deallocate memory for susc-related variables
!!
  subroutine cat_free_susc()
     implicit none

!! [body

     if ( allocated(susc_c) ) deallocate(susc_c)
     if ( allocated(susc_s) ) deallocate(susc_s)

!! body]

     return
  end subroutine cat_free_susc

!!
!! @sub cat_free_vert
!!
!! deallocate memory for vert-related variables
!!
  subroutine cat_free_vert()
     implicit none

!! [body

     if ( allocated(vert_d) ) deallocate(vert_d)
     if ( allocated(vert_m) ) deallocate(vert_m)

!! body]

     return
  end subroutine cat_free_vert

  end module context
