!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_setup_param
!!!           dt_setup_model <<<---
!!!           dt_input_mesh_
!!!           dt_input_dmft_
!!!           dt_input_latt_
!!!           dt_input_dual_
!!!           dt_input_vert_ <<<---
!!!           dt_alloc_array
!!!           dt_reset_array
!!!           dt_final_array <<<---
!!! source  : dt_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           01/03/2018 by li huang (last modified)
!!! purpose : initialize and finalize the dual fermion engine.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> config dual fermion engine                                       <<<
!!========================================================================

!!
!! @sub dt_setup_param
!!
!! setup key parameters for dual fermion engine
!!
  subroutine dt_setup_param()
     use parser, only : p_create
     use parser, only : p_parse
     use parser, only : p_get
     use parser, only : p_destroy

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control ! ALL

     implicit none

! local variables
! used to check whether the input file (dt.config.in) exists
     logical :: exists

! setup general control flags
!-------------------------------------------------------------------------
     isdia = 2       ! self-consistent scheme
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! setup common variables for interacting lattice model
!-------------------------------------------------------------------------
     nband = 1       ! number of correlated bands
     nspin = 2       ! number of spin projections
     norbs = 2       ! number of correlated orbitals
!-------------------------------------------------------------------------
     nkpts = 64      ! number of k-points
     nkp_x = 8       ! number of k-points (x_axis)
     nkp_y = 8       ! number of k-points (y_axis)
     nkp_z = 8       ! number of k-points (z_axis)
!-------------------------------------------------------------------------
     mune  = 2.00_dp ! chemical potential or fermi level
     beta  = 1.00_dp ! inversion of temperature
     part  = 1.00_dp ! hopping parameter t for Hubbard model
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! setup common variables for dual fermion engine
!-------------------------------------------------------------------------
     nffrq = 16      ! number of fermionic frequencies
     nbfrq = 7       ! number of bosonic frequncies
!-------------------------------------------------------------------------
     ndfit = 3       ! number of dual fermion iteration
     nbsit = 10      ! number of BSE iteration
!-------------------------------------------------------------------------
     dfmix = 0.70_dp ! mixing parameter (dual fermion iteration)
     bsmix = 0.70_dp ! mixing parameter (BSE solver)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: dt.config.in
         inquire (file = 'dt.config.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
! create the file parser
             call p_create()

! parse the config file
             call p_parse('dt.config.in')

! extract parameters
             call p_get('isdia' , isdia )

             call p_get('nband' , nband )
             call p_get('nspin' , nspin )
             call p_get('norbs' , norbs )

             call p_get('nkpts' , nkpts )
             call p_get('nkp_x' , nkp_x )
             call p_get('nkp_y' , nkp_y )
             call p_get('nkp_z' , nkp_z )

             call p_get('mune'  , mune  )
             call p_get('beta'  , beta  )
             call p_get('part'  , part  )

             call p_get('nffrq' , nffrq )
             call p_get('nbfrq' , nbfrq )

             call p_get('ndfit' , ndfit )
             call p_get('nbsit' , nbsit )

             call p_get('dfmix' , dfmix )
             call p_get('bsmix' , bsmix )

! destroy the parser
             call p_destroy()
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( isdia , master )
     call mp_barrier()

     call mp_bcast( nband , master )
     call mp_bcast( nspin , master )
     call mp_bcast( norbs , master )
     call mp_barrier()

     call mp_bcast( nkpts , master )
     call mp_bcast( nkp_x , master )
     call mp_bcast( nkp_y , master )
     call mp_bcast( nkp_z , master )
     call mp_barrier()

     call mp_bcast( mune  , master )
     call mp_bcast( beta  , master )
     call mp_bcast( part  , master )
     call mp_barrier()

     call mp_bcast( nffrq , master )
     call mp_bcast( nbfrq , master )
     call mp_barrier()

     call mp_bcast( ndfit , master )
     call mp_bcast( nbsit , master )
     call mp_barrier()

     call mp_bcast( dfmix , master )
     call mp_bcast( bsmix , master )
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine dt_setup_param

!!
!! @sub dt_setup_model
!!
!! prepare the model, including frequency mesh, local variables, lattice
!! variables, dual variables, vertex functions, etc
!!
  subroutine dt_setup_model()
     implicit none

! setup frequency mesh
     call dt_input_mesh_()

! setup local variables (from quantum impurity solver)
     call dt_input_dmft_()

! setup lattice variables (from scratch)
     call dt_input_latt_()

! setup dual variables (from scratch)
     call dt_input_dual_()

! setup vertex functions (from quantum impurity solver)
     call dt_input_vert_()

     return
  end subroutine dt_setup_model

!!========================================================================
!!>>> config dual fermion model                                        <<<
!!========================================================================

!!
!! @sub dt_input_mesh_
!!
!! prepare some essential meshes
!!
  subroutine dt_input_mesh_()
     use constants, only : dp
     use constants, only : one, two, pi

     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y, nkp_z
     use control, only : beta
     use control, only : part

     use context, only : kx, ky, kz, ek
     use context, only : fmesh, bmesh

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! setup k-mesh
     do i=1,nkp_x
         kx(i) = (two * pi) * float(i - 1)/ float(nkp_x)
     enddo ! over i={1,nkp_x} loop

     do i=1,nkp_y
         ky(i) = (two * pi) * float(i - 1)/ float(nkp_y)
     enddo ! over i={1,nkp_y} loop

     do i=1,nkp_z
         kz(i) = (two * pi) * float(i - 1)/ float(nkp_z)
     enddo ! over i={1,nkp_z} loop

! build a 2d lattice (band dispersion)
     k = 0
     do i=1,nkp_x
         do j=1,nkp_y
             k = k + 1
             ek(k) = -two * part * ( cos( kx(i) ) + cos( ky(j) ) )
         enddo ! over j={1,nkp_y} loop
     enddo ! over i={1,nkp_x} loop     
     call s_assert(k == nkpts) ! we have to make sure this

! setup fermionic mesh
     do i=1,nffrq
         fmesh(i) = (two * i - one - nffrq) * pi / beta
     enddo ! over i={1,nffrq} loop

! setup bosonic mesh
     do i=1,nbfrq
         bmesh(i) = (two * i - one - nbfrq) * pi / beta
     enddo ! over i={1,nbfrq} loop

     return
  end subroutine dt_input_mesh_

!!
!! @sub dt_input_dmft_
!!
!! prepare some local variables, they are from the output of quantum
!! impurity solver
!!
  subroutine dt_input_dmft_()
     use constants, only : dp
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nffrq
     use control, only : myid, master

     use context, only : dmft_g, dmft_h

     implicit none

! local variables
! loop index
     integer  :: i

! used to check whether the input file (dt.dmft_g.in) exists
     logical  :: exists

! dummy real(dp) variables
     real(dp) :: r1, r2
     real(dp) :: c1, c2

! read in impurity green's function if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'dt.dmft_g.in', exist = exists)

! find input file: dt.dmft_g.in, read it
         if ( exists .eqv. .true. ) then

! read in impurity green's function from dt.dmft_g.in
             open(mytmp, file = 'dt.dmft_g.in', form = 'formatted', status = 'unknown')
             do i=1,nffrq
                 read(mytmp,*) r1, r2, c1, c2
                 dmft_g(i,1) = dcmplx(c1, c2)
                 dmft_g(i,2) = dcmplx(c1, c2)
             enddo ! over i={1,nffrq} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in hybridization function if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'dt.dmft_h.in', exist = exists)

! find input file: dt.dmft_h.in, read it
         if ( exists .eqv. .true. ) then

! read in hybridization function from dt.dmft_h.in
             open(mytmp, file = 'dt.dmft_h.in', form = 'formatted', status = 'unknown')
             do i=1,nffrq
                 read(mytmp,*) r1, r2, c1, c2
                 dmft_h(i,1) = dcmplx(c1, c2)
                 dmft_h(i,2) = dcmplx(c1, c2)
             enddo ! over i={1,nffrq} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! since the data/arrays may be updated in master node, it is important to
! broadcast them from root to all children processes
# if defined (MPI)

! broadcast data
     call mp_bcast(dmft_g, master)
     call mp_bcast(dmft_h, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine dt_input_dmft_

!!
!! @sub dt_input_latt_
!!
!! prepare some lattice variables from scratch
!!
  subroutine dt_input_latt_()
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : ek
     use context, only : dmft_g, dmft_h
     use context, only : latt_g

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! calculate lattice green's function
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 latt_g(i,j,k) = one / ( one / dmft_g(i,j) + dmft_h(i,j) - ek(k) ) 
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine dt_input_latt_

!!
!! @sub dt_input_dual_
!!
!! prepare some dual variables from scratch
!!
  subroutine dt_input_dual_()
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     use context, only : dmft_g
     use context, only : latt_g
     use context, only : dual_g, dual_s, dual_b

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! calculate dual green's functions, dual self-energy functions, and dual
! bath's functions
     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 dual_b(i,j,k) = latt_g(i,j,k) - dmft_g(i,j) 
                 dual_g(i,j,k) = dual_b(i,j,k)
                 dual_s(i,j,k) = czero
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine dt_input_dual_

!!
!! @sub dt_input_vert_
!!
!! prepare vertex functions, they are from the output of quantum impurity
!! solver as well
!!
  subroutine dt_input_vert_()
     use constants, only : dp
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nffrq, nbfrq
     use control, only : myid, master

     use context, only : vert_d, vert_m

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: if1, if2

! used to check whether the input file (dt.vert_d.in) exists
     logical  :: exists

! dummy real(dp) variables
     real(dp) :: r1, r2
     real(dp) :: c1, c2
     real(dp) :: d1, d2
     real(dp) :: v1, v2

! read in vertex function (density channel) if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'dt.vert_d.in', exist = exists)

! find input file: dt.vert_d.in, read it
         if ( exists .eqv. .true. ) then

! read in vertex function (density channel) from dt.vert_d.in
             open(mytmp, file = 'dt.vert_d.in', form = 'formatted', status = 'unknown')
             do i=1,nbfrq
                 do if1=1,nffrq
                     do if2=1,nffrq
                         read(mytmp,*) r1, r2, c1, c2, d1, d2, v1, v2
                         vert_d(if2,if1,i) = dcmplx(v1, v2)
                     enddo ! over if2={1,nffrq} loop
                     read(mytmp,*) ! skip one line
                 enddo ! over if1={1,nffrq} loop
             enddo ! over i={1,nbfrq} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in vertex function (magentic channel) if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'dt.vert_m.in', exist = exists)

! find input file: dt.vert_m.in, read it
         if ( exists .eqv. .true. ) then

! read in vertex function (magnetic channel) from dt.vert_m.in
             open(mytmp, file = 'dt.vert_m.in', form = 'formatted', status = 'unknown')
             do i=1,nbfrq
                 do if1=1,nffrq
                     do if2=1,nffrq
                         read(mytmp,*) r1, r2, c1, c2, d1, d2, v1, v2
                         vert_m(if2,if1,i) = dcmplx(v1, v2)
                     enddo ! over if2={1,nffrq} loop
                     read(mytmp,*) ! skip one line
                 enddo ! over if1={1,nffrq} loop
             enddo ! over i={1,nbfrq} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! since the data/arrays may be updated in master node, it is important to
! broadcast them from root to all children processes
# if defined (MPI)

! broadcast data
     call mp_bcast(vert_d, master)
     call mp_bcast(vert_m, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine dt_input_vert_

!!========================================================================
!!>>> manage memory for dual fermion engine                            <<<
!!========================================================================

!!
!! @sub dt_alloc_array
!!
!! allocate memory for global variables and then initialize them
!!
  subroutine dt_alloc_array()
     use context ! ALL

     implicit none

! allocate memory for context module
     call cat_alloc_mesh()
     call cat_alloc_dmft()
     call cat_alloc_dual()
     call cat_alloc_latt()
     call cat_alloc_vert()
     
     return
  end subroutine dt_alloc_array

!!
!! @sub dt_reset_array
!!
!! reset the key variables for dual fermion engine 
!!
  subroutine dt_reset_array()
     implicit none

     return
  end subroutine dt_reset_array

!!
!! @sub dt_final_array
!!
!! garbage collection for this code, please refer to dt_alloc_array
!!
  subroutine dt_final_array()
     use context ! ALL

     implicit none

! deallocate memory for context module
     call cat_free_mesh()
     call cat_free_dmft()
     call cat_free_dual()
     call cat_free_latt()
     call cat_free_vert()

     return
  end subroutine dt_final_array
