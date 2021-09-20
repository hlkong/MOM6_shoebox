module channel44_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, initialize_sponge, set_up_sponge_field, Apply_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public channel44_initialize_sponges

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> Sets up the the inverse restoration time (Idamp), and
! the values towards which the interface heights and an arbitrary
! number of tracers should be restored within each sponge.
!
! differences from channel6's sponge:
! the target density profile is based on sb2G4tp2, which uses sb2
! topography (Scotia Arc + continental slopes on all 4 walls) and
! shoebox9 sfc boundary conditions: symmetric buoyancy & asymmetric wind
! also, an artificially deepened surface mixed layer = 100m is adopted

subroutine channel44_initialize_sponges(G, GV, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
                                             ! so as to be used as input of set_coord_from_file
  logical, intent(in) :: use_temperature    !< Switch for temperature.
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure containing pointers
                                            !! to any available thermodynamic
                                            !! fields, potential temperature and
                                            !! salinity or mixed layer density.
                                            !! Absent fields have NULL ptrs.
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                            !! open file to parse for model
                                            !! parameter values.
  type(sponge_CS),   pointer    :: CSp      !< A pointer that is set to point to
                                            !! the control structure for the
                                            !! sponge module.
  real, intent(in), dimension(SZI_(G),SZJ_(G), SZK_(G)) :: h !< Thickness field. - may not be used!
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta, m.
                                         ! eta only varies in z
  real :: Idamp(SZI_(G),SZJ_(G))         ! The inverse damping rate, in s-1.
                                         ! Idamp only varies in y; /= 0 only within the sponge layer
  real, dimension(SZK_(G)+1) :: eta0    ! target interface heights for density class in the sponge layer
  real :: damp_rate, damp, spongelen = 1.0,  min_depth, nlat, dx
  ! spongelen: thickness of sponge layer in dimensional degree
  ! dx is non-dimensional longitudinal grid increment
  integer :: i, j, k, is, ie, js, je, nz
  logical, save :: first_call = .true.
  character(len=40)  :: mdl = "channel44_initialize_sponges" ! This subroutine's name

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; eta0(:) = 0.0; 
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/G%len_lon

  ! target interface heights: all negative values
  ! corresponds to rho8a (diagnosed from sb8sG, zonal mean along 30S)
  
  eta0 = (/0.0,  0.0, 0.0, -1.3, -48.4, -92.2, -142.0, -210.6, -297.2, -413.7, &
        -544.9, -663.0, -775.2, -895.3, -1042.1, -1254.0, -1602.0, -1969.0, &
        -2287.5, -2566.3, -2810.3, -3027.1, -3225.0, -3413.9, -3606.3, -3774.4, &
        -3864.9, -3893.9, -3906.7, -3912.6, -3912.6/)

  call get_param(param_file, mdl, "SPONGE_RATE", damp_rate, &
                 "The rate at which the zonal-mean sponges damp.", units="s-1", &
                 default = 1.0/(10.0*86400.0))
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)


  nlat = G%south_lat + G%len_lat        ! should be -30.0 degree
  
  ! initialize the damping rate so it is 0 outside of the sponge layer &
  ! and increases linearly with latitude within the sponge layer
  do i = is, ie; do j = js, je
    if (G%geoLatT(i,j) >= nlat-spongelen) then
      damp = damp_rate/spongelen * (G%geoLatT(i,j)-nlat+spongelen)
    else
      damp = 0.0   ! outside of the sponge
    endif
    
    do k = 1,nz+1; eta(i,j,k) = eta0(k); enddo    ! initialize target heights for each (lat,lon) grid

    if (G%bathyT(i,j) > min_depth) then ! bathT: Ocean bottom depth (positive) at tracer points, in m.
      Idamp(i,j) = damp                 ! no need to divide this by 86400!!!
    else
      Idamp(i,j) = 0.0
    endif     ! so that at the side walls, Idamp = 0
  enddo; enddo

  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!From MOM_sponge.F90
!subroutine initialize_sponge(Idamp, eta, G, param_file, CSp, &
!                             Iresttime_i_mean, int_height_i_mean)
!
! Arguments: Idamp - The inverse of the restoring time, in s-1.
!  (in)      eta - The interface heights to damp back toward, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                 for this module


end subroutine channel44_initialize_sponges


end module channel44_initialization
