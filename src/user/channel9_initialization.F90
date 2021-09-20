module channel9_initialization

! This file is part of MOM6. See LICENSE.md for the license.

!use MOM_sponge, only : sponge_CS, initialize_sponge, set_up_sponge_field, Apply_sponge
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

public channel9_initialize_topography

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> This subroutine sets up the channel9 test case topography.
!> channel9 is similar to channel6 but has removed all topography;
!> also different from channel 7: 
!> cn7 still keeps Patagonia and Antarctic Peninsula

subroutine channel9_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m
  
  character(len=40)  :: mod = "channel9_initialize_topography" ! This subroutine's name.
  real :: y, i, j, is, ie, js, je, latext, ll=6.0

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  latext = G%len_lat

  call MOM_mesg("  channel9_initialization.F90, &
    channel9_initialize_topography: setting topography", 5)
  call log_version(param_file, mod, version, "")

  !  Calculate the depth of the bottom
  do j = js,je                                     ! meridional grid points
    do i = is,ie                                   ! zonal grid points
      y = (G%geoLatT(i,j)-G%south_lat)/latext      ! non-dimensional latitude
      D(i,j) = (1.0-spike(y, ll/latext))*max_depth ! Antarctica slope
    enddo
  enddo



end subroutine channel9_initialize_topography

! -----------------------------------------------------------------------------
! define functions used in the above subroutines

!> Returns the value of a sinusoidal bell function 
  real function spike(x,L)

    real, intent(in) :: x
    real, intent(in) :: L
    real             :: PI = 4.0*atan(1.0)

    spike = 1-sin(PI*min(abs(x)/L, 0.5))

  end function spike

! -----------------------------------------------------------------------------

end module channel9_initialization
