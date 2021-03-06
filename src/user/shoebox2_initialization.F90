module shoebox2_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type, read_param
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public shoebox2_initialize_topography
public shoebox2_initialize_thickness

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the shoebox2 test case topography.
!> Differences from shoebox: added continental slopes along all 4 side walls

subroutine shoebox2_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the shoebox2 test case topography
  real :: PI = 4.0*atan(1.0)   ! 3.1415926... calculated as 4*atan(1)
  real :: latext, lonext       ! latitude extent of the model area
  real :: ep = epsilon(1.)     ! an infinitesimally small quantity
  real :: x, y, sa_dim = 1500.0! dimensional height of Scotia Arc
  real :: sa                   ! the non-dimensional height of Scotia Arc top;
                               ! default value is 1500/4000=0.375
  real :: dx                   ! non-dimensional longitudinal grid scale
  real :: reentrants, reentrantn  ! the non-dimensional latitudes of the southern and northern
                                  ! boundary of the reentrant channel
  real :: sdp = 6.0, ll=6.0, ssp = 2.0    ! the width of the slope in Drake Passage, the half width of continental slope, and sponge width
  real :: min_depth
 
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "shoebox2_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

    call read_param(param_file, "MINIMUM_DEPTH", min_depth)       ! read in min_depth value

 
  sa = sa_dim / max_depth
  latext = G%len_lat
  lonext = G%len_lon
  reentrants = 6.0/latext          ! non-dimensional southern bound of the reentrant zone
  reentrantn = 10.0/latext         ! non-dimensional northern bound of the reentrant zone
  sdp = sdp/latext
  ssp = ssp/latext
  D = 0.0
  dx = (G%geoLonT(is+1,js) - G%geoLonT(is, js)) / lonext

  call MOM_mesg("  shoebox2_initialization.F90, shoebox2_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")

  !  Calculate the depth of the bottom.
  do j=js,je                ! meridional grid points
  do i=is,ie                ! zonal grid points
    x=(G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    y=(G%geoLatT(i,j)-G%south_lat) / latext     ! non-dimensional latitude


    D(i,j) = 1.0 - spike(x-dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &                 ! Patagonia, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &               ! Patagonia, east, original
                -sa * spike(x-1.0+dx/2, ll/lonext)*cosbell(y-12.0/latext, 2.5/latext) &                  ! Patagonia, east, extra slope
                -spike(x-dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &                   ! Antarctic Peninsula, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &               ! Antarctic Peninsula, east, original
                -sa * spike(x-1.0+dx/2, ll/lonext)*cosbell(y-4.0/latext, 2.5/latext) &                   ! Antarctic Peninsula, east, extra slope
                -spike(y, ll/latext) &                                                                   ! Antarctica
                -spike(1.0-y, ll/latext) &                                                               ! Northern Wall
                - sa *cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * homo(y-8.0/latext, 2.0/latext) &              !Scotia Arc East, center
                - sa * cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * cosbellh(y-10.0/latext-ep, 3.0/latext, 1.) & !Scotia Arc East, north slope
                - sa * cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * cosbellh(y-6.0/latext+ep, 3.0/latext, -1.) & !Scotia Arc East, south slope
                - sa * cosbell(y-12.0/latext, 2.5/latext) * cosbellh(x-dx/2-18.0/lonext-ep, 2.5/lonext, 1.) & !Scotia Arc North, east half (slope side)
                - sa * cosbell(y-12.0/latext, 2.5/latext) * homo(x-dx/2-9.0/lonext, 9.0/lonext) &             !Scotia Arc North, west half
                - sa * cosbell(y-4.0/latext, 2.5/latext) * cosbellh(x-dx/2-18.0/lonext-ep, 2.5/lonext, 1.) &  !Scotia Arc South, east half (slope side)
                - sa * cosbell(y-4.0/latext, 2.5/latext) * homo(x-dx/2-9.0/lonext, 9.0/lonext)                !Scotia Arc South, west half

    ! make sure no deeper than max depth and no shallower than Scotia Arc top IN the ocean interior
    if (D(i,j)<1.0 - sa .and. x>=dx/1.5+ll/2/lonext .and. x<=1.0-ll/2/lonext-dx/1.5 &
        .and. y>=ll/2/latext .and. y<=1.0-ll/2/latext) then
      D(i,j) = 1 - sa
    else if (D(i,j) > 1.0) then
      D(i,j) = 1.0
    endif

      ! make sure the model is not zonally reentrant outside of Drake Passage
      if (((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x<dx/1.5) &
                .or. ((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x>1.0-dx/1.5)) then
        D(i,j) = min_depth/max_depth
      endif

    D(i,j) = D(i,j) * max_depth
  enddo
  enddo


end subroutine shoebox2_initialize_topography


! -----------------------------------------------------------------------------
!> This subroutine initializes layer thicknesses for the shoebox test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine shoebox2_initialize_thickness(h, G, GV, param_file, eqn_of_state, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being
                                                              !! initialized.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse
                                                              !for model
                                                              !! parameter
                                                              !values.
  type(EOS_type),          pointer    :: eqn_of_state         !< integer that selects the
                                                              !! equation of
                                                              !state.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference
                                                              !pressure in Pa.
  ! Local variables
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                            ! negative because it is positive upward.      !
  real, dimension(SZK_(G)) :: h_profile ! Vector of initial thickness profile (m)
  real :: e_interface ! Current interface positoin (m)
  character(len=40)  :: mod = "shoebox2_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  shoebox2_initialization.F90, shoebox2_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mod, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", fail_if_missing=.true.)

! e0 is the notional position of interfaces
  e0(1) = 0. ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie
    e_interface = -G%bathyT(i,j)
    do k=nz,1,-1
      h(i,j,k) = max( GV%Angstrom_z, e0(k) - e_interface )
      e_interface = max( e0(k), e_interface - h(i,j,k) )
    enddo

  enddo ; enddo

end subroutine shoebox2_initialize_thickness


! -----------------------------------------------------------------------------
! define functions used in the above subroutines

!> Returns the value of a sinusoidal bell function 
  real function spike(x,L)

    real, intent(in) :: x
    real, intent(in) :: L
    real             :: PI = 4.0*atan(1.0)

    spike = 1-sin(PI*min(abs(x)/L, 0.5))

  end function spike

!> Returns the value of a cosine-bell function evaluated at x/L
 real function cosbell(x,L)

   real , intent(in) :: x                         !< non-dimensional position
   real , intent(in) :: L                         !< non-dimensional width
   real              :: PI = 4.0 * atan(1.0)      !< 3.1415926... calculated as 4*atan(1)

   cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
 end function cosbell

!< Return the value of a half-cosine-bell function evaluated at x/L;
!< i.e. from peak to trough only on one side of the bell
 real function cosbellh(x, L, dir)

  real, intent(in) :: x       !< non-dimensional position
  real, intent(in) :: L       !< non-dimensional width
  real             :: PI, xx  !< 3.1415926... calculated as 4*atan(1)
  real, intent(in) :: dir     !< direction flag; 1 for east/north; -1 for west/south
  PI        = 4.0*atan(1.0)

    !< if the grid falls on the opposite side of the bell, override x to be so big that x/L > 1
    if (x*dir .lt. 0.0) then
      xx = L+1
    else
      xx = x
    endif

    cosbellh  = cos(PI/2.0*MIN(abs(xx)/L, 1.0))
  end function cosbellh


 !< make sure the depth within L is homogeneous
   real function homo(x, L)

     real, intent(in) :: x       !< non-dimensional position
     real, intent(in) :: L       !< non-dimensional width

     !< if x falls within -L ~ L, assign 1 to the non-dimensional depth
    if (abs(x) .le. L) then
      homo = 1.0
    else
      homo = 0.0
    endif
   end function homo

end module shoebox2_initialization
