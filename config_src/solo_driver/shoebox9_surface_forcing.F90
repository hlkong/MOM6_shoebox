!> Wind and buoyancy forcing for the shoebox9 configurations
!> difference from shoebox8: 
!> 1) corrected rho(65N): from 1037.5 to 1037.2
!> 2) re-designed wind profile so it's symmetric about the equator
!> 3) introduced WINDPEAK_SO to change the peak wind stress over the Southern
!>    Ocean more flexibly

module shoebox9_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_variables, only : surface

implicit none ; private

public shoebox9_buoyancy_forcing
public shoebox9_wind_forcing
public shoebox9_surface_forcing_init

!> This control structure should be used to store any run-time variables
!! associated with the shoebox9 forcing.  It can be readily modified
!! for a specific case, and because it is private there will be no changes
!! needed in other code (although they will have to be recompiled).
type, public :: shoebox9_surface_forcing_CS ; private

  logical :: use_temperature !< If true, use temperature and salinity.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq
                             !! approximation, in kg m-3.
  real :: G_Earth            !< The gravitational acceleration in m s-2.
  real :: flux_const         !<  The restoring rate at the surface, in m s-1.
  real :: windpeak_SO        !< wind stress peak in the Southern Ocean
  real, dimension(:,:), pointer :: &
    buoy_restore(:,:) => NULL() !< The pattern to restore buoyancy to.
  character(len=200) :: inputdir !< The directory where NetCDF input files are.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  logical :: first_call = .true. !< True until shoebox9_buoyancy_forcing has been called
end type shoebox9_surface_forcing_CS

contains


subroutine shoebox9_buoyancy_forcing(sfc_state, fluxes, day, dt, G, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< Forcing fields.
  type(time_type),               intent(in)    :: day !< Time used for determining the fluxes.
  real,                          intent(in)    :: dt !< Forcing time step (s).
  type(ocean_grid_type),         intent(inout) :: G !< Grid structure.
  type(shoebox9_surface_forcing_CS), pointer  :: CS !< Control structure for this module.
  ! Local variables
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux, in m5 s-3 kg-1.
  real :: y, latext, sslat, eqlat ! dimensional: 65S; non-dimensional: 30S, Equator, 30N
  real :: den, rho1 = 1037.5, rho2 = 1031.3, rho3 = 1037.2 ! south, equator, north
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  latext = G%len_lat
  sslat = G%south_lat
  eqlat = -sslat/latext

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call MOM_error(FATAL, "shoebox9_buoyancy_forcing: " // &
        "Temperature and salinity mode not coded!" )
  else
    ! This is the buoyancy only mode.
    call alloc_if_needed(fluxes%buoy, isd, ied, jsd, jed)
  endif


! set up target surface density profile
  if (CS%restorebuoy .and. CS%first_call) then
    call alloc_if_needed(CS%buoy_restore, isd, ied, jsd, jed)
    CS%first_call = .false.

    do j = js, je
      y = (G%geoLatT(is,j)-G%south_lat)/latext                          ! non-dimensional latitude
      if (y <= eqlat) then
        CS%buoy_restore(is,j) = (rho2-rho1)/eqlat*y+rho1                ! southern hemisphere
      else
        CS%buoy_restore(is,j) = (rho3-rho2)/(1.0-eqlat)*(y-eqlat)+rho2        ! northern hemisphere
      endif
    enddo

    ! set up the profile for all the other longitudes
    do i = is+1, ie
      CS%buoy_restore(i, :)=CS%buoy_restore(is, :)
    enddo

  endif


  ! compute the buoyancy flux needed to restore current density to target density
  if (CS%restorebuoy) then
    if (CS%use_temperature) then      ! if temperature is used to restore buoyancy
      call MOM_error(FATAL, "shoebox9_buoyancy_surface_forcing: " // &
        "Temperature/salinity restoring not coded!" )
    else      ! density is used to restore buoyancy

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%flux_const) / CS%Rho0
      do j=js,je
        do i=is,ie
        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          (CS%buoy_restore(i, j) - sfc_state%sfc_density(i,j))
        enddo
      enddo

    endif
  endif

end subroutine shoebox9_buoyancy_forcing

!> Sets the surface wind stresses, forces%taux and forces%tauy for the
!! shoebox9 forcing configuration.
subroutine shoebox9_wind_forcing(sfc_state, forces, day, G, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface
                                                    !state of the ocean.
  type(mech_forcing),            intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),               intent(in)    :: day !< Time used for determining the fluxes.
  type(ocean_grid_type),         intent(inout) :: G !< Grid structure.
  type(shoebox9_surface_forcing_CS), pointer  :: CS !< Control structure for this module.
  ! Local variable
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  real :: x, y
  real :: PI = 4.0*atan(1.0), wp2, wp3, wp4, wp5, wp6, latext, lonext

  wp2 = -0.08  ! easterly peak in Southern subtropics; 11.75S
  wp3 = 0.1  ! easterly trough near the equator; 4.75N
  forces%taux = 0.0
  latext = G%len_lat
  lonext = G%len_lon

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec


  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true.)

  !  Set the surface wind stresses, in units of Pa.  A positive taux
  !  accelerates the ocean to the (pseudo-)east.
  do j = js, je
      y=(G%geoLatT(is, j)-G%south_lat)/latext                                 ! non-dimensional latitudes
      forces%taux(is, j) = CS%windpeak_SO*cosbellh(y-12.0/latext, 12.0/latext, -1.0) & ! south of the ACC peak
          + wp2*cosbell(y-51.0/latext, 19.0/latext)                         & ! profile in Southern subtropics
          + wp2*cosbell(y-79.0/latext, 19.0/latext)                         & ! profile in Northern subtropics
          + wp3*cosbellh(y-118.0/latext, 12.0/latext, 1.0)                  ! profile north of 53N

      if (y > 12.0/latext) then              ! north of the ACC peak
        forces%taux(is, j) = forces%taux(is, j)+ CS%windpeak_SO*cosbell(y-12.0/latext, 25.0/latext)
      endif
      if (y < 118.0/latext) then             ! south of northern hemisphere extratropics peak wind latitude
        forces%taux(is, j) = forces%taux(is, j)+ wp3*cosbell(y-118.0/latext, 25.0/latext)
      endif

      do i = is+1, ie
        forces%taux(i, j) = forces%taux(is, j)                ! same wind stress for other longitudes
      enddo
  enddo

end subroutine shoebox9_wind_forcing


!> If ptr is not associated, this routine allocates it with the given size
!! and zeros out its contents.  This is equivalent to safe_alloc_ptr in
!! MOM_diag_mediator, but is here so as to be completely transparent.
subroutine alloc_if_needed(ptr, isd, ied, jsd, jed)
  real, pointer :: ptr(:,:)
  integer :: isd, ied, jsd, jed
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(isd:ied,jsd:jed))
    ptr(:,:) = 0.0
  endif
end subroutine alloc_if_needed

!> Initializes the shoebox9 control structure.
subroutine shoebox9_surface_forcing_init(Time, G, param_file, diag, CS)
  type(time_type),         intent(in) :: Time       !< The current model time.
  type(ocean_grid_type),   intent(in) :: G          !< The ocean's grid structure.
  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file to parse for
                                                    !! model parameter values.
  type(diag_ctrl), target, intent(in) :: diag       !< A structure that is used to regulate diagnostic output.
  type(shoebox9_surface_forcing_CS), pointer :: CS !< A pointer that is set to point to the control structure
                                                    !! for this module
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  ! Local variables
  character(len=40)  :: mod = "shoebox9_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "shoebox9_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)

  call get_param(param_file, mod, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
! call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
!                "The background gustiness in the winds.", units="Pa", &
!                default=0.02)

  call get_param(param_file, mod, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back \n"//&
                 "toward some specified surface state with a rate \n"//&
                 "given by FLUXCONST.", default= .false.)
  call get_param(param_file, mod, "WINDPEAK_SO", CS%windpeak_SO, &
                 "Wind stress peak in the Southern Ocean. ", &
                 units="N m-2", default= 0.2)

  if (CS%restorebuoy) then
    call get_param(param_file, mod, "FLUXCONST", CS%flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    ! Convert CS%flux_const from m day-1 to m s-1.
    CS%flux_const = CS%flux_const / 86400.0
  endif

end subroutine shoebox9_surface_forcing_init

!------------------------------ defines the functions to used from above ----------------
  !> Returns the value of a cosine-bell function evaluated at x/L
   real function cosbell(x,L)

     real , intent(in) :: x       !< non-dimensional position
     real , intent(in) :: L       !< non-dimensional width
     real              :: PI      !< 3.1415926... calculated as 4*atan(1)

     PI      = 4.0*atan(1.0)
     cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
   end function cosbell

   !< return the value of a half-cosine-bell function evaluated at x/L;
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

     !< similar to cosbellh but takes a different shape of bell
       real function cosbellhnew(x, L, dir)

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

         cosbellhnew  = 0.5*(1+cos(PI*MIN(xx/L, 1.0)))
       end function cosbellhnew

end module shoebox9_surface_forcing
