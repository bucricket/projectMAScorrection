! Description:
!> @file
!!   Defines all derived types for RTTOV
!
!> @brief
!!   Defines all derived types for RTTOV
!!
!! @details
!!   This contains types that users will make use of in their code
!!   as well as types that RTTOV only uses internally.
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_types

  USE rttov_const, ONLY : &
      fastem_sp,          &
      ncldtyp,            &
      interp_rochon,      &
      gas_unit_compatibility

  USE parkind1, ONLY : jpim, jprb, jplm

  IMPLICIT NONE

  ! ---------------------------------------------------------------------------
  ! User-level RTTOV structures
  ! ---------------------------------------------------------------------------

  !> Specify channels and profiles to simulate, declare as array of size
  !! nchanprof which is the total number of channels to simulate
  TYPE rttov_chanprof
    INTEGER(KIND=jpim) :: chan              !< Channel index
    INTEGER(KIND=jpim) :: prof              !< Profile index
  END TYPE

  !> Input/output surface emissivities, declare as array of size nchanprof
  TYPE rttov_emissivity
    REAL(KIND=jprb)    :: emis_in           !< Input emissivity
    REAL(KIND=jprb)    :: emis_out          !< Output emissivity, value used by RTTOV
  END TYPE

  !> Input/output surface BRDFs, declare as array of size nchanprof,
  !! can also be used to specify cloud top BRDF for simple cloud scheme in VIS/NIR channels
  TYPE rttov_reflectance
    REAL(KIND=jprb)    :: refl_in           !< Input BRDF
    REAL(KIND=jprb)    :: refl_out          !< Output BRDF, value used by RTTOV
    REAL(KIND=jprb)    :: refl_cloud_top    !< Optional, cloud top BRDF for simple cloud
  END TYPE

  !> Surface skin variables
  TYPE sskin_type
    INTEGER(KIND=jpim) :: surftype          !< Surface type: 0=land, 1=sea, 2=sea-ice
    INTEGER(KIND=jpim) :: watertype         !< Water type: 0=fresh water, 1=ocean water
    REAL(KIND=jprb)    :: t                 !< Radiative skin temperature (K)
    REAL(KIND=jprb)    :: salinity          !< Practical ocean salinity unit (\%o) - FASTEM-4/5 only
    REAL(KIND=jprb)    :: foam_fraction     !< ocean foam fraction (0-1; only used if supply_foam_fraction is true)
    REAL(KIND=jprb)    :: fastem(fastem_sp) !< FASTEM land/sea-ice surface parameters
  END TYPE sskin_type

  !> Surface 2m variables
  TYPE s2m_type
    REAL(KIND=jprb) :: t                    !< Temperature (K)
    REAL(KIND=jprb) :: q                    !< Water vapour (ppmv or kg/kg)
    REAL(KIND=jprb) :: o                    !< Ozone (ppmv or kg/kg) - not currently used
    REAL(KIND=jprb) :: p                    !< Surface pressure (hPa)
    REAL(KIND=jprb) :: u                    !< U 10m wind component (m/s)
    REAL(KIND=jprb) :: v                    !< V 10m wind component (m/s)
    REAL(KIND=jprb) :: wfetc                !< Wind fetch (metres)
  END TYPE s2m_type

  !> Atmospheric profiles on model pressure levels
  TYPE profile_type
    CHARACTER(LEN=128) :: id                !< Optional profile ID string
    INTEGER(KIND=jpim) :: date(3)           !< Year, Month, Day
    INTEGER(KIND=jpim) :: time(3)           !< Hour, Minute, Second

    INTEGER(KIND=jpim) :: nlevels           !< Number of atmospheric levels
    INTEGER(KIND=jpim) :: nlayers           !< Number of atmospheric layers

    INTEGER(KIND=jpim) :: gas_units         !< Units of gas profiles (0=>compatibility mode, 1=>kg/kg, 2=>ppmv)

    REAL(KIND=jprb), POINTER :: p(:)        !< Pressure (hPa)
    REAL(KIND=jprb), POINTER :: t(:)        !< Temperature (K)
    REAL(KIND=jprb), POINTER :: q(:)        !< Water vapour (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: o3(:)       !< o3 (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: co2(:)      !< co2 (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: n2o(:)      !< n2o (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: co(:)       !< co (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: ch4(:)      !< ch4 (ppmv or kg/kg)
    REAL(KIND=jprb), POINTER :: clw(:)      !< Cloud liquid water (kg/kg)

    REAL(KIND=jprb), POINTER :: aerosols(:,:) !< Aerosol layer number densities (cm-3)
    REAL(KIND=jprb), POINTER :: cloud(:,:)    !< Cloud water ice layer densities (g/m3)
    REAL(KIND=jprb), POINTER :: cfrac(:)      !< Layer cloud fraction (0-1)
    REAL(KIND=jprb), POINTER :: icede(:)      !< Ice particle effective diameter (microns)
    INTEGER(KIND=jpim)       :: idg           !< Ice particle effective diameter parameterisation (1-4)
    INTEGER(KIND=jpim)       :: ish           !< Select ice particle shape or Baran parameterisation (1-4)

    TYPE(sskin_type) :: skin                !< Surface skin variables
    TYPE(s2m_type)   :: s2m                 !< Surface 2m variables

    REAL(KIND=jprb) :: zenangle             !< Satellite zenith angle (degrees)
    REAL(KIND=jprb) :: azangle              !< Satellite azimuth angle (degrees)
    REAL(KIND=jprb) :: sunzenangle          !< Solar azimuth angle (degrees)
    REAL(KIND=jprb) :: sunazangle           !< Solar azimuth angle (degrees)
    REAL(KIND=jprb) :: elevation            !< Surface elevation (km)
    REAL(KIND=jprb) :: latitude             !< Latitude (degrees)
    REAL(KIND=jprb) :: longitude            !< Longitude (degrees)
    REAL(KIND=jprb) :: snow_frac            !< Snow coverage fraction for IR emissivity atlas (0-1)
    REAL(KIND=jprb) :: soil_moisture        !< Soil moisture (m^3/m^3) - not currently used

    REAL(KIND=jprb) :: Be                   !< Earth magnetic field strength (Gauss)
    REAL(KIND=jprb) :: cosbk                !< Cosine of the angle between the Earth magnetic
                                            !! field and wave propagation direction

    REAL(KIND=jprb) :: ctp                  !< Black body (simple) cloud top pressure (hPa)
    REAL(KIND=jprb) :: cfraction            !< Black body (simple) cloud fraction (0-1)
  END TYPE profile_type

  !> Additional atmospheric cloud/hydrometeor profile input for RTTOV-SCATT
  TYPE profile_cloud_type
    INTEGER(KIND=jpim) :: nlevels      !< Number of atmospheric levels (same as in profile_type)
    LOGICAL(KIND=jplm) :: use_totalice !< False => separate ice and snow; True => total ice
    LOGICAL(KIND=jplm) :: mmr_snowrain !< Snow and rain input units are: False => kg/m2/s; True => kg/kg
    REAL(KIND=jprb)    :: cfrac        !< Average cloud fraction (only used if lusercfrac = TRUE)

    REAL(KIND=jprb), POINTER :: ph(:)        !< nlevels+1 of half-level model pressures (hPa)
    REAL(KIND=jprb), POINTER :: cc(:)        !< nlevels of cloud cover
    REAL(KIND=jprb), POINTER :: clw(:)       !< nlevels of cloud liquid water (kg/kg)
    REAL(KIND=jprb), POINTER :: ciw(:)       !< nlevels of cloud ice water (kg/kg)
    REAL(KIND=jprb), POINTER :: totalice(:)  !< nlevels of total ice (kg/kg)
    REAL(KIND=jprb), POINTER :: rain(:)      !< nlevels of rain (units: see mmr_snowrain)
    REAL(KIND=jprb), POINTER :: sp(:)        !< nlevels of solid precipitation (units: see mmr_snowrain)
  END TYPE profile_cloud_type

  !> Explicit optical parameters for IR scattering
  TYPE rttov_opt_param
    REAL(KIND=jprb), POINTER :: abs(:,:)     !< Absorption coef (nchannels,nlayers) (km-1)
    REAL(KIND=jprb), POINTER :: sca(:,:)     !< Scattering coef (nchannels,nlayers) (km-1)
    REAL(KIND=jprb), POINTER :: bpr(:,:)     !< b parameter (nchannels,nlayers) (no units)
    REAL(KIND=jprb), POINTER :: pha(:,:,:)   !< Phase function (nchannels,nlayers,nphangle), should be
                                             !! normalised such that integral over all angles is 4*pi
    REAL(KIND=jprb), POINTER :: phangle(:)   !< Angles over which phase fns defined (nphangle) (degrees)

    ! The following are for RTTOV internal purposes
    REAL(KIND=jprb)             :: minphadiff     !< Minimum difference between phase angles - for internal use only
    REAL(KIND=jprb),    POINTER :: cosphangle(:)  !< Cosine of phase angles (nphangle) - for internal use only
    INTEGER(KIND=jpim), POINTER :: iphangle(:)    !< Array indexes (size depends on minphadiff) - for internal use only
  END TYPE rttov_opt_param

  !> Output transmittances
  TYPE transmission_type
    REAL(KIND=jprb), POINTER  :: tau_total(:)             !< Surface-satellite transmittance (channels)
    REAL(KIND=jprb), POINTER  :: tau_levels(:,:)          !< Level-satellite transmittance (levels,channels)
    REAL(KIND=jprb), POINTER  :: tausun_total_path2(:)    !< Sun-surface-satellite solar transmittance
    REAL(KIND=jprb), POINTER  :: tausun_levels_path2(:,:) !< Sun-level-satellite solar transmittance for each level
    REAL(KIND=jprb), POINTER  :: tausun_total_path1(:)    !< Surface-satellite solar transmittance
    REAL(KIND=jprb), POINTER  :: tausun_levels_path1(:,:) !< Level-satellite solar transmittance for each level
  END TYPE

  !> Output radiances and corresponding brightness temperatures and reflectances (BRFs)
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K; BRFs are unitless.
  TYPE radiance_type
     ! Array size is (nchannels) or (nlayers,nchannels)
    REAL(KIND=jprb), POINTER  :: clear(:)       !< Clear sky radiance
    REAL(KIND=jprb), POINTER  :: total(:)       !< Cloudy radiance for given cloud
    REAL(KIND=jprb), POINTER  :: bt_clear(:)    !< Brightness temp equivalent to clear radiance
    REAL(KIND=jprb), POINTER  :: bt(:)          !< Brightness temp equivalent to total radiance
    REAL(KIND=jprb), POINTER  :: refl_clear(:)  !< Reflectance calculated from clear radiance
    REAL(KIND=jprb), POINTER  :: refl(:)        !< Reflectance calculated from total radiance
    REAL(KIND=jprb), POINTER  :: overcast(:,:)  !< Overcast radiance for opaque cloud at level bounding
                                                !!   bottom of each layer
    REAL(KIND=jprb), POINTER  :: cloudy(:)      !< 100% cloudy radiance for given cloud (simple cloud scheme)
                                                !!   or same as total (addclouds/addaerosl true)
  END TYPE radiance_type

  !> Secondary radiances optionally calculated in direct model only for clear-sky with no solar contribution
  TYPE radiance2_type
    ! Array size is (nchannels) or (nlayers,nchannels)
    REAL(KIND=jprb), POINTER  :: upclear(:)     !< Clear sky upwelling radiance without reflection term
    REAL(KIND=jprb), POINTER  :: dnclear(:)     !< Clear sky downwelling radiance
    REAL(KIND=jprb), POINTER  :: refldnclear(:) !< Reflected clear sky downwelling radiance
    REAL(KIND=jprb), POINTER  :: up(:,:)        !< Sum( B * dT ) above cloud upwelling radiance from each layer
    REAL(KIND=jprb), POINTER  :: down(:,:)      !< Sum( B / T**2 dT ) above cloud downwelling radiance from
                                                !!   each layer
    REAL(KIND=jprb), POINTER  :: surf(:,:)      !< Radiance at surface emitted from a black cloud
  END TYPE radiance2_type

  !> Output PC scores and reconstructed radiances from PC-RTTOV
  TYPE rttov_pccomp
    REAL(KIND=jprb), POINTER  :: pcscores(:)    !< Principal component scores
    REAL(KIND=jprb), POINTER  :: total_pccomp(:)!< Radiances reconstructed using principal components
    REAL(KIND=jprb), POINTER  :: bt_pccomp(:)   !< Brightness temp equivalent to radiances
                                                !! reconstructed using principal components
  END TYPE rttov_pccomp


  !> Configuration options that apply to all flavours of RTTOV
  TYPE rttov_opts_config
    LOGICAL(KIND=jplm) :: apply_reg_limits = .FALSE.                 !< Switch to restrict input profiles to coef training limits
    LOGICAL(KIND=jplm) :: verbose          = .TRUE.                  !< Switch for verbose output
    LOGICAL(KIND=jplm) :: do_checkinput    = .TRUE.                  !< Switch to apply internal profile checking
  END TYPE

  !> Options for PC-RTTOV
  TYPE rttov_opts_pc
    LOGICAL(KIND=jplm) :: addpc     = .FALSE.   !< Switch to enable PC-RTTOV
    INTEGER(KIND=jpim) :: ipcbnd    = -1_jpim   !< PC spectral band
    INTEGER(KIND=jpim) :: ipcreg    = -1_jpim   !< PC predictor channel set
    LOGICAL(KIND=jplm) :: addradrec = .FALSE.   !< Switch for calculation of reconstructed radiances
  END TYPE

  !> General radiative transfer options
  TYPE rttov_opts_rt_all
    LOGICAL(KIND=jplm) :: addrefrac     = .FALSE. !< Switch to enable atmospheric refraction
    LOGICAL(KIND=jplm) :: switchrad     = .FALSE. !< Switch for input units in AD/K models
    LOGICAL(KIND=jplm) :: use_q2m       = .TRUE.  !< Switch to enable use of 2m q variable
    LOGICAL(KIND=jplm) :: do_lambertian = .FALSE. !< Switch for setting Lambertian reflection (IR and MW)
  END TYPE

  !> VIS/IR-only radiative transfer options
  TYPE rttov_opts_rt_ir
    TYPE(rttov_opts_pc) :: pc                            !< PC-RTTOV options

    LOGICAL(KIND=jplm) :: addsolar           = .FALSE.   !< Switch to enable solar simulations
    LOGICAL(KIND=jplm) :: do_nlte_correction = .FALSE.   !< Switch to enable NLTE bias correction
    LOGICAL(KIND=jplm) :: addaerosl          = .FALSE.   !< Switch to enable IR aerosol calculations
    LOGICAL(KIND=jplm) :: addclouds          = .FALSE.   !< Switch to enable IR cloudy calculations
    LOGICAL(KIND=jplm) :: user_aer_opt_param = .FALSE.   !< Switch to supply aerosol optical properties explicitly per channel
    LOGICAL(KIND=jplm) :: user_cld_opt_param = .FALSE.   !< Switch to supply cloud optical properties explicitly per channel
    REAL(KIND=jprb)    :: cldstr_threshold   = -1.0_jprb !< Ignore cloud streams with weights lower than this
    LOGICAL(KIND=jplm) :: cldstr_simple      = .FALSE.   !< Switch for simplified cloud stream option (USE WITH CAUTION)
    LOGICAL(KIND=jplm) :: do_lambertian      = .FALSE.   !< Switch for setting Lambertian reflection (IR only)

    LOGICAL(KIND=jplm) :: ozone_data = .FALSE.           !< Switch to enable input of o3 profile
    LOGICAL(KIND=jplm) :: co2_data   = .FALSE.           !< Switch to enable input of co2 profile
    LOGICAL(KIND=jplm) :: n2o_data   = .FALSE.           !< Switch to enable input of n2o profile
    LOGICAL(KIND=jplm) :: co_data    = .FALSE.           !< Switch to enable input of co profile
    LOGICAL(KIND=jplm) :: ch4_data   = .FALSE.           !< Switch to enable input of ch4 profile
  END TYPE

  !> MW-only radiative transfer options
  TYPE rttov_opts_rt_mw
    INTEGER(KIND=jpim) :: fastem_version       = 5_jpim  !< FASTEM version (1-6)
    LOGICAL(KIND=jplm) :: clw_data             = .FALSE. !< Switch to enable input of cloud liquid water profile
    LOGICAL(KIND=jplm) :: do_lambertian        = .FALSE. !< Switch for setting Lambertian reflection (MW only)
    LOGICAL(KIND=jplm) :: supply_foam_fraction = .FALSE. !< Supply a foam fraction to FASTEM
  END TYPE

  !> Options for internal vertical interpolation and vertical grid setup
  TYPE rttov_opts_interp
    LOGICAL(KIND=jplm) :: addinterp        = .FALSE.        !< Switch to enable RTTOV interpolator
    INTEGER(KIND=jpim) :: interp_mode      = interp_rochon  !< Interpolation mode (1-5, see user guide)
    LOGICAL(KIND=jplm) :: lgradp           = .FALSE.        !< Switch to make pressure an active variable in TL/AD/K models
    LOGICAL(KIND=jplm) :: spacetop         = .TRUE.         !< Switch to assume space boundary at top-most input pressure level
    LOGICAL(KIND=jplm) :: reg_limit_extrap = .FALSE.        !< Switch to extrapolate input profiles using regression limits
  END TYPE

  !> RTTOV options structure
  TYPE rttov_options
    TYPE(rttov_opts_config)  :: config          !< General configuration options
    TYPE(rttov_opts_rt_all)  :: rt_all          !< General RT options
    TYPE(rttov_opts_rt_ir)   :: rt_ir           !< VIS/IR RT options
    TYPE(rttov_opts_rt_mw)   :: rt_mw           !< MW RT options
    TYPE(rttov_opts_interp)  :: interpolation   !< Interpolation options
  END TYPE

  !> RTTOV-SCATT options structure: RTTOV-SCATT deliberately does
  !! not give user control over certain core RT options.
  TYPE rttov_options_scatt
    TYPE(rttov_opts_config) :: config                               !< General configuration options
    LOGICAL(KIND=jplm)      :: lusercfrac           = .FALSE.       !< Switch to enable user-specified effective cloud fraction
    LOGICAL(KIND=jplm)      :: use_q2m              = .TRUE.        !< Switch to enable use of 2m q variable
    INTEGER(KIND=jpim)      :: fastem_version       = 5_jpim        !< FASTEM version (1-6)
    LOGICAL(KIND=jplm)      :: supply_foam_fraction = .FALSE.       !< Supply a foam fraction to FASTEM
    INTEGER(KIND=jpim)      :: interp_mode          = interp_rochon !< Interpolation mode (1-5, see user guide)
    LOGICAL(KIND=jplm)      :: reg_limit_extrap     = .FALSE.       !< Switch to extrapolate input profiles using regression limits
  END TYPE



  ! ---------------------------------------------------------------------------
  ! Internal RTTOV structures
  ! ---------------------------------------------------------------------------

  !> @internal Satellite geometry
  TYPE geometry_type
    REAL(KIND=jprb) :: sinzen
    REAL(KIND=jprb) :: sinzen_sq
    REAL(KIND=jprb) :: coszen
    REAL(KIND=jprb) :: coszen_sq
    REAL(KIND=jprb) :: seczen
    REAL(KIND=jprb) :: seczen_sq
    REAL(KIND=jprb) :: seczen_sqrt
    REAL(KIND=jprb) :: seczen_minus1
    REAL(KIND=jprb) :: seczen_minus1_sq
    REAL(KIND=jprb) :: sinview
    REAL(KIND=jprb) :: sinview_sq
    REAL(KIND=jprb) :: cosview_sq
    REAL(KIND=jprb) :: normzen
    REAL(KIND=jprb) :: viewang

    REAL(KIND=jprb) :: sinzen_sun
    REAL(KIND=jprb) :: sinlat
    REAL(KIND=jprb) :: coslat
  END TYPE geometry_type

  !> @internal The actual predictor arrays
  TYPE rttov_path_pred
    REAL(KIND=jprb), POINTER     :: mixedgas(:,:,:)   ! (nmixed,  nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: watervapour(:,:,:)! (nwater,  nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: ozone(:,:,:)      ! (nozone,  nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: wvcont(:,:,:)     ! (nwvcont, nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: co2(:,:,:)        ! (nco2,    nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: n2o(:,:,:)        ! (nn2o,    nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: co(:,:,:)         ! (nco,     nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: ch4(:,:,:)        ! (nch4,    nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: clw(:,:)          ! (         nlayers, nprofiles)
    REAL(KIND=jprb), POINTER     :: pmc(:,:,:,:)      ! pressure modulated cell (npmc,nlevels,nprofiles,nchannels)
  END TYPE rttov_path_pred

  !> @internal Predictors
  TYPE predictors_type
    ! the nxxxx could be set to 0 to indicate the abscence
    ! of the predictor, in that case there is no need to
    ! allocate the corresponding predictor
    INTEGER(KIND=jpim) :: nlevels   ! number of levels for predictors (all same)
    INTEGER(KIND=jpim) :: nmixed    ! number of variables for Mixed Gases
    INTEGER(KIND=jpim) :: nwater    ! number of variables for Water Vapour
    INTEGER(KIND=jpim) :: nozone    ! number of variables for Ozone
    INTEGER(KIND=jpim) :: nwvcont   ! number of variables for WV Continuum
    INTEGER(KIND=jpim) :: nco2      ! number of variables for CO2
    INTEGER(KIND=jpim) :: nn2o      ! number of variables for N2O
    INTEGER(KIND=jpim) :: nco       ! number of variables for CO
    INTEGER(KIND=jpim) :: nch4      ! number of variables for CH4
    INTEGER(KIND=jpim) :: ncloud    ! number of variables for MW Cloud
    INTEGER(KIND=jpim) :: npmc      ! number of variables for pressure modulated cell correction

    TYPE(rttov_path_pred) :: path1  ! Predictors for surface-satellite path (always required)
    TYPE(rttov_path_pred) :: path2  ! Predictors for sun-surface-satellite path (only required for solar)
  END TYPE predictors_type

  !> @internal Fast coefficients
  !  Separate structure to allow for non-PW and PW if necessary
  TYPE rttov_fast_coef
    ! separate arrays to allow different number of variables for each gas
    REAL(KIND=jprb), POINTER      :: mixedgas(:,:,:)     ! Mixed gases coefs  (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: watervapour(:,:,:)  ! Water vapour coefs (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: ozone(:,:,:)        ! Ozone coefs        (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: wvcont(:,:,:)       ! WV Cont coefs      (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: co2(:,:,:)          ! CO2 coefs          (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: n2o(:,:,:)          ! N2O coefs          (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: co(:,:,:)           ! CO coefs           (levels, channels, variables)
    REAL(KIND=jprb), POINTER      :: ch4(:,:,:)          ! CH4 coefs          (levels, channels, variables)
  END TYPE rttov_fast_coef

  !> @internal NLTE coefficients
  TYPE rttov_nlte_coef
    REAL(jprb), POINTER :: coef(:,:,:,:) ! ncoef x nsat x nsol x nchan
    REAL(jprb), POINTER :: sol_zen_angle(:), sat_zen_angle(:)
    REAL(jprb), POINTER :: cos_sol(:), sec_sat(:)
    INTEGER(jpim) :: ncoef, nsol, nsat, nchan
    INTEGER(jpim) :: start_chan, end_chan
    REAL(jprb)    :: max_sat_angle
  END TYPE rttov_nlte_coef

  !> @internal Optical depth coefs
  TYPE rttov_coef
     ! Structure for the storage of RTTOV coefficients
     ! this may differ from what is stored in the coefficient files especially
     ! for the units (ie kg/kg to ppmv)
     ! Gases are separated in MxG WV O3
     ! Number of levels is the same for all gases (taken from MxG).
     !
    INTEGER(KIND=jpim) :: id_platform  ! platform   (see documentation or MOD_CPARAM)
    INTEGER(KIND=jpim) :: id_sat    ! satellite  (.....)
    INTEGER(KIND=jpim) :: id_inst    ! instrument (.....)
    INTEGER(KIND=jpim) :: id_sensor  ! sensor
       !  1 = Infrared
       !  2 = Micro Wave
       !  3 = High resolution
    INTEGER(KIND=jpim) :: id_comp_lvl  ! RTTOV coefficient file version number
    INTEGER(KIND=jpim) :: id_comp_pc   ! Principal component coefficient file version number
    INTEGER(KIND=jpim) :: id_creation_date(3)  ! YYYY MM DD
    CHARACTER (LEN=80)    :: id_creation    ! Creation comment
    CHARACTER (LEN=32)    :: id_common_name  ! usual name of the satellite
    CHARACTER (LEN=132)   :: line_by_line(100)
    CHARACTER (LEN=132)   :: readme_srf(100) ! readme Spectral Response Function


       !FAST_MODEL_VARIABLES section
    CHARACTER (LEN=32)    :: fmv_model_def  ! FMV definition (RTTOV6 OPTRAN RTTOV7)
    INTEGER(KIND=jpim)               :: fmv_model_ver  ! fast model version compatibility level
    INTEGER(KIND=jpim)               :: fmv_ori_nchn   ! number of channels in original file
    INTEGER(KIND=jpim)               :: fmv_chn        ! number of channels read from file into coef structure
    INTEGER(KIND=jpim)               :: fmv_gas        ! number of gases in file
    INTEGER(KIND=jpim), POINTER      :: fmv_gas_id(:)    ! gas id. number i gas_id list (fmv_gas)
    INTEGER(KIND=jpim), POINTER      :: fmv_gas_pos(:)   ! respective position of each gas of gas_id list (ngases_max)
    INTEGER(KIND=jpim), POINTER      :: fmv_var(:)       ! number of variables/predictors by gas (fmv_gas)
    INTEGER(KIND=jpim), POINTER      :: fmv_coe(:)       ! number of coefficients by gas (fmv_gas)
    INTEGER(KIND=jpim), POINTER      :: fmv_int(:)       ! number of spectral intervals by gas (fmv_gas)
    INTEGER(KIND=jpim), POINTER      :: fmv_lvl(:)       ! number of levels(pres/absorber) by gas (fmv_gas)

    INTEGER(KIND=jpim)               :: nmixed         ! number of variables/predictors for Mixed Gases
    INTEGER(KIND=jpim)               :: nwater         ! number of variables/predictors for Water Vapour
    INTEGER(KIND=jpim)               :: nozone         ! number of variables/predictors for Ozone
    INTEGER(KIND=jpim)               :: nwvcont        ! number of variables/predictors for WV continuum
    INTEGER(KIND=jpim)               :: nco2           ! number of variables/predictors for CO2
    INTEGER(KIND=jpim)               :: nn2o           ! number of variables/predictors for N2O
    INTEGER(KIND=jpim)               :: nco            ! number of variables/predictors for CO
    INTEGER(KIND=jpim)               :: nch4           ! number of variables/predictors for CH4
    INTEGER(KIND=jpim)               :: nlevels        ! number of levels(pres/absorber) same for all gases
    INTEGER(KIND=jpim)               :: nlayers        ! number of layers(pres/absorber) nlevels-1
    LOGICAL(KIND=jplm)               :: IncZeeman      ! Flag to include Zeeman effect for this sensor
    INTEGER(KIND=jpim)               :: ncmixed        ! number of coefficients for Mixed Gases
    INTEGER(KIND=jpim)               :: ncwater        ! number of coefficients for Water Vapour
    INTEGER(KIND=jpim)               :: ncozone        ! number of coefficients for Ozone
    INTEGER(KIND=jpim)               :: ncwvcont       ! number of coefficients for WV continuum
    INTEGER(KIND=jpim)               :: ncco2          ! number of coefficients for CO2
    INTEGER(KIND=jpim)               :: ncn2o          ! number of coefficients for N2O
    INTEGER(KIND=jpim)               :: ncco           ! number of coefficients for CO
    INTEGER(KIND=jpim)               :: ncch4          ! number of coefficients for CH4

    ! JAH - to remove?
    INTEGER(KIND=jpim)               :: nintmixed
    INTEGER(KIND=jpim)               :: nintwater
    INTEGER(KIND=jpim)               :: nintozone
    INTEGER(KIND=jpim)               :: nintwvcont
    INTEGER(KIND=jpim)               :: nintco2
    INTEGER(KIND=jpim)               :: nintn2o
    INTEGER(KIND=jpim)               :: nintco
    INTEGER(KIND=jpim)               :: nintch4

       !GAZ_UNITS section
       ! gases are in the order of gas id codes
    INTEGER(KIND=jpim), POINTER      :: gaz_units(:)  ! unit of gas concentration for each gas
                                                      ! ppmv over dry air used in RTTOV optical depth calculations
       !FILTER_FUNCTIONS section  array size is fmv_chn
    LOGICAL(KIND=jplm) :: ff_val_bc    ! are any band corrections to be applied?
    LOGICAL(KIND=jplm) :: ff_val_gam   ! any gamma corrections?
    INTEGER(KIND=jpim) ,POINTER :: ff_ori_chn(:)   ! original chan number
    INTEGER(KIND=jpim) ,POINTER :: ff_val_chn(:)   ! validity of the channel (1=OK)
    REAL(KIND=jprb) ,POINTER :: ff_cwn (:)         ! cental wave number (cm-1)
    REAL(KIND=jprb) ,POINTER :: ff_bco (:)         ! band correction offset (K)
    REAL(KIND=jprb) ,POINTER :: ff_bcs (:)         ! band correction slope (K/K)
    REAL(KIND=jprb) ,POINTER :: ff_gam (:)         ! gamma factor transm. correction

       !TRANSMITTANCE_TRESHOLD section  array size is fmv_chn
    INTEGER(KIND=jpim) ,POINTER :: tt_chn(:)
    INTEGER(KIND=jpim) ,POINTER :: tt_val_chn(:)
    REAL(KIND=jprb)    ,POINTER :: tt_cwn (:)
    REAL(KIND=jprb)    ,POINTER :: tt_a0(:)
    REAL(KIND=jprb)    ,POINTER :: tt_a1(:)

       !PLANCK_WEIGHTED section array size if fmv_chn
    INTEGER(KIND=jpim) ,POINTER :: pw_chn(:)
    INTEGER(KIND=jpim) ,POINTER :: pw_val_chn(:)   ! 0 => non-PW thermal coefs, 1 => PW thermal coefs

       !SOLAR_SPECTRUM section array size is fmv_chn
    INTEGER(KIND=jpim) ,POINTER :: ss_chn(:)
    INTEGER(KIND=jpim) ,POINTER :: ss_val_chn(:)
    REAL(KIND=jprb)    ,POINTER :: ss_cwn (:)
    REAL(KIND=jprb)    ,POINTER :: ss_solar_spectrum(:)

       !WATER_OPTICAL_CONSTANT section array size is fmv_chn
    INTEGER(KIND=jpim) ,POINTER :: woc_chn(:)
    REAL(KIND=jprb)    ,POINTER :: woc_cwn (:)
    Complex(KIND=jprb) ,POINTER :: woc_waopc_ow(:)
    Complex(KIND=jprb) ,POINTER :: woc_waopc_fw(:)

       !WAVE_SPECTRUM section array size is ws_nomega
       !Data used to compute the frequency spectrum of the JONSWAP
       !wave model surface wave.
    INTEGER(KIND=jpim)          :: ws_nomega
    REAL(KIND=jprb)  ,POINTER   :: ws_npoint(:)
    REAL(KIND=jprb)  ,POINTER   :: ws_k_omega(:)

       !FUNDAMENTAL_CONSTANTS section
    REAL(KIND=jprb) :: fc_planck_c1      ! first radiation constant (mW/(m2*sr*cm-4))
    REAL(KIND=jprb) :: fc_planck_c2      ! second radiation constant (cm*K)
    REAL(KIND=jprb) :: fc_sat_height     ! satellite nominal altitude (km)

       !FASTEM section
    INTEGER(KIND=jpim) :: fastem_ver      ! fastem version number
    INTEGER(KIND=jpim), POINTER :: fastem_polar(:)  ! polarisation of each channel
       ! 0 = 0.5 V+H
       ! 1 = 90 - incident angle
       ! 2 = incident angle
       ! 3 = vertical
       ! 4 = horizontal
       ! 5 = V+H
       ! Full stokes vector

       !SSIREM section     array size is fmv_chn
       ! ems =   ssirem_a0
       !       - ssirem_a1*(zen**ssirem_xzn1)
       !       - ssirem_a2*(zen**ssirem_xzn2)
       ! where zen is satellite zenith angle in degrees, divided by 60.
    INTEGER(KIND=jpim) :: ssirem_ver                ! version number
    INTEGER(KIND=jpim),  POINTER  :: ssirem_chn(:)   ! original chan number
    REAL(KIND=jprb),  POINTER     :: ssirem_a0(:)    ! constant coef
    REAL(KIND=jprb),  POINTER     :: ssirem_a1(:)    ! first coef
    REAL(KIND=jprb),  POINTER     :: ssirem_a2(:)    ! second coef
    REAL(KIND=jprb),  POINTER     :: ssirem_xzn1(:)  ! 1st exponent on zenith angle
    REAL(KIND=jprb),  POINTER     :: ssirem_xzn2(:)  ! 2nd exponent on zenith angle

       !REFERENCE_PROFILE section  defined on Mixed gases pressure levels
       ! Not working for OPTRAN gas absorber levels
       ! gases are in the order of gas id codes
       ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(KIND=jprb), POINTER      :: ref_prfl_p(:)     ! pressure  (hPa)       (levels)
    REAL(KIND=jprb), POINTER      :: ref_prfl_t(:,:)   ! temperature (K)       (levels, gases)
    REAL(KIND=jprb), POINTER      :: ref_prfl_mr(:,:)  ! mixing ratio (ppmv)   (levels, gases)
       !PROFILE_LIMITS section
       ! gases are in the order of gas id codes
       ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(KIND=jprb), POINTER      :: lim_prfl_p(:)       ! pressure  (hPa)       (levels)
    REAL(KIND=jprb), POINTER      :: lim_prfl_tmax(:)    ! max temperature (K)   (levels)
    REAL(KIND=jprb), POINTER      :: lim_prfl_tmin(:)    ! min temperature (K)   (levels)
    REAL(KIND=jprb), POINTER      :: lim_prfl_gmax(:,:)  ! max mixing r (ppmv) (levels, gases)
    REAL(KIND=jprb), POINTER      :: lim_prfl_gmin(:,:)  ! min mixing r (ppmv) (levels, gases)


       !FAST_COEFFICIENTS/SOLAR_FAST_COEFFICIENTS section
       ! For non-PW instruments, "solar" will point to "thermal" coefs
       ! For instruments with solar-affected PW channels, both thermal and solar
       ! structures will be populated from coef file
    TYPE(rttov_fast_coef), POINTER :: thermal             ! FAST_COEFFICIENTS
    TYPE(rttov_fast_coef), POINTER :: solar               ! SOLAR_FAST_COEFFICIENTS when present
    LOGICAL(KIND=jplm)             :: solarcoef           ! .TRUE. if solar fast coefs present in file: need
                                                          ! to know whether solar points to thermal or is
                                                          ! allocated separately.
    LOGICAL(KIND=jplm)             :: nltecoef = .FALSE.  ! .TRUE. if nlte coefs present
    TYPE(rttov_nlte_coef), POINTER :: nlte_coef           ! nlte_coef

         !PRESSURE_MODULATED_CELL section
    LOGICAL(KIND=jplm)       :: pmc_shift = .FALSE.  ! .TRUE. if pmc shift coefs present
    REAL(KIND=jprb)          :: pmc_lengthcell       ! cell length (cm)
    REAL(KIND=jprb), POINTER :: pmc_pnominal(:)      ! nominal cell pressure (hPa) - nchannels
    REAL(KIND=jprb)          :: pmc_tempcell          ! cell temperature (K)
    REAL(KIND=jprb)          :: pmc_betaplus1        ! co2 band-average: self-HW/air-HW
    INTEGER(KIND=jpim)       :: pmc_nlay             ! number of layers used
    INTEGER(KIND=jpim)       :: pmc_nvar             ! number of variables used
    REAL(KIND=jprb), POINTER :: pmc_coef(:,:,:)     ! pressure moodulated cell corrections - nlevels, nchannels, nvariables
    REAL(KIND=jprb), POINTER :: pmc_ppmc(:)          ! actual cell pressure (hPa) - nchannels

     ! JAH - to remove?
    REAL(KIND=jprb), POINTER      :: mixedgasint(:,:)
    REAL(KIND=jprb), POINTER      :: watervapourint(:,:)
    REAL(KIND=jprb), POINTER      :: ozoneint(:,:)
    REAL(KIND=jprb), POINTER      :: wvcontint(:,:)
    REAL(KIND=jprb), POINTER      :: co2int(:,:)
    REAL(KIND=jprb), POINTER      :: n2oint(:,:)
    REAL(KIND=jprb), POINTER      :: coint(:,:)
    REAL(KIND=jprb), POINTER      :: ch4int(:,:)

       ! Auxillary variables
    REAL(KIND=jprb)               :: ratoe       ! ratio (H+R)/R  H=sat height, R=Earth radius
    INTEGER(KIND=jpim)            :: mwcldtop    ! Upper layer for MW LWP calcs
    REAL(KIND=jprb), POINTER      :: planck1(:)        ! C1 * Nu**3
    REAL(KIND=jprb), POINTER      :: planck2(:)        ! C2 * Nu
    REAL(KIND=jprb), POINTER      :: frequency_ghz(:)  ! frequency in GHz

       ! other predictor variables see Science and Validation report
    REAL(KIND=jprb), POINTER      :: dp(:)        ! interval between standard p levels (hPa)
    REAL(KIND=jprb), POINTER      :: dpp(:)       ! pressure based variable (hPa**2)
    REAL(KIND=jprb), POINTER      :: tstar(:)     ! layer temp (K)
    REAL(KIND=jprb), POINTER      :: to3star(:)   ! layer temp for O3 calculations (K)
    REAL(KIND=jprb), POINTER      :: wstar(:)     ! layer WV  (ppmv)
    REAL(KIND=jprb), POINTER      :: ostar(:)     ! layer O3  (ppmv)
    REAL(KIND=jprb), POINTER      :: co2star(:)   ! layer co2 (ppmv)
    REAL(KIND=jprb), POINTER      :: n2ostar(:)   ! layer n2o (ppmv)
    REAL(KIND=jprb), POINTER      :: costar(:)    ! layer co  (ppmv)
    REAL(KIND=jprb), POINTER      :: ch4star(:)   ! layer ch4 (ppmv)
  END TYPE rttov_coef

  !> @internal RTTOV-SCATT coefs
  TYPE rttov_scatt_coef
    ! Structure for the storage of RTTOV_SCATT coefficients
    INTEGER(KIND=jpim) :: nhydro ! Number of hydrometeors in computation
    INTEGER(KIND=jpim) :: mtype  ! Number of hydrometeors     in Mie tables
    INTEGER(KIND=jpim) :: mfreqm ! Number of frequencies      in Mie tables
    INTEGER(KIND=jpim) :: mtemp  ! Number of temperature bins in Mie tables
    INTEGER(KIND=jpim) :: mwc    ! Number of water bins       in Mie tables
    REAL(KIND=jprb)    :: offset_temp_rain       ! temperature offset in table for rain type
    REAL(KIND=jprb)    :: offset_temp_sp         ! temperature offset in table for solid prec. type
    REAL(KIND=jprb)    :: offset_temp_liq        ! temperature offset in table for cloud water type
    REAL(KIND=jprb)    :: offset_temp_ice        ! temperature offset in table for cloud ice type
    REAL(KIND=jprb)    :: offset_temp_totalice   ! temperature offset in table for total ice type
    REAL(KIND=jprb)    :: offset_water           ! liquid/ice water offset in table
    REAL(KIND=jprb)    :: scale_water            ! log10(liquid/ice water) scaling factor in table
    REAL(KIND=jprb)    :: from_scale_water       ! 10**(1._JPRB/scale_water)
    REAL(KIND=jprb)    :: conv_rain(2)           ! coefficients for rain unit conversion (mm.h-1 to g.m-3)
    REAL(KIND=jprb)    :: conv_sp  (2)           ! coefficients for solid prec. unit conversion (mm.h-1 to g.m-3)
    REAL(KIND=jprb)    :: conv_liq (2)           ! coefficients for cloud water conversion (not used)
    REAL(KIND=jprb)    :: conv_ice (2)           ! coefficients for cloud ice conversion   (not used)
    REAL(KIND=jprb)    :: conv_totalice (2)      ! coefficients for total ice conversion   (not used)
    REAL(KIND=jprb), POINTER :: mie_freq(:)      ! list of frequencies in Mie table
    REAL(KIND=jprb), POINTER :: ext(:,:,:,:)     ! extinction coefficent table
    REAL(KIND=jprb), POINTER :: ssa(:,:,:,:)     ! single scattering albedo table
    REAL(KIND=jprb), POINTER :: asp(:,:,:,:)     ! assymetry parameter table

  END TYPE rttov_scatt_coef

  !> @internal Surface and cloud fraction
  TYPE profile_aux_s
    INTEGER(KIND=jpim) :: nearestlev_surf ! nearest model level above surface
    REAL(KIND=jprb)    :: pfraction_surf  ! pressure fraction of surface in model layer (hPa)
    INTEGER(KIND=jpim) :: nearestlev_ctp  ! nearest model level above cloud top
    REAL(KIND=jprb) :: pfraction_ctp   ! pressure fraction of cloud top pressure in layer (hPa)
    REAL(KIND=jprb) :: cfraction       ! cloud fraction (0 - 1) 1 for 100% cloud cover
  END TYPE profile_aux_s

  !> @internal Auxillary profile variables
  TYPE profile_aux
    LOGICAL(KIND=jplm) :: on_coef_levels
    TYPE(profile_aux_s), POINTER :: s(:)
    REAL(KIND=jprb), POINTER :: debye_prof(:,:,:)  ! Debye terms
    REAL(KIND=jprb), POINTER :: relhum(:,:)        !Relative humidity
    REAL(KIND=jprb), POINTER :: relhumref(:,:)
    REAL(KIND=jprb), POINTER :: dg(:,:)            !Generalized effective diameter
    REAL(KIND=jprb), POINTER :: fac1_dg(:,:)       !Intermediate variables used to compute the
    REAL(KIND=jprb), POINTER :: fac2_dg(:,:)       !generalized diameter
    REAL(KIND=jprb), POINTER :: fac3_dg(:,:)
    INTEGER(KIND=jpim), POINTER :: iaertyp(:,:,:)
    INTEGER(KIND=jpim), POINTER :: iaernum(:,:)

    REAL(KIND=jprb), POINTER :: t_layer(:,:)       ! avg layer temperature
    REAL(KIND=jprb), POINTER :: w_layer(:,:)       ! avg layer humidity
    REAL(KIND=jprb), POINTER :: o3_layer(:,:)      ! avg layer humidity
    REAL(KIND=jprb), POINTER :: dt(:,:)            ! deviation from ref prof
    REAL(KIND=jprb), POINTER :: dto(:,:)           ! deviation from ref prof
    REAL(KIND=jprb), POINTER :: tr(:,:), tr_r(:,:)  ! ratio t / ref_t
    REAL(KIND=jprb), POINTER :: wr(:,:), wr_sqrt(:,:), wr_rsqrt(:,:)           !
    REAL(KIND=jprb), POINTER :: or(:,:), or_sqrt(:,:)
    REAL(KIND=jprb), POINTER :: tw(:,:), tw_sqrt(:,:), tw_4rt(:,:)            ! ratio t / ref_t
    REAL(KIND=jprb), POINTER :: ww(:,:), ww_r(:,:)            !
    REAL(KIND=jprb), POINTER :: ow(:,:), ow_r(:,:), ow_sqrt(:,:), ow_rsqrt(:,:)           !
    REAL(KIND=jprb), POINTER :: SUM(:,:)
  END TYPE profile_aux

  !> @internal Auxillary profile variables for RTTOV_SCATT
  TYPE profile_scatt_aux
    REAL(KIND=jprb), POINTER :: cfrac(:)        ! horizontal cloud fraction (one value used for all layers)
    REAL(KIND=jprb), POINTER :: ems_bnd(:)      ! surface emissivity for boundary conditions
    REAL(KIND=jprb), POINTER :: ref_bnd(:)      ! surface emissivity for boundary conditions
    REAL(KIND=jprb), POINTER :: ems_cld(:)      ! surface emissivity taking into account cloud/rain impact on od
    REAL(KIND=jprb), POINTER :: ref_cld(:)      ! surface reflectivity taking into account cloud/rain impact on od
    REAL(KIND=jprb), POINTER :: dz(:,:)         ! layer depth   [km]
    REAL(KIND=jprb), POINTER :: tbd(:,:)        ! temperature at layer boundary [K]
    REAL(KIND=jprb), POINTER :: clw(:,:)        ! cloud liquid water (g/m3)
    REAL(KIND=jprb), POINTER :: ciw(:,:)        ! cloud ice water (g/m3)
    REAL(KIND=jprb), POINTER :: totalice(:,:)   ! total ice (g/m3)
    REAL(KIND=jprb), POINTER :: rain(:,:)       ! rain (g/m3)
    REAL(KIND=jprb), POINTER :: sp(:,:)         ! solid precipitation (g/m3)
    !RWS  REAL(KIND=jprb), POINTER :: mclayer(:)  ! upper level cloud layer
    INTEGER(KIND=jpim), POINTER :: mclayer(:)   ! upper level cloud layer
    REAL(KIND=jprb), POINTER :: delta(:,:)      ! (= ext*dz/coszen)
    REAL(KIND=jprb), POINTER :: tau(:,:)        ! optical depths (= exp(-delta))
    REAL(KIND=jprb), POINTER :: ext(:,:)        ! extinction coefficient integreated over hydrometeor types
    REAL(KIND=jprb), POINTER :: ssa(:,:)        ! single scattering albedo integreated over hydrometeor types
    REAL(KIND=jprb), POINTER :: asm(:,:)        ! asymetry parameter integreated over hydrometeor types [-1,1]
    REAL(KIND=jprb), POINTER :: lambda(:,:)     ! eddington approx. variable
                                    ! (= sqrt( 3*ext*ext*(1-ssa)*(1-ssa*asm) )
    REAL(KIND=jprb), POINTER :: h (:,:)         ! boundary condition variable (= 1.5_JPRB*ext(1-ssa*asm))
    REAL(KIND=jprb), POINTER :: b0(:,:)         ! lower level temperature
    REAL(KIND=jprb), POINTER :: b1(:,:)         ! temperature gradient
    REAL(KIND=jprb), POINTER :: bn(:,:)         ! upper level temperature
  END TYPE profile_scatt_aux

  !> @internal Path optical depths as predicted or interpolated (unitless)
  TYPE opdp_path_type
    REAL(KIND=jprb), POINTER :: atm_level(:,:) ! neg optical depth for thermal radiation (levels to space),
                                               ! size (levels, channels)
    REAL(KIND=jprb), POINTER :: sun_level_path2(:,:) ! neg optical depth for solar radiation (levels to space) for
                                                     ! combined sun-surface-satellite path, size (levels, channels)
  END TYPE opdp_path_type

  !> @internal Transmissions and optical depths (unitless)
  TYPE rttov_path_transmission
    REAL(KIND=jprb), POINTER  :: tau_surf_p(:,:)         ! transmittance from surface (streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_surf_p_r(:,:)       ! reciprocal transmittance from surface (streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_surf(:,:)           ! transmittance from surface (streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_surf_r(:,:)         ! reciprocal transmittance from surface (streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_level(:,:,:)        ! transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_level_r(:,:,:)      ! reciprocal transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_level_p(:,:,:)      ! transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    REAL(KIND=jprb), POINTER  :: tau_level_p_r(:,:,:)    ! reciprocal transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    REAL(KIND=jprb), POINTER  :: od_singlelayer(:,:,:)   ! single-layer optical depth
    REAL(KIND=jprb), POINTER  :: od_singlelayer_r(:,:,:) ! reciprocal single-layer optical depth
    REAL(KIND=jprb), POINTER  :: od_sfrac(:,:)
    REAL(KIND=jprb), POINTER  :: od_sfrac_r(:,:)
    REAL(KIND=jprb), POINTER  :: od_frac_ac(:,:)
    REAL(KIND=jprb), POINTER  :: tau_surf_ac(:,:)

    REAL(KIND=jprb), POINTER  :: fac2(:,:,:)    ! Mask for integration calculation: thermal and solar path1
  END TYPE

  !> @internal Auxillary transmittance variables
  TYPE transmission_type_aux
    REAL(KIND=jprb), POINTER  :: fac1(:,:,:)    ! Mask for integration calculation
    REAL(KIND=jprb), POINTER  :: surf_fac(:,:)
    REAL(KIND=jprb) :: anynegtau ! used to store information about presence of any negative transmittances

    TYPE(rttov_path_transmission), POINTER :: thermal_path1
    TYPE(rttov_path_transmission), POINTER :: solar_path2
    TYPE(rttov_path_transmission), POINTER :: solar_path1
  END TYPE transmission_type_aux

  !> @internal Auxillary radiance variables
  TYPE radiance_aux
    ! auxillary calculation arrays for RTE integration
    ! Direct model arrays need to be passed to TL AD and K codes
    ! array size is of (nchannels) or (nlevels, nchannels)
    REAL(KIND=jprb), POINTER :: air(:,:)
    REAL(KIND=jprb), POINTER :: surfair(:)
    REAL(KIND=jprb), POINTER :: skin(:)
    REAL(KIND=jprb), POINTER :: cosmic(:)
    REAL(KIND=jprb), POINTER :: air_t_eff(:,:)
    REAL(KIND=jprb), POINTER :: surf_t_eff(:)
    REAL(KIND=jprb), POINTER :: skin_t_eff(:)
    REAL(KIND=jprb), POINTER :: cosmic_t_eff(:)

    REAL(KIND=jprb), POINTER :: up(:,:,:)               ! sum( B * dT )
    REAL(KIND=jprb), POINTER :: down(:,:,:)             ! sum ( B / T**2 dT )
    REAL(KIND=jprb), POINTER :: up_solar(:,:,:)         ! sum( B * dT )
    REAL(KIND=jprb), POINTER :: down_solar(:,:,:)       ! sum ( B / T**2 dT )
    REAL(KIND=jprb), POINTER :: meanrad_up(:,:)
    REAL(KIND=jprb), POINTER :: meanrad_down(:,:)
    REAL(KIND=jprb), POINTER :: meanrad_up_solar(:,:)
    REAL(KIND=jprb), POINTER :: meanrad_down_solar(:,:)
    REAL(KIND=jprb), POINTER :: down_ref(:,:,:)
    REAL(KIND=jprb), POINTER :: down_ref_solar(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC1_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC2_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC3_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC4_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC5_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC6_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC7_2(:,:,:)
    REAL(KIND=jprb), POINTER :: FAC1_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC2_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC3_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC4_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC5_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC6_3(:,:)
    REAL(KIND=jprb), POINTER :: FAC7_3(:,:)
    REAL(KIND=jprb), POINTER :: cloudy(:,:)
  END TYPE radiance_aux

  !> @internal Raytracing variables
  TYPE raytracing_type
    REAL(KIND=jprb), POINTER :: LTICK   (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: HGPL    (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: DMAIR_R (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: REFRACTIVITY(:,:)   !(levels,profiles)
    REAL(KIND=jprb), POINTER :: R       (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: R_R     (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: Z_R     (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: RATOESUN(:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: RATOESAT(:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: ZASUN   (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: ZASAT   (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: INT     (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: ZTEMP   (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: PPW     (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: DISPCO2 (:,:)       !(levels,profiles)
    REAL(KIND=jprb), POINTER :: PATHSAT (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: PATHSAT_rsqrt (:,:) !(layers,profiles)
    REAL(KIND=jprb), POINTER :: PATHSAT_sqrt (:,:)  !(layers,profiles)
    REAL(KIND=jprb), POINTER :: PATHSUN (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: PATHEFF (:,:)       !(layers,profiles)
    REAL(KIND=jprb), POINTER :: CO2_CM  (:)         !(profiles)
  END TYPE raytracing_type

  !> @internal Sea-surface BRDF model variables
  TYPE sunglint_type_s
    REAL(KIND=jprb)          :: CSI
    REAL(KIND=jprb)          :: ALFA
    REAL(KIND=jprb)          :: C_SHAD
    REAL(KIND=jprb)          :: P_PRIME
    REAL(KIND=jprb)          :: PXY_GAMMAXY
    REAL(KIND=jprb)          :: GAMMA_O
    REAL(KIND=jprb)          :: GAMMA_P
    REAL(KIND=jprb)          :: G_SHAD
    REAL(KIND=jprb)          :: GAMMAX
    REAL(KIND=jprb)          :: Q_SHAD
    REAL(KIND=jprb)          :: ZENSAT
    REAL(KIND=jprb)          :: ZENSUN
    REAL(KIND=jprb)          :: DAZNG
    REAL(KIND=jprb)          :: FAC1
    REAL(KIND=jprb)          :: A_SHAD
    REAL(KIND=jprb)          :: B_SHAD
    REAL(KIND=jprb)          :: LAMBDA_A
    REAL(KIND=jprb)          :: LAMBDA_B
    REAL(KIND=jprb)          :: X_U
    REAL(KIND=jprb)          :: ALFA1
    REAL(KIND=jprb)          :: OMEGA_M
    REAL(KIND=jprb)          :: WINDSP
    REAL(KIND=jprb)          :: WANGL
    REAL(KIND=jprb)          :: GAMMA_SQ
    REAL(KIND=jprb)          :: GLINT
    REAL(KIND=jprb)          :: OMEGA
  END TYPE sunglint_type_s

  !> @internal Sea-surface BRDF model variables
  TYPE sunglint_type
    TYPE(sunglint_type_s), POINTER :: s(:)
    REAL(KIND=jprb),POINTER  :: BETA (:,:)
    REAL(KIND=jprb),POINTER  :: PSI  (:,:)
  END TYPE sunglint_type

  !> @internal IR scattering variables
  TYPE transmission_scatt_ir_type
    REAL(KIND=jprb)   , POINTER  :: opdps         (:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpa         (:,:)
    REAL(KIND=jprb)   , POINTER  :: gpar          (:,:)
    REAL(KIND=jprb)   , POINTER  :: gpartot       (:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpscls      (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpacls      (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: gparcls       (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpaerla     (:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpcldla     (:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpsaer      (:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpaaer      (:,:)
    REAL(KIND=jprb)   , POINTER  :: gparaera      (:,:)
    REAL(KIND=jprb)   , POINTER  :: gparaer       (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphup        (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphdo        (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphupcls     (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: azphdocls     (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: azphuptot     (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphdotot     (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphaerup     (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphaerdo     (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphaerupa    (:,:)
    REAL(KIND=jprb)   , POINTER  :: azphaerdoa    (:,:)
    REAL(KIND=jprb)   , POINTER  :: phasintupref  (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: phasintdoref  (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpabs       (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpsca       (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpac        (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpacl       (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpacsun     (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpaclsun    (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: azphacup      (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: azphacdo      (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: opdpext       (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: ssa           (:,:,:)
  END TYPE transmission_scatt_ir_type

  !> @internal PC coefs for each predictor channel set
  TYPE rttov_coef_pccomp1
    INTEGER(KIND=jpim)           :: fmv_pc_npred          ! Number of predictors in the regression set
    INTEGER(KIND=jpim), POINTER  :: predictindex  (:)     ! Predictors channel indices
    REAL(KIND=jprb)   , POINTER  :: coefficients  (:,:)   ! Regression coefficients
    REAL(KIND=jprb)   , POINTER  :: coefficients_t  (:,:) ! Regression coefficients transposed
  END TYPE rttov_coef_pccomp1

  !> @internal PC eigenvectors
  TYPE rttov_coef_pccomp2
    REAL(KIND=jprb)   , POINTER  :: eigenvectors  (:,:)   ! Eigenvectors
    REAL(KIND=jprb)   , POINTER  :: eigenvectors_t(:,:)   ! Transposed Eigenvectors
  END TYPE rttov_coef_pccomp2

  !> @internal PC-RTTOV coefs
  TYPE rttov_coef_pccomp
    INTEGER(KIND=jpim)           :: fmv_pc_comp_pc
    INTEGER(KIND=jpim)           :: fmv_pc_cld
    INTEGER(KIND=jpim)           :: fmv_pc_msets          ! Maximum number of regression sets
    INTEGER(KIND=jpim)           :: fmv_pc_bands          ! Number of bands
    INTEGER(KIND=jpim)           :: fmv_pc_mnum           ! Maximum number of eigenvectors
    INTEGER(KIND=jpim)           :: fmv_pc_mchn           ! Maximum number of channels
    INTEGER(KIND=jpim)           :: fmv_pc_nchn           ! Number of channels
    INTEGER(KIND=jpim)           :: fmv_pc_nchn_noise     ! Number of channels for which instrument noise is available
    INTEGER(KIND=jpim)           :: fmv_pc_nche           ! Number of channels for which emissisity coefs are available
    INTEGER(KIND=jpim)           :: fmv_pc_gas            ! Number of gases for which a reference profile is given
    INTEGER(KIND=jpim), POINTER  :: fmv_pc_sets   (:)     ! Number of regression sets in each band
    INTEGER(KIND=jpim), POINTER  :: emiss_chn     (:)     ! Number of channels for which emissivity coefficients are
    REAL   (KIND=jprb), POINTER  :: emiss_c1      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c2      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c3      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c4      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c5      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c6      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c7      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c8      (:)     ! Emissivity coefficient
    REAL   (KIND=jprb), POINTER  :: emiss_c9      (:)     ! Emissivity coefficient
    INTEGER(KIND=jpim)           :: fmv_pc_nlev           ! Number of reference profile levels
    REAL(KIND=jprb), POINTER     :: ref_pc_prfl_p (:)     ! pressure  (hPa)       (levels)
    REAL(KIND=jprb), POINTER     :: ref_pc_prfl_mr(:,:)   ! mixing ratio (ppmv)   (levels)
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_tmin(:)   ! Profile limit :temperature
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_tmax(:)   ! Profile limit :temperature
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_qmin(:)   ! Profile limit :water vapour
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_qmax(:)   ! Profile limit :water vapour
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_ozmin(:)  ! Profile limit :ozone
    REAL(KIND=jprb), POINTER     :: lim_pc_prfl_ozmax(:)  ! Profile limit :ozone
    REAL(KIND=jprb)              :: lim_pc_prfl_pmin      ! Surface pressure
    REAL(KIND=jprb)              :: lim_pc_prfl_pmax      ! Surface pressure
    REAL(KIND=jprb)              :: lim_pc_prfl_tsmin     ! Surface temperature
    REAL(KIND=jprb)              :: lim_pc_prfl_tsmax     ! Surface temperature
    REAL(KIND=jprb)              :: lim_pc_prfl_skmin     ! Skin temperature
    REAL(KIND=jprb)              :: lim_pc_prfl_skmax     ! Skin temperature
    REAL(KIND=jprb)              :: lim_pc_prfl_wsmin     ! 10m wind speed
    REAL(KIND=jprb)              :: lim_pc_prfl_wsmax     ! 10m wind speed
    REAL(KIND=jprb), POINTER     :: co2_pc_ref    (:)     ! Fixed co2 profile to be used in the computation of PC's
    REAL(KIND=jprb), POINTER     :: n2o_pc_ref    (:)     ! Fixed n2o profile to be used in the computation of PC's
    REAL(KIND=jprb), POINTER     :: co_pc_ref     (:)     ! Fixed co  profile to be used in the computation of PC's
    REAL(KIND=jprb), POINTER     :: ch4_pc_ref    (:)     ! Fixed ch4 profile to be used in the computation of PC's
    REAL(KIND=jprb), POINTER     :: noise_in      (:)     ! Noise values for the channels whose radiances are
                                                          ! reconstrucetd using principal components
    REAL(KIND=jprb), POINTER     :: noise         (:)     ! Noise values for the channels whose radiances are
                                                          ! used as predictors in the computation of principal components
    REAL(KIND=jprb), POINTER     :: noise_r       (:)     ! Reciprocal noise
    INTEGER(KIND=jpim), POINTER  :: ff_ori_chn_in(:)
    REAL(KIND=jprb),    POINTER  :: ff_cwn_in(:)          ! central wave number of reconstructed radiances
    REAL(KIND=jprb),    POINTER  :: ff_bco_in (:)         ! band correction offset (K)
    REAL(KIND=jprb),    POINTER  :: ff_bcs_in (:)         ! band correction slope (K/K)
    REAL(KIND=jprb),    POINTER  :: planck1_in(:)         ! C1 * Nu**3
    REAL(KIND=jprb),    POINTER  :: planck2_in(:)         ! C2 * Nu
    TYPE(rttov_coef_pccomp1), POINTER:: pcreg     (:,:)
    TYPE(rttov_coef_pccomp2), POINTER:: eigen     (:)
  END TYPE rttov_coef_pccomp

  !> @internal IR scattering coefs
  TYPE rttov_coef_scatt_ir
    INTEGER(KIND=jpim)           :: fmv_aer_chn    ! number of channels for which aerosol optical parameters are stored
    INTEGER(KIND=jpim)           :: fmv_wcl_chn
    INTEGER(KIND=jpim)           :: fmv_icl_chn
    INTEGER(KIND=jpim)           :: fmv_aer_pha_chn
    INTEGER(KIND=jpim)           :: fmv_wcl_pha_chn
    INTEGER(KIND=jpim)           :: fmv_icl_pha_chn
    INTEGER(KIND=jpim)           :: fmv_aer_comp
    INTEGER(KIND=jpim)           :: fmv_wcl_comp
    INTEGER(KIND=jpim)           :: fmv_icl_comp
    INTEGER(KIND=jpim)           :: fmv_icl_ishp
    INTEGER(KIND=jpim)           :: fmv_aer_pha_ioff
    INTEGER(KIND=jpim)           :: fmv_wcl_pha_ioff
    INTEGER(KIND=jpim)           :: fmv_icl_pha_ioff
    INTEGER(KIND=jpim)           :: fmv_aer_ph
    INTEGER(KIND=jpim)           :: fmv_wcl_ph
    INTEGER(KIND=jpim)           :: fmv_icl_ph
    INTEGER(KIND=jpim)           :: icl_nabs
    INTEGER(KIND=jpim)           :: icl_nsca
    INTEGER(KIND=jpim)           :: icl_nbpr
    CHARACTER(LEN=4),   POINTER  :: fmv_aer_comp_name(:)
    CHARACTER(LEN=4),   POINTER  :: fmv_wcl_comp_name(:)
    CHARACTER(LEN=16),  POINTER  :: fmv_icl_comp_name(:,:)
    INTEGER(KIND=jpim), POINTER  :: fmv_aer_rh    (:)
    INTEGER(KIND=jpim), POINTER  :: fmv_wcl_rh    (:)
    REAL   (KIND=jprb), POINTER  :: fmv_aer_rh_val(:)
    REAL   (KIND=jprb), POINTER  :: fmv_wcl_rh_val(:)
    INTEGER(KIND=jpim), POINTER  :: ifmv_aer_ph_val(:)
    INTEGER(KIND=jpim), POINTER  :: ifmv_wcl_ph_val(:)
    INTEGER(KIND=jpim), POINTER  :: ifmv_icl_ph_val(:)
    REAL   (KIND=jprb), POINTER  :: fmv_aer_ph_val(:)
    REAL   (KIND=jprb), POINTER  :: fmv_wcl_ph_val(:)
    REAL   (KIND=jprb), POINTER  :: fmv_icl_ph_val(:)
    REAL   (KIND=jprb), POINTER  :: fmv_aer_ph_val_cos(:)
    REAL   (KIND=jprb), POINTER  :: fmv_wcl_ph_val_cos(:)
    REAL   (KIND=jprb), POINTER  :: fmv_icl_ph_val_cos(:)
    REAL   (KIND=jprb)           :: fmv_aer_ph_val_min
    REAL   (KIND=jprb)           :: fmv_wcl_ph_val_min
    REAL   (KIND=jprb)           :: fmv_icl_ph_val_min
    REAL   (KIND=jprb), POINTER  :: fmv_icl_dg    (:,:)
    INTEGER(KIND=jpim), POINTER  :: aer_pha_chanlist(:)
    INTEGER(KIND=jpim), POINTER  :: wcl_pha_chanlist(:)
    INTEGER(KIND=jpim), POINTER  :: icl_pha_chanlist(:)
    INTEGER(KIND=jpim), POINTER  :: aer_pha_index(:)
    INTEGER(KIND=jpim), POINTER  :: wcl_pha_index(:)
    INTEGER(KIND=jpim), POINTER  :: icl_pha_index(:)
    REAL(KIND=jprb)   , POINTER  :: abs           (:,:)
    REAL(KIND=jprb)   , POINTER  :: sca           (:,:)
    REAL(KIND=jprb)   , POINTER  :: bpr           (:,:)
    REAL(KIND=jprb)   , POINTER  :: pha           (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: confac        (:)
  END TYPE rttov_coef_scatt_ir

  !> @internal Baran scheme variables
  TYPE rttov_coef_optpiclb
    ! interpolation factors for frequency
    INTEGER(KIND=jpim), POINTER  :: iwn2013(:),iwn2014(:) ! channel
    INTEGER(KIND=jpim), POINTER  :: jwn2013(:),jwn2014(:)
    REAL(KIND=jprb)   , POINTER  :: dx_dwn2013(:),dx_dwn2014(:)
  END TYPE rttov_coef_optpiclb

  !> @internal IR scattering optical properties
  TYPE rttov_optpar_ir
    TYPE(rttov_coef_scatt_ir), POINTER :: optpaer(:)
    TYPE(rttov_coef_scatt_ir), POINTER :: optpwcl(:)
    TYPE(rttov_coef_scatt_ir), POINTER :: optpicl(:)
    TYPE(rttov_coef_optpiclb), POINTER :: optpiclb
  END TYPE rttov_optpar_ir

  !> @internal IR cloud stream variables
  TYPE ircld_type
    INTEGER(KIND=jpim), POINTER  :: nstream(:)
    INTEGER(KIND=jpim), POINTER  :: nstreamref(:)
    INTEGER(KIND=jpim), POINTER  :: iloop(:)
    INTEGER(KIND=jpim), POINTER  :: icount(:)
    INTEGER(KIND=jpim), POINTER  :: icounstr(:)
    INTEGER(KIND=jpim), POINTER  :: icount1(:)
    REAL(KIND=jprb)   , POINTER  :: xstrclr(:)
    INTEGER(KIND=jpim), POINTER  :: icldarr   (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrref1  (:,:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrref2  (:,:,:)
    INTEGER(KIND=jpim), POINTER  :: cldtyp    (:,:,:)
    INTEGER(KIND=jpim), POINTER  :: indexstr  (:,:)
    INTEGER(KIND=jpim), POINTER  :: icount1ref(:,:)
    INTEGER(KIND=jpim), POINTER  :: iloopin   (:,:)
    INTEGER(KIND=jpim), POINTER  :: iflag     (:,:)
    REAL(KIND=jprb)   , POINTER  :: xstr      (:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrminref(:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrref   (:,:)
    REAL(KIND=jprb)   , POINTER  :: cldcfr    (:,:)
    REAL(KIND=jprb)   , POINTER  :: maxcov    (:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrmax   (:,:)
    REAL(KIND=jprb)   , POINTER  :: xstrmin   (:,:)
    REAL(KIND=jprb)   , POINTER  :: a         (:,:)
    REAL(KIND=jprb)   , POINTER  :: ntotref   (:,:)

    REAL(KIND=jprb)   , POINTER  :: tave      (:,:)
    REAL(KIND=jprb)   , POINTER  :: wmixave   (:,:)
    REAL(KIND=jprb)   , POINTER  :: xpresave  (:,:)
    REAL(KIND=jprb)   , POINTER  :: ppv       (:,:)
    REAL(KIND=jprb)   , POINTER  :: esw       (:,:)
    REAL(KIND=jprb)   , POINTER  :: esi       (:,:)
    LOGICAL(KIND=jplm), POINTER  :: flag      (:,:)
  END TYPE ircld_type

  !> @internal RTTOV coefs
  TYPE rttov_coefs
    LOGICAL(KIND=jplm)         :: initialised = .FALSE.
    TYPE(rttov_coef)           :: coef
    TYPE(rttov_coef_scatt_ir)  :: coef_scatt_ir
    TYPE(rttov_optpar_ir)      :: optp
    TYPE(rttov_coef_pccomp)    :: coef_pccomp
  END TYPE

  !> @internal RTTOV internal state
  TYPE rttov_traj
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! and their dimensions are known before running RTTOV (nlevels, nprofiles, nchannels)
! it is possible to allocate these variables from outside RTTOV
!
    TYPE(profile_type), POINTER :: profiles_coef(:)
    TYPE(profile_type), POINTER :: profiles_dry(:)
    TYPE(predictors_type)       :: predictors
    TYPE(raytracing_type)       :: raytracing
    TYPE(raytracing_type)       :: raytracing_coef
    TYPE(ircld_type)            :: ircld
    TYPE(opdp_path_type)        :: opdp_path
    TYPE(opdp_path_type)        :: opdp_path_coef

    REAL(KIND=jprb), POINTER :: thermrefl(:)    ! Surface refl for thermal calcs (nchanprof)
    REAL(KIND=jprb), POINTER :: fresnrefl(:)    ! Fresnel reflection coefficients (nchanprof)
    TYPE(profile_aux)  :: aux_prof
    TYPE(profile_aux)  :: aux_prof_coef
    TYPE(transmission_scatt_ir_type)  :: transmission_scatt_ir
    TYPE(sunglint_type):: sunglint

    TYPE(rttov_coefs),         POINTER :: coefs
    INTEGER(KIND=jpim)   :: nchannels
    INTEGER(KIND=jpim)   :: nlevels
    INTEGER(KIND=jpim)   :: nlayers
    TYPE(rttov_options)  :: opts
  END TYPE

  !> @internal RTTOV internal state
  TYPE rttov_traj_dyn
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! but their dimensions are known when RTTOV starts running (nstreams)
!
    INTEGER(KIND=jpim)               :: nstreams = -1  ! This initialisation used to determine alloc status
    TYPE(radiance_aux              ) :: auxrad_stream
    TYPE(transmission_scatt_ir_type) :: transmission_scatt_ir_stream
    TYPE(transmission_type_aux     ) :: transmission_aux
  END TYPE

  !> @internal RTTOV internal state (optical depths, transmittances)
  TYPE rttov_path_traj_trans
    ! Structure to hold optical depth and transmittance data
    ! within static trajectory.
    REAL(KIND=jprb),     POINTER :: tau_ref       (:,:)
    REAL(KIND=jprb),     POINTER :: tau_ref_surf  (:)
    REAL(KIND=jprb),     POINTER :: tau_level     (:,:)! sat to level transmittance
    REAL(KIND=jprb),     POINTER :: tau_surf      (:)
    REAL(KIND=jprb),     POINTER :: od_level      (:,:)! sat to level optical depth
    REAL(KIND=jprb),     POINTER :: opdp_ref_coef (:,:)! layer optical depth before threshold

    REAL(KIND=jprb),     POINTER :: od_singlelayer(:,:)! single layer optical depth
    REAL(KIND=jprb),     POINTER :: od_frac       (:)
  END TYPE

  !> @internal RTTOV internal state (direct model only)
  TYPE rttov_traj_sta
!
! Hold RTTOV trajectory; these variables do not have counterparts in TL, AD, K
!
    LOGICAL(KIND=jplm),  POINTER         :: thermal(:)        ! switch for thermal calculations (nchanprof)
    LOGICAL(KIND=jplm),  POINTER         :: solar(:)          ! switch for solar calculations (nchanprof)
    LOGICAL(KIND=jplm)                   :: dothermal         ! flag to indicate thermal calculations required
    LOGICAL(KIND=jplm)                   :: dosolar           ! flag to indicate solar calculations required

    LOGICAL(KIND=jplm), POINTER          :: do_lambertian(:)
    REAL(KIND=jprb), POINTER             :: solar_spec_esd(:) ! Solar spectrum adjusted for esd (nchanprof)
    REAL(KIND=jprb), POINTER             :: refl_norm(:)      ! Normalisation factor for solar surface reflectance
    TYPE(rttov_path_traj_trans), POINTER :: thermal_path1
    TYPE(rttov_path_traj_trans), POINTER :: solar_path2
    TYPE(rttov_path_traj_trans), POINTER :: solar_path1
    TYPE(geometry_type), POINTER         :: angles(:)         ! geometry angles
    TYPE(geometry_type), POINTER         :: angles_coef(:)    ! geometry angles
    TYPE(profile_type),  POINTER         :: profiles_coef_ref(:)
    TYPE(radiance_aux)                   :: auxrad
    TYPE(rttov_chanprof), POINTER        :: chanprof_in(:)
    TYPE(rttov_chanprof), POINTER        :: chanprof_pc(:)
  END TYPE

  !> @internal Used for coef testing
  TYPE rttov_lbl_check
    REAL(KIND=jprb), POINTER :: atm_layer(:,:)
    REAL(KIND=jprb), POINTER :: atm_layer_path2(:,:)
    LOGICAL(KIND=jplm) :: plane_geometry
  END TYPE

END MODULE rttov_types
