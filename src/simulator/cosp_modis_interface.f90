MODULE MOD_COSP_Modis_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF
  use mod_modis_sim,   ONLY: num_trial_res,min_OpticalThickness,CO2Slicing_PressureLimit,&
                             CO2Slicing_TauLimit,phase_TauLimit,size_TauLimit,re_fill,   &
                             get_ssa_nir,phaseDiscrimination_Threshold,re_water_min,     &
                             re_water_max,get_g_nir,re_ice_min,re_ice_max,               &
                             highCloudPressureLimit,lowCloudPressureLimit,phaseIsNone,   &
                             phaseIsLiquid,phaseIsIce,phaseIsUndetermined,trial_re_w,    &
                             trial_re_i,g_w,g_i,w0_w,w0_i
  implicit none

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  TYPE modis_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type modis_IN
     integer,pointer :: &
          Npoints,        & ! Number of horizontal gridpoints
          Ncolumns,       & ! Number of subcolumns
          Nlevels           ! Number of vertical levels
     integer :: &
          Nsunlit           ! Number of sunlit lit pixels
     real(wp),allocatable,dimension(:) :: &
          sunlit,         & ! Sunlit scenes
          notSunlit         ! Dark scenes
     real(wp),allocatable,dimension(:,:) :: &
          pres              ! Gridmean pressure at layer edges (Pa) 
     real(wp),pointer ::  &
          tau(:,:,:),     & ! Subcolumn optical thickness @ 0.67 microns.
          liqFrac(:,:,:), & ! Liquid water fraction
          g(:,:,:),       & ! Subcolumn assymetry parameter  
          w0(:,:,:)         ! Subcolumn single-scattering albedo
  end type modis_IN
contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROTUINE cosp_modis_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_MODIS_INIT()
    integer :: k, l,i 
    
    ! Retrieval parameters
    min_OpticalThickness          = 0.3_wp     ! Minimum detectable optical thickness
    CO2Slicing_PressureLimit      = 70000._wp  ! Cloud with higher pressures use thermal 
                                               ! methods, units Pa
    CO2Slicing_TauLimit           = 1._wp      ! How deep into the cloud does CO2 slicing 
                                               ! see? 
    phase_TauLimit                = 1._wp      ! How deep into the cloud does the phase 
                                               ! detection see?
    size_TauLimit                 = 2._wp      ! Depth of the re retreivals
    phaseDiscrimination_Threshold = 0.7_wp     ! What fraction of total extincton needs to 
                                               ! be in a single category to make phase 
                                               ! discrim. work? 
    re_fill                       = -999._wp   ! Fill value
    re_water_min                  = 4._wp      ! Minimum effective radius (liquid)
    re_water_max                  = 30._wp     ! Maximum effective radius (liquid)
    re_ice_min                    = 5._wp      ! Minimum effective radius (ice)
    re_ice_max                    = 90._wp     ! Minimum effective radius (ice)
    highCloudPressureLimit        = 44000._wp  ! High cloud pressure limit (Pa)
    lowCloudPressureLimit         = 68000._wp  ! Low cloud pressure limit (Pa)
    phaseIsNone                   = 0          ! No cloud
    phaseIsLiquid                 = 1          ! Liquid cloud
    phaseIsIce                    = 2          ! Ice cloud
    phaseIsUndetermined           = 3          ! Undetermined cloud

    ! Precompute near-IR optical params vs size for retrieval scheme    
    trial_re_w(1:num_trial_res) = re_water_min + (re_water_max - re_water_min) /         &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
    trial_re_i(1:num_trial_res) = re_ice_min   + (re_ice_max -   re_ice_min) /           &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)

    ! Initialize estimates for size retrieval
    g_w(1:num_trial_res)  = get_g_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    w0_w(1:num_trial_res) = get_ssa_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    g_i(1:num_trial_res)  = get_g_nir(phaseIsIce,trial_re_i(1:num_trial_res))
    w0_i(1:num_trial_res) = get_ssa_nir(phaseIsIce,trial_re_i(1:num_trial_res))

  END SUBROUTINE COSP_MODIS_INIT
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE MOD_COSP_Modis_INTERFACE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_Modis_INTERFACE
