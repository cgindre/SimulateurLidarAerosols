! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Centre National de la Recherche Scientifique
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2007: ActSim code of M. Chiriaco and H. Chepfer rewritten by S. Bony
!
! May 2008, H. Chepfer:
! - Units of pressure inputs: Pa 
! - Non Spherical particles : LS Ice NS coefficients, CONV Ice NS coefficients
! - New input: ice_type (0=ice-spheres ; 1=ice-non-spherical)
!
! June 2008, A. Bodas-Salcedo:
! - Ported to Fortran 90 and optimisation changes
!
! August 2008, J-L Dufresne:
! - Optimisation changes (sum instructions suppressed)
!
! October 2008, S. Bony,  H. Chepfer and J-L. Dufresne :  
! - Interface with COSP v2.0:
!      cloud fraction removed from inputs
!      in-cloud condensed water now in input (instead of grid-averaged value)
!      depolarisation diagnostic removed
!      parasol (polder) reflectances (for 5 different solar zenith angles) added
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - Modification of the integration of the lidar equation.
! - change the cloud detection threshold
!
! April 2008, A. Bodas-Salcedo:
! - Bug fix in computation of pmol and pnorm of upper layer
!
! April 2008, J-L. Dufresne
! - Bug fix in computation of pmol and pnorm, thanks to Masaki Satoh: a factor 2 
! was missing. This affects the ATB values but not the cloud fraction. 
!
! May 2015 - D. Swales - Modified for COSPv2.0
!
!---------------------------------------------------------------------------------
!
! Purpose   : Compute lidar signal from model-simulated profiles of
!             extinction and backscatter in each model gridbox.
! Author    : Po-Lun Ma (Po-Lun.Ma@pnnl.gov), Helene Chepfer, David M. Winker
! Reference : Ma et al., 2018
! History   : March 2015, Initial code development from COSPv1.4
!
! August 2019 : R. Guzman adapted the code for COSPv2
!
! May 2020 : R. Guzman introduced linear interpolation routine COSP_INTERP_NEW_GRID
!            to get rid of the striping features appearing when model vertical
!            resolution is significantly coarser than the output grid vertical
!            resolution (issue #34 GitHub).
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_lidar_aerosols_simulator
  USE COSP_KINDS,         ONLY: wp
  USE MOD_COSP_CONFIG,    ONLY: R_UNDEF, atbmin, c0, det1, &
                                use_vgrid_aerosols, vgrid_zl_aerosols, &
                                vgrid_zu_aerosols, vgrid_z_aerosols, LIDAR_AEROSOLS_FLAGS
  USE MOD_COSP_STATS,     ONLY: COSP_INTERP_NEW_GRID, hist1d
  USE MOD_LIDAR_SIMULATOR,ONLY: cmp_backsignal
  implicit none

contains
  ! ######################################################################################
  ! SUBROUTINE lidar_aerosols_column
  ! Inputs with a vertical dimensions (nlevels) should ordered in along the vertical 
  ! dimension from TOA-2-SFC, for example: varIN(nlevels) is varIN @ SFC. 
  ! ######################################################################################
  subroutine lidar_aerosols_column(npoints, nlevels, llm, llm_aerosols, alpha_aer, &
                            beta_mol, tau_mol, betatot, tautot, &
                            pplay, zlev, zlev_half, vgrid_z_aerosols, cloud_fraction, &
                            pnorm_aero, ext_aero, sr_aero)

    ! INPUTS
    INTEGER,intent(in) :: & 
         npoints,      & ! Number of horizontal gridpoints
         nlevels,      & ! Number of vertical layers (model grid)
         llm,          & ! Number of vertical layers (clouds output grid)
         llm_aerosols    ! Number of vertical layers (aerosols output grid)
    REAL(wp),intent(in),dimension(npoints,nlevels) :: &
         beta_mol,     & ! Molecular backscatter coefficient
         tau_mol,      & ! Molecular optical depth
         betatot,      & ! Molecular + aerosols backscatter coefficient
         tautot,       & ! Molecular + aerosols Optical thickess integrated from top
         alpha_aer,    & ! Aerosol extinction
         pplay,        & ! Pressure on model levels (Pa)
         zlev            ! Model full levels
    real(wp),intent(in),dimension(npoints,nlevels+1) :: &
         zlev_half       ! Model half levels
    real(wp),intent(in),dimension(llm_aerosols) :: & 
         vgrid_z_aerosols  ! mid-level altitude of the aerosols output vertical grid
    real(wp),intent(in),dimension(npoints,llm) :: &
         cloud_fraction  ! 3D "lidar" cloud fraction

    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,llm_aerosols,LIDAR_AEROSOLS_FLAGS) :: &
         pnorm_aero,  & ! Molecular + aerosols backscatter signal power (m^-1.sr^-1)
         ext_aero,    & ! extinction (m^-1)
         sr_aero        ! Scattering Ratio (no unit)


    ! LOCAL VARIABLES
    INTEGER :: i,j,k,zi,zf,zinc,nb_aero_sublay,k_llm
    integer,dimension(npoints,nlevels)  :: layxx ! temporary arrays
    REAL(wp),dimension(npoints,nlevels) :: &
         pmol,         & ! Molecular attenuated backscatter lidar signal power(m^-1.sr^-1)
         pnorm           ! Molecular + aerosols backscatter signal power (m^-1.sr^-1)
    real(wp),dimension(npoints,1,nlevels)      :: ph_in, beta_mol_in, pmol_in, &
                                                  pnorm_in, zlev_in, alpha_in
    real(wp),dimension(npoints,1,llm_aerosols) :: pplayFlip, beta_molFlip, &
                                                  pmolFlip, pnormFlip, zlevFlip, &
                                                  alphaFlip 
    real(wp),dimension(npoints,llm_aerosols)   :: pmolm, sr, pnorm_c, cld_aerosols
    real(wp) :: cloud_fraction_threshold = 99.

    ! COSPv2 convention for spaceborne lidar (default)
       zi   = 1
       zf   = nlevels
       zinc = 1
    ! Number of aerosols sub layers per cloud layer
       nb_aero_sublay = llm_aerosols/llm

!------------------------------------------------------------
!---- 1. Initialization
!------------------------------------------------------------

    pnorm_aero(:,:,:) = R_UNDEF
!    dp_aero1(:,:)    = R_UNDEF
    ext_aero(:,:,:)   = R_UNDEF
    sr_aero(:,:,:)    = R_UNDEF
    sr(:,:)          = R_UNDEF

    ! ####################################################################################
    ! *) Molecular signal
    ! ####################################################################################
    call cmp_backsignal(nlevels,npoints,beta_mol(1:npoints,zi:zf:zinc),&
                        tau_mol(1:npoints,zi:zf:zinc),pmol(1:npoints,zi:zf:zinc))
                        

    ! #################################################################################
    ! *) Total Backscatter signal (molecular + aerosols)
    ! #################################################################################
    call cmp_backsignal(nlevels,npoints,betatot(1:npoints,zi:zf:zinc),&
            tautot(1:npoints,zi:zf:zinc),pnorm(1:npoints,zi:zf:zinc))

    ! Vertically regrid input data
    if (use_vgrid_aerosols) then
       !!! All the vertical regridding here can be performed with the linear interpolation
       !!! routine COSP_INTERP_NEW_GRID because the output grid for these aerosols diagnotics
       !!! should always be as fine or finer than the model vertical grid.
       ph_in(:,1,:) = pplay(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            ph_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),pplayFlip(:,1,llm_aerosols:1:-1))

       beta_mol_in(:,1,:) = beta_mol(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            beta_mol_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),beta_molFlip(:,1,llm_aerosols:1:-1))

       pmol_in(:,1,:) = pmol(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            pmol_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),pmolFlip(:,1,llm_aerosols:1:-1))

       pnorm_in(:,1,:) = pnorm(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            pnorm_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),pnormFlip(:,1,llm_aerosols:1:-1))

       zlev_in(:,1,:) = zlev(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            zlev_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),zlevFlip(:,1,llm_aerosols:1:-1))

       alpha_in(:,1,:) = alpha_aer(:,nlevels:1:-1)
       call cosp_interp_new_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            alpha_in,llm_aerosols,vgrid_z_aerosols(llm_aerosols:1:-1),vgrid_zu_aerosols(llm_aerosols:1:-1),alphaFlip(:,1,llm_aerosols:1:-1))

    endif

    ! Compute LIDAR scattering ratio
    if (use_vgrid_aerosols) then
      do k=1,llm_aerosols
      do i=1,npoints
          pnorm_c(i,k) = max(pnormFlip(i,1,k),pmolFlip(i,1,k))
        if (pnorm_c(i,k) .gt. 0. .and. pmolFlip(i,1,k) .gt. 0.) then
          sr(i,k) = pnorm_c(i,k)/pmolFlip(i,1,k)
        endif
      enddo 
      enddo
	
!          pnorm_c = pnormFlip(i,1,k)
!          where ((pnorm_c .lt. xmax) .and. (betamolFlip(:,1,:) .lt. xmax) .and.          &
!                (betamolFlip(:,1,:) .gt. 0.0 ))
!             sr_aero1 = pnorm_c/betamolFlip(:,1,:)
!          elsewhere
!             sr_aero1 = R_UNDEF
!          end where
    else ! use_vgrid_aerosols is always .true., otherwise the code would crash here...
         ! because sr as dimensions (npoints,llm_aerosols) and not (npoints,nlevels)
      do k=1,nlevels
      do i=1,npoints
          pnorm(i,k) = max(pnorm(i,k),pmol(i,k))
        if (pnorm(i,k) .gt. 0. .and. pmol(i,k) .gt. 0.) then
          sr(i,k) = pnorm(i,k)/pmol(i,k)
        endif
      enddo 
      enddo
    endif

!          pnorm_c = pnorm(:,1,:)
!          where ((pnorm_c .lt. xmax) .and. (pmol .lt. xmax) .and. (pmol .gt. 0.0 ))
!             sr_aero1 = pnorm_c/pmol
!          elsewhere
!             sr_aero1 = R_UNDEF
!          end where


 do i=1,npoints
 do k=1,llm_aerosols
   pmolm(i,k)=1.+det1*c0*(sqrt(1.e-7/beta_molFlip(i,1,k))+(1.e-7/pmolFlip(i,1,k)))
   if ( zlevFlip(i,1,k) .gt. 1.2e4     .or. &
        sr(i,k)        .le. pmolm(i,k)  .or. &
        pnorm_c(i,k)   .le. atbmin     ) then
   layxx(i,k)= 0
   else
   layxx(i,k)= 1
   endif
! Transposing cloud fraction profile values (llm) to aerosols vertical grid (llm_aerosols)
   k_llm = 1 + int((k-1)/nb_aero_sublay)
   cld_aerosols(i,k) = cloud_fraction(i,k_llm)
 enddo
 enddo


! Filling in the output variables

! Flag 0 : No masking at all
   DO i = 1,npoints
      DO k=1,llm_aerosols

               ext_aero(i,k,1)   = alphaFlip(i,1,k)
               pnorm_aero(i,k,1) = pnorm_c(i,k)
!               pnorm_aero(i,k,1) = pmolFlip(i,1,k) !RG: testing interpolation routine
               sr_aero(i,k,1)    = sr(i,k)

      ENDDO
   ENDDO

! Flag 1 : Cloud masking only
   DO i = 1,npoints
      DO k=1,llm_aerosols

         if ( cld_aerosols(i,k) .lt. cloud_fraction_threshold .and. cld_aerosols(i,k) .ge. 0.) then 
!               aerotype1(i,k)   = aertype(i,k)
               ext_aero(i,k,2)   = ext_aero(i,k,1)
               pnorm_aero(i,k,2) = pnorm_aero(i,k,1)
               sr_aero(i,k,2)    = sr_aero(i,k,1)
!               dp_aero1(i,k)    = dp_all(i,k)
!               ext_tmp(i,k)     = alpha_aer(i,k)
!
!               aerotype2(i,k)   = aertype(i,k)
!               ext_aero2(i,k)   = ext_aer(i,k)
!               pnorm_aero2(i,k) = pnorm_c(i,k)
!               sr_aero2(i,k)    = sr(i,k)
!               dp_aero2(i,k)    = dp_all(i,k)
!               ext_aero0(i,k)   = alpha_aer(i,k)
         else
           exit
         endif
      ENDDO
   ENDDO

! Flag 2 : Detection Threshold masking only
   DO i = 1, npoints
      DO k=1,llm_aerosols

               ext_aero(i,k,3)   = ext_aero(i,k,1)
               pnorm_aero(i,k,3) = pnorm_aero(i,k,1)
               sr_aero(i,k,3)    = sr_aero(i,k,1)

!         if (aerotype1(i,k) .lt. 10. ) then
         if ( (sr(i,k)    .lt. pmolm(i,k)  .or.  &
              pnorm_c(i,k) .lt. atbmin)    .and. &
              (pnorm_aero(i,k,3) .ge. 0.) ) then 
               ext_aero(i,k,3)   = 0.
               pnorm_aero(i,k,3) = pmolFlip(i,1,k)
               sr_aero(i,k,3)    = 1.
!               aerotype2(i,k)   = 0.
!               ext_aero2(i,k)   = 0.
!               pnorm_aero2(i,k) = pmol(i,k)
!               sr_aero2(i,k)    = 1.
!               dp_aero2(i,k)    = 0.
         endif
!         endif
      ENDDO
   ENDDO

! Flag 3 : Cloud + Detection Threshold masking
    DO i = 1,npoints
      DO k=1,llm_aerosols

         ! Cloud mask applied to Detection Threshold masked fields
         if ( cld_aerosols(i,k) .lt. cloud_fraction_threshold .and. cld_aerosols(i,k) .ge. 0.) then 
               ext_aero(i,k,4)   = ext_aero(i,k,3)
               pnorm_aero(i,k,4) = pnorm_aero(i,k,3)
               sr_aero(i,k,4)    = sr_aero(i,k,3)
         else
           exit
         endif

      ENDDO
   ENDDO


  end subroutine lidar_aerosols_column

  ! ####################################################################################
  ! SUBROUTINE calipso_aerosols_typing
  ! Conventions: 
  ! ####################################################################################
!  subroutine calipso_aerosols_typing(Npoints, Ncolumns, Nlevels, Ncat, Nphase, tmp, x, &
!                          ATB, ATBperp, pplay, S_att, S_cld, S_cld_att, undef, lidarcld, &
!                          cldlayer, lidarcldphase, cldlayerphase, lidarcldtemp)




!  end subroutine calipso_aerosols_typing


  ! ####################################################################################
  ! SUBROUTINE read_calipso_observations
  ! Conventions:
  ! ####################################################################################
!  subroutine read_calipso_observations(Npoints, Ncolumns, Nlevels, Ncat, Nphase, tmp, x, &
!                          ATB, ATBperp, pplay, S_att, S_cld, S_cld_att, undef, lidarcld, &
!                          cldlayer, lidarcldphase, cldlayerphase, lidarcldtemp)




!  end subroutine read_calipso_observations



end module mod_lidar_aerosols_simulator
