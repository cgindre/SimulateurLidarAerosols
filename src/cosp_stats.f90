MODULE MOD_COSP_STATS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF,R_GROUND
  IMPLICIT NONE
CONTAINS
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,M,zl,zu,r,log_units)
    ! Input arguments
    integer,intent(in) :: Npoints  !# of grid points
    integer,intent(in) :: Nlevels  !# of levels
    integer,intent(in) :: Ncolumns !# of columns
    real(wp),dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
    real(wp),dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
    integer,intent(in) :: M  !# levels in the new grid
    real(wp),dimension(M),intent(in) :: zl ! Lower boundary of new levels  [m]
    real(wp),dimension(M),intent(in) :: zu ! Upper boundary of new levels  [m]
    logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
    ! Output
    real(wp),dimension(Npoints,Ncolumns,M),intent(out) :: r ! Variable on new grid
    
    ! Local variables
    integer :: i,j,k
    logical :: lunits
    real(wp) :: ws
    real(wp),dimension(Nlevels) :: xl,xu ! Lower and upper boundaries of model grid
    real(wp),dimension(M) :: dz          ! Layer depth
    real(wp),dimension(Nlevels,M) :: w   ! Weights to do the mean at each point
    real(wp),dimension(Ncolumns,Nlevels) :: yp  ! Variable to be changed to a different grid.
                                            ! Local copy at a particular point.
                                            ! This allows for change of units.
   
    lunits=.false.
    if (present(log_units)) lunits=log_units
    
    r = R_UNDEF
    do i=1,Npoints
       ! Vertical grid at that point
       xl = zhalf(i,:)
       xu(1) = zfull(i,1)+zfull(i,1)-zhalf(i,1)
       xu(2:Nlevels) = xl(1:nlevels-1)
       dz = zu - zl
       yp = y(i,:,:) ! Temporary variable to regrid
       
       
       ! Find weights
       w = 0.0
       do k=1,M
          do j=1,Nlevels
             if ((xl(j) < zl(k)).and.(xu(j) > zl(k)).and.(xu(j) <= zu(k))) then
                !xl(j)-----------------xu(j)
                !      zl(k)------------------------------zu(k)
                w(j,k) = xu(j) - zl(k)
             else if ((xl(j) >= zl(k)).and.(xu(j) <= zu(k))) then
                !           xl(j)-----------------xu(j)
                !      zl(k)------------------------------zu(k)
                w(j,k) = xu(j) - xl(j)
             else if ((xl(j) >= zl(k)).and.(xl(j) < zu(k)).and.(xu(j) >= zu(k))) then
                !                           xl(j)-----------------xu(j)
                !      zl(k)------------------------------zu(k)
                w(j,k) = zu(k) - xl(j)
             else if ((xl(j) <= zl(k)).and.(xu(j) >= zu(k))) then
                !  xl(j)---------------------------xu(j)
                !        zl(k)--------------zu(k)
                w(j,k) = dz(j)
             endif
          enddo
       enddo
       ! Check for dBZ and change if necessary
       if (lunits) then
          where (yp /= R_UNDEF)
             yp = 10.0**(yp/10.0)
          elsewhere
             yp = 0.0
          end where
       endif
       ! Do the weighted mean
       do j=1,Ncolumns
          do k=1,M
             ws = sum(w(:,k))
             if (ws > 0.0) r(i,j,k) = sum(w(:,k)*yp(j,:))/ws
          enddo
       enddo
       ! Check for dBZ and change if necessary
       if (lunits) then
          where (r(i,:,:) <= 0.0)
             r(i,:,:) = R_UNDEF
          elsewhere
             r(i,:,:) = 10.0*log10(r(i,:,:))
          end where
       endif
    enddo
  END SUBROUTINE COSP_CHANGE_VERTICAL_GRID
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
  ! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
  ! All rights reserved.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_LIDAR_ONLY_CLOUD(Npoints,Ncolumns,Nlevels,temp_tot,beta_tot, &
                                   betaperp_tot,beta_mol,Ze_tot, &
                                   lidar_only_freq_cloud,tcc)
    ! Inputs
    integer,intent(in) :: &
         Npoints,       & ! Number of horizontal gridpoints
         Ncolumns,      & ! Number of subcolumns
         Nlevels          ! Number of vertical layers
    real(wp),dimension(Npoints,Nlevels),intent(in) :: &
         beta_mol         ! Molecular backscatter
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         beta_tot,      & ! Total backscattered signal
         temp_tot,      & ! Total backscattered signal
         betaperp_tot,  & ! perpendicular Total backscattered signal
         Ze_tot           ! Radar reflectivity
    ! Outputs
    real(wp),dimension(Npoints,Nlevels),intent(out) :: &
         lidar_only_freq_cloud
    real(wp),dimension(Npoints),intent(out) ::&
         tcc
    
    ! local variables
    real(wp) :: sc_ratio
    real(wp),parameter :: &
         s_cld=5.0, &
         s_att=0.01
    integer :: flag_sat,flag_cld,pr,i,j
    
    lidar_only_freq_cloud = 0._wp
    tcc = 0._wp
    do pr=1,Npoints
       do i=1,Ncolumns
          flag_sat = 0
          flag_cld = 0
          do j=1,Nlevels
             sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
             if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
             if (Ze_tot(pr,i,j) .lt. -30.) then  !radar can't detect cloud
                if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
                   lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1. !top->surf
                   flag_cld=1
                endif
             else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
                flag_cld=1
             endif
          enddo !levels
          if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1._wp
       enddo !columns
    enddo !points
    lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
    tcc=tcc/Ncolumns
    
    ! Unit conversion
    where(lidar_only_freq_cloud /= R_UNDEF) &
            lidar_only_freq_cloud = lidar_only_freq_cloud*100._wp
    where(tcc /= R_UNDEF) tcc = tcc*100._wp
    
  END SUBROUTINE COSP_LIDAR_ONLY_CLOUD
  
  ! ######################################################################################
  ! FUNCTION hist1D
  ! ######################################################################################
  function hist1d(Npoints,var,nbins,bins)
    ! Inputs
    integer,intent(in) :: &
         Npoints, & ! Number of points in input array
         Nbins      ! Number of bins for sorting
    real(wp),intent(in),dimension(Npoints) :: &
         var        ! Input variable to be sorted
    real(wp),intent(in),dimension(Nbins+1) :: &
         bins       ! Histogram bins [lowest,binTops]  
    ! Outputs
    real(wp),dimension(Nbins) :: &
         hist1d     ! Output histogram      
    ! Local variables
    integer :: ij
    
    do ij=2,Nbins+1  
       hist1D(ij-1) = count(var .ge. bins(ij-1) .and. var .lt. bins(ij))
    enddo
    
  end function hist1D
  
  ! ######################################################################################
  ! SUBROUTINE hist2D
  ! ######################################################################################
  subroutine hist2D(var1,var2,npts,bin1,nbin1,bin2,nbin2,jointHist)
    implicit none
    
    ! INPUTS
    integer, intent(in) :: &
         npts,  & ! Number of data points to be sorted
         nbin1, & ! Number of bins in histogram direction 1 
         nbin2    ! Number of bins in histogram direction 2
    real(wp),intent(in),dimension(npts) :: &
         var1,  & ! Variable 1 to be sorted into bins
         var2     ! variable 2 to be sorted into bins
    real(wp),intent(in),dimension(nbin1+1) :: &
         bin1     ! Histogram bin 1 boundaries
    real(wp),intent(in),dimension(nbin2+1) :: &
         bin2     ! Histogram bin 2 boundaries
    ! OUTPUTS
    real(wp),intent(out),dimension(nbin1,nbin2) :: &
         jointHist
    
    ! LOCAL VARIABLES
    integer :: ij,ik
    
    do ij=2,nbin1+1
       do ik=2,nbin2+1
          jointHist(ij-1,ik-1)=count(var1 .ge. bin1(ij-1) .and. var1 .lt. bin1(ij) .and. &
               var2 .ge. bin2(ik-1) .and. var2 .lt. bin2(ik))        
       enddo
    enddo
  end subroutine hist2D
END MODULE MOD_COSP_STATS
