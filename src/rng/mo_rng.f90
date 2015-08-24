! *****************************COPYRIGHT****************************
! (c) British Crown Copyright 2009, the Met Office.
! All rights reserved.
! $Revision: 23 $, $Date: 2011-03-31 07:41:37 -0600 (Thu, 31 Mar 2011) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/branches/dustinswales/icarus-scops-4.1-bsd/congvec.f $
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************


MODULE mod_rng
  USE cosp_kinds, ONLY: dp, sp, wp 
  IMPLICIT NONE 
  
  INTEGER, parameter :: huge32 = 2147483647
  INTEGER, parameter :: i2_16  = 65536
  INTEGER, parameter :: ki9    = selected_int_kind(R=9)
  
  TYPE rng_state
     INTEGER(ki9) :: seed ! 32 bit integer
  END TYPE rng_state
  
  INTERFACE init_rng
     MODULE PROCEDURE init_rng_1, init_rng_n
  END INTERFACE init_rng
  
  INTERFACE get_rng
     MODULE PROCEDURE get_rng_1, get_rng_n, get_rng_v
  END INTERFACE get_rng
  
CONTAINS 
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Set single seed
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE init_rng_1(s, seed_in)
    TYPE(rng_state), INTENT(INOUT) :: s
    INTEGER,         INTENT(IN   ) :: seed_in
    s%seed = seed_in     
  END SUBROUTINE init_rng_1
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Set vector of seeds
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE init_rng_n(s, seed_in)
    TYPE(rng_state), DIMENSION(:), INTENT(INOUT) :: s
    INTEGER,         DIMENSION(:), INTENT(IN   ) :: seed_in
    
    INTEGER :: i 
    DO i = 1, SIZE(seed_in)
       s(i)%seed = seed_in(i) 
    END DO
  END SUBROUTINE init_rng_n
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Create single random number from seed
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION get_rng_1(s)  
    TYPE(rng_state), INTENT(INOUT) :: s
    REAL(WP)                       :: get_rng_1
    
    REAL(SP)                       :: r
    ! Return the next random numbers 
    
    ! Marsaglia CONG algorithm
    s%seed=69069*s%seed+1234567
    ! mod 32 bit overflow
    s%seed=mod(s%seed,2_ki9**30_ki9)   
    r = s%seed*0.931322574615479E-09
    
    ! convert to range 0-1 (32 bit only)
    if ( i2_16*i2_16 .le. huge32 ) then
       r=r+1
       r=r-int(r)
    endif
    get_rng_1 = REAL(r, KIND = WP) 
    
  END FUNCTION get_rng_1
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Create single random number N times
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION get_rng_n(s, n) RESULT (r) 
    integer,intent(inout) :: n
    TYPE(rng_state),INTENT(INOUT) :: s
    ! Return the next N random numbers 
    REAL(WP), DIMENSION (n)        :: r
    
    INTEGER :: i 
    
    DO i = 1, N
       r(i) = get_rng_1(s)
    END DO
  END FUNCTION get_rng_n
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Create a vector of random numbers from a vector of input seeds
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION get_rng_v(s) RESULT (r) 
    ! Return the next N random numbers 
    TYPE(rng_state), DIMENSION(:), INTENT(INOUT) :: s
    REAL(WP),        DIMENSION (SIZE(S))         :: r
    
    INTEGER :: i
    
    DO i = 1, size(s)
       r(i) = get_rng_1(s(i))
    END DO
  END FUNCTION get_rng_v
  
END MODULE mod_rng
