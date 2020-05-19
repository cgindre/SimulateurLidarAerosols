! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
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
! History:
! May 2020 - R. Guzman        - Introducing new linear interpolation routine to perform
!                               interpolations for the lidar and radar simulators. Used
!                               along the full vertical extent when regridding low
!                               vertical variability fields compared to the output grid,
!                               and mixing it with the legacy vertical regridding
!                               routine COSP_CHANGE_VERTICAL_GRID in order to avoid
!                               striping features to appear in mid-high atmopsheric
!                               levels of high variability fields (ATB, Reflectivities). 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_INTERP
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF

  IMPLICIT NONE
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------------- SUBROUTINE INTERP_LINEAR --------------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE INTERP_LINEAR(m, data_num, t_data, p_data, interp_num, &
                           t_interp, p_interp, lunits)

  implicit none

  integer :: data_num
  integer :: m
  integer :: interp_num

  integer :: interp
  integer :: left
  real :: p_data(m,data_num)
  real :: p_interp(m,interp_num)
  logical :: ascends_ok
  integer :: right
  logical :: lunits
  integer :: i
  real :: t
  real :: t_data(data_num)
  real :: t_interp(interp_num)

  call vec_ascends_strictly ( data_num, t_data, ascends_ok )

  if ( .not. ascends_ok ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
    write ( *, '(a)' ) &
      '  Independent variable array T_DATA is not strictly increasing.'
    stop 1
  end if

  do interp = 1, interp_num

    t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
    call vec_bracket (data_num, t_data, t, left, right)

    if (lunits) then
      do i = 1, m
        if ( (p_data(i,left) /= R_UNDEF) .and. &
             (p_data(i,right) /= R_UNDEF) ) then
      p_interp(i,interp) = &
        ( ( t_data(right) - t                ) * 10._wp**(p_data(i,left)/10._wp)   &
        + (                 t - t_data(left) ) * 10._wp**(p_data(i,right)/10._wp) ) &
        / ( t_data(right)     - t_data(left) )
        else
      p_interp(i,interp) = 0._wp
        endif
      enddo
    else
      p_interp(1:m,interp) = &
        ( ( t_data(right) - t                ) * p_data(1:m,left)   &
        + (                 t - t_data(left) ) * p_data(1:m,right) ) &
        / ( t_data(right)     - t_data(left) )
    endif

  end do

  return
END SUBROUTINE INTERP_LINEAR


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------------- SUBROUTINE VEC_BRACKET ----------------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE VEC_BRACKET(n, x, xval, left, right)

  implicit none

  integer :: n

  integer :: i
  integer :: left
  integer :: right
  real :: x(n)
  real :: xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
END SUBROUTINE VEC_BRACKET


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !--------------- SUBROUTINE VEC_ASCENDS_STRICTLY -----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE VEC_ASCENDS_STRICTLY(n, x, ascends_ok)

  implicit none

  integer :: n

  integer :: i
  logical :: ascends_ok
  real :: x(n)

  do i = 1, n - 1
    if ( x(i+1) <= x(i) ) then
      ascends_ok = .false.
      return
    end if
  end do

  ascends_ok = .true.

  return
END SUBROUTINE VEC_ASCENDS_STRICTLY

END MODULE MOD_INTERP
