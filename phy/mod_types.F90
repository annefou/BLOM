! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_types
! ------------------------------------------------------------------------------
! This module defines numeric data types.
! ------------------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env, only: &
      int8, int16, int32, int64, real32, real64

   implicit none

   private

   integer, parameter :: &
      i1 = int8, &
      i2 = int16, &
      i4 = int32, &
      i8 = int64, &
      r4 = real32, &
      r8 = real64

   public :: i1, i2, i4, i8, r4, r8

end module mod_types
