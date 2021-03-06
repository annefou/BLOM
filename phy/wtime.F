! ------------------------------------------------------------------------------
! Copyright (C) 2005 HYCOM Consortium and contributors
! Copyright (C) 2006 Mats Bentsen
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

#if defined(AIX)
      REAL*8 FUNCTION WTIME()
      IMPLICIT NONE
C
C     USE THE FUNCTION  RTC  TO RETURN WALL TIME.
C
      REAL*8 RTC
      
      WTIME = RTC()
      RETURN
C     END OF WTIME.
      END
#elif defined(MPI)
      REAL*8 FUNCTION WTIME()
      IMPLICIT NONE
C
C     USE THE MPI FUNCTION  MPI_WTIME  TO RETURN WALL TIME.
C
      DOUBLE PRECISION MPI_WTIME
C
      WTIME = MPI_WTIME()
      RETURN
C     END OF WTIME.
      END
#else
      REAL*8 FUNCTION WTIME()
      IMPLICIT NONE
C
C     USE THE F90 INTRINSIC  SYSTEM_CLOCK  TO RETURN WALL TIME.
C
C     WILL FAIL IF THE COUNT IS EVER NEGATIVE, BUT THE STANDARD
C     SAYS THAT IT IS AWAYS NON-NEGATIVE IF A CLOCK EXISTS.
C     NOT THREAD-SAFE, UNLESS LCOUNT AND IOVER ARE THREADPRIVATE.
C
      REAL*8     ZERO,ONE
      PARAMETER (ZERO=0.0, ONE=1.0)
C
      INTEGER COUNT, MCOUNT, RATE
C
      REAL*8  OFFSEC, OFFSET, PERSEC
      INTEGER ICOUNT, IOVER,  LCOUNT, NCOUNT
      SAVE    OFFSEC, OFFSET, PERSEC
      SAVE    ICOUNT, IOVER,  LCOUNT, NCOUNT
C
      DATA IOVER, LCOUNT / -1, -1 /
C
      CALL SYSTEM_CLOCK(COUNT)
C
      IF     (COUNT.LT.LCOUNT) THEN
C
C        COUNT IS SUPPOSED TO BE NON-DECREASING EXCEPT WHEN IT WRAPS,
C        BUT SOME IMPLEMENTATIONS DON''T DO THIS.  SO IGNORE ANY
C        DECREASE OF LESS THAN ONE PERCENT OF THE RANGE.
C
         IF     (LCOUNT-COUNT.LT.NCOUNT) THEN
           COUNT  = LCOUNT
         ELSE
           IOVER  = IOVER + 1
           OFFSET = OFFSET + OFFSEC
         ENDIF
      ENDIF
      LCOUNT = COUNT
C
      IF     (IOVER.EQ.0) THEN
C
C       FIRST CYCLE, FOR ACCURACY WITH 64-BIT COUNTS.
C
         WTIME = (COUNT - ICOUNT) * PERSEC
      ELSEIF (IOVER.GT.0) THEN
C
C        ALL OTHER CYCLES.
C
         WTIME = COUNT * PERSEC + OFFSET
      ELSE
C
C        INITIALIZATION.
C
         CALL SYSTEM_CLOCK(ICOUNT, RATE, MCOUNT)
         NCOUNT =  MCOUNT/100
         PERSEC =  ONE/RATE
         OFFSEC =  MCOUNT * PERSEC
         OFFSET = -ICOUNT * PERSEC
         IOVER  =  0
         WTIME  =  ZERO
      ENDIF
      RETURN
C     END OF WTIME.
      END
#endif  /* MPI:else */
