MODULE init


CONTAINS
  !--------------------------------------------------------------
  SUBROUTINE initializeRandomNumberGenerator
  !
  ! Markus Miettinen, Tuesday, January 15, 2013
  !--------------------------------------------------------------
    USE datastructures, ONLY : seed
    IMPLICIT NONE

    CALL rmarin(seed,0,0) !-- initialize random number generator

  END SUBROUTINE initializeRandomNumberGenerator


  !--------------------------------------------------------------
  SUBROUTINE randomInitPositions
  !
  !
  !   A subroutine that builds a PBC configuration of 'nDonors'
  !   donors and 'nAcceptors' acceptors of radius 'rMin' randomly
  !   placed in a simulation box of size 'boxLenX' x 'boxLenY'
  !   x 'boxLenZ'.
  !
  !   Markus Miettinen, Tuesday, January 15, 2013
  !--------------------------------------------------------------
    USE datastructures, ONLY : & !-- constants
         nDonors, nAcceptors, &
         boxLenX, boxLenY, boxLenZ 
    USE datastructures, ONLY : &  !-- variables
         xDonors, yDonors, zDonors, &
         xAcceptors, yAcceptors, zAcceptors
    IMPLICIT NONE
  
    REAL, PARAMETER :: rMin2 = 1.0**2 !-- squared minimum distance between particles (nm2)

    INTEGER         :: i, j
    LOGICAL         :: overlap

    REAL    :: dx, dy, dz, r2

    REAL, DIMENSION(nDonors+nAcceptors) :: x,y,z
    REAL, DIMENSION(3)  :: randValue

    !-- set the donor positions ----------------------------------------
    i = 1
    DO WHILE (i <= nDonors+nAcceptors)
       overlap = .FALSE.
       !-- choose a random position inside the box
       CALL ranmar(randValue,3) !-- 3vector of uniform values [0 to 1]
       x(i) = boxLenX*randValue(1)
       y(i) = boxLenY*randValue(2)
       z(i) = boxLenZ*randValue(3)
     
       !-- check for overlaps with previously set particles
       j = 1
       DO WHILE (.NOT.overlap .AND. j<i)
          dx = x(i) - x(j)
          dy = y(i) - y(j)
          dz = z(i) - z(j)

          dx = dx - ANINT(dx/boxLenX)*boxLenX
          dy = dy - ANINT(dy/boxLenY)*boxLenY
          dz = dz - ANINT(dz/boxLenZ)*boxLenZ
          
          r2 = dx**2 + dy**2 + dz**2
          overlap = (r2 < rMin2) !-- overlap if distance less than sqrt(rMin2)
          j = j+1 !-- take the next particle
       END DO

       IF (.NOT.overlap) i = i+1
    END DO

    xDonors = x(1:nDonors)
    yDonors = y(1:nDonors)
    zDonors = z(1:nDonors)

    xAcceptors = x(nDonors+1:nDonors+nAcceptors)
    yAcceptors = y(nDonors+1:nDonors+nAcceptors)
    zAcceptors = z(nDonors+1:nDonors+nAcceptors)

    WRITE(*,*) 'Random initiation done!'

  END SUBROUTINE randomInitPositions

END MODULE init
