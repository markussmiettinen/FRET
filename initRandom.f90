PROGRAM initRandom
  !
  !   A program that builds a PBC configuration of 'nDonors'
  !   donors and 'nAcceptors' acceptors of radius 'rMin' randomly
  !   placed in a simulation box of size 'boxLenX' x 'boxLenY'
  !   x 'boxLenZ', and writes it to the screen.
  !
  !   Markus Miettinen, Thursday, February 21, 2013
  !--------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, PARAMETER :: nDonors = 1500
  INTEGER, PARAMETER :: nAcceptors = 1500
  INTEGER, PARAMETER :: seed = 5

  REAL, PARAMETER :: rMin2 = 1.0**2 !-- squared minimum distance between particles (nm2)

  REAL, PARAMETER :: boxLenX = 232.5
  REAL, PARAMETER :: boxLenY = 232.5
  REAL, PARAMETER :: boxLenZ = 232.5

  REAL, DIMENSION(nDonors)    :: xDonors, yDonors, zDonors
  REAL, DIMENSION(nAcceptors) :: xAcceptors, yAcceptors, zAcceptors

  INTEGER         :: i, j
  LOGICAL         :: overlap
  
  REAL    :: dx, dy, dz, r2

  REAL, DIMENSION(nDonors+nAcceptors) :: x,y,z
  REAL, DIMENSION(3)  :: randValue


  CALL rmarin(seed,0,0) !-- initialize random number generator

  
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
  
  WRITE(*,*) boxLenX, boxLenY, boxLenZ
  WRITE(*,*) 1, 1, 1
  WRITE(*,*) nDonors
  DO i=1,nDonors
     WRITE(*,*) xDonors(i), yDonors(i), zDonors(i)
  END DO
  WRITE(*,*) nAcceptors
  DO i=1,nAcceptors
     WRITE(*,*) xAcceptors(i), yAcceptors(i), zAcceptors(i)
  END DO

END PROGRAM initRandom
