!-------------------------------------------------
PROGRAM doFRET
!
! For simulating FRET in a complex geometry.
!
! Markus Miettinen, Wednesday, July  2, 2014
!-------------------------------------------------
USE init, ONLY : &
     initializeRandomNumberGenerator, randomInitPositions
USE datastructures, ONLY : &
     firstDonor, &
     t, maxTime, maxNphotons, photonsFromDonors, &
     assignAcceptorsToDonors, &
     deExciteAcceptors, exciteDonors, &
     moveDyesAndUpdateFRETefficiencies, &
     checkDonorsForFRET, &
     writeOutResults, readInStructure, &
     reExciteRightAfterAllDonorsDeexcited
IMPLICIT NONE
INTEGER :: percentReady = 0

!-- initialize the random number generator
CALL initializeRandomNumberGenerator

!-- readInStructure from file 'structure.xyz'
CALL readInStructure

!-- Assign acceptors to their donors:
CALL assignAcceptorsToDonors

collectPhotonsFromDonors: DO WHILE (MAXVAL(photonsFromDonors).LT.maxNphotons)
   !-- report current status of simulation:
   writeProgressIntoStandardError: IF (percentReady < 100*MAXVAL(photonsFromDonors)/maxNphotons) THEN
      percentReady = 100*MAXVAL(photonsFromDonors)/maxNphotons
      WRITE(UNIT=0,FMT="(I3,'% done',A,$)",ADVANCE='NO') percentReady,CHAR(13)
   END IF writeProgressIntoStandardError

   !-- Excite donors (at least one will be excited)
   !-- build a list of excited donors, and de-excite
   !-- all acceptors:
   CALL exciteDonors

   !-- move dyes and calculate FRET efficiencies
!   CALL moveDyesAndUpdateFRETefficiencies

   !-- Check if donors FRET, fluoresce or quench:
   DO t=1,maxTime
      CALL deExciteAcceptors
      CALL checkDonorsForFRET
      noMoreExcitedDonors: IF (firstDonor.EQ.0) THEN
         IF (reExciteRightAfterAllDonorsDeexcited) THEN
            ! WRITE(*,*) 'No more excited donors, will excite again. -> Acceptor photon count will be inaccurate!'
            CYCLE collectPhotonsFromDonors
         END IF
      END IF noMoreExcitedDonors
   END DO
END DO collectPhotonsFromDonors

CALL writeOutResults

END PROGRAM doFRET

