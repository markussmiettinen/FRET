!==================================================================
MODULE datastructures
!
! Module containing the constants and variables needed for FRET
! simulation. And the functions to operate on them.
!
! Markus Miettinen, Wednesday, July  2, 2014
!==================================================================
  IMPLICIT NONE

  !-- all structures by definition visible outside this module
  PUBLIC

  !== Simulation constants =========================================
  !-- control parameters:
  LOGICAL, PARAMETER :: reExciteRightAfterAllDonorsDeexcited = .FALSE.

  !-- numerical parameters:
  INTEGER, PARAMETER  :: seed = 5 !-- the seed for random numbers
  REAL, PARAMETER :: pi = 3.14159265359 !-- pi

  !== Physical parameters ==========================================
  !REAL, PARAMETER :: boxLenX = 232.5 !-- size of the unit cell (nm)
  !REAL, PARAMETER :: boxLenY = 232.5
  !REAL, PARAMETER :: boxLenZ = 232.5
  !REAL, PARAMETER :: A = boxLenX*boxLenY
  !
  INTEGER, PARAMETER :: nDonors    = 1497 !-- total number of donors in unit cell
  INTEGER, PARAMETER :: nAcceptors = 129 !-- total numer of acceptors

  !------- read these in >>>
  REAL :: boxLenX !-- size of the unit cell (nm)
  REAL :: boxLenY
  REAL :: boxLenZ
  REAL :: A       !-- area of a side of the unit cell
  INTEGER :: periodicX, periodicY, periodicZ !-- 1 (0) if (not) periodic
  !<<< in eseht dear -------

  REAL, DIMENSION(1:3)    :: boxLens    !-- boxlengths in a vector 
  INTEGER, DIMENSION(1:3) :: periodics  !-- periodicities in a vector

  REAL, DIMENSION(nDonors) :: &
       xDonors, yDonors, zDonors !-- positions of donors
  REAL, DIMENSION(nAcceptors) ::  &
       xAcceptors, yAcceptors, zAcceptors !-- positions of acceptors

  INTEGER, DIMENSION(nDonors) :: &
       nextDonor !-- the index of the donor to be assessed next
  INTEGER, DIMENSION(nAcceptors) :: &
       excitableAcceptors = 1 !-- 0 if acceptor excited, 1 if still excitable

  REAL, PARAMETER :: experimentDuration = 50 !-- in nanoseconds
  INTEGER, PARAMETER :: maxTime = 1024 !-- experiment broken in this many time intervals
  REAL, PARAMETER :: dt = experimentDuration / maxTime !-- time step length (ns)

  INTEGER, PARAMETER :: maxNphotons = 10000 !-- max number of photons per time interval per experiment

  REAL, PARAMETER :: r0  = 5.0643 !-- Foerster radius (nm)
  REAL, PARAMETER :: r02 = r0**2

  REAL, PARAMETER :: &
       rDonor = 0.3, & !-- radius within which the donor dye moiety can move (nm)
       rAcceptor = 0.3 !-- radius within which the acceptor dye moiety can move (nm)

  REAL, PARAMETER :: &
       donorAcceptorCutoff  &
       = 2*r0 + rDonor + rAcceptor !-- beyond this donor will not excite acceptor (nm)
  REAL, PARAMETER :: &
       donorAcceptorCutoff2 = donorAcceptorCutoff**2

!!  REAL, PARAMETER :: irradiance = 25e-18 !-- Laser irradiance (J/s.nm2)
!!  REAL, PARAMETER :: hc = 1.98644568e-16 !-- Planck times light (J.nm)
!  REAL, PARAMETER :: irradiancePerHC = 0.25/1.98644568 !-- 1/s.nm3
!  REAL, PARAMETER :: laserWavelength = 450 !-- nm
!  REAL, PARAMETER :: photonDensity = irradiancePerHC*laserWavelength !-- 1/s.nm2
!  REAL, PARAMETER :: nPhotonsPerSecond = photonDensity*A !-- 1/s
!  REAL, PARAMETER :: excitationTime = 60e-12 !-- s
!!  REAL, PARAMETER :: nPhotonsPerExcitation = REALLY REALLY FEW???

!!  REAL, PARAMETER :: excitementProbability = &
!!       (irradiance*dt*A/(Ephoton*nDonors))*(1-10**(extinctCoeff*nDonors/A))
  REAL, PARAMETER :: excitementProbability = 0.51 !==== PUT HERE THE RIGHT FORMULA!!! ==========

  REAL, PARAMETER :: donorFluorescenceRate    = 0.1738  !-- events/ns 
  REAL, PARAMETER :: donorQuenchRate          = 0.20004 !-- events/ns 
  REAL, PARAMETER :: acceptorDeExcitationRate = 0.2     !-- events/ns
  REAL, PARAMETER :: acceptorDeExcitationProb = 1-EXP(-acceptorDeExcitationRate*dt)

  REAL, PARAMETER :: &
       maxDeexcitationRateWithoutFRET = donorQuenchRate + donorFluorescenceRate, &
       onePerTauD = maxDeexcitationRateWithoutFRET

  !-- Simulation variables ------------------------------------------
  !-- numerical variables:
  INTEGER :: t !-- current time step
  INTEGER :: firstDonor !-- the index of the donor to be assessed first

  TYPE donortype
     !-- The fundamental data structure, containing 2 vectors:
     !-- (1) the identities of the possible acceptors for this donor,
     !-- and
     !-- (2) the rates of FRET to these acceptors,
     !-- as well
     !-- the total number of possible acceptors for this donor,
     !-- and the maximal possible de-excitation rate and -probability.
     INTEGER, DIMENSION(:), ALLOCATABLE :: acceptorID
     REAL, DIMENSION(:), ALLOCATABLE :: FRETrate
     INTEGER :: nDonorsAcceptors
     REAL :: maxDeexcitationRate, maxDeexcitationProb
  END TYPE donortype

  TYPE(donortype), DIMENSION(nDonors) :: donors

  !-- Measurables ---------------------------------------------------
  INTEGER, DIMENSION(maxTime) :: &
       photonsFromAcceptors = 0, photonsFromDonors = 0
  INTEGER, DIMENSION(maxTime) :: &
       FRETevents = 0, donorQuenchEvents = 0

CONTAINS
  !-----------------------------------------------------
  SUBROUTINE readInStructure
  !
  ! Read in (a unit cell of) the structure, given
  ! in the file 'structure.xyz'.
  !
  ! Markus Miettinen, Thursday, February 21, 2013
  !-----------------------------------------------------

    IMPLICIT NONE

    INTEGER :: allocStat, k
    INTEGER :: nDonorsIN, nAcceptorsIN

    OPEN(UNIT=101,FILE='structure.xyz',STATUS='OLD',ACTION='READ')

    READ(101,*) boxLenX,  boxLenY,  boxLenZ
    boxLens(1:3) = (/ boxLenX,  boxLenY,  boxLenZ /)
    READ(101,*) periodicX, periodicY, periodicZ
    periodics(1:3) = (/ periodicX, periodicY, periodicZ /)

    READ(101,*) nDonorsIN
    IF (nDonorsIN.NE.nDonors) THEN
       WRITE(*,*) 'Number of donors = ',nDonorsIN,', not what expected = ',nDonors
       STOP 
    END IF
    readInDonorPositions: DO k=1,nDonors
       READ(101,*) xDonors(k), yDonors(k), zDonors(k)
    END DO readInDonorPositions

    READ(101,*) nAcceptorsIN
    IF (nAcceptorsIN.NE.nAcceptors) THEN
       WRITE(*,*) 'Number of acceptors = ',nAcceptorsIN,', not what expected = ',nAcceptors
       STOP 
    END IF
    readInAcceptorPositions: DO k=1,nAcceptors
       READ(101,*) xAcceptors(k), yAcceptors(k), zAcceptors(k)
    END DO readInAcceptorPositions

  END SUBROUTINE readInStructure


  !----------------------------------------------------------------
  SUBROUTINE assignAcceptorsToDonors
  !
  ! Acceptors within 'donorAcceptorCutoff' from a given donor are
  ! assigned to it and their rates for FRET are calculated using
  ! the formula:
  !
  !             1    (R_0)^6
  ! k_{FRET} = --- * -------
  !            t_D     r^6
  !
  ! i.e., Eq. (1) of Biophys. J. 89 3822 (2005).
  !
  ! Markus Miettinen, Saturday, November 2, 2013
  !----------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: allocStat
    INTEGER :: i, j, nAcceptorsOfThisDonor
    REAL    :: dx, dy, dz, r2
    REAL    :: FRETrate, maxDeexcitationRate
    INTEGER, DIMENSION(nAcceptors) :: acceptorIndecesOfThisDonor
    REAL, DIMENSION(nAcceptors)    :: FRETratesOfThisDonor
    CHARACTER(LEN=72) :: outputform 
    CHARACTER(LEN=5) :: cutoffINtext

    WRITE (cutoffINtext,'(F5.2)') donorAcceptorCutoff
    outputform = &
         "('These donors have no acceptors within the cutoff (" &
         // cutoffINtext // &
         " nm):',/,I6)"
    overDonors: DO i=1,nDonors
       nAcceptorsOfThisDonor = 0
       maxDeexcitationRate = donorQuenchRate + donorFluorescenceRate
       overAcceptors: DO j=1,nAcceptors
          dx = xDonors(i) - xAcceptors(j)
          dy = yDonors(i) - yAcceptors(j)
          dz = zDonors(i) - zAcceptors(j)

          dx = dx - ANINT(dx/boxLenX)*boxLenX*periodicX
          dy = dy - ANINT(dy/boxLenY)*boxLenY*periodicY
          dz = dz - ANINT(dz/boxLenZ)*boxLenZ*periodicZ
          
          r2 = dx**2 + dy**2 + dz**2

          withinCutoff: IF (r2 < donorAcceptorCutoff2) THEN
             nAcceptorsOfThisDonor = nAcceptorsOfThisDonor + 1
             acceptorIndecesOfThisDonor(nAcceptorsOfThisDonor) = j
             FRETrate = onePerTauD*(r02/r2)**3
             FRETratesOfThisDonor(nAcceptorsOfThisDonor) = FRETrate
             maxDeexcitationRate = maxDeexcitationRate + FRETrate
          END IF withinCutoff
       END DO overAcceptors

       donorHasAcceptors: IF (nAcceptorsOfThisDonor > 0) THEN
          !-- allocate space for the acceptor indeces and probabilities:
          ALLOCATE(donors(i)%acceptorID(nAcceptorsOfThisDonor), STAT = allocstat)
          IF ( allocstat /= 0 ) STOP 'Memory allocation error!'
          ALLOCATE(donors(i)%FRETrate(nAcceptorsOfThisDonor), STAT = allocstat)
          IF ( allocstat /= 0 ) STOP 'Memory allocation error!'
          !-- save the acceptor indeces and probabilities
          donors(i)%acceptorID(:) = acceptorIndecesOfThisDonor(1:nAcceptorsOfThisDonor)
          donors(i)%FRETrate(:)   = FRETratesOfThisDonor(1:nAcceptorsOfThisDonor)
          donors(i)%nDonorsAcceptors    = nAcceptorsOfThisDonor
          donors(i)%maxDeexcitationRate = maxDeexcitationRate
          donors(i)%maxDeexcitationProb = 1-EXP(-maxDeexcitationRate*dt)
       ELSE !-- donor has no acceptors, however...
          !---- ... allocate single-space, and put the rate of FRET to zero.
          ALLOCATE(donors(i)%acceptorID(1), STAT = allocstat)
          IF ( allocstat /= 0 ) STOP 'Memory allocation error!'
          ALLOCATE(donors(i)%FRETrate(1), STAT = allocstat)
          IF ( allocstat /= 0 ) STOP 'Memory allocation error!'
          donors(i)%acceptorID(1) = 1   !-- points to the first acceptor...
          donors(i)%FRETrate(1)   = 0.0 !-- ...but will never see FRET.
          donors(i)%nDonorsAcceptors    = nAcceptorsOfThisDonor
          donors(i)%maxDeexcitationRate = maxDeexcitationRate
          donors(i)%maxDeexcitationProb = 1-EXP(-maxDeexcitationRate*dt)
          WRITE(*,FMT=outputform,ADVANCE='NO') i
          outputform="(I6)"
       END IF donorHasAcceptors
    END DO overDonors
    WRITE(*,"(/)")
  END SUBROUTINE assignAcceptorsToDonors


  !----------------------------------------------------------------
  SUBROUTINE exciteDonors
  !
  ! Excite the donors with a light pulse. In other words, make them
  ! excited with a probability
  !
  ! P(excitement) = 
  !
  ! The excited donors are put in a randomly ordered list, following
  ! which they will be later assessed.
  !
  ! NOTE. Also all the acceptors will be made excitable. This is
  ! expected to correspond best to the experimental setup where
  ! each donor excitation only takes place after such a long time
  ! time after the previous excitation that no more excited donors
  ! or acceptors should be present in the sample.
  !
  ! NOTE. At least one donor will always be excited. This is made
  ! sure by trying to excite until at least one is excited.
  !
  ! Markus Miettinen, Wednesday, July  2, 2014
  !----------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: k
    INTEGER :: prev,this,next
    REAL, DIMENSION(nDonors)  :: randValue

    !-- make all the acceptors again excitable:
    excitableAcceptors = 1

!    WRITE(*,*) '**** WARNING! Wrong formula used for Donor excitation probabilities! WARNING! ****'
    firstDonor = 0 !-- the first donor to be assessed so far undetermined

    toExciteAtLeastOneDonor: DO WHILE (firstDonor.EQ.0)
       CALL ranmar(randValue,nDonors) !-- vector of uniform values [0 to 1]
       overDonors: DO this = 1,nDonors
          excited: IF (randValue(this).LT.excitementProbability) THEN
             listExists: IF (firstDonor.GT.0) THEN
                next = firstDonor
                newFirst: IF (randValue(next).GT.randValue(this)) THEN
                   nextDonor(this) = firstDonor
                   firstDonor = this
                ELSE !-- this is larger than the first, so comes later in the list 
                   whileNotLastOfTheList: DO WHILE (nextDonor(next).GT.0)
                      prev = next
                      next = nextDonor(next)
                      thisGoesBeforeNext: IF (randValue(next).GT.randValue(this)) THEN
                         nextDonor(this) = next
                         nextDonor(prev) = this
                         CYCLE overDonors
                      END IF thisGoesBeforeNext
                   END DO whileNotLastOfTheList
                   !-- that was last of the list, make this the new last
                   nextDonor(next) = this 
                   nextDonor(this) = 0
                END IF newFirst
             ELSE !-- no list exists, so create it
                firstDonor = this
                nextDonor(this) = 0
             END IF listExists
          END IF excited
       END DO overDonors
    END DO toExciteAtLeastOneDonor

  END SUBROUTINE exciteDonors


  !-----------------------------------------------------------------
  SUBROUTINE moveDyesAndUpdateFRETefficiencies
  !
  ! Move the relevant dye moieties (excited donors and their
  ! acceptors) around their fixed points and calculate the
  ! FRET efficiencies corresponding to the new positions.
  !
  ! Markus Miettinen, Saturday, November  2, 2013
  !-----------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: i, thisDonor, thisAcceptor
    REAL :: r2, FRETrate, sumOfFRETrates, maxDeexcitationRate

    REAL, DIMENSION(1:3) :: &
         xyzDonorNew, dxdydz
    LOGICAL, DIMENSION(nAcceptors) :: &
         thisAcceptorWasNotYetMoved
    REAL, DIMENSION(1:3,nAcceptors) :: &
         xyzAcceptorNew

    thisAcceptorWasNotYetMoved = .TRUE.

    thisDonor = firstDonor
    forAllExcitedDonors: DO
       !-- move thisDonor
       xyzDonorNew = moveDonor(thisDonor)
       sumOfFRETrates = 0.0
       !-- move the acceptors of thisDonor
       overAcceptors: DO i=1,donors(thisDonor)%nDonorsAcceptors
          thisAcceptor = donors(thisDonor)%acceptorID(i)
          needToMove: IF (thisAcceptorWasNotYetMoved(thisAcceptor)) THEN
             xyzAcceptorNew(1:3,thisAcceptor) = moveAcceptor(thisAcceptor)
             thisAcceptorWasNotYetMoved(thisAcceptor) = .FALSE.
          END IF needToMove
          !-- calculate distance to donor
          dxdydz = xyzDonorNew(1:3) - xyzAcceptorNew(1:3,thisAcceptor)
          dxdydz = dxdydz - ANINT(dxdydz/boxLens)*boxLens*periodics
          r2 = SUM(dxdydz**2)
          !-- update values for FRET
          FRETrate = onePerTauD*(r02/r2)**3
          donors(thisDonor)%FRETrate(i) = FRETrate
          sumOfFRETrates = sumOfFRETrates + FRETrate
       END DO overAcceptors
       maxDeexcitationRate = &
            maxDeexcitationRateWithoutFRET + sumOfFRETrates
       donors(thisDonor)%maxDeexcitationRate = maxDeexcitationRate
       donors(thisDonor)%maxDeexcitationProb = 1-EXP(-maxDeexcitationRate*dt)
       thisDonor = nextDonor(thisDonor)
       thatWasTheLastExcitedDonor: IF (thisDonor.EQ.0) THEN
          EXIT forAllExcitedDonors
       END IF thatWasTheLastExcitedDonor
    END DO forAllExcitedDonors
    
  CONTAINS
    !--------------------------------------------------
    FUNCTION moveDonor(thisDye) RESULT(here)
    !
    ! Returns the position at which the donor will be
    ! during this round of donor excitations.
    !
    ! Markus Miettinen,Tuesday, November 5, 2013
    !--------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: thisDye
      REAL, DIMENSION(1:3) :: here
      
      REAL, DIMENSION(1:2) :: &
           psiPhi, sinPsiPhi, cosPsiPhi

      CALL ranmar(psiPhi,2) !-- between [0,1]
      psiPhi = pi*psiPhi !-- [0,pi]
      psiPhi(2) = 2*psiPhi(2) !-- [0,2pi]
      sinPsiPhi = SIN(psiPhi)
      cosPsiPhi = COS(psiPhi)
      here(1) = xDonors(thisDye) + rDonor*sinPsiPhi(1)*cosPsiPhi(2)
      here(2) = yDonors(thisDye) + rDonor*sinPsiPhi(1)*sinPsiPhi(2)
      here(3) = zDonors(thisDye) + rDonor*cosPsiPhi(1)
    END FUNCTION moveDonor

    !----------------------------------------------------
    FUNCTION moveAcceptor(thisDye) RESULT(here)
    !
    ! Returns the position at which the acceptort will be
    ! during this round of donor excitations.
    !
    ! Markus Miettinen, Tuesday, November 5, 2013
    !----------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: thisDye
      REAL, DIMENSION(1:3) :: here
      
      REAL, DIMENSION(1:2) :: &
           psiPhi, sinPsiPhi, cosPsiPhi

      CALL ranmar(psiPhi,2) !-- between [0,1]
      psiPhi = pi*psiPhi !-- [0,pi]
      psiPhi(2) = 2*psiPhi(2) !-- [0,2pi]
      sinPsiPhi = SIN(psiPhi)
      cosPsiPhi = COS(psiPhi)
      here(1) = xAcceptors(thisDye) + rAcceptor*sinPsiPhi(1)*cosPsiPhi(2)
      here(2) = yAcceptors(thisDye) + rAcceptor*sinPsiPhi(1)*sinPsiPhi(2)
      here(3) = zAcceptors(thisDye) + rAcceptor*cosPsiPhi(1)
    END FUNCTION moveAcceptor


  END SUBROUTINE moveDyesAndUpdateFRETefficiencies



  !--------------------------------------------------------
  SUBROUTINE deExciteAcceptors
  !
  ! Try de-exciting the excited acceptors. 
  !
  ! Markus Miettinen, Thursday, January 17, 2013
  !--------------------------------------------------------

    IMPLICIT NONE
    INTEGER :: i
    REAL :: randValue

    DO i=1,nAcceptors
       acceptorExcited: IF (excitableAcceptors(i).EQ.0) THEN
          CALL ranmar(randValue,1)
          deExciteAcceptor: IF (randValue.LT.acceptorDeExcitationProb) THEN
             excitableAcceptors(i) = 1
             photonsFromAcceptors(t) = photonsFromAcceptors(t) + 1             
          END IF deExciteAcceptor
       END IF acceptorExcited
    END DO

  END SUBROUTINE deExciteAcceptors



  !-----------------------------------------------------------------
  SUBROUTINE checkDonorsForFRET
  !
  ! Go through the excited donors in the order defined in subroutine
  ! exciteDonors, and check if they perform FRET to their acceptors,
  ! or otherwise de-excite (either via fluorescing or quenching).
  !
  ! Markus Miettinen, Wednesday, January 16, 2013
  !-----------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: k
    INTEGER :: lastExcited, prev, this

    startingFromFirstDonor: DO
       thereStillAreExcitedDonors: IF (firstDonor.GT.0) THEN
          firstDonorDeExcites: IF (donorDeExcites(firstDonor)) THEN
             firstDonor = nextDonor(firstDonor)
             CYCLE startingFromFirstDonor
          ELSE !-- first donor does not de-excite
             lastExcited = firstDonor
             this = firstDonor
             whileNotTheLastInList: DO WHILE (nextDonor(this).GT.0)
                prev = this
                this = nextDonor(this)
                thisDonorStaysExcited: IF (.NOT.donorDeExcites(this)) THEN
                   previousWasDeExcited: IF (prev.NE.lastExcited) THEN
                      nextDonor(lastExcited) = this
                   END IF previousWasDeExcited
                   lastExcited = this
                ELSE !-- this donor de-excites!
                   nextDonor(lastExcited) = nextDonor(this)
                END IF thisDonorStaysExcited
             END DO whileNotTheLastInList
             !-- end of the list has been reached, exit for making next time step
             EXIT startingFromFirstDonor
          END IF firstDonorDeExcites
       ELSE !-- all donors have been de-excited
          EXIT startingFromFirstDonor
       END IF thereStillAreExcitedDonors
    END DO startingFromFirstDonor

  CONTAINS
    !---------------------------------------------------
    FUNCTION donorDeExcites(thisDonor) RESULT(deExcited)
    !
    ! Returns true if the donor de-excited. Checks for
    ! Quenching, Fluorescence and FRET. And updates then
    ! corresponding measurables if needed.
    !
    ! Markus Miettinen, Thursday, February 21, 2013
    !---------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: thisDonor
      LOGICAL :: deExcited
      
      INTEGER :: i
      REAL :: randValue
      REAL :: currentDeexcitationRate
      REAL :: rateLimit, rateSum

      CALL ranmar(randValue,1)
      
      thisDonorPossiblyDeexcites: IF (randValue.LT.donors(thisDonor)%maxDeexcitationProb) THEN
         currentDeexcitationRate = donors(thisDonor)%maxDeexcitationRate
         getCurrentDeexcitationRate: DO i=1,donors(thisDonor)%nDonorsAcceptors 
            acceptorUnavailable: IF (excitableAcceptors(donors(thisDonor)%acceptorID(i)).EQ.0) THEN
               currentDeexcitationRate = &
                    currentDeexcitationRate - donors(thisDonor)%FRETrate(i)
            END IF acceptorUnavailable
         END DO getCurrentDeexcitationRate
         thisDonorDeexcites: IF (randValue.LT.(1-EXP(-currentDeexcitationRate*dt))) THEN
            deExcited = .TRUE.
            !-- so this donor excites, but how, exactly?
            CALL ranmar(randValue,1) !-- new random number [0,1]
            rateLimit = randValue*currentDeexcitationRate !-- [0,currentDeexcitationRate]
            donorQuenches: IF (rateLimit.LT.donorQuenchRate) THEN
               donorQuenchEvents(t) = donorQuenchEvents(t) + 1
            ELSE !-- no quenching, so maybe fluorescence?
               rateSum = donorQuenchRate + donorFluorescenceRate
               donorFluoresces: IF (rateLimit.LT.rateSum) THEN
                  photonsFromDonors(t) = photonsFromDonors(t) + 1
               ELSE !-- no fluorescence either, so must be FRET
                  FRETevents(t) = FRETevents(t) + 1
                  !-- so this donor FRETs, but to which acceptor?
                  overThisDonorsAcceptors: DO i=1,donors(thisDonor)%nDonorsAcceptors
                     rateSum = &
                          rateSum + &
                          donors(thisDonor)%FRETrate(i) * &
                          excitableAcceptors(donors(thisDonor)%acceptorID(i)) !-- 0 if excited, 1 if not
                     FRETtoThisAcceptor: IF (rateLimit.LT.rateSum) THEN
                        !-- this (the i:th) acceptor excites:
                        excitableAcceptors(donors(thisDonor)%acceptorID(i)) = 0
                     END IF FRETtoThisAcceptor
                  END DO overThisDonorsAcceptors
               END IF donorFluoresces
            END IF donorQuenches
         ELSE !-- this donor stays excited with the current de-excitation rate
            deExcited = .FALSE.
         END IF thisDonorDeexcites
      ELSE !-- this donor will stay excited even if all its acceptors are free
         deExcited = .FALSE.
      END IF thisDonorPossiblyDeexcites

    END FUNCTION donorDeExcites

  END SUBROUTINE checkDonorsForFRET



  !-----------------------------------------------------
  SUBROUTINE writeOutResults
  !
  ! Write out the results.
  !
  ! Markus Miettinen, Thursday, January 17, 2013
  !-----------------------------------------------------

    IMPLICIT NONE

    INTEGER :: k
!    REAL :: time

!    WRITE(*,*) '#  time (ns)   A-photons   D-photons       FRETs    Quenches'
    WRITE(*,*) '#  timestep   A-photons   D-photons       FRETs    Quenches'
    DO k=1,maxTime
!       time = k*dt
!       WRITE(*,*) time, &
       WRITE(*,*) k, &
            photonsFromAcceptors(k), &
            photonsFromDonors(k), &
            FRETevents(k), &
            donorQuenchEvents(k)
    END DO
  END SUBROUTINE writeOutResults



END MODULE datastructures
