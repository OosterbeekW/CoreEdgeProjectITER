      MODULE PARAM
      PARAMETER (NRMAX = 401,  NPMAX = 401)
      PARAMETER (MAXNODE=NRMAX*NPMAX,MBMAX=1026)
      PARAMETER (NRMMAX = 801, NPMMAX =4098)
      PARAMETER (MAXMNODE = NRMMAX*NPMMAX)
      PARAMETER (NPTSMAX = 1601)
      END
C-----------------------------------------------------------------------
      MODULE COMMAX
      PARAMETER (NPSIMAX=801, NCHIMAX=4098, NMAX=NPSIMAX)
      PARAMETER (NPNC=NPSIMAX*NCHIMAX)
      END
C-----------------------------------------------------------------------
      MODULE COMPIO
      PARAMETER (NIN=10, NOUT=20, NMAP=12)
      END
C-----------------------------------------------------------------------
      MODULE CORNERS
      USE PARAM
      REAL    RS(4,2)
      INTEGER IJ(4,2), NODENO(MAXMNODE,4)
      END
C-----------------------------------------------------------------------
      MODULE GAUSINT
      REAL XGAUSS(4),WGAUSS(4)
      REAL H(4,4,4,4),HR(4,4,4,4),HS(4,4,4,4),HRS(4,4,4,4)
      END
C-----------------------------------------------------------------------
      MODULE COMDAT
      REAL    ELLIP,TRIA,QUAD,PAR1,PAR2,PAR3,PAR4,
     >        AGA,BGA,CGA,DGA,EGA,FGA,GGA,HGA,
     >        API,BPI,CPI,DPI,EPI,FPI,GPI,HPI,PEDGE,
     >        ACUR,BCUR,CCUR,DCUR,ECUR,FCUR,GCUR,HCUR, BETAP_PED,
     >        ERRIT,ERRCUR,EPS,ALFA,B,C,XIAB,Q95,BETAP,AMIX,BMIX,
     >        ABB, BBB, AMPL, RVAC,BVAC,ZEFF,ZN0,RPE,ETAEI, COREP,
     >	      PSIPED !cWO
      INTEGER IAS,IAV,ICUR,NRCUR,NPCUR,NMESH,NBB,NQB,
     >        MHARM,ISHAPE,ISOL,IGAM,IPAI,NR,NP,NRMAP,NPMAP,NITER
      END
C-----------------------------------------------------------------------
      MODULE COMMAP
      USE PARAM
      REAL CS(NRMMAX),QS(NRMMAX),DQS(NRMMAX),CURJ(NRMMAX),CHI(NPMMAX),
     >     GEM11(MAXMNODE),GEM12(MAXMNODE),GEM33(MAXMNODE),
     >     CHIKN(NPMMAX),P0(NRMMAX),RBPHI(NRMMAX),
     >     DP(NRMMAX),DRBPHI(NRMMAX),
     >     DQEC,DJ0,DJE,CPSURF,RADIUS,RAXIS,DP0,DPE,DRBPHI0,DRBPHIE
      INTEGER JS0,NCHI,NPSI
      END
C-----------------------------------------------------------------------
      MODULE COMPRI
      INTEGER NPR1,NPR2,NROUT,NDIAG
      END
C-----------------------------------------------------------------------
      MODULE COMPLO
      INTEGER       NPL1
      CHARACTER*100 TXTOUT(40)
      END
C-----------------------------------------------------------------------
      MODULE COMANG
      REAL ANGLE
      END
C-----------------------------------------------------------------------
      MODULE MESH
      USE PARAM
      REAL    XXOLD(4,MAXNODE),YYOLD(4,MAXNODE),PSIOLD(4*MAXNODE)
      END
C-----------------------------------------------------------------------
      MODULE MESHAC
      USE PARAM
      REAL    AMESH,BMESH,CMESH
      REAL    XR1,SIG1,XR1DONE,SIG1DONE
      REAL    SG(NRMMAX),DSG(NRMMAX),DDSG(NRMMAX)
      INTEGER IMESH,IARC,NRDONE
      END
C-----------------------------------------------------------------------
      MODULE NODES
      USE PARAM
      REAL    PSIKN(NRMMAX),THTKN(NPMMAX),RADPSI(NRMMAX)
      REAL    DPSIKN(NRMMAX),DDPSIKN(NRMMAX)
      END
C-----------------------------------------------------------------------
      MODULE TOLERA
      REAL    PSITOL,THTTOL,TOL
      END
C-----------------------------------------------------------------------
      MODULE FF
      REAL    CPSI,CTHT,CXAXIS,CYAXIS
      INTEGER N1,N2,N3,N4
      END
C-----------------------------------------------------------------------
      MODULE FAXIS
      USE PARAM
      REAL    PSI(4*MAXMNODE)
      INTEGER NAXIS
      END
C-----------------------------------------------------------------------
      MODULE COMOUT
      REAL    BETAPL,BETA, BETAP_F, BETAP_PEDF
      END
C-----------------------------------------------------------------------
      MODULE COMPROF
      USE PARAM
      REAL    DPR(NPTSMAX),DF2(NPTSMAX),ZJZ(NPTSMAX),QIN(NPTSMAX),
     >        DPRES(1001),DGAM(1001),PINT(1001),GINT(1001)
      INTEGER NPTS
      END
C-----------------------------------------------------------------------
      MODULE COMSOLV
c      USE PARAM
c      REAL KKBIG(KKLDA,4*MAXNODE)
      REAL, ALLOCATABLE :: KKBIG(:,:)
      INTEGER KKLDA
      END
C-----------------------------------------------------------------------
      MODULE COMB02
      USE COMMAX
      REAL    B02(NPNC), DTB02(NPNC), DSB02(NPNC)
      END
C-----------------------------------------------------------------------
      MODULE COMPQ
      USE COMMAX
      REAL    CP0(NPNC),CP1(NPNC),CP2(NPNC),CQ0(NPNC),CQ1(NPNC)
      INTEGER NCPQ
      END
C-----------------------------------------------------------------------
      MODULE COMSPL
      USE COMMAX
      REAL     Q1(NPSIMAX), Q2(NPSIMAX), Q3(NPSIMAX), Q4(NPSIMAX),
     R         P1(NPSIMAX), P2(NPSIMAX), P3(NPSIMAX), P4(NPSIMAX),
     R         RBP1(NPSIMAX), RBP2(NPSIMAX), RBP3(NPSIMAX),
     R         RBP4(NPSIMAX)
      END
C-----------------------------------------------------------------------
      MODULE COMNAM
      REAL     QAXIS,TBB,TBF
      END
C-----------------------------------------------------------------------
      MODULE COMPIE
      REAL PI
      END
C
************************************************************************
*DECK HELENA
      PROGRAM HELENA
C-----------------------------------------------------------------------
C
C MAIN PROGRAM HELENA :        (VERSION 9  DATE 28-09-95)
C ---------------------
C
C      - SOLVES THE 2D GRAD-SHAFRANOV EQUATION FOR  ARBITRARY UP/DOWN
C        SYMMETRIC CONTINUOUS PLASMA BOUNDARIES AND EQUILIBRIUM PRESSURE
C        AND GAMMA PROFILES.
C      - 2D CUBIC ISOPARAMETRIC FINITE ELEMENTS ARE USED FOR AN ACCURATE
C        REPRESENTATION OF THE SOLUTION.
C      - THE FINAL SOLUTION IS OBTAINED AN A FLUXSURFACE GRID.
C
C      - HELENA RUNNING ON CRAY
C
C
C PROGRAM ORGANIZATION :
C ----------------------
C
C        BLOCKDATA                       : INITIALIZE NAMELIST VAR.
C        HELENA                          : MAIN PROGRAM
C          (XUFLOW)                      : IBM ERROR HANDLING
C          INIVAL                        : INITIALIZE VALUES
C          OUT                           : PRINT OUTPUT UNIT 20
C          GAUSS                         : INITIALIZE GAUSSIAN POINTS
C            CUBICH                      : DEF. CUBIC ELEMENTS
C          SOLSHP                        : BOUNDARY SHAPE SOLOVIEV EQ.
C            (ZERO)                      : HGOLIB ROUTINE
C            FSOL
C              (RFT2)                    : HGOLIB ROUTINE
C          FSHAPE
C            (GRIDINV)                   : HGOLIB ROUTINE
C            (RFT2)
C          ELMNO                         : INITIALIZE ELEMENT NUMBERING
C          INIGRID                       : DEF. INITIAL GRID
C            RADB                        : RADIAL PROFILE OF ELEMENTS
C          (BEGPLT)                      : PPPLIB ROUTINE
C          (LBLTOP)                      : PPPLIB ROUTINE
C          (LBLBOT)                      : PPPLIB ROUTINE
C          PLOTGR                        : PLOT A ISOPARAMETRIC GRID
C            (NFRAME)                    : PPPLIB ROUTINE
C            PLOTCU                      : PLOT CUBIC LINE IN X,Y PLANE
C              CUB1D                     : 1D CUBIC INTERPOLATION
C              (LPLOT6)                  : PPPLIB ROUTINE
C          INITKQ                        : INIT. MATRIX KK AND VECTOR Q
C          FORMKQ                        : CALC. MATIX KK AND VECTOR Q
C            DPDPSI                      : DERIVATIVE OF PRESSURE
C            DGDPSI                      :     ,,     OF GAMMA
C            INTERP                      : CUBIC ELM. INTERPOLATION
C          BOUNKQ                        : INSERT BOUNDARY COND.
C            INTERP
C            DELRC                       : DELETE ROWS AND COLUMNS
C          SOLVE                         : SOLVE MATRIX EQUATIONS
C            CONJGR                      : CONJUGATE GRADIENTS
C              SHRINK                    : REMOVE ROWS/COLUMNS B.C.
C              SCALE                     : SCALE MATRIX PROBLEM
C              ASUB                      : MATRIX VECTOR INPRODUCT
C                EXPAND                  : INSERT ROWS/COLUMNS B.C.
C                SHRINK
C           RESTORE                      : INSERT B.C. IN SOLUTION
C           FINDAXIS                     : FIND MAGNETIC AXIS
C             ROOT                       : SOLVE QUADRATIC EQUATION
C             CUB1D
C           NORMAL                       : NORMALIZE PSI SOLUTION
C           REMESH                       : CALC. NEW ISOPARAMETRIC MESH
C             RADMESH
C               RPACK
C               DRPACK
C               DDRPACK
C             PSIMIMA                    : FIND MINIMUM PSI IN ELM.
C               CUB1D
C             THTMIMA                    : FIND MINIMUM THETA IN ELM.
C             POSITION                   : FIND PSI/THETA POINT IN ELM.
C               (C05ZAF)                 : NAG ROUTINE
C               FZERO
C                 INTERP
C               (C05PBF)                 : NAG ROUTINE
C             POSBND                     : FIND PSI/THETA POINT BOUNDARY
C               SOLVP3                   : SOLVE CUBIC EQUATION
C             INTERP
C             RADMESH
C           MAPPING                      : CALC. METRIC FLUX COORD.
C             PROFILES                   : EQ. PROFILES
C               PRES                     : PRESSURE
C               XGAMMA                   : GAMMA
C             RADMESH
C             INTERP
C             PLOTM
C             SOLVP3
C             (LPLOT6)                   : PPPLIB ROUTINE
C           (FINPLT)                     : PPPLIB ROUTINE
C           DIAGNO                       : CALC. BETA, BETAPL
C
C
C THE NUMBERING USED FOR THE INTERPOLATING FUNCTIONS H :
C
C                R0,S0
C         H    H
C         H    H        WITH I,J = 0 OR 1 INDICATES A DERIVATIVE
C         HHHHHH        WITH RESPECT TO R,S. R0 AND S0 REFER TO THE
C         H    H        4 NODES OF ONE ELEMENT AND CAN BE +1. OR -1.
C         H    H
C                I,J
C                                S
C                                .
C                             +1 .
C                    N2----------.-----------N3
C                     |          .           |
C                     |          .           |
C                     |          .           |        --> R
C           .......-1.|......................|.+1.......
C                     |          .           |
C                     |          .           |
C                     |          .           |
C                     |          .           |
C                    N1----------.-----------N4
C                                . -1
C                                .
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMMAP
      USE COMPRI
      USE COMPLO
      USE CORNERS
      USE COMPROF
      USE TOLERA
      USE MESH
      USE MESHAC
      USE FAXIS
      USE COMNAM
      USE COMSOLV
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),FR(2*MBMAX+2)
      REAL QQ(4*MAXNODE),DIAG(4*MAXNODE)
      REAL PSPLOT(101),DPPLOT(101),DGPLOT(101),ZJPLOT(101),QPLOT(101)
      REAL DF2OLD(NPTSMAX)
      REAL FM(MBMAX)
      REAL FCIRC(NRMMAX),B02AV(NRMMAX),B0MAX(NRMMAX),RAV(NRMMAX)
      REAL DRMERC(NRMMAX),DIMERC(NRMMAX),HH(NRMMAX),QPROF(NRMMAX)
      REAL DQPROF(NRMMAX),GEONC(NRMMAX),ZJPAR(NRMMAX), ZJAR(NRMMAX)
      REAL DUMMY1(NRMMAX),DUMMY2(NRMMAX),ZVOL(NRMMAX),ZVOLP(NRMMAX)

C------------------------------------------- DEFINE INPUT NAMELISTS ----
      NAMELIST/SHAPE/  ELLIP,TRIA,QUAD,MHARM,ISHAPE,ISOL,FM,MFM,
     >                 XR1,SIG1,
     >                 PAR1,PAR2,PAR3,PAR4,AMESH,BMESH,CMESH,IMESH,
     >                 IAS,IARC
      NAMELIST/PROFILE/AGA,BGA,CGA,DGA,EGA,FGA,GGA,HGA,
     >                 API,BPI,CPI,DPI,EPI,FPI,GPI,HPI, PEDGE,
     >                 ACUR,BCUR,CCUR,DCUR,ECUR,FCUR,GCUR,HCUR,
     >                 DPR,DF2,ZJZ,QIN,
     >                 ICUR,IGAM,IPAI,NPTS,IAV, COREP, PSIPED
      NAMELIST/PHYS/   EPS,ALFA,B,C,XIAB,Q95,BETAP,
     >                 RVAC,BVAC,ZEFF,ZN0,RPE,ETAEI, BETAP_PED
      NAMELIST/NUM/    NR,NP,NRMAP,NPMAP,NCHI,NITER,NMESH,AMIX,BMIX,
     >                 ERRIT,ERRCUR,NRCUR,NPCUR,
     >                 NBB,ABB,BBB,AMPL
      NAMELIST/PRI/    NPR1,NPR2,NROUT,NDIAG
      NAMELIST/PLOT/   NPL1
      NAMELIST/BALL/   QB1,QB2,NQB

C----------------------------------------- READ INPUT PARAMETERS -------
      CALL INIVAL
C----------------------------------------- REMOVE OPEN STAT. ON IBM ----
      READ(10,SHAPE)
      READ(10,PROFILE)
cWO
C      write(*,*) "READ THE PROFILE !!!!!", ipai
      READ(10,PHYS)
      READ(10,NUM)
      READ(10,PRI)
      READ(10,PLOT)
      READ(10,BALL)
      IF (ICUR.EQ.0) NMESH = 1
      IF ((IAS.EQ.1) .AND. (MOD(NCHI,2).NE.0)) THEN
        STOP ' NCHI MUST BE 2^n'
      ENDIF
      IF ((IAS.EQ.0) .AND. (MOD(NCHI,2).NE.1)) THEN
        STOP ' NCHI MUSt BE 2^n + 1'
      ENDIF
C------------------------------------------ALLOCATE LARGE MATRIX
      IF ((NR.GT.NRMAX).OR.(NP.GT.NPMAX)) THEN
        WRITE(*,*) ' NR or NP too large, NRMAX =',NRMAX,' NPMAX=',NPMAX
        STOP
      ENDIF
      KKLDA = 4*NP+9

      ALLOCATE(KKBIG(KKLDA,4*(NR+2)*NP))
C----------------------------------------- INITIALIZE PLOTTING
      IF (NPL1.EQ.1) THEN
        CALL BEGPLT('PCUBEQ')
        CALL LBLTOP(' ',1)
        CALL LBLBOT('HELENA EQUILIBRIUM VERSION 12',34)
      ENDIF

C-------------------------------- INITIALIZE INTERPOLATING FUNCTIONS ---
      CALL GAUSS
C-------------------------------- INITIAL GUESS FOR A ------------------
      A = 4. * B/ABS(B)
C------------------------------------ CALCULATE SHAPE OF SOLOVIEV ------
      IF (ISOL.EQ.1) THEN
        A = 2*(1. + (1.- EPS**2 /4.)/ELLIP**2)
        B = 2*EPS + TRIA/ELLIP**2  /(1.+ (1.- EPS**2 /4.)/ELLIP**2)
        CALL SOLSHP(FR,MHARM)
        CPSURF = .5*EPS**2/(1.+EPS**2)**2 * ELLIP / SQRT(1.-EPS**2/4.)
        RADIUS = SQRT(EPS**2/(1.+EPS**2))
        B0 = SQRT(1.+EPS**2)
        ALFA = RADIUS**2 * B0 / CPSURF
        NITER = 1
        AGA = 0.
	BGA = 0.
	CGA = 0.
	DGA = 0.
	EGA = 0.
	FGA = 0.
	GGA = 0.
	HGA = 0.
	API = 0.
	BPI = 0.
	CPI = 0.
	DPI = 0.
	EPI = 0.
	FPI = 0.
	GPI = 0.
	HPI = 0.
	PEDGE=0.
	IPAI = 1
	IGAM = 1
      ELSE
C------------------------------------ PLASMA BOUNDARY IN FOURIER SERIES
        IF (ABS(ISHAPE).NE.2) THEN
          CALL FSHAPE(FR,MHARM)
        ELSE
          IF (IAS.EQ.1) THEN
            DO M=1,MFM
              FR(M) = FM(M)
            ENDDO
          ELSEIF (ISHAPE.GT.0) THEN
            DO M=1,MFM
              FR(2*M-1) = FM(M)
              FR(2*M) = 0.
            ENDDO
          ELSE
            DO M=1,MFM
              FR(2*M-1) = FM(2*M-1)
              FR(2*M)   = 0.
            ENDDO
          ENDIF
        ENDIF
	TXTOUT(1) = ' '
      ENDIF
C------------------------------------ ROTATE BOUNDARY

      CALL ROTATE_BND(FR,MFM,PAR1,IAS)

C------------------------------------ NORMALIZE PROFILES ---------------
      IF ((IGAM.EQ.7).OR.(IGAM.EQ.8).OR.
     >    (IGAM.EQ.4).OR.(IGAM.EQ.2).OR.(IGAM.EQ.11)) THEN
         FSCALE = DF2(1)
         DO I=1,NPTS
            DF2(I) = DF2(I) / FSCALE
         ENDDO
      ENDIF
      IF ((IPAI.EQ.7).OR.(IPAI.EQ.8)
     >     .OR.(IPAI.EQ.9).OR.(IPAI.EQ.11)  ) THEN
         PSCALE = DPR(1)
         DO I=1,NPTS
            DPR(I) = DPR(I) !/ PSCALE
         ENDDO
      ENDIF
      IF ((ICUR.EQ.2).OR.(ICUR.EQ.11).OR.(ICUR.EQ.13)) THEN
         CSCALE = ZJZ(1)
         DO I=1,NPTS
            ZJZ(I) = ZJZ(I) / CSCALE
         ENDDO
      ENDIF
C--------------------------------- INITIALIZE PRES AND GAM PROFILES ----
      CALL INIPRES
      CALL INIGAM
C------------------------------------------- LOOP OVER CURRENT PROFILE
      NRTMP = NR
      NPTMP = NP
      IF (ICUR.NE.0) THEN
         DO I = 1, NPTS
	    IF ((IGAM.GE.11).OR.(IGAM.EQ.2)) THEN
	       PSTEMP = (FLOAT(I-1)/FLOAT(NPTS-1))**2
            ELSE
	       PSTEMP = FLOAT(I-1)/FLOAT(NPTS-1)
	    ENDIF
            DF2(I) = DGDPSI(PSTEMP)
         ENDDO
         IF (IGAM.LE.4) THEN
            IF (IGAM.NE.2) IGAM = 4
         ELSEIF (IGAM.LT.11) THEN
            IGAM = 7
         ELSE
            IGAM = 11
         ENDIF
         NRTMP = NR
         NPTMP = NP
         NR = NRCUR
         NP = NPCUR
         NRM = NR
         NPM = NP
         IMESHTMP = IMESH
c         IMESH = 0
         ICURTMP = ICUR
         IARCTMP = IARC
         IARC = 0
         NMESHTMP = NMESH
      ENDIF
C---------------------------------------------- PROFILES -------------
      DO I=1,101
         SS = (I-1)*0.01
         PS = SS*SS
	 PSPLOT(I) = SS
	 DPPLOT(I) = DPDPSI(PS)
	 DGPLOT(I) = DGDPSI(PS)
	 IF (ICUR.NE.0)  ZJPLOT(I) = CURPHI(PS)
	 IF (ICUR.GE.90) QPLOT(I)  = QPRFL(PS)
      ENDDO
C------------------------------------ NAMELIST TO (PLOT) OUTPUT
      IF (NPL1.NE.0) THEN
        CALL LPLOT(4,3,1,PSPLOT,DPPLOT,101,1,'dP/d$y$',-5,
     >             's',1,' ',1)
        IF (IGAM.GT.4) THEN
	   CALL LPLOT(5,3,1,PSPLOT,DGPLOT,101,1,'FdF/d$y$',-6,
     >             's',1,' ',1)
        ELSE
           CALL LPLOT(5,3,1,PSPLOT,DGPLOT,101,1,'Gamma',5,
     >             's',1,' ',1)
        ENDIF
        IF (ICUR.EQ.99) THEN
           CALL LPLOT(6,3,1,PSPLOT,QPLOT,101,1,'q-profile',9,
     >             's',1,' ',1)
        ELSEIF (ICUR.NE.0) THEN
           CALL LPLOT(6,3,1,PSPLOT,ZJPLOT,101,1,'<J_phi>',-7,
     >             's',1,' ',1)
	ENDIF
      ENDIF
      CALL OUT(1)
      IF (ISOL.EQ.1) THEN
        WRITE(20,*)
        WRITE(20,*) '               SOLOVIEV EQUILIBRIUM '
        WRITE(20,2) A
        WRITE(20,3) B
        WRITE(20,4) ALFA
        WRITE(20,5) RADIUS
        WRITE(20,6) B0
        WRITE(20,7) CPSURF
	TXTOUT(1) = 'SOLOVIEV EQUILIBRIUM'
      ENDIF
      IF (NPR1.NE.0) THEN
        WRITE(20,*) '***********************************************'
        WRITE(20,*) '*        INPUT PROFILES :                     *'
	WRITE(20,*) '***********************************************'
	IF ((ICUR.EQ.0).AND.(IGAM.GT.4)) THEN
          WRITE(20,*) '*  PSI,       dP/dPSI,    FdF/dPSI            *'
	ELSEIF ((ICUR.EQ.0).AND.(IGAM.LE.4)) THEN
	  WRITE(20,*) '*  PSI,       dP/dPSI,     Gamma              *'
	ELSEIF ((ICUR.EQ.99).AND.(IGAM.LE.4)) THEN
          WRITE(20,*) '*  PSI,       dP/dPSI,     Gamma,         q   *'
	ELSEIF ((ICUR.EQ.99).AND.(IGAM.GT.4)) THEN
          WRITE(20,*) '*  PSI,       dP/dPSI,     FdF/dPSI,      q   *'
	ELSEIF ((ICUR.NE.0).AND.(IGAM.GT.4)) THEN
          WRITE(20,*) '*  PSI,       dP/dPSI,    FdF/dPSI,     J_phi *'
	ELSEIF ((ICUR.NE.0).AND.(IGAM.LE.4)) THEN
          WRITE(20,*) '*  PSI,       dP/dPSI,     Gamma,       J_phi *'
	ENDIF
	WRITE(20,*) '***********************************************'
        NR_PR = 21
        DO I=1,NR_PR
          PS = FLOAT(I-1)/FLOAT(NR_PR-1)
          PSPLOT(I) = PS
          DPPLOT(I) = DPDPSI(PS)
          DGPLOT(I) = DGDPSI(PS)
          ZJPLOT(I) = CURPHI(PS)
          IF (ICUR.EQ.0) THEN
            WRITE(20,622) PS,DPDPSI(PS),DGDPSI(PS)
          ELSEIF (ICUR.EQ.99) THEN
	    WRITE(20,622) PS,DPDPSI(PS),DGDPSI(PS),QPRFL(PS)
	  ELSE
            WRITE(20,622) PS,DPDPSI(PS),DGDPSI(PS),CURPHI(PS)
          ENDIF
        ENDDO
        WRITE(20,*)
      ENDIF
 622  FORMAT(4E12.4)

C----------------------------------- START OF MAIN LOOP 555 -----------

C------------------------------------ DEFINE NUMBERING OF ELEMENTS ----
 555  CALL ELMNO(NR,NP,NODENO)
C------------------------------------ INITIALIZE PSI ------------------
C------------------------------------ DEFINE INITIAL GRID X,Y ---------
      CALL INIGRID(XX,YY,PSI,NR,NP,FR,MHARM,IAS)

c      IF (ICUR.EQ.0) THEN
c         CALL DBLELM(XX,YY,PSI,NODENO,NR,NP,1)
c      ENDIF

C------------------------------------ PLOT THE INITIAL GRID
      IF ((NPL1.NE.0).AND.(ICUR.EQ.0)) THEN
        CALL LBLTOP(' ',1)
        CALL LBLBOT('HELENA EQUILIBRIUM Version 12',34)
        CALL PLOTGR(0,XX,YY,NR,NP,IAS)
      ENDIF


C----------------------------------- START OF LOOP OVER CURRENT/BETA---
      DO 888 NMG = 1,NMESH
         IF (ICUR.NE.0) THEN
            WRITE(20,*) '***************************************'
            WRITE(20,21) NMG
            WRITE(20,*) '***************************************'
   21       FORMAT(' * ITERATION CURRENT PROFILE,  NMG=',I3,'  *')
            NIT = 50
         ELSE
            WRITE(20,*) '***************************************'
            WRITE(20,22) NR,NP
            WRITE(20,*) '***************************************'
   22       FORMAT(' * FINAL ITERATION, NR=',I3,' NP=',I3,'      *')
            NIT = NITER
         ENDIF
      DO 10 NI = 1, NIT
         DO J=1,4*NR*NP
           QQ(J) = 0.
         ENDDO
C------------------------------------ INITIALIZE KK AND QQ TO ZERO -----
         DO I=1,4*NR*NP
            PSIOLD(I) = PSI(I)
         ENDDO
C------------------------------------ FORM MATRIX, NO CONDITIONS -------
        IF (ICUR.EQ.0) THEN
          CALL FORMKQ(XX,YY,PSI,NR,NP,QQ,A,B,C,EPS,IGAM,ISOL,NI,IAS)
        ELSE
          CALL FORMKQ(XX,YY,PSI,NR,NP,QQ,A,B,C,EPS,IGAM,ISOL,NI*NMG,IAS)
        ENDIF
C------------------------------------ SOLVE SET OF EQUATIONS -----------
        CALL SOLVE2(QQ,NR,NP,PSI,NI*NMG,IAS)
C------------------------------------ FIND MAGN. AXIS AND PSI AT AXIS --
        CALL FINDAXIS(XX,YY,NR,NP,PSAXIS,XAXIS,YAXIS,NAX,RAX,SAX,IAS)
        A = A / (1.-PSAXIS)
        IF (NDIAG.EQ.1)  WRITE(*,13) XAXIS,YAXIS
C------------------------------------ NORMALIZE FLUX TO ZERO ON AXIS ---
        CALL NORMAL(PSI,NR,NP,PSAXIS)
        IF (AMIX.NE.0.) THEN
          DO 14 I=1,4*NR*NP
            PSI(I) = (1.- AMIX)*PSI(I) + AMIX*PSIOLD(I)
   14     CONTINUE
          CALL FINDAXIS(XX,YY,NR,NP,PSAXIS,XAXIS,YAXIS,NAX,RAX,SAX,IAS)
          AOLD = A
          A = A / (1.-PSAXIS)
          CALL NORMAL(PSI,NR,NP,PSAXIS)
        ENDIF
        ERR1 = 0.
	ERR2 = 0.
	ERR3 = 0.
	ERR4 = 0.
        DO 16 I=1,NR*NP
          NBASE = 4*(I-1)+1
          ERR1 = ERR1 + ABS(PSI(NBASE)  - PSIOLD(NBASE))
          NBASE2 = 4*(I-1)+2
          ERR2 = ERR2 + ABS(PSI(NBASE2)  - PSIOLD(NBASE2))
          NBASE3 = 4*(I-1)+3
          ERR3 = ERR3 + ABS(PSI(NBASE3)  - PSIOLD(NBASE3))
          NBASE4 = 4*(I-1)+4
          ERR4 = ERR4 + ABS(PSI(NBASE4)  - PSIOLD(NBASE4))
   16   CONTINUE
        ERR1 = ERR1 / FLOAT(NR*NP)
        ERR2 = ERR2 / FLOAT(NR*NP)
        ERR3 = ERR3 / FLOAT(NR*NP)
        ERR4 = ERR4 / FLOAT(NR*NP)
        IF (ISOL.EQ.1) THEN
C         CALL ERRSOL(XX,YY,PSI,NR,NP)
        ENDIF
      IF (ERR1+ERR2+ERR3+ERR4.LT.ERRIT) GOTO 777
   10 CONTINUE
  777 CONTINUE
      WRITE(20,19) ERR1+ERR2+ERR3+ERR4,NI
      IF (NMG.NE.NMESH) THEN
        IF (ICUR.EQ.0) THEN
           NRM = 2*(NR-1) + 1
           NPM = 2*(NP-1) + 1
        ENDIF
        CALL REMESH(XX,YY,PSI,A,B,C,EPS,NR,NP,NRM,NPM,MESHNO,
     >              CX,CY,XAXIS,YAXIS,NAX,RAX,SAX,IGAM,IAS)
        IF (ICUR.NE.0) THEN
           DO 730 I = 1, NPTS
	      DF2OLD(I) = DF2(I)
 730       CONTINUE
C-------------------------------- update FF' based on current profile
           IF (ICUR.LT.90) THEN
             CALL FLXINT(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,
     >                   IAS,IAV)
           ELSE
C-------------------------------- update FF' based on q-profile
             CALL FLXINT2(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,
     >                   ISOL,CX,CY,IAS,AMPL,SUMDQ)
           ENDIF


           CALL BETAPOL(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                 IGAM,ISOL,BPL,IAS)
           BETAP_F=BPL
C       WRITE(*,'(A18,e12.4)') 'BETAP = ', BPL

           CALL BETAPOL_PED(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                 IGAM,ISOL,BPL,IAS,PSIPED)
C       WRITE(*,'(A18,e12.4)') 'BETAP_PED = ', BPL

           BETAP_PEDF =  BPL
           IF ((BETAP_PED.GE.0).AND.(ABS(BPL-BETAP_PED).GT.1e-3))  THEN
             IF (NMG.EQ.NBB) THEN
               B1 = B
               BP1 = BPL
               IF (B.LT.0.) THEN
                  IF (BPL .GT. BETAP_PED) THEN
                     B = 1.02*B1
                  ELSE
                     B = 0.98*B1
                  ENDIF
ccc                  B =  BPL /(BBB * BETAP_PED) * B1
               ELSE
                  B =   B1 / BPL * BETAP_PED
               ENDIF
               WRITE(*,'(A18,2e12.4,A4,e12.4)') 'first change of B',
     >                               B1,BPL,' >> ',B
             ELSE
                IF (MOD(NMG,NBB).EQ.0) THEN
                 BOLD = B
                 DELTAB = (B1-BOLD)/(BP1-BPL) * (BETAP_PED-BPL)
                 IF (ABS(DELTAB/BOLD) .gt. BBB) THEN
                   B = BOLD + BBB*DELTAB/ABS(DELTAB)*ABS(BOLD)
                 ELSE
                   B = BOLD + ABB*(B1-BOLD)/(BP1-BPL) * (BETAP_PED-BPL)
                 ENDIF

                 WRITE(*,'(A17,4f12.4,A4,f12.4)') 'next change of B',
     >                                BOLD,BPL,B1,BP1,' >> ',B
                 B1 = BOLD
                 BP1 = BPL

                 IF (B .GT. 200.) THEN
                   B =  - B
                   B1 = - BOLD
                 ENDIF

               ENDIF
             ENDIF
           ENDIF
           CALL INIGAM
           DIFDF2 = 0.
           DO I = 1, NPTS
              DIFDF2 = DIFDF2 + ABS(DF2(I) - DF2OLD(I))
           ENDDO
           DO I = 1, NR*NP
              DO K = 1, 4
                 XX(K,I) = XXOLD(K,I)
                 YY(K,I) = YYOLD(K,I)
                 PSI(4*(I-1)+K) = PSIOLD(4*(I-1)+K)
              ENDDO
           ENDDO
	   IF (ICUR.LT.90) THEN
C             WRITE(*,711) NMG,DIFDF2,B,BPL
	   ELSE
C             WRITE(*,711) NMG,DIFDF2,B,BPL,SUMDQ
           ENDIF

 711       FORMAT(' ITERATION : ',I3,' DIFF : ',1PE10.2,0P2f9.3,1PE10.2)
           IF (((DIFDF2.LT.ERRCUR).AND.(ABS(BPL-BETAP_PED).LT.0.001))
     >          .OR.(NMG.EQ.NMESH-1)) THEN
		WRITE(*,711) NMG,DIFDF2,B,BPL		
		GOTO 666
	   ENDIF         
	ENDIF
C------------------------------------ DEFINE NUMBERING OF ELEMENTS ----
        IF (ICUR.EQ.0) THEN
           CALL ELMNO(NR,NP,NODENO)
        ENDIF
      ENDIF
 888  CONTINUE
      GOTO 999
 666  NMESH = 1
      NR = NRTMP
      NP = NPTMP
      IMESH = IMESHTMP
      IARC = IARCTMP
      ICUR = 0

      GOTO 555
 999  CONTINUE

      DEALLOCATE(KKBIG)
C
C-------------------------------------- REMESH ON FLUX COORDINATES -----
C
      IF (NDIAG.EQ.1) write(*,*) ' final remesh'
      CALL REMESH(XX,YY,PSI,A,B,C,EPS,NR,NP,NRMAP,NPMAP,MESHNO,
     >            CX,CY,XAXIS,YAXIS,NAX,RAX,SAX,IGAM,IAS)
      IF (NDIAG.EQ.1) write(*,*) ' done final remesh'

C
C-------------------------------------- PLOT FINAL GRID ----------------
C
      IF (NPL1.EQ.1) THEN
        CALL LBLBOT('HELENA EQUILIBRIUM Version 12',34)
        CALL PLOTGR(1,XX,YY,NR,NP,IAS)
      ENDIF
C
C-------------------------------------- SCALE ALFA TO INPUT CURRENT ----
C

      IF (Q95.GT.0.) THEN
        IF (NDIAG.EQ.1) WRITE(*,*) '  profq',Q95
        CALL PROFQ(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,
     >                   CX,CY,DUMMY1,DUMMY2,Q95OUT,Q1OUT,IAS)
        IF (NDIAG.EQ.1) WRITE(*,*) Q1OUT,Q95,ALFA
        ALFA = ALFA * Q95/Q1OUT
        CALL CURRENT(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,CUR,IAS)
        XIAB = CUR
      ELSEIF (XIAB.GT.0.) THEN
        IF (NDIAG.EQ.1) WRITE(*,*) '  current',XIAB
        CALL CURRENT(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,CUR,IAS)
        ALFA = ABS(CUR / XIAB)
        IF (NDIAG.EQ.1) WRITE(*,*) '   alfa : ',ALFA
      ENDIF
C      CALL MOMENTS(XX,YY,PSI,NR,NP,XIAB,ALFA,EPS,IAS)
C--------------------------------------- CALCULATE SOME QUANTITIES ------
      WRITE(TXTOUT(13),13) XAXIS,YAXIS

      ICUR = ICURTMP
      CALL DIAGNO(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,IAS,XAXIS,
     >            YAXIS,IAV,ZVOL,ZVOLP,ZJAR)

!--------------------------------------- WRITE RESTART FILE ------------
      IF (NPR2.EQ.1) THEN
        CALL INOUT(XAXIS,CX,CY,XX,YY,PSI,FM,MFM,A,NRTMP,NPTMP,
     >             ICURTMP,NMESHTMP)
      ENDIF
C
C-------------------------------------- MAPPING OF FLUX COORDINATES ----
C
      CALL MAPPING(XX,YY,PSI,CX,CY,XAXIS,A)


c      DO K=1,80,2
c        node =K*NP+1
c        WRITE(20,*) NP,PSI(4*(NODE-1)+1)
c        WRITE(20,81) (XX(1,node+i-1),YY(1,node+i-1),I=1,NP)
c      ENDDO
c   81 format(2e16.8)

C
C--------------------------------------- CALC. BALLOONING STABILITY -----
C

      CALL MERCIER(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,
     >             DIMERC,DRMERC,HH,QPROF,DQPROF,GEONC,ZJPAR)
      IF (NQB.NE.0) THEN

        CALL HELBAL(ZVOL,ZVOLP,XAXIS)
      ENDIF

      IF (NPR1.NE.0) THEN
        CALL CIRCUL(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,
     >            FCIRC,B02AV,B0MAX,RAV)

        CALL WORLD(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,IAS,
     >             RVAC,BVAC,ZN0,ZEFF,RPE,ETAEI,FCIRC,QPROF,DQPROF,
     >             RAV,GEONC,ZJPAR,ZJAR,DIMERC,DRMERC,HH)

      ENDIF

      IF (NPR2 .EQ.1) THEN
!--------------------------------------- WRITE EQDSK FILE --------------
      CALL EQDSK(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,XIAB,RVAC,BVAC,
     >           QPROF)
      ENDIF
C--------------------------------------- CLOSE PLOTFILE ----------------
      IF (NPL1.EQ.1) THEN
        CALL FINPLT
      ENDIF


C
C--------------------------------------- FORMAT STATEMENTS --------------
    2 FORMAT('A   : ',E12.4)
    3 FORMAT('B   : ',E14.6)
    4 FORMAT('ALFA   : ',E14.6)
    5 FORMAT('RADIUS : ',E14.6)
    6 FORMAT('B0     : ',E14.6)
    7 FORMAT('CPSURF : ',E14.6)
   11 FORMAT('ITERATION     : NO  = ',I3)
   13 FORMAT('MAGNETIC AXIS : X   = ',E14.6,'  Y   = ',E14.6)
   15 FORMAT('PSI AT AXIS   : PSI = ',E14.6)
   17 FORMAT('AMPLITUDE     : A   = ',E14.6)
   19 FORMAT('    REST = ',1PE10.3,' ITERATIONS : ',I3)
      STOP 'FINISHED HELENA '
      END

***********************************************************************
*DECK INIVAL
      SUBROUTINE INIVAL
C-----------------------------------------------------------------------
C SUBROUTINE TO INITIALIZE THE NAMELIST INPUT VARIABLES
C-----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMDAT
      USE COMPRI
      USE COMPLO
      USE TOLERA
      USE COMMAP
      USE COMPIE
      USE COMPROF
      USE MESHAC

      RS(1,1) = -1.
      RS(1,2) = -1.
      RS(2,1) = -1.
      RS(2,2) = +1.
      RS(3,1) = +1.
      RS(3,2) = +1.
      RS(4,1) = +1.
      RS(4,2) = -1.
      IJ(1,1) = 0
      IJ(1,2) = 0
      IJ(2,1) = 1
      IJ(2,2) = 0
      IJ(3,1) = 0
      IJ(3,2) = 1
      IJ(4,1) = 1
      IJ(4,2) = 1
      XGAUSS(1) = -0.86113 63115 94053
      XGAUSS(2) = -0.33998 10435 84856
      XGAUSS(3) =  0.33998 10435 84856
      XGAUSS(4) =  0.86113 63115 94053
      WGAUSS(1) =  0.34785 48451 37454
      WGAUSS(2) =  0.65214 51548 62546
      WGAUSS(3) =  0.65214 51548 62546
      WGAUSS(4) =  0.34785 48451 37454

      PI = 2.*ASIN(1.)
      IAS = 0
      IAV = 0
      ELLIP = 1.0
      TRIA  = 0.0
      QUAD  = 0.0
      MHARM = 128
      ISHAPE = 1
      IMESH = 0
      AMESH = 0.
      BMESH = 0.
      ISOL = 0
      IARC = 0
      XR1 = 99.
      SIG1 = 99.
      AGA = -1.
      BGA = 0.
      CGA = 0.
      DGA = 0.
      EGA = 1.
      API = -1.
      BPI = 0.
      CPI = 0.
      DPI = 0.
      EPI = 1.
      ACUR=-1.
      BCUR=0.
      CCUR=0.
      DCUR=0.
      ECUR=1.
      ICUR = 0
      IPAI = 1
      IGAM = 1
      EPS = 0.1
      ALFA = 1.
      XIAB = -1.
      Q95 = -1.
      BETAP = -1.
      BETAP_PED = -1.
      COREP = 1.0
      PSIPED = 0.96
      B = 0.1
      C = 1.
      AMIX = 0.
      BMIX = 0.
      NR=5
      NP=9
      NRMAP = 5
      NPMAP = 9
      NCHI  = 8
      NITER = 1
      NMESH = 1
      NPTS = 100

      NBB = 4
      ABB = 1.
      BBB = 1.
      AMPL = 1.
      NPR1 = 1
      NPR2 = 0
      NPL1 = 1
      ERRCUR = 1.e-5
      ERRIT  = 1.e-8
      PSITOL = 1.E-8
      THTTOL = 1.E-8
      TOL    = 1.E-8
      NQB = 1

      RVAC = 1.0
      BVAC = 1.0
      ZEFF = 1.0
      ZN0  = 1.0
      RPE  = 0.5
      ETAEI = 1.0
      NROUT = 21

      RETURN
      END

************************************************************************
*DECK OUT
      SUBROUTINE OUT(I)
C-----------------------------------------------------------------------
C SUBROUTINE TO PRINT THE OUTPUT ON UNIT 20
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPRI
      USE COMPLO
      USE COMMAP
      IF (I.EQ.1) THEN
        WRITE(20,*) '*****************************************'
        WRITE(20,10)
        WRITE(20,*) '*****************************************'
        WRITE(20,*) '*          VERSION 12 (30-7-1999)       *'
        WRITE(20,*) '*****************************************'

	WRITE(20,*)
	WRITE(20,*)
        WRITE(20,20) ELLIP,TRIA,QUAD
        WRITE(20,21) MHARM,ISHAPE,ISOL,IAS
	WRITE(20,*)
        WRITE(20,30) IGAM,AGA,BGA,CGA,DGA
	WRITE(20,31) EGA,FGA,GGA,HGA
        WRITE(20,32) IPAI,API,BPI,CPI,DPI
	WRITE(20,33) EPI,FPI,GPI,HPI
	WRITE(20,37) COREP, PSIPED
	IF (ICUR.NE.0) THEN
	  WRITE(20,34) ICUR,ACUR,BCUR,CCUR,DCUR
	  WRITE(20,35) ECUR,FCUR,GCUR,HCUR
	  WRITE(20,36) IAV
	ENDIF
	WRITE(20,*)
        WRITE(20,40) EPS,ALFA,B,C
	WRITE(20,41) XIAB,BETAP,BETAP_PED
	WRITE(20,*)
        WRITE(20,50) NR,NP,NRMAP,NPMAP,NCHI,NITER
	WRITE(20,*)
	IF (ICUR.NE.0) THEN
	  WRITE(20,51) NRCUR,NPCUR,NMESH,NPTS,ERRCUR
	ENDIF
	WRITE(20,*)
        WRITE(20,60) NPR1
	WRITE(20,*)
        WRITE(20,70) NPL1
	WRITE(20,*)
	WRITE(20,*)
	IF (NPL1.NE.0) THEN
c	  WRITE(TXTOUT(1),10)
          WRITE(TXTOUT(2),20) ELLIP,TRIA,QUAD
          WRITE(TXTOUT(3),21) MHARM,ISHAPE,ISOL,IAS
          WRITE(TXTOUT(4),30) IGAM,AGA,BGA,CGA,DGA
	  WRITE(TXTOUT(5),31) EGA,FGA,GGA,HGA
cWO output to screen
          WRITE(TXTOUT(6),32) IPAI,API,BPI,CPI,DPI
		  WRITE(*,32) IPAI,API,BPI,CPI,DPI
	  WRITE(TXTOUT(7),33) EPI,FPI,GPI,HPI
		WRITE(*,33) EPI,FPI,GPI,HPI 
	  WRITE(TXTOUT(8),37) COREP, PSIPED
		WRITE(*,37) COREP, PSIPED
	  IF (ICUR.NE.0) THEN
            WRITE(TXTOUT(9),34) ICUR,ACUR,BCUR,CCUR,DCUR
		WRITE(*,34) ICUR,ACUR,BCUR,CCUR,DCUR
	    WRITE(TXTOUT(10),35) ECUR,FCUR,GCUR,HCUR
		WRITE(*,35) ECUR,FCUR,GCUR,HCUR
	  ENDIF
          WRITE(TXTOUT(11),40) EPS,ALFA,B,C
		WRITE(*,40) EPS,ALFA,B,C
	  WRITE(TXTOUT(12),41) XIAB,BETAP,BETAP_PED
		WRITE(*,41) XIAB,BETAP,BETAP_PED
          WRITE(TXTOUT(13),50) NR,NP,NRMAP,NPMAP,NCHI,NITER
		  WRITE(*,50) NR,NP,NRMAP,NPMAP,NCHI,NITER
	  IF (ICUR.NE.0) THEN
	    WRITE(TXTOUT(14),51) NRCUR,NPCUR,NMESH,NPTS,ERRCUR
	  ENDIF
          WRITE(TXTOUT(15),60) NPR1
          WRITE(TXTOUT(16),70) NPL1
	  IX=100
	  IY=722
	  TXTOUT(16) = ' HELENA NAMELIST INPUT PARAMETERS : '
	  CALL DLCH(IX,IY,TXTOUT(16),35,-2)
	  IY = IY-30
	  DO J=1,16
	    CALL DLCH(IX,IY,TXTOUT(J),100,-1)
	    IY=IY-18
     	  ENDDO
	ENDIF
      ENDIF
   10 FORMAT(' *  HELENA : GRAD-SHAFRANOV EQUILIBRIUM  *')
   20 FORMAT(' $SHAPE    ELLIP =',F7.3,', TRIA  =',F7.3,
     >           ', QUAD  =',F7.3)
   21 FORMAT('           MHARM =',I3,', ISHAPE =',I2,', ISOL =',I2,
     >       ', IAS =',I2)
   30 FORMAT(' $PROFILE  IGAM =',I2,' AGA =',F7.3,', BGA =',F7.3,
     >       ' CGA =',F7.3,' DGA =',F7.3,',')
   31 FORMAT(9X,'           EGA =',F7.3,', FGA =',F7.3,', GGA =',F7.3,
     >       ', HGA =',F7.3)
   32 FORMAT('           IPAI =',I2,' API =',F7.3,', BPI =',F7.3,
     >       ' CPI =',F7.3,' DPI =',F7.3,',')
   33 FORMAT(9X'           EPI =',F7.3,', FPI =',F7.3,', GPI =',F7.3,
     >       ', HPI =',F7.3)
   34 FORMAT('           ICUR =',I2,' ACUR=',F7.3,', BCUR=',F7.3,
     >       ' CCUR=',F7.3,' DCUR=',F7.3,',')
   35 FORMAT('           ECUR=',F7.3,', FCUR=',F7.3,', GCUR=',F7.3,
     >       ', HCUR=',F7.3,',')
   36 FORMAT('           AVERAGING : IAV =',I3,',')
   37 FORMAT('           COREP =',F7.3", PSIPED =",F7.3)
   40 FORMAT(' $PHYS     EPS  =',F7.3,', ALFA =',F7.3,
     >           ', B =',F7.3,', C =',F7.3,',')
   41 FORMAT('           XIAB =',F7.3,', BETAP = ',F7.3,
     >           ', BETAP_PED =', F7.3)
   50 FORMAT(' $NUM      NR =',I3,', NP =',I3,', NRMAP =',I3,
     >       ', NPMAP =',I3,', NCHI = ',I4,', NITER = ',I3)
   51 FORMAT('           NRCUR =',I3,', NPCUR =',I3,', NMESH =',I3,
     >       ', NPTS =',I3,', ERRCUR =',1PE10.2)
   60 FORMAT(' $PRI      NPR1  =',I2)
   70 FORMAT(' $PLOT     NPL1  =',I2)
      RETURN
      END

************************************************************************
*DECK INOUT
      SUBROUTINE INOUT(XAXIS,CX,CY,XX,YY,PSI,FM,MFM,A,NRTMP,NPTMP,
     >                 ICURTMP,NMESHTMP)
C-----------------------------------------------------------------------
C SUBROUTINE TO PRINT THE NAMELIST VALUES UNIT 30 USED FOR RESTARTING
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMDAT
      USE MESHAC
      USE COMPROF
      USE COMPRI
      USE COMPLO
      USE COMMAP
      REAL XX(4,*),YY(4,*),PSI(*),QPRF(NRMMAX),FM(*)

C------------------------------------------- DEFINE INPUT NAMELISTS ----
C      NAMELIST/SHAPE/  ELLIP,TRIA,QUAD,MHARM,ISHAPE,ISOL,FM,MFM,
C     >                 XR1,SIG1,
C     >                 PAR1,PAR2,PAR3,PAR4,AMESH,BMESH,CMESH,IMESH,
C     >                 IAS,IARC
C      NAMELIST/PROFILE/AGA,BGA,CGA,DGA,EGA,FGA,GGA,HGA,
C     >                 API,BPI,CPI,DPI,EPI,FPI,GPI,HPI,
C     >                 ACUR,BCUR,CCUR,DCUR,ECUR,FCUR,GCUR,HCUR,
C     >                 DPR,DF2,ZJZ,QIN,
C     >                 ICUR,IGAM,IPAI,NPTS,IAV
C      NAMELIST/PHYS/   EPS,ALFA,B,C,XIAB,Q95,BETAP,
C     >                 RVAC,BVAC,ZEFF,ZN0,RPE,ETAEI
C      NAMELIST/NUM/    NR,NP,NRMAP,NPMAP,NCHI,NITER,NMESH,AMIX,BMIX,
C     >                 ERRIT,ERRCUR,NRCUR,NPCUR,
C     >                 NBB,ABB,BBB,AMPL
C      NAMELIST/PRI/    NPR1,NPR2,NROUT,NDIAG
C      NAMELIST/PLOT/   NPL1
C      NAMELIST/BALL/   QB1,QB2,NQB


      CALL PROFQ(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,
     >                   CX,CY,QPRF,ZJZ,Q95OUT,Q1OUT,IAS)

      CSCALE = ZJZ(1)
      DO I=1,NR
        ZJZ(I) = ZJZ(I) / CSCALE
cWO Output ZJZ values
	    WRITE(78,*) ZJZ(I)
		
      ENDDO

	  
	  
      ID=30
      WRITE(ID,*)  ' HELENA INPUT-OUTPUT FILE'
      WRITE(ID,*)
      WRITE(ID,20) ELLIP,TRIA,QUAD
      WRITE(ID,21) MFM,MHARM,ISHAPE,ISOL
      WRITE(ID,22) IMESH,AMESH,BMESH
      WRITE(ID,23) XR1,SIG1,IAS,IARC
      WRITE(ID,24) PAR1,PAR2,PAR3,PAR4
      IF (ISHAPE.EQ.2) THEN
        DO M=1,MFM
          WRITE(ID,25) M,FM(M)
	ENDDO
      ENDIF

      IF (IPAI .EQ. 7) THEN
        IPAI=11
	EPI = SQRT(EPI)
      ENDIF
      IF (IGAM .EQ. 7) THEN
        IGAM = 11
        EGA = SQRT(EGA)
      ENDIF

      WRITE(ID,*)  '$END'
      WRITE(ID,30) IGAM,AGA,BGA,CGA
      WRITE(ID,31) DGA,EGA,FGA
      WRITE(ID,32) GGA,HGA
      WRITE(ID,33) IPAI,API,BPI,CPI
      WRITE(ID,34) DPI,EPI,FPI
      WRITE(ID,35) GPI,HPI
      WRITE(ID,36) ICURTMP,ACUR,BCUR,CCUR
      WRITE(ID,37) DCUR,ECUR,FCUR
      WRITE(ID,38) GCUR,HCUR
cWO
      WRITE(ID,39) NROUT,COREP,PSIPED

      ICUR = 11

      DO I=1,NROUT
        SS   = FLOAT(I-1)/FLOAT(NROUT-1)
        ZPSI = SS * SS
        DPRI = DPDPSI(ZPSI)
        DF2I = DGDPSI(ZPSI)
        ZJZI = CURPHI(ZPSI)
	QOUT = QPRFL(ZPSI)
!         WRITE(ID,392)  I,DPRI,I,DF2I,I,ZJZI,I,QIN(I)
!         WRITE(ID,392)  I,DPRI,I,DF2I,I,ZJZI,I,QOUT
         WRITE(ID,392)  I,DPRI,I,DF2I,I,ZJZI,I,QPRF(I)
         WRITE(67,392)  I,DPRI,I,DF2I,I,ZJZI,I,QPRF(I)

      ENDDO

      WRITE(ID,*)  '$END'

      WRITE(ID,40) EPS,ALFA,B,C
      WRITE(ID,41) XIAB,BETAP,BETAP_PED
      WRITE(ID,*)  '$END'

      WRITE(ID,50) NRTMP,NPTMP,NRMAP,NPMAP,NCHI,NITER
      WRITE(ID,55) NRCUR,NPCUR,NMESHTMP
      WRITE(ID,56) ERRCUR,AMIX
      WRITE(ID,57) ABB,BBB,NBB,AMPL
      WRITE(ID,*)  '$END'

      WRITE(ID,60) NPR1
      WRITE(ID,70) NPL1

      WRITE(ID,71) NQB

   20 FORMAT(/1X,'$SHAPE    ELLIP =',E12.4,', TRIA  =',E12.4,
     >           ', QUAD  =',E12.4,',')
   21 FORMAT(11X,'MFM = ',I4,', MHARM =',I4,', ISHAPE =',I2,
     >           ', ISOL =',I2,',')
   22 FORMAT(11X,'IMESH =',I3,', AMESH=',E12.4,', BMESH=',E12.4,',')
   23 FORMAT(11X,'XR1 = ',F8.4,', SIG1 = ',F8.4,', IAS = ',I3
     >          ', IARC = ',I3)
   24 FORMAT(11X,'PAR1 = ',F8.4,', PAR2 = ',F8.4,', PAR3 = ',F8.4
     >          ', PAR4 = ',F8.4)
   25 FORMAT(11X,' FM(',I3,') = ',1PE12.4,',')
   30 FORMAT(/1X,'$PROFILE  IGAM =',I2,', AGA =',1PE12.4,', BGA =',
     >        E12.4,', CGA =',E12.4,',')
   31 FORMAT(11X,'DGA =',1PE12.4,', EGA ='E12.4,', FGA =',E12.4,',')
   32 FORMAT(11X,'GGA =',1PE12.4,', HGA =',E12.4,',')
   33 FORMAT(1X,'          IPAI =',I2,', API =',1PE12.4,', BPI =',
     >            E12.4,', CPI =',E12.4,',')
   34 FORMAT(11X,'DPI =',1PE12.4,', EPI =',E12.4,', FPI =',E12.4,',')
   35 FORMAT(11X,'GPI =',1PE12.4,', HPI =',E12.4,',')
   36 FORMAT(1X,'          ICUR =',I2,', ACUR =',1PE12.4,', BCUR =',
     >            E12.4,', CCUR =',E12.4,',')
   37 FORMAT(11X,'DCUR =',1PE12.4,', ECUR =',E12.4,', FCUR =',E12.4,',')
   38 FORMAT(11X,'GCUR =',1PE12.4,', HCUR =',E12.4,',')
cWO
   39 FORMAT(11X,'NPTS = ',I4,', COREP =',E12.4,', PSIPED =',E12.4,",")
  391 FORMAT(11x,'DPR(',I3,')=',1PE13.5,', DF2(',I3,')=',E13.5,
     >          ', QIN(',I3,')=',E13.5,',')
  392 FORMAT(3x,'DPR(',I3,')=',1PE13.5,', DF2(',I3,')=',E13.5,
     >          ', ZJZ(',I3,')=',E13.5,', QIN(',I3,')=',E13.5,',')
   40 FORMAT(/1X,'$PHYS     EPS  =',F10.6,', ALFA =',F10.6,
     >           ', B =',E12.4,', C =',F10.6,',')
   41 FORMAT(11x,'XIAB =',F10.6,', BETAP =',F10.6,
     > 		', BETAP_PED =',F10.6,',')
   50 FORMAT(/1X,'$NUM      NR =',I5,', NP =',I5,', NRMAP =',I5,
     >           ', NPMAP =',I5,', NCHI =',I5,', NITER = ',I3,',')
   55 FORMAT(11X,'NRCUR =',I5,', NPCUR =',I5,', NMESH = ',I4,',')
   56 FORMAT(11X,'ERRCUR =',1PE12.4,', AMIX = ',0PF8.4',')
   57 FORMAT(11X,'ABB = ',F8.4,', BBB = ',F8.4,', NBB = ',I4,', AMPL =',
     >       F8.4)
   60 FORMAT(/1X,'$PRI      NPR1  =',I2,' $END ')
   70 FORMAT(/1X,'$PLOT     NPL1  =',I2,' $END ')
   71 FORMAT(/1X,'$BALL     NQB = ',I4,' $END')
      RETURN
      END

************************************************************************
*DECK GAUSS
      SUBROUTINE GAUSS
C-----------------------------------------------------------------------
C CALCULATE THE SIXTHEEN FUNCTIONS AT THE SIXTEEN GAUSSIAN POINTS
C-----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT

      DO 10 IRS=1,4
        R0 = RS(IRS,1)
        S0 = RS(IRS,2)
        DO 20 IM=1,4
          MI = IJ(IM,1)
          MJ = IJ(IM,2)
          DO 30 NGR=1,4
          R = XGAUSS(NGR)
            DO 40 NGS=1,4
              S  = XGAUSS(NGS)
              CALL CUBICH(MI,MJ,R0,S0,R,S,H(IM,IRS,NGR,NGS),
     >                    HR(IM,IRS,NGR,NGS),HS(IM,IRS,NGR,NGS),
     >                    HRS(IM,IRS,NGR,NGS),DUMRR,DUMSS)
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      RETURN
      END

************************************************************************
*DECK CUBICH
      SUBROUTINE CUBICH(I,J,R0,S0,R,S,H,HR,HS,HRS,HRR,HSS)

C------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE VALUE OF THE CUBIC POLYNOMIALS AND
C THE DERIVATIVES OF THE CORNER MARKED BY R0,S0 AT THE POSITION R,S
C------------------------------------------------------------------
      REAL H,HR,HS,HRS,HI,HRI,HJ,HSJ,HRR,HSS,HRRI,HSSJ

      IF (I.EQ.0) THEN
        HI = - (R+R0)**2 * (R*R0-2.) / 4.
        HRI = - (R+R0)*(R*R0-2.)/2. - R0*(R+R0)**2 / 4.
        HRRI = - 1.5 * R * R0
      ELSE
        HI = + R0 * (R+R0)**2 * (R*R0 - 1.) / 4.
        HRI = + R0*(R+R0)*(R*R0-1.)/2. + (R+R0)**2 /4.
        HRRI = 1.5*R + .5*R0
      ENDIF
      IF (J.EQ.0) THEN
        HJ = - (S+S0)**2 * (S*S0-2.) / 4.
        HSJ = - (S+S0)*(S*S0-2.)/2. - S0*(S+S0)**2 / 4.
        HSSJ = - 1.5 * S * S0
      ELSE
        HJ = + S0 * (S+S0)**2 * (S*S0 - 1.) / 4.
        HSJ = + S0*(S+S0)*(S*S0-1.)/2. + (S+S0)**2 / 4.
        HSSJ = 1.5*S + .5*S0
      ENDIF
      H = HI * HJ
      HR = HRI * HJ
      HS = HI * HSJ
      HRS = HRI * HSJ
      HRR = HRRI * HJ
      HSS = HI   * HSSJ
      RETURN
      END

************************************************************************
*DECK CUBICH1
      SUBROUTINE CUBICH1(I,J,R0,S0,R,S,H)
C------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE VALUE OF THE CUBIC POLYNOMIALS AND
C THE DERIVATIVES OF THE CORNER MARKED BY R0,S0 AT THE POSITION R,S
C------------------------------------------------------------------
      REAL H,HI,HJ

      IF (I.EQ.0) THEN
        HI = - (R+R0)**2 * (R*R0-2.) / 4.
      ELSE
        HI = + R0 * (R+R0)**2 * (R*R0 - 1.) / 4.
      ENDIF
      IF (J.EQ.0) THEN
        HJ = - (S+S0)**2 * (S*S0-2.) / 4.
      ELSE
        HJ = + S0 * (S+S0)**2 * (S*S0 - 1.) / 4.
      ENDIF
      H = HI * HJ
      RETURN
      END

************************************************************************
*DECK CUBICH2
      SUBROUTINE CUBICH2(I,J,R0,S0,R,S,H,HR,HS)
C------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE VALUE OF THE CUBIC POLYNOMIALS AND
C THE DERIVATIVES OF THE CORNER MARKED BY R0,S0 AT THE POSITION R,S
C------------------------------------------------------------------
      REAL H,HR,HS,HI,HRI,HJ,HSJ

      IF (I.EQ.0) THEN
        HI = - (R+R0)**2 * (R*R0-2.) / 4.
        HRI = - (R+R0)*(R*R0-2.)/2. - R0*(R+R0)**2 / 4.
      ELSE
        HI = + R0 * (R+R0)**2 * (R*R0 - 1.) / 4.
        HRI = + R0*(R+R0)*(R*R0-1.)/2. + (R+R0)**2 /4.
      ENDIF
      IF (J.EQ.0) THEN
        HJ = - (S+S0)**2 * (S*S0-2.) / 4.
        HSJ = - (S+S0)*(S*S0-2.)/2. - S0*(S+S0)**2 / 4.
      ELSE
        HJ = + S0 * (S+S0)**2 * (S*S0 - 1.) / 4.
        HSJ = + S0*(S+S0)*(S*S0-1.)/2. + (S+S0)**2 / 4.
      ENDIF
      H = HI * HJ
      HR = HRI * HJ
      HS = HI * HSJ
      RETURN
      END

************************************************************************
*DECK SOLSHP
      SUBROUTINE SOLSHP(FR,MF)
C-----------------------------------------------------------------------
C THE SHAPE R=SQRT(X**2+Y**2)=FR(THETA) OF THE PLASMA CROSS-
C SECTION IN THE Z-PLANE IS DERIVED FROM A SHAPE CORRESPONDING
C TO THE PSI=1 SURFACE OF THE SOLOVIEV EQUILIBRIUM.
C-----------------------------------------------------------------------
      USE COMDAT
      USE COMANG
      USE COMPIE

      REAL FR(*)
      EXTERNAL FSOL

      DO 10 J=1,MF
        ANGLE= 2*PI*(J-1.)/ MF
        F1=.5*(ELLIP+(1.-ELLIP)*COS(2.*ANGLE))
        FUN1=FSOL(F1)
        DO 20 N=1,20
          F2=F1+.05*N
          FUN2=FSOL(F2)
          IF(FUN2.GT.0.) GOTO 30
   20   CONTINUE
   30   CALL ZERO(F1,FUN1,F2,FUN2,FSOL,1.E-8,FR(J),FUN,IZERO,0)
   10 CONTINUE
CCC   CALL PRARR1('FR  :',FR,MF,203)
      CALL RFT2(FR,MF,1)
      DO 50 M=1,MF
        FR(M) = 2 * FR(M) / REAL(MF)
   50 CONTINUE
      RETURN
      END

***********************************************************************
*DECK FSOL
      FUNCTION FSOL(F)
C-----------------------------------------------------------------------
C EXTERNAL FUNCTION FSOL=PSI(F)-1 NEEDED FOR THE COMPUTATION
C OF THE RADIUS FR(ANGLE) OF THE PSI=1 SURFACE OF THE SOLOVIEV
C EQUILIBRIUM.
C-----------------------------------------------------------------------
      USE COMDAT
      USE COMANG

      X=F*COS(ANGLE)
      Y=F*SIN(ANGLE)
      PAR4 = QUAD
      PAR3 = TRIA
      PAR2 = EPS
      PAR1 = ELLIP
      FSOL= (1. -.25*PAR2**2) * (1.+PAR2*X)**2 * (Y/PAR1)**2
     >    + (X - .5*PAR2*(1. - X*X))**2
     >    + (4.-PAR2**2) * (X -.5*PAR2*(1-X*X))*PAR4*Y/(2*PAR2) - 1.
      RETURN
      END

************************************************************************
*DECK FSHAPE
      SUBROUTINE FSHAPE(FR,MF)
C-----------------------------------------------------------------------
C  THE SHAPE R=SQRT(X**2+Y**2)=FR(THETA) OF THE PLASMA CROSS-SECTION IN
C  Z-PLANE IS COMPUTED FROM A FORMULA WHICH GIVES ELLIPSES, D-SHAPES ETC

C     X=A*COS(GAMMA+C*SIN(GAMMA)+D*SIN(2.*GAMMA)),
C     Y=B*SIN(GAMMA),
C  WHERE THETA IS FOUND FROM GAMMA BY INVERSION.
C  THE FOURIER COEFFICIENTS FRFNUL,FRF(M) OF FR(J) ARE ALSO CALCULATED.
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPIE
      REAL FR(*), THETA(2*MBMAX+2), GAMMA(2*MBMAX+2)
      REAL XV(2*MBMAX+2),YV(2*MBMAX+2)
C
      IF (ISHAPE.EQ.3) THEN
        WRITE(20,*)
        WRITE(20,*) ' SHAPE OF BOUNDARY : '
        WRITE(20,*)
        DANG = PI*PAR3/2.
        DO 200 J=1,MF
          ANGLE = 2*PI*(J-1.)/REAL(MF)
          FR(J) = 1.
          IF ((ANGLE.GT.(PAR2-DANG)).AND.(ANGLE.LT.(PAR2+DANG))) THEN
            FR(J) = FR(J) + PAR1 * COS((ANGLE-PAR2)/PAR3)**PAR4
          ENDIF
          WRITE(20,*) ANGLE,FR(J)
  200   CONTINUE
        WRITE(20,*)
      ELSEIF (ISHAPE.EQ.4) THEN
        WRITE(20,*)
        WRITE(20,*) ' SHAPE OF BOUNDARY : '
        WRITE(20,*)
        DANG = PAR3 /SQRT(3.)
        DO 210 J=1,MF
          ANGLE = 2*PI*(J-1.)/REAL(MF)
          FR(J) = 1.
          IF ((ANGLE.GT.(PAR2-DANG)).AND.(ANGLE.LT.(PAR2+DANG))) THEN
            FR(J) = FR(J) + PAR1 * (2*PAR3/3**1.5 * PAR2**3 -
     >                         ABS((ANGLE-PAR2)*(ANGLE-PAR2-PAR3)*
     >                                          (ANGLE-PAR2+PAR3)))
          ENDIF
          WRITE(20,*) ANGLE,FR(J)
  210   CONTINUE
        WRITE(20,*)
      ELSE
C------------------------------------------ THETA(GAMMA(J)) ------------
        DO 10 J=1,MF
          GA = 2*PI*(J-1.)/REAL(MF)
          XJO= COS(GA + TRIA*SIN(GA) + QUAD*SIN(2*GA))
          YJO= ELLIP * SIN(GA + PAR4*COS(GA))
          XJ = XJO * COS(PAR1) - YJO*SIN(PAR1)
          YJ = XJO * SIN(PAR1) + YJO*COS(PAR1) + PAR2
          THETA(J) = ATAN2(YJ,XJ)
   10   CONTINUE
C-------------- INVERSION OF THETA(GAMMA(J)) TO GAMMA(THETA(J)) ------  --
        CALL GRID2NV(THETA,GAMMA,MF,1.E-8,IGRINV,0)
        DO 30 J=1,MF
          GAMM=GAMMA(J)
          XJO= COS(GAMM + TRIA*SIN(GAMM) + QUAD*SIN(2*GAMM))
          YJO= ELLIP * SIN(GAMM + PAR4*COS(GAMM))
          XV(J) = XJO * COS(PAR1) - YJO*SIN(PAR1)
          YV(J) = XJO * SIN(PAR1) + YJO*COS(PAR1) + PAR2
   30   CONTINUE
        XEAST = XV(1)
        XWEST = XV(MF/2+1)
        X0 = (XEAST + XWEST)/2.
        XL = (XEAST - XWEST)/2.
        DO 35 J=1,MF
          XN = (XV(J) - X0) / XL
          YN =  YV(J) / XL
          FR(J) = SQRT(XN**2 + YN**2)
   35   CONTINUE
      ENDIF
C------------------- FOURIER COEFFICIENTS FRFNUL AND FRF(M) OF FR(J).
      CALL RFT2(FR,MF,1)
      DO 50 M=1,MF
        FR(M) = 2. * FR(M) / REAL(MF)
   50 CONTINUE
ccc      CALL PRARR1('FR : ',FR,MF,203)
      RETURN
      END

************************************************************************
*DECK ROTATE_BND
      SUBROUTINE ROTATE_BND(FR,MF,ANGLE,IAS)
C-------------------------------------------------------------------
C subroutine rotates the plasma boundary by the angle angle
C-------------------------------------------------------------------
      REAL FR(*),ANGLE
      INTEGER MF,IAS
      REAL RBND(MF+2)

      PI = 2.*asin(1.)

      IF (IAS .EQ. 1) THEN

        DO J=1,MF

          THTJ = 2.*PI * FLOAT(J-1)/FLOAT(MF) + ANGLE
          XP   = FR(1) * COS(THTJ) / 2.
          YP   = FR(1) * SIN(THTJ) / 2.

          DO M=2,MF/2
            RM = FR(2*M-1) * COS((M-1)*THTJ) - FR(2*M) * SIN((M-1)*THTJ)
            XP = XP + RM * COS(THTJ)
            YP = YP + RM * SIN(THTJ)
          ENDDO

          RBND(J) = SQRT(XP*XP + YP*YP)

        ENDDO

        CALL RFT2(RBND,MF,1)

        DO M=1,MF
          FR(M) = 2. * RBND(M) / REAL(MF)
        ENDDO

      ENDIF

      RETURN

      END


************************************************************************
*DECK ELMNO
      SUBROUTINE ELMNO(NR,NP,NODENO)
C-------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE FOUR NODENUMBERS OF EVERY ELEMENT
C-------------------------------------------------------------------
      USE PARAM
      INTEGER NODENO(MAXMNODE,4)
      NO = 0
      DO 10 N=1,NR-1
        DO 20 M=1,NP-1
          NO = NO + 1
          NODENO(NO,1) = (N-1)*NP + M
          NODENO(NO,2) = NODENO(NO,1) + 1
          NODENO(NO,3) = NODENO(NO,2) + NP
          NODENO(NO,4) = NODENO(NO,1) + NP
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
************************************************************************
*DECK MESHAC2
      SUBROUTINE MESHAC2(NR,SG,DSG,DDSG,XR1,SIG1)
C-----------------------------------------------------------------------
C subroutine to construct non-equidistant radial mesh in Helena
C-----------------------------------------------------------------------
      PARAMETER (NMAX=1001)
      REAL SG(*),DSG(*),DDSG(*)
      REAL S1(NMAX),F1(NMAX),F2(NMAX),F3(NMAX),F4(NMAX),FSUM(NMAX)
      REAL ABLTG(3)

C--------------------------------------- set parameters of gaussians
      BGF  = 0.3
      XR2  = 9999.
      SIG2 = 1.
      FACT = 1.
C--------------------------------------- integrate gaussian
      DSI = 1. / FLOAT(NMAX-1)
      S1(1) = 0.
      FSUM(1)  = 0.
      SUM = 0.
      FINT2 = FGAUS(S1(1),BGF,XR1,XR2,SIG1,SIG2,FACT,DFG)
      DO I=2,NMAX
        S1(I) = FLOAT(I-1) * DSI
        FINT1 = FINT2
        FINT2 = FGAUS(S1(I),BGF,XR1,XR2,SIG1,SIG2,FACT,DFG)
        SUM = SUM + (FINT1+FINT2)/2. * DSI
        FSUM(I) = SUM
      ENDDO

      DO I=1,NMAX-1
        FSUM(I) = FSUM(I)/FSUM(NMAX)
      ENDDO
      FSUM(NMAX) = 1.
      ALFA = 0.
      BETA = 0.
      TYP = 2
      CALL SPLINE(NMAX,FSUM,S1,ALFA,BETA,TYP,F1,F2,F3,F4)

      SG(1)   = 0.
      DSG(1)  = F2(1)
      DDSG(1) = F3(1)
      DO I=2,NR-1
        FI = FLOAT(I-1)/FLOAT(NR-1)
        SG(I)   = SPWERT(NMAX,FI,F1,F2,F3,F4,FSUM,ABLTG)
        DSG(I)  = ABLTG(1)
        DDSG(I) = ABLTG(2)
      ENDDO
      SG(NR)   = 1.
      DSG(NR)  = F2(NMAX)
      DDSG(NR) = F3(NMAX)
      RETURN
      END

************************************************************************
*DECK FGAUS
      FUNCTION FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT,DFGAUSS)
C-----------------------------------------------------------------------
C     BGF + (1 - BGF) * (GAUSS1 + FACT * GAUSS2) / FACT
C-----------------------------------------------------------------------
C
      ZNORM1 = 0.39894 / SIG1
      ZNORM2 = 0.39894 / SIG2
      ZEX1   = -0.5 * (ZS - XR1)**2 / SIG1**2
      ZEX2   = -0.5 * (ZS - XR2)**2 / SIG2**2
      DEX1   = -(ZS-XR1)/SIG1**2
      DEX2   = -(ZS-XR2)/SIG2**2

      F1     = ZNORM1 * EXP(ZEX1)
      F2     = ZNORM2 * EXP(ZEX2)
      DF1    = ZNORM1 * DEX1 * EXP(ZEX1)
      DF2    = ZNORM2 * DEX2 * EXP(ZEX2)
C
      FGAUS  = BGF + (1.0 - BGF) * (F1 + FACT * F2) / FACT
      DFGAUSS = (1.0-BGF) * (DF1 + FACT * DF2) / FACT
C
      RETURN
      END

************************************************************************
*DECK ARCLENGTH
      SUBROUTINE ARCLENGTH(FR,MF,THETA,DTC,NP,IAS,WR,WS)
C-----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE ARCLENGTH OF THE PLASMA BOUNDARY
C ROUTINE RETURNS THE THETA VALUES AT EQUIDISTANT ARCLENGTH
C PARAMETERS :
C                FR    : FOURIER SERIES OF BOUNDARY (INPUT)
C                MF    : NUMBER OF HARMONICS IN FR  (INPUT)
C                THETA : RESULTING VALUEs OF THETA
C                DTC   : THE DERIVATIVE OF THETA TO EQUIDISTANT ANGLE
C                NP    : NUMBER OF POINTS IN THETA,DTC
C                IAS   : CONTROLS SYM/ASYM (INPUT)
C-----------------------------------------------------------------------
      PARAMETER (NPTS=512)
      REAL FR(*),THETA(*),DTC(*),WR(4),WS(4),XL(NPTS),DXL(NPTS)
      INTEGER NF,NP,IAS

      PI = 2.*ASIN(1.)
      XL(1) = 0.
C---------------------------------- calculate lenght(theta)
      DO J=1,NPTS-1
        TI1 = FLOAT(J-1)/FLOAT(NPTS-1) * FLOAT(IAS+1)*PI
        TI2 = FLOAT(J)  /FLOAT(NPTS-1) * FLOAT(IAS+1)*PI
C---------------------------------- 4 point Gaussian integration
        DL = 0.
        DO K=1,4
          TK = TI1 + (TI2-TI1)*(WR(K)+1.)/2.
          RR = FR(1)/2.
          DR = 0.
          DO M=2,MF/2
            RR = RR + FR(2*M-1)*COS((M-1)*TK)
     >              + FR(2*M)  *SIN((M-1)*TK)
            DR = DR - FLOAT(M-1)*FR(2*M-1)*SIN((M-1)*TK)
     >              + FLOAT(M-1)*FR(2*M)  *COS((M-1)*TK)
          ENDDO
          DL = DL + SQRT(RR*RR + DR*DR) * WS(K)
        ENDDO
        XL(J+1) = XL(J) + DL * FLOAT(IAS+1)/2 *PI/FLOAT(NPTS-1)
        RR = FR(1)/2.
        DR = 0.
        DO M=2,MF/2
           RR = RR + FR(2*M-1)*COS((M-1)*TI2)
     >             + FR(2*M)  *SIN((M-1)*TI2)
           DR = DR - FLOAT(M-1)*FR(2*M-1)*SIN((M-1)*TI2)
     >             + FLOAT(M-1)*FR(2*M)  *COS((M-1)*TI2)
        ENDDO
      DXL(J+1) = SQRT(RR*RR + DR*DR)
      ENDDO
      TI2 = 0.
      RR = FR(1)/2.
      DR = 0.

      DO M=2,MF/2
        RR = RR + FR(2*M-1)*COS((M-1)*TI2)
     >          + FR(2*M)  *SIN((M-1)*TI2)
        DR = DR - FLOAT(M-1)*FR(2*M-1)*SIN((M-1)*TI2)
     >          + FLOAT(M-1)*FR(2*M)  *COS((M-1)*TI2)
      ENDDO
      DXL(1) = SQRT(RR*RR + DR*DR)
c      WRITE(*,2) (I,XL(I),DXL(I),I=1,NPTS)
    2 format(i4,2f10.6)

      THETA(1) = 0.
      THETA(NP) = FLOAT(IAS+1)*PI
      DTC(1)  = 1./(DXL(1)    * FLOAT(IAS+1)*PI/XL(NPTS))
      DTC(NP) = 1./(DXL(NPTS) * FLOAT(IAS+1)*PI/XL(NPTS))
      I=1
      DO J=2,NP-1
        CHIL = FLOAT(J-1)/FLOAT(NP-1) * XL(NPTS)
        DO WHILE ((CHIL.GT.XL(I)).AND.(I.LT.NPTS))
          I = I+1
         ENDDO
c--------------------------- interval located, use linear interpolation
        THETA(J) = (FLOAT(I-1) + (CHIL-XL(I))/(XL(I+1)-XL(I)) )
     >           / FLOAT(NPTS-1) * FLOAT(IAS+1)*PI
        DCT   = ( (CHIL-XL(I))/(XL(I+1)-XL(I)) * (DXL(I+1)-DXL(I))
     >           +  DXL(I) ) * FLOAT(IAS+1)*PI/XL(NPTS)
        DTC(J) =  1./ DCT
      ENDDO
      RETURN
      END
************************************************************************
*DECK INIGRID
      SUBROUTINE INIGRID(XX,YY,PSI,NR,NP,FR,MHARM,IAS)
C--------------------------------------------------------------------
C ON EXIT XX AND YY ARE FILLED WITH THE VALUES OF X,XR,XS,XRS AND
C Y,YR,YS,YRS OF EVERY NODE :
C XX(1,NODE) = X, XX(2,NODE) = XR, XX(3,NODE) = XS, XX(4,NODE) = XRS
C THE SHAPE OF THE BOUNDARY IS GIVEN BY FR(M)
C--------------------------------------------------------------------
      USE PARAM
      USE GAUSINT
      USE MESHAC
      REAL XX(4,*),YY(4,*),FR(*),PSI(*),THETA(NPMMAX),DTC(NPMMAX)

      PI = 2.*ASIN(1.)
      DT = (1.+FLOAT(IAS))*PI/REAL(NP-1)
      DR = 1./REAL(NR-1)

C--------------------------- change theta grid to constant arclenght
      IF (IARC.EQ.1) THEN
        MF = 2.*MHARM
        CALL ARCLENGTH(FR,MF,THETA,DTC,NP,IAS,XGAUSS,WGAUSS)
      ELSE
        DO J=1,NP
          THETA(J) = DT * REAL(J-1)
          DTC(J) = 1.
        ENDDO
      ENDIF

      IF ((IMESH.EQ.2).AND.(XR1.LE.1.).AND.(SIG1.LT.1.)) THEN
        CALL MESHAC2(NR,SG,DSG,DDSG,XR1,SIG1)
        NRDONE = NR
        XR1DONE = XR1
	SIG1DONE = SIG1
      ELSE
        DO I=1,NR
          SG(I)  = FLOAT(I-1)/FLOAT(NR-1)
          DSG(I)  = 1.
          DDSG(I) = 0.
        ENDDO
      ENDIF

      DO 10 I=1,NR
        DO 20 J=1,NP
          NODE = NP*(I-1) + J
          THTJ = THETA(J)
          RADIUS = SG(NR-I+1)
          XX(1,NODE) = RADIUS * FR(1) * COS(THTJ) / 2.
          XX(2,NODE) = FR(1) * COS(THTJ)          / 2.
          XX(3,NODE) = - RADIUS * FR(1) * SIN(THTJ) / 2.
          XX(4,NODE) = - FR(1) * SIN(THTJ)          / 2.
          YY(1,NODE) = RADIUS * FR(1) * SIN(THTJ)   / 2.
          YY(2,NODE) = FR(1) * SIN(THTJ)            / 2.
          YY(3,NODE) = RADIUS * FR(1) * COS(THTJ)   / 2.
          YY(4,NODE) = FR(1) * COS(THTJ)            / 2.
C---------------------------- KEEP ELLIPTICITY ON AXIS -----------
          DO 30 M = 2,MHARM
            IF (M.EQ.2) THEN
              RM = RADIUS * ( FR(2*M-1) * COS((M-1)*THTJ)
     >                      + FR(2*M)   * SIN((M-1)*THTJ) )
              DRM = ( FR(2*M-1) * COS((M-1)*THTJ)
     >              + FR(2*M)   * SIN((M-1)*THTJ))
              DRMT = RADIUS*(- FR(2*M-1) *REAL(M-1)*SIN((M-1)*THTJ)
     >                       + FR(2*M)   *REAL(M-1)*COS((M-1)*THTJ))
              DRMTR= (-FR(2*M-1)*REAL(M-1)*SIN((M-1)*THTJ)
     >                +FR(2*M) *REAL(M-1)*COS((M-1)*THTJ))
            ELSE
              RM = RADIUS**(M-1) * ( FR(2*M-1) * COS((M-1)*THTJ)
     >                           + FR(2*M)   * SIN((M-1)*THTJ) )
              DRM =(M-1)*RADIUS**(M-2) * ( FR(2*M-1) * COS((M-1)*THTJ)
     >                                  + FR(2*M)   * SIN((M-1)*THTJ))
              DRMT = RADIUS**(M-1)*(- FR(2*M-1) *(M-1)*SIN((M-1)*THTJ)
     >                              + FR(2*M)   *(M-1)*COS((M-1)*THTJ))
            DRMTR=(M-1)*RADIUS**(M-2)*(-FR(2*M-1)*(M-1)*SIN((M-1)*THTJ)
     >                                 +FR(2*M) *(M-1)*COS((M-1)*THTJ))
            ENDIF
            XX(1,NODE) = XX(1,NODE) + RM * COS(THTJ)
            YY(1,NODE) = YY(1,NODE) + RM * SIN(THTJ)
            XX(2,NODE) = XX(2,NODE) + DRM * COS(THTJ)
            YY(2,NODE) = YY(2,NODE) + DRM * SIN(THTJ)
            XX(3,NODE) = XX(3,NODE) - RM * SIN(THTJ) + DRMT * COS(THTJ)
            YY(3,NODE) = YY(3,NODE) + RM * COS(THTJ) + DRMT * SIN(THTJ)
            XX(4,NODE) = XX(4,NODE) - DRM * SIN(THTJ)+DRMTR * COS(THTJ)
            YY(4,NODE) = YY(4,NODE) + DRM * COS(THTJ)+DRMTR * SIN(THTJ)
   30     CONTINUE
          XX(2,NODE) = - XX(2,NODE) * DR/2. * DSG(NR-I+1)
          XX(3,NODE) =   XX(3,NODE) * DT/2. * DTC(J)
          XX(4,NODE) = - XX(4,NODE) * DR/2. * DT/2. * DTC(J)*DSG(NR-I+1)
          YY(2,NODE) = - YY(2,NODE) * DR/2. * DSG(NR-I+1)
          YY(3,NODE) =   YY(3,NODE) * DT/2. * DTC(J)
          YY(4,NODE) = - YY(4,NODE) * DR/2. * DT/2. * DTC(J)*DSG(NR-I+1)
          PSI(4*(NODE-1)+1) = RADIUS **2
          PSI(4*(NODE-1)+2) = - 2.* RADIUS * DR / 2. * DSG(NR-I+1)
          PSI(4*(NODE-1)+3) = 0.
          PSI(4*(NODE-1)+4) = 0.
   20   CONTINUE
   10 CONTINUE
      RETURN
      END


************************************************************************
*DECK PLOTGR
      SUBROUTINE PLOTGR(MESHNO,XX,YY,NR,NP,IAS)
C-------------------------------------------------------------------
C THE X,Y GRID IS PLOTTED USING THE ISOPARAMETRIC REPRESENTATION
C-------------------------------------------------------------------
      REAL         XX(4,*), YY(4,*)
      CHARACTER*3  STRING
      CHARACTER*13 LABEL
      INTEGER      IPLT, MESHNO, IAS, NR, NP

      WRITE(STRING,'(I3)') MESHNO
      LABEL= 'GRID  NO: '//STRING

      IF (IAS.EQ.1) THEN
        CALL NFRAME(22,11,1,-1.2,1.2,-2.5,2.5,LABEL,13,'X',1,'Y',1)
      ELSE
        CALL NFRAME(22,11,1,-1.2,1.2,0.,4.0,LABEL,13,'X',1,'Y',1)
      ENDIF
      DO 10 IR=1,NR-1
        DO 20 IP=1,NP-1
          NBASE = IP + (IR-1)*NP
          CALL PLOTCU(XX(1,NBASE),  XX(3,NBASE),
     >                YY(1,NBASE),  YY(3,NBASE),
     >                XX(1,NBASE+1),XX(3,NBASE+1),
     >                YY(1,NBASE+1),YY(3,NBASE+1))
          CALL PLOTCU(XX(1,NBASE),   XX(2,NBASE),
     >                YY(1,NBASE),   YY(2,NBASE),
     >                XX(1,NBASE+NP),XX(2,NBASE+NP),
     >                YY(1,NBASE+NP),YY(2,NBASE+NP))
   20   CONTINUE
   10 CONTINUE
      RETURN
      END

************************************************************************
*DECK PLOTCU
      SUBROUTINE PLOTCU(X1,X1S,Y1,Y1S,X2,X2S,Y2,Y2S)
C-----------------------------------------------------------------------
C PLOTS A CUBIC LINE FROM X1,Y1 TO X2,Y2 GIVEN BY THE ARRAYS XI(1..4)
C AND YI(1..4)
C-----------------------------------------------------------------------
      REAL XP(5),YP(5)
      DO 10 I=1,5
        S = -1.+2*(I-1)/4.
        CALL CUB1D(X1,X1S,X2,X2S,S,XP(I),DUMMY)
        CALL CUB1D(Y1,Y1S,Y2,Y2S,S,YP(I),DUMMY)
   10 CONTINUE
      CALL LPLOT6(2,1,XP,YP,-5,' ')
      RETURN

      END

************************************************************************
*DECK CUB1D
      SUBROUTINE CUB1D(X1,X1S,X2,X2S,S,X,XS)
C-----------------------------------------------------------------------
C CUBIC HERMITE INTERPOLATION IN ONE DIMENSION
C-----------------------------------------------------------------------
      REAL H0M,H0P,H1M,H1P,H0MS,H0PS,H1MS,H1PS

      H0M  =  (S-1.)**2 *(S+2.) * 0.25
      H0MS =  (S-1.)*(S+2.)/2. + (S-1.)**2 * 0.25
      H0P  = -(S+1.)**2 *(S-2.) * 0.25
      H0PS = -(S+1.)*(S-2.)/2. - (S+1.)**2 * 0.25
      H1M  =  (S-1.)**2 *(S+1.) * 0.25
      H1MS =  (S-1.)*(S+1.)/2. + (S-1.)**2 * 0.25
      H1P  =  (S+1.)**2 *(S-1.) * 0.25
      H1PS =  (S+1.)*(S-1.)/2. + (S+1.)**2 * 0.25

      X  = X1*H0M  + X1S*H1M +  X2*H0P  + X2S*H1P
      XS = X1*H0MS + X1S*H1MS + X2*H0PS + X2S*H1PS
      RETURN
      END


************************************************************************
*DECK FORMKQ
      SUBROUTINE FORMKQ(XX,YY,PSI,NR,NP,QQ,A,B,C,EPS,IGAM,ISOL,ITER,IAS)
C----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE MATRIX KK AND THE RIGHT HAND SIDE ARRAY QQ
C NO BOUNDARY CONDITION ARE YET USED.
C NUMBER OF ROWS AND COLUMNS : 4*NR*NP
C NR : NUMBER OF RADIAL NODES
C NP : NUMBER OF POLOIDAL NODES
C A  : THE TOTAL AMPLITUDE OF THE RHS (HBT DEFINITION)
C B  : MEASURE OF THE TOTAL PRESSURE   (HBT DEFINITION)
C EPS : THE INVERS ASPECT RATIO
C IGAM=1-4 : HBT GAMMA PROFILE INPUT, 5- : FF' AS INPUT PROFILE.
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMSOLV
      REAL QQ(*)
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL GPX(4,4,MAXNODE),GPJAC(4,4,MAXNODE)
      INTEGER INDEX(MAXNODE)


      SAVE GPX,GPJAC,INDEX

      NELM = (NP-1)*(NR-1)

      IF (ITER.EQ.1) THEN

C-----------------------------------------------------------------------
C CALCULATE THE INDEX IN THE MATRIX. IT IS REORDERED FROM
C CLOCKWISE ORDERING IN THE REST OF HELENA TO SAVE REDUCE THE MATRIX SIZE
C BY A FACTOR OF 2.
C THE POLOIDAL INDEX J CHANGES TO:
C          JN = 1           (J=1, NP)
C          JN = 2*(J-1),    (2 .LE. J .GE. (NP-1)/2+1)
C          JN = 2*(NP-J)+1, ((NP-1)/2 + 2 .LE. J .GE. NP-1)
C-----------------------------------------------------------------------
      IF (IAS.EQ.1) THEN
        DO I=1,NR
          J  = 1
          JN = 1
          IJ1 = (I-1)*NP     + J
          IJN = (I-1)*(NP-1) + JN
          INDEX(IJ1) = IJN
          DO J=2,(NP-1)/2+1
            JN = 2*(J-1)
            IJ1 = (I-1)*NP     + J
	    IJN = (I-1)*(NP-1) + JN
	    INDEX(IJ1) = IJN
          ENDDO
	  DO J=(NP-1)/2+2,NP
	    JN = 2*(NP-J)+1
	    IJ1 = (I-1)*NP     + J
	    IJN = (I-1)*(NP-1) + JN
	    INDEX(IJ1) = IJN
          ENDDO
        ENDDO
      ELSE
        DO I=1,NR
	  DO J=1,NP
	    IJ1 = (I-1)*NP + J
	    INDEX(IJ1) = IJ1
	  ENDDO
	ENDDO
      ENDIF
C------------------------------------- NELM ELEMENTS ------------------
      DO 10 I=1,4*NR*NP
        DO 20 J=1,KKLDA
           KKBIG(J,I) = 0.
   20   CONTINUE
   10 CONTINUE
      DO 50 N=1,NELM
        N1 = NODENO(N,1)
        N2 = NODENO(N,2)
        N3 = NODENO(N,3)
        N4 = NODENO(N,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 60 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 70 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            WRS = WR * WS
!            CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
!     >                  R,S,X,XR,XS)
!            CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
!     >                  R,S,Y,YR,YS)
!            CALL INTERP1(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1)
!     >                  ,PSI(4*(N4-1)+1),R,S,PS)
            CALL INTERP3(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                   YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),R,S,X,XR,XS,YR,YS,PS)

            XJAC =  XR*YS - XS*YR
            GPX(NGR,NGS,N) = X
            GPJAC(NGR,NGS,N) = XJAC
C----------------------------------------- CALCULATE RIGHT HAND SIDE ---
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS = A * ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = A * (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
C------------------------------------- 4 NODES PER ELEMENT ------------
            DO 80 I=1,4
C------------------------------------- 4 FUNCTIONS V PER NODE ---------
              DO 90 J=1,4
                NROW = 4*(INDEX(NODENO(N,I)) - 1) + J
                SUMQ = - ARHS * H(J,I,NGR,NGS) * XJAC
                VX =  YS*HR(J,I,NGR,NGS)  - YR*HS(J,I,NGR,NGS)
                VY = -XS*HR(J,I,NGR,NGS)  + XR*HS(J,I,NGR,NGS)
                QQ(NROW) = QQ(NROW) - WRS * SUMQ
C------------------------------------- 4 NODES OF FUNCTION PSI --------
                DO 100 K=1,4
C------------------------------------- 4 FUNCTIONS H IN PSI -----------
                  DO 110 L=1,4
                    NCOL = 4*(INDEX(NODENO(N,K)) - 1) + L
                    NOFF = NROW - NCOL
                    IF (NOFF.GE.0) THEN
                      PSIX =  YS*HR(L,K,NGR,NGS)  - YR*HS(L,K,NGR,NGS)
                      PSIY = -XS*HR(L,K,NGR,NGS)  + XR*HS(L,K,NGR,NGS)
                      SUMK = - 1./(1.+EPS*X) * (PSIX*VX+PSIY*VY) / XJAC
                      KKBIG(NOFF+1,NCOL) = KKBIG(NOFF+1,NCOL)
     >                                   + WRS * SUMK
                    ENDIF
  110             CONTINUE
  100           CONTINUE
   90         CONTINUE
   80       CONTINUE
   70     CONTINUE
   60   CONTINUE
   50 CONTINUE
C------------------------------------------- REMOVE EMPTY COLUMNS (BND.)
      NEND = NP - 1
      IF (IAS.EQ.0) NEND=NP
      DO J=1,NEND
          KKBIG(1,4*J-3) = 1.E20
          KKBIG(1,4*J-1) = 1.E20
      ENDDO
C------------------------------------------- IF MATRIX EXISTS THEN
      ELSE
      DO 150 N=1,NELM
        N1 = NODENO(N,1)
        N2 = NODENO(N,2)
        N3 = NODENO(N,3)
        N4 = NODENO(N,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 160 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----

          DO 170 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            X = GPX(NGR,NGS,N)
            XJAC = GPJAC(NGR,NGS,N)
            CALL INTERP1(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1)
     >                  ,PSI(4*(N4-1)+1),R,S,PS)
C----------------------------------------- CALCULATE RIGHT HAND SIDE ---
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS = A * ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = A * (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
C------------------------------------- 4 NODES PER ELEMENT ------------
            DO 180 I=1,4
C------------------------------------- 4 FUNCTIONS V PER NODE ---------
              DO 190 J=1,4
                NROW = 4*(INDEX(NODENO(N,I)) - 1) + J
                SUMQ = - ARHS * H(J,I,NGR,NGS) * XJAC
                QQ(NROW) = QQ(NROW) - WR * WS * SUMQ
  190         CONTINUE
  180       CONTINUE
  170     CONTINUE
  160   CONTINUE
  150 CONTINUE
      ENDIF
      RETURN
      END

************************************************************************
*DECK INIPRES
      SUBROUTINE INIPRES
C------------------------------------------------------------------
C INITIALIZING THE ARRAY DPRES(1001) USED BY FUNCTION DPDPSI
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF

      NPT=1001
      DO I=1,NPT
        PSI = REAL(I-1)/REAL(NPT-1)
        IF (IPAI.EQ.1) THEN
          DPDPSI = 1.0 + API*PSI + BPI*PSI**2 + CPI*PSI**3
     >         + DPI*PSI**4 + EPI*PSI**5  + FPI*PSI**6
     >         + GPI*PSI**7 + HPI*PSI**8
        ELSEIF (IPAI.EQ.2) THEN
          DPDPSI=FPI+
     >   (1.-FPI)*(1.+API*PSI+BPI*PSI**2+CPI*PSI**3+DPI*PSI**4)**EPI
        ELSEIF ((IPAI.EQ.7).OR.(IPAI.EQ.8).OR.(IPAI.EQ.9)) THEN
          DPS = 1. / (NPTS-1.)
          NINT = INT((NPTS-1)*PSI) + 1
          IF (PSI.GE.1.) NINT=NPTS-1
          DPDPSI = DPR(NINT) + (PSI-DPS*(NINT-1))/DPS *
     >                         (DPR(NINT+1)-DPR(NINT))

          IF (IPAI.EQ.7) THEN
cWO
C	    FACTOR = 1.d0 + (FPI-1.d0)*(1.d0+TANH((PSI-EPI)/GPI))/2.d0
C            DPDPSI = FACTOR * DPDPSI
		WRITE(50,*) PSI, DPDPSI
C	SQRT(PSI) because IPAI=7 uses array input with points equidistant in psi
C	This way: cutoff at s=psiped
	    DPDPSI = DPDPSI*((-0.5*(TANH((SQRT(PSI)-PSIPED)/0.02)-1.0))
     >                     *(COREP-1)+1)
	  ENDIF

          IF (IPAI.EQ.8) THEN
            IF (PSI.GT.EPI) THEN
              DPDPSI = DPDPSI + FPI*(
     >          (PSI-EPI)**2 * (3.-2*PSI-EPI)/(1.-EPI)**3 )**GPI
            ENDIF
          ENDIF
          IF (IPAI.EQ.9) DPDPSI = DPDPSI**API
        ELSEIF (IPAI.EQ.4) THEN
          DP = 1. + API * PSI + BPI * PSI**2 + CPI * PSI**3
     >          + DPI * PSI**4
          IF (PSI.GT.EPI) THEN
            DP = DP + FPI*(
     >          (PSI-EPI)**2 * (3.-2*PSI-EPI)/(1.-EPI)**3 )**GPI
          ENDIF
          DPDPSI = DP
        ELSEIF (IPAI.EQ.5) THEN
          DP = 1. + API * PSI + BPI * PSI**2 + CPI * PSI**3
     >          + DPI * PSI**4
          IF (HPI.LE.EPI) HPI=1.0
	  EHPI = (EPI+HPI)/2.
          SS = SQRT(PSI)
          IF ((SS.GT.EPI).AND.(SS.LE.EHPI)) THEN
	    CALL CUBFCT(SS,EPI,EHPI,DP1,DP2)
            DP = DP + FPI*DP1**GPI
	  ELSEIF ((SS.GT.EHPI).AND.(SS.LT.HPI)) THEN
	    CALL CUBFCT(SS,EHPI,HPI,DP1,DP2)
            DP = DP + FPI*DP2**GPI
          ENDIF
          DPDPSI = DP
	ELSEIF (IPAI .EQ. 6) THEN
          Z  = (PSI-API)/ BPI
          Z0 =     -API / BPI
	  DPDPSI = ((DPI+2.d0*EPI*PSI+ 3.d0*FPI*PSI**2)* exp(-Z)
     >           - ((CPI+DPI*PSI+EPI*PSI**2+FPI*PSI**3)/BPI * EXP(-Z))
     >           - EXP(Z) / BPI ) / ( EXP(Z) + EXP(-Z))
     >           - ((CPI + DPI * PSI + EPI*PSI**2 + FPI*PSI**3)
     >		 * EXP(-Z) - EXP(Z))
     >           * (EXP(Z)/BPI - EXP(-Z)/BPI) / (EXP(Z) + EXP(-Z))**2

          DP0 = ( (DPI - CPI/BPI)*EXP(-Z0) - EXP(Z0)/BPI )
     >        / (EXP(Z0) + EXP(-Z0))
     >        - ( CPI * EXP(-Z0) - EXP(Z0))/BPI
     >        * (EXP(Z0) - EXP(-Z0)) / (EXP(Z0) + EXP(-Z0))**2

	  DPDPSI = DPDPSI / DP0

        ELSEIF (IPAI .EQ. 11) THEN
          DPS = 1. / (NPTS-1.)
          NINT = INT((NPTS-1)*SQRT(PSI)) + 1
          IF (PSI.GE.1.) NINT=NPTS-1
          DPDPSI = DPR(NINT) + (SQRT(PSI)-DPS*(NINT-1))/DPS *
     >                         (DPR(NINT+1)-DPR(NINT))
C          FACTOR = 1.d0 + (FPI-1.d0)*(1.d0+TANH((PSI-EPI)/GPI))/2.d0
C          DPDPSI = FACTOR * DPDPSI
C!          FACTOR = 1.d0 + COREP*0.5d0*(1.d0-TANH((PSI)/0.3d0))
C!          DPDPSI = DPDPSI + COREP*EXP(-(PSI/0.38d0)**2)
C          DPDPSI = DPDPSI + COREP*(1.d0-PSI*PSI)*0.5d0*
C     >     (1.d0-TANH((PSI-0.7d0)/0.2d0))

cWO
	  DPDPSI = DPDPSI*((-0.5*(TANH((PSI-PSIPED+0.03)/0.02)-1.0))
     >                   *(COREP-1)+1)
	WRITE(50,*) PSI, DPDPSI, FACTOR

        ELSEIF (IPAI.EQ.12) THEN
	  SS = SQRT(PSI)
          DPDPSI = 1.0 + API*SS + BPI*SS**2 + CPI*SS**3
     >         + DPI*SS**4 + EPI*SS**5  + FPI*SS**6
     >         + GPI*SS**7 + HPI*SS**8
        ELSEIF (IPAI.EQ.20) THEN
	  PSI1 = MAX(0.d0,PSI - 1.d-3)
	  PSI2 = MIN(PSI + 1.d-3, 1.d0)
	  PRES1 = (1.d0+API*PSI1    + BPI*PSI1**2
     >                 +CPI*PSI1**3 + DPI*PSI1**4)
     >           *(1.d0-TANH((PSI1-EPI)/GPI))/2.d0
	  PRES2 = (1.d0+API*PSI2    + BPI*PSI2**2
     >                 +CPI*PSI2**3 + DPI*PSI2**4)
     >           *(1.d0-TANH((PSI2-EPI)/GPI))/2.d0
          DPDPSI = (PRES2 - PRES1)/ (PSI2 - PSI1)
        ENDIF
        DPRES(I) = DPDPSI
      ENDDO
      DO I=2,NPT
        DPRES(I) = DPRES(I) / DPRES(1)
      ENDDO
      DPRES(1) = 1.d0

      DPS = 1./REAL(NPT-1)
      SUM = PEDGE
      PINT(NPT) = PEDGE
      DO I=1,NPT-1
        SUM = SUM + (DPRES(NPT-I+1)+DPRES(NPT-I))*DPS/2.
        PINT(NPT-I) = SUM
      ENDDO
      RETURN
      END

      SUBROUTINE CUBFCT(S,SL,SU,H1,H2)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      REAL  H(4)
        DS= SU-SL
        Q1= (S-SL)/DS
        Q2= (SU-S)/DS
        H1= 3.*Q1**2 - 2.*Q1**3
        H2= 3.*Q2**2 - 2.*Q2**3
c        H(3)= (S-SU)*Q1**2
c        H(4)= (S-SL)*Q2**2
      RETURN
      END

************************************************************************
*DECK DPDPSI
      FUNCTION DPDPSI(PSI)
C------------------------------------------------------------------
C THE NORMALIZED PROFILE OF THE PRESSURE GRADIENT VERSUS FLUX
C THIS ROUTINE MUST BE INITIALIZED BY A CALL TO INIPRES
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF
      NPT = 1001
      DPS = 1. /REAL(NPT-1)
      NINT = MAX(INT((NPT-1)*(PSI)) + 1,1)
      IF (PSI.GE.1.) NINT=NPT-1
      DPDPSI = DPRES(NINT) + ((PSI)-DPS*(NINT-1))/DPS *
     >                         (DPRES(NINT+1)-DPRES(NINT))
      RETURN
      END
************************************************************************
*DECK INIGAM
      SUBROUTINE INIGAM
C---------------------------------------------------------------------
C SUBROUTINE TO INITIALIZE THE ARRAY DGAM USED IN THE FUNCTION DGDPSI
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF

      NPT=1001
      DO 10 I=1,NPT
        PSI = REAL(I-1)/REAL(NPT-1)
        IF ((IGAM.EQ.1).OR.(IGAM.EQ.5)) THEN
          DGDPSI = 1.0 + AGA*PSI + BGA*PSI**2 + CGA*PSI**3
     >             + DGA*PSI**4
     >             + EGA*PSI**5 + FGA*PSI**6 + GGA*PSI**7
     >             + HGA*PSI**8
        ELSEIF (IGAM .EQ. -1) THEN
          DGDPSI=FGA+
     >   (1.0-FGA)*(1.+AGA*PSI+BGA*PSI**2+CGA*PSI**3+DGA*PSI**4)**EGA
        ELSEIF (IGAM.EQ.3) THEN
          SS = SQRT(PSI)
          DGDPSI = 1.0 + AGA*SS + BGA*SS**2 + CGA*SS**3
     >             + DGA*SS**4
     >             + EGA*SS**5 + FGA*SS**6 + GGA*SS**7
     >             + HGA*SS**8
        ELSEIF ((IGAM.EQ.7).OR.(IGAM.EQ.8).OR.(IGAM.EQ.4)) THEN
          DPS = 1. / (NPTS-1.)
          NINT = INT((NPTS-1)*PSI) + 1
          IF (PSI.GE.1.) NINT=NPTS-1
          DGDPSI = DF2(NINT) + (PSI-DPS*(NINT-1))/DPS *
     >                           (DF2(NINT+1)-DF2(NINT))
        ELSEIF ((IGAM.EQ.11).OR.(IGAM.EQ.2)) THEN
          DPS = 1. / (NPTS-1.)
          NINT = INT((NPTS-1)*SQRT(PSI)) + 1
          IF (PSI.GE.1.) NINT=NPTS-1
          DGDPSI = DF2(NINT) + (SQRT(PSI)-DPS*(NINT-1))/DPS *
     >                           (DF2(NINT+1)-DF2(NINT))
        ELSEIF (IGAM.EQ.12) THEN
          SS = SQRT(PSI)
          DGDPSI = 1.0 + AGA*SS + BGA*SS**2 + CGA*SS**3
     >             + DGA*SS**4
     >             + EGA*SS**5 + FGA*SS**6 + GGA*SS**7
     >             + HGA*SS**8
        ENDIF
        DGAM(I) = DGDPSI
   10 CONTINUE
      DPS = 1./REAL(NPT-1)
      SUM = 0.
      GINT(NPT) = 0.
      DO I=1,NPT-1
        SUM = SUM + (DGAM(NPT-I+1)+DGAM(NPT-I))*DPS/2.
        GINT(NPT-I) = SUM
      ENDDO

      RETURN
      END
************************************************************************
*DECK DGDPSI
      FUNCTION DGDPSI(PSI)
C------------------------------------------------------------------
C THE NORMALIZED PROFILE OF GRADIENT OF GAMMA VERSUS FLUX
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF
      NPT=1001
      DPS = 1./REAL(NPT-1)
      NINT = MAX(INT((NPT-1)*(PSI)) + 1,1)
      IF (PSI.GE.1.) NINT=NPT-1
      DGDPSI = DGAM(NINT) + ((PSI)-DPS*(NINT-1))/DPS *
     >                         (DGAM(NINT+1)-DGAM(NINT))
      RETURN
      END
************************************************************************
*DECK CURPHI
      FUNCTION CURPHI(PSI)
C------------------------------------------------------------------
C THE NORMALIZED PROFILE OF THE CURRENT DENSITY VERSUS FLUX
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF
      IF (ICUR.EQ.1) THEN
        CURPHI = (1.0 + ACUR*PSI + BCUR*PSI**2 + CCUR*PSI**3)**HCUR
     >         + DCUR*(PSI-PSI**2)**GCUR
     >         * EXP(-((PSI-ECUR)/FCUR)**2)
      ELSEIF (ICUR.EQ.3) THEN
        CURPHI = 1.0+ACUR*PSI   +BCUR*PSI**2+CCUR*PSI**3+DCUR*PSI**4
     >              +ECUR*PSI**5+FCUR*PSI**6+GCUR*PSI**7+HCUR*PSI**8
      ELSEIF (ICUR.EQ.2) THEN
        DPS = 1./REAL(NPTS-1)
        NINT = INT((NPTS-1)*PSI) + 1
        IF (PSI.GE.1.) NINT=NPTS-1
        CURPHI = ZJZ(NINT) + (PSI-DPS*(NINT-1))/DPS *
     >                         (ZJZ(NINT+1)-ZJZ(NINT))
c        IF (PSI.GT.ECUR) THEN
c            CURPHI = CURPHI + FCUR*(
c     >          (PSI-ECUR)**2 * (3.-2*PSI-ECUR)/(1.-ECUR)**3 )**GCUR
        IF (PSI.GT.ECUR) THEN
            CURPHI = CURPHI * FCUR
        ENDIF
      ELSEIF (ICUR.EQ.4) THEN
        CUR = (1.+ACUR*PSI + BCUR * PSI**2 + CCUR * PSI**3)**INT(DCUR)
        IF (PSI.GT.ECUR) THEN
cWO original: (PSI-ECUR)**2 * (3.-2*PSI-ECUR)/(1.-ECUR)**3 )**GCUR
c		gives zero derivative at psi=1
            CUR = CUR + FCUR*(
     >          (PSI-ECUR)**2 * (2.955-2*PSI-ECUR)/(1.-ECUR)**3 )**GCUR
        ENDIF
        CURPHI = CUR
      ELSEIF (ICUR .EQ.5) THEN
        PSIN = PSI / ECUR
        CUR = (1.+ACUR*PSIN + BCUR * PSIN**2 + CCUR * PSIN**3)
        IF (PSIN .GE. 1.) CUR = 0.
	CURPHI = CUR
      ELSEIF (ICUR .EQ.6) THEN
        CUR0 = (0.5 - 0.5*TANH(( - ACUR)/ BCUR))
        CUR  = (0.5 - 0.5*TANH((PSI- ACUR)/ BCUR))
        CDCUR = (CCUR+DCUR)/2.
        IF ((PSI.GT.CCUR).AND.(PSI.LE.CDCUR)) THEN
	  CALL CUBFCT(PSI,CCUR,CDCUR,DC1,DC2)
          CUR = CUR + ECUR* (DC1**FCUR)
	ELSEIF ((PSI.GT.CDCUR).AND.(PSI.LT.DCUR)) THEN
	  CALL CUBFCT(PSI,CDCUR,DCUR,DC1,DC2)
          CUR = CUR + ECUR* (DC2**FCUR)
        ENDIF
        CUR0   = CUR0 *(0.5 - 0.5*TANH((-DCUR)/GCUR))
        CURPHI = CUR * (0.5 - 0.5*TANH((PSI-DCUR)/GCUR)) / CUR0
      ELSEIF (ICUR.EQ.11) THEN
        DPS = 1./REAL(NPTS-1)
        NINT = INT((NPTS-1)*SQRT(PSI)) + 1
        IF (PSI.GE.1.) NINT=NPTS-1
        CUR = ZJZ(NINT) + (SQRT(PSI)-DPS*(NINT-1))/DPS *
     >                         (ZJZ(NINT+1)-ZJZ(NINT))
cWO  Now change the bootstrap current (relative to corep)
C    Basically same formula as icur = 4 (different ECUR, FCUR)
        IF (SQRT(PSI).GE.ECUR) THEN
CC  This looks like good shape but current doesnt go to zero			
C            CUR = CUR + FCUR*( (SQRT(PSI)-ECUR)**2 *
C     >			(2.973-2*SQRT(PSI)-ECUR)/(1.-ECUR)**3)**GCUR
C Still not to zero, but slightly higher dependent of FCUR
C 
			CURSHAPE = (ECUR-1.)/FCUR+ECUR+2.0035
			CUR = CUR * (FCUR*( (SQRT(PSI)-ECUR)**2 *
     >			(CURSHAPE-2*SQRT(PSI)-ECUR)/(1.-ECUR)**3)**GCUR + 1.)
        ENDIF
  
cWO  If FCUR=1.0, Factor doesnt do anything
c					  ! previously used for bootstrap current
        FACTOR = 1.d0 !+ (FCUR-1.d0)*(1.d0+TANH((PSI-ECUR)/GCUR))/2.d0
        CURPHI = FACTOR * CUR

		
c	 Using tanh function
C 		CURMUL = 2.5
C 		MAXMUL = 3.0
C 		IF ((PSI.GE.PSIPED).AND.(CURMUL.NE.1.)) THEN
C 			 CURPHI = CURPHI*CURMUL*FCUR*
C      >		 ((0.5*((tanh((PSI-0.992)/0.01)-1.)+2.))*(MAXMUL-1.)+1.)*
C      >				(0.5*(TANH((PSI-PSIPED)/0.001)-1.0)+1.)
C 	    ENDIF
		
      ELSEIF (ICUR.EQ.13) THEN
        DPS = 1./REAL(NPTS-1)
        NINT = INT((NPTS-1)*SQRT(PSI)) + 1

        IF (PSI.GE.1.) NINT=NPTS-1
        DC = ZJZ(NINT) + (SQRT(PSI)-DPS*(NINT-1))/DPS *
     >                         (ZJZ(NINT+1)-ZJZ(NINT))
        IF (HCUR.LE.ECUR) HCUR=1.0
	EHCUR = (ECUR+HCUR)/2.
        SS = SQRT(PSI)
        IF ((SS.GT.ECUR).AND.(SS.LE.EHCUR)) THEN
	  CALL CUBFCT(SS,ECUR,EHCUR,DC1,DC2)
          DC = DC + FCUR*DC1**GCUR
		ELSEIF ((SS.GT.EHCUR).AND.(SS.LT.HCUR)) THEN
	  CALL CUBFCT(SS,EHCUR,HCUR,DC1,DC2)
          DC = DC + FCUR*DC2**GCUR
        ENDIF
         CURPHI=DC
      ELSEIF (ICUR.EQ.12) THEN
        SS = SQRT(PSI)
        CURPHI = 1.0+ACUR*SS   +BCUR*SS**2+CCUR*SS**3+DCUR*SS**4
     >              +ECUR*SS**5+FCUR*SS**6+GCUR*SS**7+HCUR*SS**8
      ELSEIF (ICUR.EQ.14) THEN
        DC = 1.0+ACUR*PSI   +BCUR*PSI**2+CCUR*PSI**3+DCUR*PSI**4
        IF (HCUR.LE.ECUR) HCUR=1.0
	EHCUR = (ECUR+HCUR)/2.
        SS = SQRT(PSI)
        IF ((SS.GT.ECUR).AND.(SS.LE.EHCUR)) THEN
	  CALL CUBFCT(SS,ECUR,EHCUR,DC1,DC2)
          DC = DC + FCUR*DC1**GCUR
		ELSEIF ((SS.GT.EHCUR).AND.(SS.LT.HCUR)) THEN
	  CALL CUBFCT(SS,EHCUR,HCUR,DC1,DC2)
          DC = DC + FCUR*DC2**GCUR
        ENDIF
        CURPHI=DC
      ELSEIF (ICUR.EQ.15) THEN
        DC = 1.0+ACUR*PSI   +BCUR*PSI**2+CCUR*PSI**3+DCUR*PSI**4
        IF (HCUR.LE.ECUR) HCUR=1.0
	EHCUR = (0.2*ECUR+0.8*HCUR)
        IF ((PSI.GT.ECUR).AND.(PSI.LE.EHCUR)) THEN
	  CALL CUBFCT(PSI,ECUR,EHCUR,DC1,DC2)
          DC = DC + FCUR*DC1**GCUR
		ELSEIF ((PSI.GT.EHCUR).AND.(PSI.LT.HCUR)) THEN
	  CALL CUBFCT(PSI,EHCUR,HCUR,DC1,DC2)
          DC = DC + FCUR*DC2**GCUR
        ENDIF
        CURPHI=DC
      ELSE
        CURPHI = 0.
      ENDIF
      RETURN
      END

************************************************************************
*DECK QPRFL
      FUNCTION QPRFL(PSI)
C------------------------------------------------------------------
C THE NORMALIZED PROFILE OF THE CURRENT DENSITY VERSUS FLUX
C------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF

      DPS = 1./REAL(NPTS-1)
      NINT = INT((NPTS-1)*SQRT(PSI)) + 1
      IF (PSI.GE.1.) NINT=NPTS-1
      QPRFL = QIN(NINT) + (SQRT(PSI)-DPS*(NINT-1))/DPS *
     >                         (QIN(NINT+1)-QIN(NINT))
      RETURN
      END


************************************************************************
*DECK INTERP
      SUBROUTINE INTERP(XN1,XN2,XN3,XN4,R,S,X,XR,XS,XRS,XRR,XSS)
C----------------------------------------------------------------------
C SUBROUTINE CALCULATES THE INTERPOLATED VALUE OF THE FUNCTION X GIVEN
C BY XI(1..4) AT THE FOUR NODES USING BI-CUBIC HERMITE ELEMENTS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      REAL XN1(4),XN2(4),XN3(4),XN4(4)

      HI0M   = - (R-1.)**2 * (-R-2.) * 0.25
      HRI0M  = - (R-1.)*(-R-2.)*0.5 +(R-1.)**2 * 0.25
      HRRI0M = + 1.5 * R
      HI1M   = - (R-1.)**2 * (-R-1.) * 0.25
      HRI1M  = - (R-1.)*(-R-1.)*0.5 + (R-1.)**2 *0.25
      HRRI1M = + 1.5 * R - .5

      HJ0M   = - (S-1.)**2 * (-S-2.) * 0.25
      HSJ0M  = - (S-1.)*(-S-2.)*0.5 +(S-1.)**2 * 0.25
      HSSJ0M = + 1.5 * S
      HJ1M   = - (S-1.)**2 * (-S-1.) * 0.25
      HSJ1M  = - (S-1.)*(-S-1.)*0.5 + (S-1.)**2 * 0.25
      HSSJ1M = + 1.5 * S - .5

      HI0P   = - (R+1.)**2 * (R-2.) * 0.25
      HRI0P  = - (R+1.)*(R-2.)*0.5 - (R+1.)**2 * 0.25
      HRRI0P = - 1.5 * R
      HI1P   = + (R+1.)**2 * (R-1.) * 0.25
      HRI1P  = + (R+1.)*(R-1.)*0.5 + (R+1.)**2 * 0.25
      HRRI1P = + 1.5 * R + .5

      HJ0P  = - (S+1.)**2 * (S-2.) * 0.25
      HSJ0P = - (S+1.)*(S-2.)*0.5 - (S+1.)**2 * 0.25
      HSSJ0P = - 1.5 * S
      HJ1P  = + (S+1.)**2 * (S-1.) * 0.25
      HSJ1P = + (S+1.)*(S-1.)*0.5 + (S+1.)**2 * 0.25
      HSSJ1P = + 1.5 * S + .5

      X = HI0M*HJ0M * XN1(1) + HI1M*HJ0M * XN1(2)
     >  + HI0M*HJ1M * XN1(3) + HI1M*HJ1M * XN1(4)
     >  + HI0M*HJ0P * XN2(1) + HI1M*HJ0P * XN2(2)
     >  + HI0M*HJ1P * XN2(3) + HI1M*HJ1P * XN2(4)
     >  + HI0P*HJ0M * XN4(1) + HI1P*HJ0M * XN4(2)
     >  + HI0P*HJ1M * XN4(3) + HI1P*HJ1M * XN4(4)
     >  + HI0P*HJ0P * XN3(1) + HI1P*HJ0P * XN3(2)
     >  + HI0P*HJ1P * XN3(3) + HI1P*HJ1P * XN3(4)

      XR = HRI0M*HJ0M * XN1(1) + HRI1M*HJ0M * XN1(2)
     >   + HRI0M*HJ1M * XN1(3) + HRI1M*HJ1M * XN1(4)
     >   + HRI0M*HJ0P * XN2(1) + HRI1M*HJ0P * XN2(2)
     >   + HRI0M*HJ1P * XN2(3) + HRI1M*HJ1P * XN2(4)
     >   + HRI0P*HJ0M * XN4(1) + HRI1P*HJ0M * XN4(2)
     >   + HRI0P*HJ1M * XN4(3) + HRI1P*HJ1M * XN4(4)
     >   + HRI0P*HJ0P * XN3(1) + HRI1P*HJ0P * XN3(2)
     >   + HRI0P*HJ1P * XN3(3) + HRI1P*HJ1P * XN3(4)

      XS = HI0M*HSJ0M * XN1(1) + HI1M*HSJ0M * XN1(2)
     >   + HI0M*HSJ1M * XN1(3) + HI1M*HSJ1M * XN1(4)
     >   + HI0M*HSJ0P * XN2(1) + HI1M*HSJ0P * XN2(2)
     >   + HI0M*HSJ1P * XN2(3) + HI1M*HSJ1P * XN2(4)
     >   + HI0P*HSJ0M * XN4(1) + HI1P*HSJ0M * XN4(2)
     >   + HI0P*HSJ1M * XN4(3) + HI1P*HSJ1M * XN4(4)
     >   + HI0P*HSJ0P * XN3(1) + HI1P*HSJ0P * XN3(2)
     >   + HI0P*HSJ1P * XN3(3) + HI1P*HSJ1P * XN3(4)

      XRR = HRRI0M*HJ0M * XN1(1) + HRRI1M*HJ0M * XN1(2)
     >    + HRRI0M*HJ1M * XN1(3) + HRRI1M*HJ1M * XN1(4)
     >    + HRRI0M*HJ0P * XN2(1) + HRRI1M*HJ0P * XN2(2)
     >    + HRRI0M*HJ1P * XN2(3) + HRRI1M*HJ1P * XN2(4)
     >    + HRRI0P*HJ0M * XN4(1) + HRRI1P*HJ0M * XN4(2)
     >    + HRRI0P*HJ1M * XN4(3) + HRRI1P*HJ1M * XN4(4)
     >    + HRRI0P*HJ0P * XN3(1) + HRRI1P*HJ0P * XN3(2)
     >    + HRRI0P*HJ1P * XN3(3) + HRRI1P*HJ1P * XN3(4)

      XSS = HI0M*HSSJ0M * XN1(1) + HI1M*HSSJ0M * XN1(2)
     >    + HI0M*HSSJ1M * XN1(3) + HI1M*HSSJ1M * XN1(4)
     >    + HI0M*HSSJ0P * XN2(1) + HI1M*HSSJ0P * XN2(2)
     >    + HI0M*HSSJ1P * XN2(3) + HI1M*HSSJ1P * XN2(4)
     >    + HI0P*HSSJ0M * XN4(1) + HI1P*HSSJ0M * XN4(2)
     >    + HI0P*HSSJ1M * XN4(3) + HI1P*HSSJ1M * XN4(4)
     >    + HI0P*HSSJ0P * XN3(1) + HI1P*HSSJ0P * XN3(2)
     >    + HI0P*HSSJ1P * XN3(3) + HI1P*HSSJ1P * XN3(4)

      XRS = HRI0M*HSJ0M * XN1(1) + HRI1M*HSJ0M * XN1(2)
     >    + HRI0M*HSJ1M * XN1(3) + HRI1M*HSJ1M * XN1(4)
     >    + HRI0M*HSJ0P * XN2(1) + HRI1M*HSJ0P * XN2(2)
     >    + HRI0M*HSJ1P * XN2(3) + HRI1M*HSJ1P * XN2(4)
     >    + HRI0P*HSJ0M * XN4(1) + HRI1P*HSJ0M * XN4(2)
     >    + HRI0P*HSJ1M * XN4(3) + HRI1P*HSJ1M * XN4(4)
     >    + HRI0P*HSJ0P * XN3(1) + HRI1P*HSJ0P * XN3(2)
     >    + HRI0P*HSJ1P * XN3(3) + HRI1P*HSJ1P * XN3(4)

      RETURN
      END


************************************************************************
*DECK INTERP1
      SUBROUTINE INTERP1(XN1,XN2,XN3,XN4,R,S,X)
C----------------------------------------------------------------------
C SUBROUTINE CALCULATES THE INTERPOLATED VALUE OF THE FUNCTION X GIVEN
C BY XI(1..4) AT THE FOUR NODES USING BI-CUBIC HERMITE ELEMENTS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      REAL XN1(4),XN2(4),XN3(4),XN4(4)

      HI0M  = - (R-1.)**2 * (-R-2.) * 0.25
      HI1M  = - (R-1.)**2 * (-R-1.) * 0.25

      HJ0M  = - (S-1.)**2 * (-S-2.) * 0.25
      HJ1M  = - (S-1.)**2 * (-S-1.) * 0.25

      HI0P  = - (R+1.)**2 * (R-2.)  * 0.25
      HI1P  = + (R+1.)**2 * (R-1.)  * 0.25

      HJ0P  = - (S+1.)**2 * (S-2.)  * 0.25
      HJ1P  = + (S+1.)**2 * (S-1.)  * 0.25

      X = HI0M*HJ0M * XN1(1) + HI1M*HJ0M * XN1(2)
     >  + HI0M*HJ1M * XN1(3) + HI1M*HJ1M * XN1(4)
     >  + HI0M*HJ0P * XN2(1) + HI1M*HJ0P * XN2(2)
     >  + HI0M*HJ1P * XN2(3) + HI1M*HJ1P * XN2(4)
     >  + HI0P*HJ0M * XN4(1) + HI1P*HJ0M * XN4(2)
     >  + HI0P*HJ1M * XN4(3) + HI1P*HJ1M * XN4(4)
     >  + HI0P*HJ0P * XN3(1) + HI1P*HJ0P * XN3(2)
     >  + HI0P*HJ1P * XN3(3) + HI1P*HJ1P * XN3(4)

      RETURN
      END

************************************************************************
*DECK INTERP2
      SUBROUTINE INTERP2(XN1,XN2,XN3,XN4,R,S,X,XR,XS)
C----------------------------------------------------------------------
C SUBROUTINE CALCULATES THE INTERPOLATED VALUE OF THE FUNCTION X GIVEN
C BY XI(1..4) AT THE FOUR NODES USING BI-CUBIC HERMITE ELEMENTS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      REAL XN1(4),XN2(4),XN3(4),XN4(4)

      HI0M  = - (R-1.)**2 * (-R-2.) * 0.25
      HRI0M = - (R-1.)*(-R-2.)*0.5 +(R-1.)**2 * 0.25
      HI1M  = - (R-1.)**2 * (-R-1.) * 0.25
      HRI1M = - (R-1.)*(-R-1.)*0.5 + (R-1.)**2 * 0.25
      HJ0M  = - (S-1.)**2 * (-S-2.) * 0.25
      HSJ0M = - (S-1.)*(-S-2.)*0.5 +(S-1.)**2  * 0.25
      HJ1M  = - (S-1.)**2 * (-S-1.) * 0.25
      HSJ1M = - (S-1.)*(-S-1.)*0.5 + (S-1.)**2 * 0.25

      HI0P  = - (R+1.)**2 * (R-2.) * 0.25
      HRI0P = - (R+1.)*(R-2.)*0.5 - (R+1.)**2 * 0.25
      HI1P  = + (R+1.)**2 * (R-1.) * 0.25
      HRI1P = + (R+1.)*(R-1.)*0.5 + (R+1.)**2 * 0.25

      HJ0P  = - (S+1.)**2 * (S-2.) * 0.25
      HSJ0P = - (S+1.)*(S-2.)*0.5 - (S+1.)**2 * 0.25
      HJ1P  = + (S+1.)**2 * (S-1.) * 0.25
      HSJ1P = + (S+1.)*(S-1.)*0.5 + (S+1.)**2 * 0.25

      X = HI0M*HJ0M * XN1(1) + HI1M*HJ0M * XN1(2)
     >  + HI0M*HJ1M * XN1(3) + HI1M*HJ1M * XN1(4)
     >  + HI0M*HJ0P * XN2(1) + HI1M*HJ0P * XN2(2)
     >  + HI0M*HJ1P * XN2(3) + HI1M*HJ1P * XN2(4)
     >  + HI0P*HJ0M * XN4(1) + HI1P*HJ0M * XN4(2)
     >  + HI0P*HJ1M * XN4(3) + HI1P*HJ1M * XN4(4)
     >  + HI0P*HJ0P * XN3(1) + HI1P*HJ0P * XN3(2)
     >  + HI0P*HJ1P * XN3(3) + HI1P*HJ1P * XN3(4)

      XR = HRI0M*HJ0M * XN1(1) + HRI1M*HJ0M * XN1(2)
     >   + HRI0M*HJ1M * XN1(3) + HRI1M*HJ1M * XN1(4)
     >   + HRI0M*HJ0P * XN2(1) + HRI1M*HJ0P * XN2(2)
     >   + HRI0M*HJ1P * XN2(3) + HRI1M*HJ1P * XN2(4)
     >   + HRI0P*HJ0M * XN4(1) + HRI1P*HJ0M * XN4(2)
     >   + HRI0P*HJ1M * XN4(3) + HRI1P*HJ1M * XN4(4)
     >   + HRI0P*HJ0P * XN3(1) + HRI1P*HJ0P * XN3(2)
     >   + HRI0P*HJ1P * XN3(3) + HRI1P*HJ1P * XN3(4)

      XS = HI0M*HSJ0M * XN1(1) + HI1M*HSJ0M * XN1(2)
     >   + HI0M*HSJ1M * XN1(3) + HI1M*HSJ1M * XN1(4)
     >   + HI0M*HSJ0P * XN2(1) + HI1M*HSJ0P * XN2(2)
     >   + HI0M*HSJ1P * XN2(3) + HI1M*HSJ1P * XN2(4)
     >   + HI0P*HSJ0M * XN4(1) + HI1P*HSJ0M * XN4(2)
     >   + HI0P*HSJ1M * XN4(3) + HI1P*HSJ1M * XN4(4)
     >   + HI0P*HSJ0P * XN3(1) + HI1P*HSJ0P * XN3(2)
     >   + HI0P*HSJ1P * XN3(3) + HI1P*HSJ1P * XN3(4)

      RETURN
      END

      SUBROUTINE INTERP3(XN1,XN2,XN3,XN4,
     >                   YN1,YN2,YN3,YN4,
     >                   PN1,PN2,PN3,PN4,R,S,
     >                   X,XR,XS,YR,YS,PS)
C----------------------------------------------------------------------
C SUBROUTINE CALCULATES THE INTERPOLATED VALUE OF THE FUNCTION X GIVEN
C BY XI(1..4) AT THE FOUR NODES USING BI-CUBIC HERMITE ELEMENTS
C----------------------------------------------------------------------
      REAL XN1(4),XN2(4),XN3(4),XN4(4)
      REAL YN1(4),YN2(4),YN3(4),YN4(4)
      REAL PN1(4),PN2(4),PN3(4),PN4(4)

      HI0M  = - (R-1.)**2 * (-R-2.) * 0.25
      HRI0M = - (R-1.)*(-R-2.)*0.5 +(R-1.)**2 * 0.25
      HI1M  = - (R-1.)**2 * (-R-1.) * 0.25
      HRI1M = - (R-1.)*(-R-1.)*0.5 + (R-1.)**2 * 0.25
      HJ0M  = - (S-1.)**2 * (-S-2.) * 0.25
      HSJ0M = - (S-1.)*(-S-2.)*0.5 +(S-1.)**2  * 0.25
      HJ1M  = - (S-1.)**2 * (-S-1.) * 0.25
      HSJ1M = - (S-1.)*(-S-1.)*0.5 + (S-1.)**2 * 0.25

      HI0P  = - (R+1.)**2 * (R-2.) * 0.25
      HRI0P = - (R+1.)*(R-2.)*0.5 - (R+1.)**2 * 0.25
      HI1P  = + (R+1.)**2 * (R-1.) * 0.25
      HRI1P = + (R+1.)*(R-1.)*0.5 + (R+1.)**2 * 0.25

      HJ0P  = - (S+1.)**2 * (S-2.) * 0.25
      HSJ0P = - (S+1.)*(S-2.)*0.5 - (S+1.)**2 * 0.25
      HJ1P  = + (S+1.)**2 * (S-1.) * 0.25
      HSJ1P = + (S+1.)*(S-1.)*0.5 + (S+1.)**2 * 0.25

      X = HI0M*HJ0M * XN1(1) + HI1M*HJ0M * XN1(2)
     >  + HI0M*HJ1M * XN1(3) + HI1M*HJ1M * XN1(4)
     >  + HI0M*HJ0P * XN2(1) + HI1M*HJ0P * XN2(2)
     >  + HI0M*HJ1P * XN2(3) + HI1M*HJ1P * XN2(4)
     >  + HI0P*HJ0M * XN4(1) + HI1P*HJ0M * XN4(2)
     >  + HI0P*HJ1M * XN4(3) + HI1P*HJ1M * XN4(4)
     >  + HI0P*HJ0P * XN3(1) + HI1P*HJ0P * XN3(2)
     >  + HI0P*HJ1P * XN3(3) + HI1P*HJ1P * XN3(4)

      XR = HRI0M*HJ0M * XN1(1) + HRI1M*HJ0M * XN1(2)
     >   + HRI0M*HJ1M * XN1(3) + HRI1M*HJ1M * XN1(4)
     >   + HRI0M*HJ0P * XN2(1) + HRI1M*HJ0P * XN2(2)
     >   + HRI0M*HJ1P * XN2(3) + HRI1M*HJ1P * XN2(4)
     >   + HRI0P*HJ0M * XN4(1) + HRI1P*HJ0M * XN4(2)
     >   + HRI0P*HJ1M * XN4(3) + HRI1P*HJ1M * XN4(4)
     >   + HRI0P*HJ0P * XN3(1) + HRI1P*HJ0P * XN3(2)
     >   + HRI0P*HJ1P * XN3(3) + HRI1P*HJ1P * XN3(4)

      XS = HI0M*HSJ0M * XN1(1) + HI1M*HSJ0M * XN1(2)
     >   + HI0M*HSJ1M * XN1(3) + HI1M*HSJ1M * XN1(4)
     >   + HI0M*HSJ0P * XN2(1) + HI1M*HSJ0P * XN2(2)
     >   + HI0M*HSJ1P * XN2(3) + HI1M*HSJ1P * XN2(4)
     >   + HI0P*HSJ0M * XN4(1) + HI1P*HSJ0M * XN4(2)
     >   + HI0P*HSJ1M * XN4(3) + HI1P*HSJ1M * XN4(4)
     >   + HI0P*HSJ0P * XN3(1) + HI1P*HSJ0P * XN3(2)
     >   + HI0P*HSJ1P * XN3(3) + HI1P*HSJ1P * XN3(4)

      PS= HI0M*HJ0M * PN1(1) + HI1M*HJ0M * PN1(2)
     >  + HI0M*HJ1M * PN1(3) + HI1M*HJ1M * PN1(4)
     >  + HI0M*HJ0P * PN2(1) + HI1M*HJ0P * PN2(2)
     >  + HI0M*HJ1P * PN2(3) + HI1M*HJ1P * PN2(4)
     >  + HI0P*HJ0M * PN4(1) + HI1P*HJ0M * PN4(2)
     >  + HI0P*HJ1M * PN4(3) + HI1P*HJ1M * PN4(4)
     >  + HI0P*HJ0P * PN3(1) + HI1P*HJ0P * PN3(2)
     >  + HI0P*HJ1P * PN3(3) + HI1P*HJ1P * PN3(4)

      YR = HRI0M*HJ0M * YN1(1) + HRI1M*HJ0M * YN1(2)
     >   + HRI0M*HJ1M * YN1(3) + HRI1M*HJ1M * YN1(4)
     >   + HRI0M*HJ0P * YN2(1) + HRI1M*HJ0P * YN2(2)
     >   + HRI0M*HJ1P * YN2(3) + HRI1M*HJ1P * YN2(4)
     >   + HRI0P*HJ0M * YN4(1) + HRI1P*HJ0M * YN4(2)
     >   + HRI0P*HJ1M * YN4(3) + HRI1P*HJ1M * YN4(4)
     >   + HRI0P*HJ0P * YN3(1) + HRI1P*HJ0P * YN3(2)
     >   + HRI0P*HJ1P * YN3(3) + HRI1P*HJ1P * YN3(4)

      YS = HI0M*HSJ0M * YN1(1) + HI1M*HSJ0M * YN1(2)
     >   + HI0M*HSJ1M * YN1(3) + HI1M*HSJ1M * YN1(4)
     >   + HI0M*HSJ0P * YN2(1) + HI1M*HSJ0P * YN2(2)
     >   + HI0M*HSJ1P * YN2(3) + HI1M*HSJ1P * YN2(4)
     >   + HI0P*HSJ0M * YN4(1) + HI1P*HSJ0M * YN4(2)
     >   + HI0P*HSJ1M * YN4(3) + HI1P*HSJ1M * YN4(4)
     >   + HI0P*HSJ0P * YN3(1) + HI1P*HSJ0P * YN3(2)
     >   + HI0P*HSJ1P * YN3(3) + HI1P*HSJ1P * YN3(4)

      RETURN
      END


************************************************************************
*DECK SOLVE2
      SUBROUTINE SOLVE2(QQ,NR,NP,PSI,ITER,IAS)
C-----------------------------------------------------------------------
C SUBROUTINE TO SOLVE THE SYSTEM OF EQUATIONS USING GAUSSIAN ELIMINATION
C-----------------------------------------------------------------------
      USE PARAM
      USE COMSOLV
      REAL    QQ(*),PSI(*)
      INTEGER INDEX(MAXNODE)
      CHARACTER*1 UPLO
      EXTERNAL DPBTRF

      SAVE INDEX

      IF (IAS.EQ.0) THEN
         ND = 4*NR*NP
      ELSE
         ND = 4*NR*(NP-1)
      ENDIF
      NROW = 4*MAXMNODE
      NVAR = 4*NP+7

      IF (ITER.EQ.1) THEN
C------------------------------INVERSE OF INDEX IN FORMKQ TO RESTORE
      IF (IAS.EQ.1) THEN
        DO I=1,NR
          J  = 1
          JN = 1
          IJ1 = (I-1)*NP     + J
          IJN = (I-1)*(NP-1) + JN
	  INDEX(IJN) = IJ1
          DO J=2,(NP-1)/2+1
            JN = 2*(J-1)
	    IJ1 = (I-1)*NP     + J
	    IJN = (I-1)*(NP-1) + JN
	    INDEX(IJN) = IJ1
          ENDDO
	  DO J=(NP-1)/2+2,NP-1
	    JN = 2*(NP-J)+1
	    IJ1 = (I-1)*NP     + J
	    IJN = (I-1)*(NP-1) + JN
	    INDEX(IJN) = IJ1
          ENDDO
        ENDDO
      ENDIF

c---------------------------------------------------- ESSL version
c        CALL DPBF(KKBIG,KKLDA,ND,4*NP+8)
c---------------------------------------------------- lapack version
        UPLO='L'
        CALL DPBTRF(UPLO,ND,NVAR,KKBIG,KKLDA,INFO)
      ENDIF

      DO 220 I=1,4*NR*NP
        PSI(I) = QQ(I)
  220 CONTINUE
c-------------------------------------------------- ESSL version
c      CALL DPBS(KKBIG,KKLDA,ND,4*NP+8,PSI)
c-------------------------------------------------- lapack version

      CALL DPBTRS('L',ND,NVAR,1,KKBIG,KKLDA,PSI,NROW,INFO)

c-------------- restore to simple clockwise numbering
      IF (IAS.EQ.1) THEN
        DO I=1,4*NR*NP
          QQ(I) = PSI(I)
        ENDDO
        DO I=1,NR*(NP-1)
          DO K=1,4
	    IK = 4*(I-1)+K
	    IKN = 4*(INDEX(I)-1)+K
            PSI(IKN) = QQ(IK)
          ENDDO
        ENDDO
        DO I=1,NR
          DO K=1,4
	    IK  = 4*(I-1)*NP + K
	    IK2 = 4*(I-1)*NP + 4*(NP-1) + K
            PSI(IK2) = PSI(IK)
          ENDDO
        ENDDO
      ENDIF
c-------------------------------- fill in boundary conditions
      DO 225 I=1,NR*NP
        PSI(4*(I-1)+1) = PSI(4*(I-1)+1) + 1.
  225 CONTINUE
      DO 100 J=1,NP
        PSI(4*J-3) = 1.
        PSI(4*J-1) = 0.
  100 CONTINUE
      IF (IAS.EQ.1) THEN
        DO 230 I=1,NR
          DO 240 K=1,4
            NBASE = 4*(I-1)*NP + K
            NBASE2 = NBASE + 4*(NP-1)
            PSI(NBASE2) = PSI(NBASE)
  240     CONTINUE
  230   CONTINUE
      ENDIF
CCC      CALL PRARR1('PSI : ',PSI,4*NR*NP,203)

      RETURN
      END


************************************************************************
*DECK FINDAXIS
      SUBROUTINE FINDAXIS(XX,YY,NR,NP,PSAXIS,XAXIS,YAXIS,
     >                    NAX,RAX,SAX,IAS)
C-----------------------------------------------------------------------
C SUBROUTINE TO LOCALIZE THE POSITION OF THE MAGNETIC AXIS ; THE MINIMUM
C OF PSI OF ALL ELEMENTS
C-----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE FAXIS
      REAL XX(4,*),YY(4,*)
      DIMENSION X(2),FVEC(2),FJAC(2,2)

      PSMIN = 1.E20
      IF (IAS.EQ.1) THEN
c----------------------------------- asymmetric part -----------
      NELM = (NR-1)*(NP-1)
      XERR = 1E-8
      XTOL = 1E-4
      TOLL = 1. + XTOL
      NEQ2 = 2
      LDFJAC = 2
      LWA = 15
      NTRIAL=50
      TOLX = 1.E-8
      TOLF = 1.E-8
      DO I=NELM,NELM/2,-1
        NAXIS = I
        IFAIL = 1
        X(1) = 0.
        X(2) = 0.
        call mnewtax(ntrial,x,neq2,tolx,tolf,xerr,ferr)
        if ((xerr.le.tolx).or.(ferr.le.tolf)) then
          ifail = 0.
        else
c         WRITE(20,*) ' accuracy not reached : ',xerr,ferr
        endif
        R = X(1)
        S = X(2)
        IF ((IFAIL.EQ.0).AND.
     >      ((ABS(R).LE.TOLL).AND.(ABS(S).LE.TOLL))) THEN
          N1 = NODENO(NAXIS,1)
          N2 = NODENO(NAXIS,2)
          N3 = NODENO(NAXIS,3)
          N4 = NODENO(NAXIS,4)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                R,S,ZPSI,ZPSIR,ZPSIS,ZPSIRS,ZPSIRR,ZPSISS)
          IF (ZPSI.LT.PSMIN) THEN
            PSMIN = ZPSI
            NAX = NAXIS
            NNAX = NAXIS
            RAX = R
            SAX = S
          ENDIF
        ENDIF
      ENDDO
      N1 = NODENO(NNAX,1)
      N2 = NODENO(NNAX,2)
      N3 = NODENO(NNAX,3)
      N4 = NODENO(NNAX,4)
      CALL INTERP1(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),RAX,SAX,XAXIS)
      CALL INTERP1(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),RAX,SAX,YAXIS)
      PSAXIS = PSMIN

      ELSE
c**************************************** SYMMETRIC PART ***************
      DO N=1, (NR-1)*(NP-1), NP-1
        N1 = NODENO(N,1)
        N2 = NODENO(N,4)
        IF  (PSI(4*(N1-1)+2)*PSI(4*(N2-1)+2).LE.0.) THEN
C------------------------------------- QUAD. EQ FOR R VALUE AT MINIMUM -
          PSIM  = PSI(4*(N1-1)+1)
          PSIMR = PSI(4*(N1-1)+2)
          PSIP  = PSI(4*(N2-1)+1)
          PSIPR = PSI(4*(N2-1)+2)
          AA =  3. * (PSIM + PSIMR - PSIP + PSIPR ) / 4.
          BB =  ( - PSIMR + PSIPR ) / 2.
          CC =  ( - 3*PSIM - PSIMR + 3*PSIP - PSIPR) / 4.
          DET = BB*BB - 4.*AA*CC
          R  = ROOT(AA,BB,CC,DET,1.)
          IF (ABS(R).GT.1.+1.E-8) THEN
            R = ROOT(AA,BB,CC,DET,-1.)
          ENDIF
C-------- THE SIGN OF R CHANGES FOR ELEMENTS ON THE LEFT  (SEE REMESH) -
          CALL CUB1D(XX(1,N1),XX(2,N1),XX(1,N2),XX(2,N2),R,XAXIS,DUMMY)
          CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSAXIS,DUMMY)
          IF (PSAXIS.LT.PSMIN) THEN
            PSMIN = PSAXIS
            NAX  = N
            NAX1 = N1
            NAX2 = N2
            RAX = R
            SAX = -1.
c            WRITE(20,*) 'LOCAL MINIMUM AT XAXIS = ',XAXIS
          ENDIF
        ENDIF
      ENDDO
      DO  N=NP-1, (NR-1)*(NP-1), NP-1
        N1 = NODENO(N,2)
        N2 = NODENO(N,3)
        IF  (PSI(4*(N1-1)+2)*PSI(4*(N2-1)+2).LT.0.) THEN
C------------------------------------- QUAD. EQ FOR R VALUE AT MINIMUM -
          PSIM  = PSI(4*(N1-1)+1)
          PSIMR = PSI(4*(N1-1)+2)
          PSIP  = PSI(4*(N2-1)+1)
          PSIPR = PSI(4*(N2-1)+2)
          AA =  3. * (PSIM + PSIMR - PSIP + PSIPR ) / 4.
          BB =  ( - PSIMR + PSIPR ) / 2.
          CC =  ( - 3*PSIM - PSIMR + 3*PSIP - PSIPR) / 4.
          DET = BB*BB - 4.*AA*CC
          R  = ROOT(AA,BB,CC,DET,1.)
          IF (ABS(R).GT.1.+1.E-8) THEN
            R = ROOT(AA,BB,CC,DET,-1.)
          ENDIF
C-------- THE SIGN OF R CHANGES FOR ELEMENTS ON THE LEFT  (SEE REMESH) -
          CALL CUB1D(XX(1,N1),XX(2,N1),XX(1,N2),XX(2,N2),R,XAXIS,DUMMY)
          CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSAXIS,DUMMY)
          IF (PSAXIS.LT.PSMIN) THEN
            PSMIN = PSAXIS
            NAX  = N
            NAX1 = N1
            NAX2 = N2
            RAX = R
	    SAX = 1.
c            WRITE(20,*) 'LOCAL MINIMUM AT XAXIS = ',XAXIS
          ENDIF
        ENDIF
      ENDDO
      YAXIS = 0.
      ENDIF
      RETURN
      END

***********************************************************************
*DECK FZERO2
      SUBROUTINE FZERO2(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C----------------------------------------------------------------------
C SOLUTION DETERMINES THE MINIMUM OF THE FLUX IN ONE ELEMENT
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE FAXIS
      REAL X(N),FVEC(N),FJAC(LDFJAC,N)

      R = X(1)
      S = X(2)
      N1 = NODENO(NAXIS,1)
      N2 = NODENO(NAXIS,2)
      N3 = NODENO(NAXIS,3)
      N4 = NODENO(NAXIS,4)
      CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >            PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >            R,S,ZPSI,ZPSIR,ZPSIS,ZPSIRS,ZPSIRR,ZPSISS)
      IF (IFLAG.EQ.1) THEN
        FVEC(1) = ZPSIR
        FVEC(2) = ZPSIS
      ENDIF
      IF (IFLAG.EQ.2) THEN
        FJAC(1,1) = ZPSIRR
        FJAC(1,2) = ZPSIRS
        FJAC(2,1) = ZPSIRS
        FJAC(2,2) = ZPSISS
      ENDIF
      RETURN
      END

************************************************************************
*DECK ROOT
      FUNCTION ROOT(A,B,C,D,SGN)
C---------------------------------------------------------------------
C THIS FUNCTION GIVES BETTER ROOTS OF QUADRATICS BY AVOIDING
C CANCELLATION OF SMALLER ROOT
C---------------------------------------------------------------------
      IF (B*SGN .GE. 0.0) THEN
        ROOT = -2.0*C/(B+SGN*SQRT(D))
      ELSE
        ROOT = (-B + SGN*SQRT(D)) / (2.0 * A)
      ENDIF
      RETURN
      END
************************************************************************
*DECK NORMAL
      SUBROUTINE NORMAL(PSI,NR,NP,PSAXIS)
C-----------------------------------------------------------------------
C SUBROUTINE TO NORMALIZE PSI TO ONE ON THE BOUNDARY AND ZERO ON AXIS
C-----------------------------------------------------------------------
      REAL PSI(*)

      DO 60 I=1,NR*NP
        PSI(4*(I-1)+1) = 1. - (1.- PSI(4*(I-1)+1)) / (1.- PSAXIS)
        DO 70 L=2,4
          PSI(4*(I-1)+L) = PSI(4*(I-1)+L) / (1. - PSAXIS)
   70   CONTINUE
   60 CONTINUE
      RETURN
      END

************************************************************************
*DECK REMESH
       SUBROUTINE REMESH(XX,YY,PSI,A,B,C,EPS,NR,NP,NRNEW,NPNEW,MESHNO,
     >                  CX,CY,XAXIS,YAXIS,NAX,RAX,SAX,IGAM,IAS)
C---------------------------------------------------------------------
C FORM THE SYSTEM OF NEW FINITE ELEMENTS AS FLUX COORDINATES USING
C THE EXACT INTERPOLATION
C   XX,YY,PSI : ON EXIT CONTAIN THE VALUES ON THE NEW GRID
C   NR,NP     : THE NUMBER RADIAL AND POLOIDAL POINTS IN THE OLD GRID
C   NRNEW,NPNEW :         ,,                ,,         IN THE NEW GRID
C   XAXIS : POSITION OF MAGNETIC AXIS
C   NAX1,NAX2 : NODENUMBERS OF ELEMENT WITH MAGNETIC AXIS
C   RAX : R VALUE OF MAGNETIC AXIS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE MESH
      USE NODES
      USE TOLERA
      USE MESHAC
      REAL XX(4,*), YY(4,*), PSI(*)
      INTEGER NSYM(2*NR-2),SSYM(2*NR-2)
      LOGICAL FOUND,NOBRACK,CHAGR(MAXMNODE)

      PI = 2. * ASIN(1.)
      TOLTHT = 1.e-10
      TOLPSI = 1.e-10

      FACTAS = 1.
      IF (IAS.EQ.0) FACTAS = 2.

      DO I=1,NR-1
        NSYM(I) = (I-1)*(NP-1) + 1
        SSYM(I) = -1.
        NSYM(NR-1+I) = (NR-1)*(NP-1) - (I-1)*(NP-1)
        SSYM(NR-1+I) = +1.
      ENDDO

      IF ((IMESH.EQ.2).AND.(XR1.LE.1.).AND.(SIG1.LT.1.))THEN
          IF ((NRDONE.NE.NRNEW).OR.(XR1DONE.NE.XR1)
     >                        .OR.(SIG1DONE.NE.SIG1)) THEN
            CALL MESHAC2(NRNEW,SG,DSG,DDSG,XR1,SIG1)
            NRDONE = NRNEW
	    XR1DONE = XR1
	    SIG1DONE = SIG1
	  ENDIF
          DO I=1,NRNEW
            PSIKN(NRNEW-I+1)   = SG(I)**2
            DPSIKN(NRNEW-I+1)  = 2.*SG(I) * DSG(I)
            DDPSIKN(NRNEW-I+1) = 2.*SG(I) * DDSG(I) + 2.*DSG(I)**2
            RADPSI(NRNEW-I+1)  =  REAL(I-1)/REAL(NRNEW-1)
	  ENDDO
      ELSE
        DO  I=1,NRNEW
          RPSI =  REAL(I-1)/REAL(NRNEW-1)
          CALL RADMESH(RPSI,PSID,DPSID,DDPSID)
          PSIKN(NRNEW-I+1)   = PSID
	  DPSIKN(NRNEW-I+1)  = DPSID
	  DDPSIKN(NRNEW-I+1) = DDPSID
          RADPSI(NRNEW-I+1)  = RPSI
	ENDDO
      ENDIF
      RADPSI(NRNEW) = 0.
      PSIKN(NRNEW) = 0.
C      CALL PRARR1('PSI VALUES : ',PSIKN,NRNEW,203)
      DO 6 J=1,NPNEW
        THTKN(J) = (1.+FLOAT(IAS)) * PI * REAL(J-1)/REAL(NPNEW-1)
    6 CONTINUE
      MESHNO = MESHNO + 1
C------------------------- UPDATE OLD MESH FOR THE FIRST ITERATION ----
      DO 8 I= 1,NRNEW*NPNEW
        CHAGR(I) = .FALSE.
    8 CONTINUE
      DO 10 I = 1,NR*NP
        DO 20 K=1,4
          XXOLD(K,I) = XX(K,I)
          YYOLD(K,I) = YY(K,I)
          PSIOLD(4*(I-1)+K) = PSI(4*(I-1)+K)
   20   CONTINUE
   10 CONTINUE

C--------------------------- LOOP OVER ALL FLUXSURFACES -----------
      DO 30 I=1,NRNEW-1
        PSIVAL=PSIKN(I)
C        WRITE(*,*) ' FINDING SURFACE : ',PSIVAL,NAX
C--------------------------- FIND STARTING POINT OF FLUXCONTOUR
        DO 40 N=NAX, 1, (1-NP)
          N1 = NODENO(N,1)
          N2 = NODENO(N,2)
          N3 = NODENO(N,3)
          N4 = NODENO(N,4)
C------------------------------------- QUAD. EQ FOR R VALUE AT MINIMUM -
          RR = -1.
          CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                 PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                 RR,SAX,PSIM,PSIMR,DPSIS)
          RR = +1.
          CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                 PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                 RR,SAX,PSIP,PSIPR,DPSIS)
          AA =  3. * (PSIM + PSIMR - PSIP + PSIPR ) / 4.
          BB =  ( - PSIMR + PSIPR ) / 2.
          CC =  ( - 3*PSIM - PSIMR + 3*PSIP - PSIPR) / 4.
          DET = BB*BB - 4.*AA*CC
          R = 999.
          IF (DET .GT. 0) THEN
            R = ROOT(AA,BB,CC,DET,1.)
            IF (ABS(R).GT.1.+1.E-5) THEN
              R = ROOT(AA,BB,CC,DET,-1.)
            ENDIF
          ENDIF
          IF (ABS(R).GT.1.) THEN
            PSIMIN=MIN(PSIM,PSIP)
            PSIMAX=MAX(PSIM,PSIP)
          ELSE
            CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                 PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                 R,SAX,PS,DPSIR,DPSIS)
            PSIMIN=MIN(MIN(PSIM,PS),PSIP)
            PSIMAX=MAX(MAX(PSIM,PS),PSIP)
          ENDIF
          IF ((PSIVAL.GE.PSIMIN-PSITOL).AND.
     >        (PSIVAL.LE.PSIMAX+PSITOL))THEN
            A3 = (PSIM+PSIMR-PSIP+PSIPR)/4.
            A2 = (- PSIMR + PSIPR)/4.
            A1=(-3*PSIM-PSIMR+3*PSIP-PSIPR)/4.
            A0=( 2*PSIM+PSIMR+2*PSIP-PSIPR)/4.-PSIVAL
            CALL SOLVP3(A0,A1,A2,A3,RR,R2,R3,IFAIL)
            CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                   PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                   RR,SAX,PS,DPSIR,DPSIS)
            CALL INTERP2(XXOLD(1,N1),XXOLD(1,N2),
     >                   XXOLD(1,N3),XXOLD(1,N4),
     >                   RR,SAX,ZX,DXR,DXS)
            IF (ZX.LT.XAXIS) THEN
C              WRITE(*,*) ' NODE ON WRONG SIDE : ',ZX,XAXIS
              RR = R2
              CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                     PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                     RR,SAX,PS,DPSIR,DPSIS)
              CALL INTERP2(XXOLD(1,N1),XXOLD(1,N2),
     >                     XXOLD(1,N3),XXOLD(1,N4),
     >                     RR,SAX,ZX,DXR,DXS)
            ENDIF
            IF ((ABS(RR).LE.1.+1.E-5).AND.(ZX.GE.XAXIS)) GOTO 45
c------------------ special case psi=1.
            IF ((I.EQ.1).AND.(ABS(PS-PSIVAL).LT.1.E-5)) THEN
              RR = -1.
              GOTO 45
            ENDIF
          ENDIF
   40   CONTINUE
        WRITE(*,*) ' STARTING VALUE NOT FOUND : ',PSIVAL,RR,PS
   45   CONTINUE
c--------------------------- starting position found ----------------
        CALL INTERP2(YYOLD(1,N1),YYOLD(1,N2),
     >               YYOLD(1,N3),YYOLD(1,N4),RR,SAX,ZY,DYR,DYS)
        THT0 = ATAN2(ZY-YAXIS,ZX-XAXIS)
        IF (THT0.LT.0.) THT0 = THT0 + 2.*PI
C----------------------------- trace fluxsurface to find theta values
        DD = 0.25*(PSIVAL)**0.33 / FACTAS
        IPREV = 0
        JPREV = 0
        SS = SAX
        NN = N
        THT1 = THT0
        THT2 = THT0
        ITMAX = 5000
        ITTEST = 0
	NOBRACK = .TRUE.
        DO 50 J=1,NPNEW
	  If (IAS.EQ.1) THEN
            JINDEX = MOD(INT(THT0/(2.*PI)*NPNEW) + J,NPNEW)+1
	  ELSE
	    JINDEX = J
	  ENDIF
          THTVAL = THTKN(JINDEX)
          FOUND = .FALSE.

c---------------------------------- treat theta=pi as special point
          IF ((J.EQ.NPNEW).AND.(IAS.EQ.0)) THEN

            DO 140 NS=1, 2*NR-2
              N  = NSYM(NS)
              SS = SSYM(NS)
              N1 = NODENO(N,1)
              N2 = NODENO(N,2)
              N3 = NODENO(N,3)
              N4 = NODENO(N,4)
              RR = -1.
              CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                 PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                 RR,SS,PSIM,PSIMR,DPSIS)
              RR = +1.
              CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                 PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                 RR,SS,PSIP,PSIPR,DPSIS)
              A3 = (PSIM+PSIMR-PSIP+PSIPR)/4.
              A2 = (- PSIMR + PSIPR)/4.
              A1=(-3*PSIM-PSIMR+3*PSIP-PSIPR)/4.
              A0=( 2*PSIM+PSIMR+2*PSIP-PSIPR)/4.-PSIVAL
              CALL SOLVP3(A0,A1,A2,A3,RR,R2,R3,IFAIL)
              CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                     PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                     RR,SS,PS,DPSIR,DPSIS)
              CALL INTERP2(XXOLD(1,N1),XXOLD(1,N2),
     >                     XXOLD(1,N3),XXOLD(1,N4),
     >                     RR,SS,ZX,DXR,DXS)
              IF (ZX.GT.XAXIS) THEN
C                WRITE(*,*) ' NODE ON WRONG SIDE : ',ZX,XAXIS
                 RR = R2
                 CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                        PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                        RR,SS,PS,DPSIR,DPSIS)
                 CALL INTERP2(XXOLD(1,N1),XXOLD(1,N2),
     >                        XXOLD(1,N3),XXOLD(1,N4),
     >                        RR,SS,ZX,DXR,DXS)
              ENDIF
              IF ((ABS(RR).LE.1.+1.E-5).AND.(ZX.LT.XAXIS)
     >                                 .AND.(ZX.GE.-1.-1.E-8)) THEN
                FOUND = .TRUE.
                NN = N
                NOBRACK= .FALSE.
                GOTO 145
              ENDIF
  140       CONTINUE
          ENDIF
  145     CONTINUE


      	  DO ITN=1,ITMAX
            N1 = NODENO(NN,1)
            N2 = NODENO(NN,2)
            N3 = NODENO(NN,3)
            N4 = NODENO(NN,4)
            CALL INTERP2(XXOLD(1,N1),XXOLD(1,N2),
     >                   XXOLD(1,N3),XXOLD(1,N4),RR,SS,ZX,ZXR,ZXS)
            CALL INTERP2(YYOLD(1,N1),YYOLD(1,N2),
     >                   YYOLD(1,N3),YYOLD(1,N4),RR,SS,ZY,ZYR,ZYS)
            CALL INTERP2(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                   PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >                   RR,SS,PS,DPSIR,DPSIS)
            THT = ATAN2(ZY-YAXIS,ZX-XAXIS)

            IF (THT.LT.0.) THT = THT + 2*PI
            IF (NOBRACK) THEN
              THT1 = THT2
              THT2 = THT
	    ELSEIF (FOUND) THEN
	      THT1 = THT
            ENDIF
            IF (((ABS(THTVAL-THT).LT.TOLTHT) .OR.
     >           (ABS(THTVAL-THT+2*PI) .LT. TOLTHT) .OR.
     >           (ABS(THTVAL-THT-2*PI) .LT. TOLTHT)) .AND.
     >           (ABS(PSIVAL-PS).LT.TOLPSI))  THEN
c
c----------------------------- node located ---------------------
c
c              IF ((j.eq.np)) THEN
c              WRITE(*,*) ' NODE LOCATED : '
c              WRITE(*,41) J,ZX,ZY,PS,THT,PS-PSIVAL,THT-THTVAL
c              ENDIF
              FOUND = .TRUE.
              IF ((ABS(RR).GT.1.+TOLPSI).OR.
     >            (ABS(SS).GT.1.+TOLTHT)) THEN
                WRITE(*,*) ' WARNING : ',RR,SS
              ENDIF
              NODE = (I-1)*NPNEW + JINDEX
              N1 = NODENO(NN,1)
              N2 = NODENO(NN,2)
              N3 = NODENO(NN,3)
              N4 = NODENO(NN,4)
              CALL INTERP(XXOLD(1,N1),XXOLD(1,N2),XXOLD(1,N3),
     >                    XXOLD(1,N4),RR,SS,X,XR,XS,XRS,XRR,XSS)
              CALL INTERP(YYOLD(1,N1),YYOLD(1,N2),YYOLD(1,N3),
     >                    YYOLD(1,N4),RR,SS,Y,YR,YS,YRS,YRR,YSS)
              CALL INTERP(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >                    PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >               RR,SS,ZPSI,ZPSIR,ZPSIS,ZPSIRS,ZPSIRR,ZPSISS)
c-------------------------------------------------------------------
              RAD = (X-XAXIS)**2 + (Y-YAXIS)**2
              THY = (X-XAXIS) / RAD
              THX = -(Y-YAXIS) / RAD
              THXX = 2*(Y-YAXIS)*(X-XAXIS) / RAD**2
              THYY = - THXX
              THXY = ( (Y-YAXIS)**2 - (X-XAXIS)**2 ) / RAD**2
              THS = THX * XS + THY * YS
              THR = THX * XR + THY * YR
              THRR = THXX*XR*XR + 2*THXY*XR*YR + THX*XRR
     >             + THYY*YR*YR + THY*YRR
              THRS = THXX*XR*XS + THXY*XR*YS + THX*XRS
     >             + THXY*YR*XS + THYY*YR*YS + THY*YRS
              THSS = THXX*XS*XS + 2*THXY*XS*YS + THX*XSS
     >             + THYY*YS*YS + THY*YSS
              PTJAC = ZPSIR*THS - ZPSIS*THR
              RT = - ZPSIS / PTJAC
              ST =   ZPSIR / PTJAC
              PTJR = ZPSIRR*THS+ZPSIR*THRS-ZPSIRS*THR-ZPSIS*THRR
              PTJS = ZPSIRS*THS+ZPSIR*THSS-ZPSISS*THR-ZPSIS*THRS
              RPT = (-PTJR*THS/PTJAC**2 + THRS/PTJAC) * RT
     >            + (-PTJS*THS/PTJAC**2 + THSS/PTJAC) * ST
              SPT = ( PTJR*THR/PTJAC**2 - THRR/PTJAC) * RT
     >            + ( PTJS*THR/PTJAC**2 - THRS/PTJAC) * ST
              XPT = - XRR * THS*ZPSIS/PTJAC**2
     >              + XRS * (THS*ZPSIR + THR*ZPSIS)/PTJAC**2
     >              + XR  * RPT        + XS * SPT
     >              - XSS * THR*ZPSIR/PTJAC**2
              YPT = - YRR * THS*ZPSIS/PTJAC**2
     >              + YRS * (THS*ZPSIR + THR*ZPSIS)/PTJAC**2
     >              + YR  * RPT        + YS * SPT
     >              - YSS * THR*ZPSIR/PTJAC**2
              XYJAC = XR*YS - XS*YR
              PSIX = ( YS*ZPSIR - YR*ZPSIS) / XYJAC
              PSIY = (-XS*ZPSIR + XR*ZPSIS) / XYJAC
c              CALL RADMESH(RADPSI(I),DUM,DDUM,DDDUM)
              ETAP = 1./ DPSIKN(I)
              EJAC = ETAP * (PSIX*THY - PSIY*THX)
              XX(1,NODE) = X
              YY(1,NODE) = Y
C-------------------------------- WATCH MINUS SIGN FROM R ORIENTATION --
              XX(2,NODE) = - ( THY  / EJAC)   /  (2.*(NRNEW-1.))
              YY(2,NODE) = - (-THX  / EJAC)   /  (2.*(NRNEW-1.))
              XX(3,NODE) =+(-ETAP*PSIY / EJAC) / (FACTAS*(NPNEW-1)/PI)
              YY(3,NODE) =+( ETAP*PSIX / EJAC) / (FACTAS*(NPNEW-1)/PI)
              XX(4,NODE) =-(XPT/ETAP)/(2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
              YY(4,NODE) =-(YPT/ETAP)/(2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
              PSI(4*(NODE-1)+1) = ZPSI
              PSI(4*(NODE-1)+2) = - DPSIKN(I) / (2.*(NRNEW-1.))
              PSI(4*(NODE-1)+3) = 0.
              PSI(4*(NODE-1)+4) = 0.
              CHAGR(NODE) = .TRUE.
              IF (ITN.GT.ITTEST) ITTEST = ITN
              GOTO 50
            ELSEIF ( ((THT1.LE.THTVAL+1.e-5).AND.
     >                (THT2.GE.THTVAL-1.e-5))
     >        .OR.
     >         ((IAS.EQ.1).AND.(THT1.GT.THT2+1.57)
     >                    .AND.(THT1.LE.THTVAL+1.e-5)
     >                    .AND.(THT2+2*PI.GE.THTVAL-1.e-5))
     >        .OR.
     >         ((IAS.EQ.1).AND.(THT1.GT.THT2+1.57)
     >                    .AND.(THT1-2*PI.LE.THTVAL+1.e-5)
     >                    .AND.(THT2.GE.THTVAL-1.e-5)) )  THEN


   47         FORMAT(' THETA BRACK : ',3E14.6)
c
c----------------------------- theta value bracketed ------------
c
              FVEC1 = -(ZY - YAXIS) + (ZX - XAXIS) * TAN(THTVAL)
              FVEC2 = -PS + PSIVAL
              FJAC11 = ZYR - ZXR*TAN(THTVAL)
              FJAC12 = ZYS - ZXS*TAN(THTVAL)
              FJAC21 = DPSIR
              FJAC22 = DPSIS
              DIS = FJAC22*FJAC11-FJAC12*FJAC21
              DR = (FJAC22*FVEC1-FJAC12*FVEC2)/DIS
              DS = (FJAC11*FVEC2-FJAC21*FVEC1)/DIS
              DR = SIGN(MIN(ABS(DR),0.05),DR)
              DS = SIGN(MIN(ABS(DS),0.05),DS)
              RR = RR + DR
              SS = SS + DS
   46         FORMAT(' 2D NEWTON : ',1P5E12.4)
   56         FORMAT(' No brack : ',1P5E12.4)
              NOBRACK = .FALSE.
            ELSE
c
c------------------------ theta not found yet, track flux surface
c
              NOBRACK = .TRUE.
              DS = - DD * DPSIR/ABS(DPSIR)
              IF (ABS((PSIVAL-PS)/PSIVAL).GT.0.01) THEN
                DS=0.
                DR = (PSIVAL-PS) / DPSIR
                RR = RR + DR
              ELSE
                DR = (PSIVAL-PS - DPSIS * DS ) / DPSIR
                DTOT = SQRT(DR*DR+DS*DS)
                DR = DR * DD/DTOT
                DS = DS * DD/DTOT
                RR = RR + DR
                SS = SS + DS
              ENDIF
            ENDIF
            IN = (NN-1)/(NP-1)+1
            JN = MOD(NN-1,NP-1)+1
            IF (SS.LT.-1.) THEN
              IF (JPREV.EQ.-1) THEN
                NN = NN - 1
                IF (JN.EQ.1) NN = NN + (NP-1)
                SS = SS + 2.
                JPREV = 0
              ELSE
                SS = -1.
                JPREV = -1
              ENDIF
            ELSEIF (SS.GT.1) THEN
              IF (JPREV.EQ.1) THEN
                NN = NN + 1
                IF (JN.EQ.NP-1) NN = NN - (NP-1)
                SS = SS - 2.
                JPREV = 0
              ELSE
                SS = 1.
                JPREV = 1
              ENDIF
            ENDIF
            IF (RR.LT.-1) THEN
              IF ((IN.GT.1).AND.(IPREV.EQ.-1))THEN
                NN = NN - NP + 1
                RR = RR + 2.
                IPREV = 0
              ELSE
                RR = -1.
                IPREV = -1
              ENDIF
            ELSEIF (RR.GT.1.) THEN
              IF ((IN.LT.NR-1).AND.(IPREV.EQ.1)) THEN
                NN = NN + NP - 1
                RR = RR - 2.
                IPREV = 0
              ELSE
                RR = 1.
                IPREV = 1
              ENDIF
            ENDIF
          ENDDO
          WRITE(20,*) ' FATAL NODE NOT FOUND : ',I,J,PSIVAL,THTVAL
   50   CONTINUE
C        WRITE(*,*) ' MAX NUMBER ITERATIONS : ',ITTEST,DD
   30 CONTINUE
   41 FORMAT(i3,2e12.4,4e14.6)
   42 FORMAT(i3,6e12.4)

      IF (IAS.EQ.1) THEN
C------------------------------------ copy tht=0 to tht=2PI
      DO I=1,NRNEW-1
        J = NPNEW
        NODE = (I-1)*NPNEW + J
        J0 = 1
        N0 = (I-1)*NPNEW + J0
        XX(1,NODE) = XX(1,N0)
        XX(2,NODE) = XX(2,N0)
        XX(3,NODE) = XX(3,N0)
        XX(4,NODE) = XX(4,N0)
        YY(1,NODE) = YY(1,N0)
        YY(2,NODE) = YY(2,N0)
        YY(3,NODE) = YY(3,N0)
        YY(4,NODE) = YY(4,N0)
        PSI(4*(NODE-1)+1) = PSI(4*(N0-1)+1)
        PSI(4*(NODE-1)+2) = PSI(4*(N0-1)+2)
        PSI(4*(NODE-1)+3) = PSI(4*(N0-1)+3)
        PSI(4*(NODE-1)+4) = PSI(4*(N0-1)+4)
        CHAGR(NODE)=.TRUE.
      ENDDO
      ENDIF

      N1 = NODENO(NAX,1)
      N2 = NODENO(NAX,2)
      N3 = NODENO(NAX,3)
      N4 = NODENO(NAX,4)
      S = SAX
      R = RAX
      CALL INTERP(XXOLD(1,N1),XXOLD(1,N2),XXOLD(1,N3),XXOLD(1,N4),
     >            R,S,X,XR,XS,XRS,XRR,XSS)
      CALL INTERP(YYOLD(1,N1),YYOLD(1,N2),YYOLD(1,N3),YYOLD(1,N4),
     >            R,S,Y,YR,YS,YRS,YRR,YSS)
      CALL INTERP(PSIOLD(4*(N1-1)+1),PSIOLD(4*(N2-1)+1),
     >            PSIOLD(4*(N3-1)+1),PSIOLD(4*(N4-1)+1),
     >            R,S,EQPSI,PSIR,PSIS,PSRS,PSRR,PSSS)
      EJAC = XR*YS - XS*YR
      RY = - XS / EJAC
      RX =   YS / EJAC
      SY =   XR / EJAC
      SX = - YR / EJAC
      PSIXX = PSRR*RX*RX + 2.*PSRS*RX*SX + PSSS*SX*SX
      PSIYY = PSRR*RY*RY + 2.*PSRS*RY*SY + PSSS*SY*SY
      IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
        PSIYY0 = A * ( 1. + B*X*(1+EPS*X/2.) ) - PSIXX
      ELSE
        PSIYY0 = A * ( C + B*(1.+EPS*X)**2   ) - PSIXX
      ENDIF
      CX = PSIXX/2.
      CY = PSIYY0/2.
ccc      WRITE(20,3) CX,CY
      DO 60 J = 1, NPNEW
        NODE = (NRNEW-1)*NPNEW + J
        XX(1,NODE) = XAXIS
        YY(1,NODE) = YAXIS
        TN = TAN(THTKN(J))
        TN2 = TN**2
        CN = COS(THTKN(J))
        IF (THTKN(J).EQ.(PI/2.)) THEN
          XX(2,NODE) = 0.
          YY(2,NODE) = -1./(SQRT(CY)*2.*(NRNEW-1))
          XX(4,NODE) = +1./(SQRT(CY)*2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
          YY(4,NODE) = 0.
        ELSEIF (THTKN(J).EQ.(3*PI/2)) THEN
          XX(2,NODE) = 0.
          YY(2,NODE) = -1./(SQRT(CY)*2.*(NRNEW-1))
          XX(4,NODE) = +1./(SQRT(CY)*2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
          YY(4,NODE) = 0.
        ELSE
         XX(2,NODE) = - SIGN(1.,CN)/(SQRT(CX + CY*TN2) * 2.*(NRNEW-1))
         YY(2,NODE) = - ABS(TN) /(SQRT(CX + CY*TN2) * 2.*(NRNEW-1))
         XX(4,NODE) = + (CX+CY*TN2)**(-1.5) * CY * ABS(TN)
     >              / (CN**2 * 2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
         YY(4,NODE) = - CX * (CX + CY*TN2)**(-1.5) / ( CN*ABS(CN)
     >              * 2.*FACTAS*(NRNEW-1)*(NPNEW-1)/PI)
        ENDIF
        IF (THTKN(J).GT.PI) THEN
          YY(2,NODE) = - YY(2,NODE)
          XX(4,NODE) = - XX(4,NODE)
        ENDIF
        XX(3,NODE) = 0.
        YY(3,NODE) = 0.
        PSI(4*(NODE-1)+1) = 0.
        PSI(4*(NODE-1)+2) = 0.
        PSI(4*(NODE-1)+3) = 0.
        PSI(4*(NODE-1)+4) = 0.
        CHAGR(NODE) = .TRUE.
   60 CONTINUE
      DO 80 I=1,NRNEW*NPNEW
        IF (.NOT. CHAGR(I)) WRITE(20,*) 'NODE MISSED AT I = ',I
   80 CONTINUE
      NR = NRNEW
      NP = NPNEW
      CALL ELMNO(NR,NP,NODENO)
    3 FORMAT('ELONGATION ON AXIS : CX,CY = ',2E12.4)

c      CALL PRARR1('X : ',XX,4*NRNEW*NPNEW,203)
c      CALL PRARR1('Y : ',YY,4*NRNEW*NPNEW,203)
c      CALL PRARR1('PSI : ',PSI,4*NRNEW*NPNEW,203)

      RETURN
      END


************************************************************************
*DECK RADMESH
      SUBROUTINE RADMESH(RPSI,ZPSI,DZPSI,DDZPSI)
C-------------------------------------------------------------------
C FUNCTION TO DETERMINE A NON-EQUIDISTANT GRID IN SQRT(PSI). GRID
C WILL BE EQUIDISTANT IN RPSI.
C RPSI MUST BE A MONOTONIC FUNCTION OF PSI BETWEEN 0. AND 1.
C-------------------------------------------------------------------
      USE MESHAC

      ZPSI = RPSI**2
      DZPSI= 2*RPSI
      DDZPSI = 2.
      IF (IMESH.EQ.1) THEN
        A1 = AMESH
        A2 = BMESH

        ZPSI  = RPSI**2  + A1*RPSI**3 + A2*RPSI**4 -(A1+A2)*RPSI**5
        DZPSI = 2*RPSI   + 3*A1*RPSI**2 + 4*A2*RPSI**3
     >                 - 5*(A1+A2)*RPSI**4
        DDZPSI = 2.      + 6*A1*RPSI + 12*A2*RPSI**2
     >                 - 20*(A1+A2)*RPSI**3

      ENDIF
      RETURN
      END
************************************************************************
*DECK RPACK
      FUNCTION RPACK(X,X0,B,T,AMP)
C-------------------------------------------------------------------
C POLYNOMIAL TO DESCRIBE LOCAL CHANGES OF THE RADIAL COORDINATE
C       X0 : CENTRAL X VALUE
C       B  : THE POSITION OF THE MAXIMUM
C       T  : THE TOTAL HALFWIDTH
C------------------------------------------------------------------
      DX = ABS(X-X0)
      IF (DX.LT.B) THEN
        P4 = 2*DX/B - (DX/B)**2
      ELSE
        DX = DX - B
        P4 = 1. - 3*(DX/(T-B))**2 + 2*(DX/(T-B))**3
      ENDIF
      RPACK = AMP * P4
      IF (X.LT.X0) RPACK = - RPACK
      IF (ABS(X-X0).GT.T) RPACK = 0.
      RETURN
      END

************************************************************************
*DECK DRPACK
      FUNCTION DRPACK(X,X0,B,T,AMP)
C-------------------------------------------------------------------
C DERIVATIVE OF FUNCTION RPACK

C------------------------------------------------------------------
      DX = ABS(X-X0)
      IF (DX.LT.B) THEN
        P4 = 2./ B - 2*DX / B**2
      ELSE
        DX = DX - B
        P4 = - 6*DX/((T-B)**2) + 6*DX**2/(T-B)**3
      ENDIF
      DRPACK = AMP * P4
      IF (X.LT.X0) DRPACK = - DRPACK
      IF (ABS(X-X0).GT.T) DRPACK = 0.
      RETURN
      END

************************************************************************
*DECK DDRPACK
      FUNCTION DDRPACK(X,X0,B,T,AMP)
C-------------------------------------------------------------------
C SECOND DERIVATIVE OF FUNCTION RPACK
C------------------------------------------------------------------
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DDRPACK NOT CONTINUOUS !!!!!!
      DX = ABS(X-X0)
      IF (DX.LT.B) THEN
        P4 =  - 2. / B**2
      ELSE
        DX = DX - B
        P4 = - 6./((T-B)**2) + 12.*DX/(T-B)**3
      ENDIF
      DDRPACK = AMP * P4
      IF (X.LT.X0) DDRPACK = - DDRPACK
      IF (ABS(X-X0).GT.T) DDRPACK = 0.
      RETURN
      END



************************************************************************
*DECK PSIMIMA
      SUBROUTINE PSIMIMA(N,PSI,PSIMIN,PSIMAX)
C-----------------------------------------------------------------------
C SUBROUTINE TO DETERMINE THE MINIMA AND MAXIMA OF PSI AT THE ELEMENT
C BOUNDARIES FOR USE IN THE REMESH SUBROUTINE
C   N : NUMBER OF THE ELEMENT
C   PSI : VECTOR WITH PSI VALUES
C   PSIMIN,PSIMAX : THE RESULTING MINIMUM AND MAXIMUM VALUES OF PSI
C-----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      REAL PSI(*)

      PSIMIN = 1E10
      PSIMAX =-1E10
      DO 10 I=1,4
        IM = MOD(I,4) + 1
        N1 = NODENO(N,I)
        N2 = NODENO(N,IM)
        IF (ABS(N1-N2).EQ.1) THEN
           PSIM  = PSI(4*(N1-1)+1)
           PSIMR = PSI(4*(N1-1)+3)
           PSIP  = PSI(4*(N2-1)+1)
           PSIPR = PSI(4*(N2-1)+3)
        ELSE
           PSIM  = PSI(4*(N1-1)+1)
           PSIMR = PSI(4*(N1-1)+2)
           PSIP  = PSI(4*(N2-1)+1)
           PSIPR = PSI(4*(N2-1)+2)
        ENDIF
        PSMA = MAX(PSIM,PSIP)
        PSMI = MIN(PSIM,PSIP)
        AA =  3. * (PSIM + PSIMR - PSIP + PSIPR ) / 4.
        BB =  ( - PSIMR + PSIPR ) / 2.
        CC =  ( - 3*PSIM - PSIMR + 3*PSIP - PSIPR) / 4.
        DET = BB**2 - 4*AA*CC
        IF (DET.GT.0.) THEN
           R = (-BB + SQRT(BB**2-4*AA*CC) ) / (2*AA)
           IF (ABS(R).GT.1.) THEN
              R = (-BB - SQRT(BB**2-4*AA*CC) ) / (2*AA)
           ENDIF
           IF (ABS(R).LE.1) THEN
              CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSMIMA,DUMMY)
              IF (PSMIMA.GT.PSMA) THEN
                 PSMA = PSMIMA
              ELSE
                 PSMI = PSMIMA
              ENDIF
           ENDIF
        ENDIF
        IF (PSMI.LT.PSIMIN) PSIMIN = PSMI
        IF (PSMA.GT.PSIMAX) PSIMAX = PSMA
   10 CONTINUE
      RETURN
      END

************************************************************************
*DECK THTMIMA
      SUBROUTINE THTMIMA(N,NP,XX,YY,XAXIS,YAXIS,THTMIN,THTMAX)
C-----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      REAL XX(4,*), YY(4,*)

      THTMIN = 1E10
      THTMAX =-1E10
      DO 10 I=1,4
        NODE = NODENO(N,I)
        THETA = ATAN2(YY(1,NODE)-YAXIS,XX(1,NODE)-XAXIS)
        IF (THETA.LT.-1.E-3) THEN
          THETA = THETA + 6.28318530717958624
        ENDIF
        IF (THETA.LT.THTMIN) THTMIN = THETA
        IF (THETA.GT.THTMAX) THTMAX = THETA
   10 CONTINUE
      IF ((THTMIN.LT.1.57).AND.(THTMAX.GT.4.71)) THEN
        THTMIN = 1E10
        THTMAX =-1E10
        DO 20 I=1,4
          NODE = NODENO(N,I)
          THETA = ATAN2(YY(1,NODE)-YAXIS,XX(1,NODE)-XAXIS)
          IF (THETA.LT.THTMIN) THTMIN = THETA
          IF (THETA.GT.THTMAX) THTMAX = THETA
   20   CONTINUE
      ENDIF
      RETURN
      END

************************************************************************
*DECK SOLVP3
      SUBROUTINE SOLVP3(C0,C1,C2,C3,X1,X2,X3,IFAIL)
C-----------------------------------------------------------------------
C SOLVES A CUBIC EQUATION WITH A SOLUTION WITH -1.< X < 1
C CN : THE COEFFICIENT OF X**N, X : THE REAL SOLUTION WITH -1.< X < 1.
C-----------------------------------------------------------------------
      REAL C0,C1,C2,C3,X1,X2,X3

      X1 = 99.
      X2 = 999.
      X3 = 9999.
      TOL = 1E-8
      IFAIL = 0
c------------------------------------- 2nd order poly for small c3
      IF (ABS(C3)/(ABS(C1)+ABS(C2)+ABS(C3)).LT.1.e-9) THEN
        AA = C2
        BB = C1
        CC = C0
        DET = BB**2 - 4*AA*CC
        IF (DET.GT.0.) THEN
          X1 = ROOT(AA,BB,CC,DET,1.)
          IF (ABS(X1).GT.1.+TOL) THEN
            X1 = ROOT(AA,BB,CC,DET,-1.)
          ENDIF
        ELSE
           IFAIL = 1
        ENDIF

      ELSE
c------------------------------------- 3rd order poly solution
      PI = 2*ASIN(1.)
      A0 = C0 / C3
      A1 = C1 / C3
      A2 = C2 / C3
      P = - (A2**2)/3. + A1
      Q = 2./27.*(A2**3) - A2 * A1/3. + A0
      DET = (P/3.)**3 + (Q/2.)**2
      IF (DET.GE.0) THEN
        U  = SIGN(1.,-Q/2.+SQRT(DET))*ABS(-Q/2. + SQRT(DET))**(1./3.)
        V  = SIGN(1.,-Q/2.-SQRT(DET))*ABS(-Q/2. - SQRT(DET))**(1./3.)
        X1 =  U + V - A2/3.
        IF (ABS(X1).GE.(1+TOL)) IFAIL = 1
      ELSE
        P = -P
        ANGLE = SIGN(1.,P)*ACOS((Q/2.)/SQRT(ABS(P)/3.)**3)
        X1 = -2.*SQRT(ABS(P)/3)*COS(ANGLE/3.) - A2/3.
        X2 = -2.*SQRT(ABS(P)/3.)*COS(2*PI/3. - ANGLE/3.) - A2/3.
        X3 = -2.*SQRT(ABS(P)/3.)*COS(2*PI/3. + ANGLE/3.) - A2/3.
      ENDIF
      IF (ABS(X1).GT.ABS(X2)) THEN
        DUM = X1
        X1 = X2
        X2 = DUM
      ENDIF
      IF (ABS(X2).GT.ABS(X3)) THEN
        DUM = X2
        X2 = X3
        X3 = DUM
      ENDIF
      IF (ABS(X1).GT.ABS(X2)) THEN
        DUM = X1
        X1 = X2
        X2 = DUM
      ENDIF
      ENDIF
      IF (ABS(X1).GT.(1.+TOL)) IFAIL=1
      RETURN
      END

************************************************************************
*DECK MAPPING
      SUBROUTINE MAPPING(XX,YY,PSI,CX,CY,XAXIS,A)
C-----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE METRIC COEFFICIENTS NEEDED FOR CASTOR
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE NODES
      USE COMPRI
      USE COMPLO
      USE GAUSINT
      USE COMMAP
      USE COMNAM
      COMMON/COMLOCAL/CCHI,SCHI,XCHI,YCHI,XOUT,YOUT,IJCHI
      REAL XX(4,*),YY(4,*),PSI(*),MAXERR
      REAL CCHI(4,2*MAXMNODE),SCHI(2*MAXMNODE)
      REAL XCHI(2*MAXMNODE),YCHI(2*MAXMNODE)
      REAL XOUT(2*MAXMNODE),YOUT(2*MAXMNODE)
      INTEGER IJCHI(2*MAXMNODE)
      LOGICAL CHIN
      REAL XPLOT(2*NRMMAX),PPLOT(2*NRMMAX),PSIPLOT(2*NRMMAX)
      REAL QPLOT(2*NRMMAX),CUPLOT(2*NRMMAX),DQPLOT(2*NRMMAX)
C--------------------------------------------- VARIABLES FOR VACUUM --
      REAL VX(2*NPMMAX-1),VY(2*NPMMAX-1)
C---------------------------------------------------------------------
      PI = 2.*ASIN(1.)
      MAXERR = -1.e20
      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.
C--------------------------------------------- NORM FACTORS ----------
      CALL PROFILES(P0,RBPHI,DP,DRBPHI,A)

      IF (NPR1.NE.0) THEN
        WRITE(20,*)
        WRITE(20,*) '*************************************'
        WRITE(20,*) '* PRESSURE PROFILE BEFORE NORM.  :  *'
        WRITE(20,*) '*************************************'
        WRITE(20,91) (P0(I)*EPS/ALFA**2,I=1,NR)
        WRITE(20,*) '*************************************'
        WRITE(20,*)
      ENDIF
   91 FORMAT(4E16.8)
      PAX = P0(1)*EPS/(ALFA*ALFA)

      RAXIS = 1. + EPS * XAXIS
      BM = RBPHI(1) / RAXIS
      B0 = 1./ BM
      R0 = 1./ RAXIS
      RADIUS = EPS * R0
      PSCALE = B0**2 * EPS / ALFA**2
      RBSCALE = B0 * R0
      DO 90 I = 1,NR
        P0(I) = P0(I) * PSCALE
        DP(I) = DP(I) * PSCALE
        RBPHI(I) = RBPHI(I) * RBSCALE
        DRBPHI(I) = DRBPHI(I) * RBSCALE
   90 CONTINUE
      CPSURF = RADIUS**2 * B0 / ALFA
      RAXIS = 1.
      PAXIS = P0(1)
      WRITE(20,2) RADIUS
      WRITE(20,3) B0
      WRITE(20,4) CPSURF
    2 FORMAT(' RADIUS : ',E12.4)
    3 FORMAT(' B0     : ',E12.4)
    4 FORMAT(' CPSURF : ',E12.4)
C------------------------------------------------- DERIVATIVE DP/DS ----
      DPE = DP(NR)
      DP0 = 0.
      DRBPHI0 = 0.


C--------------------------------------------- DATA FOR VECTOR PLOT ----
      WRITE(21,*) NR,NCHI,EPS
C------------------------------------------ - Q PROFILE ----------------
C                                           - DQ/DS PROFILE
C                                           - CHI VALUES AT THETA NODES
      DO 10 I = 1, NR-1
        SUMQ = 0.
        SUMQR = 0.
        ZPSI = PSIKN(I)
c        CALL RADMESH(RADPSI(I),ZPSI,DZPSI,DDZPSI)
        PSIR  = - DPSIKN(I)  /(2.*REAL(NR-1))
        PSIRR =   DDPSIKN(I) /(2.*REAL(NR-1))**2
        CS(NR-I+1) = SQRT(ZPSI)
        DO 20 J = 1, NP-1
          N1 = (I-1)*NP + J
          N2 = N1 + 1
          N3 = N2 + NP
          N4 = N1 + NP
          DO 30 K = 1,4
            S = XGAUSS(K)
            CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),-1.,S,
     >                  X,XR,XS,XRS,XRR,XSS)
            CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),-1.,S,
     >                  Y,YR,YS,YRS,YRR,YSS)
            EJAC = XR*YS - XS*YR
            ER = XRR*YS + XR*YRS - XRS*YR - XS*YRR
            BIGR  = (1. + EPS * X)
            SUMQ  = SUMQ - WGAUSS(K) * EJAC / ( BIGR * ABS(PSIR))
            SUMQR = SUMQR + PSIRR * EJAC / ((PSIR**2)*BIGR)* WGAUSS(K)
            SUMQR = SUMQR - ER / (BIGR*PSIR)               * WGAUSS(K)
            SUMQR = SUMQR + EJAC*EPS*XR/((BIGR**2)*PSIR)  * WGAUSS(K)
   30     CONTINUE
          CCHI(1,(I-1)*NP+J+1) = SUMQ
          CCHI(2,(I-1)*NP+J+1) = SUMQR
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),-1.,1.,
     >                X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),-1.,1.,
     >                Y,YR,YS,YRS,YRR,YSS)
          EJAC = XR*YS - XS*YR
          ER = XRR*YS + XR*YRS - XRS*YR - XS*YRR
          BIGR = (1. + EPS * X )
          ZSUMQ  = - EJAC / ( BIGR * ABS(PSIR))
          ZSUMQR = + PSIRR * EJAC / (PSIR**2 *BIGR)
          ZSUMQR =  ZSUMQR - ER / (BIGR*PSIR)
          ZSUMQR =  ZSUMQR + EJAC*EPS*XR/((BIGR**2) * PSIR)
          CCHI(3,(I-1)*NP+J+1) = ZSUMQ
          CCHI(4,(I-1)*NP+J+1) = ZSUMQR
   20   CONTINUE
        CCHI(1,(I-1)*NP+1) = 0.
        CCHI(2,(I-1)*NP+1) = 0.
        N1 = (I-1)*NP + 1
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
        CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),-1.,-1.,
     >              X,XR,XS,XRS,XRR,XSS)
        CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),-1.,-1.,
     >              Y,YR,YS,YRS,YRR,YSS)
        EJAC = XR*YS - XS*YR
        ER = XRR*YS + XR*YRS - XRS*YR - XS*YRR
        BIGR = (1. + EPS * X )
        ZSUMQ = - EJAC / (BIGR * ABS(PSIR))
        ZSUMQR =  + PSIRR * EJAC / (PSIR**2 *BIGR)
        ZSUMQR = ZSUMQR - ER / (BIGR*PSIR)
        ZSUMQR = ZSUMQR + EJAC*EPS*XR/((BIGR**2) * PSIR)
        CCHI(3,(I-1)*NP+1) =  ZSUMQ
        CCHI(4,(I-1)*NP+1) =  ZSUMQR
        QS(NR-I+1) =  0.5*FACTAS * SUMQ * RBPHI(NR-I+1) / RBSCALE
        DQS(NR-I+1)=SUMQR*RBPHI(NR-I+1)+SUMQ*DRBPHI(NR-I+1)/(2.*(NR-1))
        DQS(NR-I+1)=  0.5*FACTAS * DQS(NR-I+1) / RBSCALE
   10 CONTINUE
      IF ((EPS.LE.0.01).AND.(P0(1).LE.0.)) THEN
        S2 = SQRT(PSIKN(NR-1))
        S3 = SQRT(PSIKN(NR-2))
        QS(1) = QS(2) - (QS(3)-QS(2))/(S3-S2) * S2
      ELSE
        QS(1) = RBPHI(1)* PI/(2.*SQRT(CX*CY)*(1.+EPS*XAXIS)*RBSCALE)
      ENDIF
      DQS(1) = 0.
      DO 40 I = 1, NR-1
        DO 50 J = 1, NP
          DUM = CCHI(1,I*NP)
          NO = (I-1)*NP+J
          CCHI(1,NO) = FLOAT(1+IAS)*PI*CCHI(1,NO) / DUM
          DUM2 = CCHI(2,I*NP)
          CCHI(2,NO) = FLOAT(1+IAS)*PI*CCHI(2,NO) / DUM
          CCHI(3,NO) = FLOAT(1+IAS)*PI*CCHI(3,NO) / DUM
          CCHI(4,NO) = FLOAT(1+IAS)*PI*CCHI(4,NO) / DUM
C------------------------------ FINAL VALUES OF CHI --------------------
          QQ  = QS(NR-I+1)
          DQQ = DQS(NR-I+1)
          RB  = RBPHI(NR-I+1)
          DRB = DRBPHI(NR-I+1)  /  (2*(NR-1))
          CCHI(2,NO)= +(DQQ/QQ - DRB/RB) * CCHI(1,NO) - CCHI(2,NO)
          CCHI(4,NO)= +(DQQ/QQ - DRB/RB) * CCHI(3,NO) - CCHI(4,NO)
   50   CONTINUE
   40 CONTINUE
c      IF (NPR1.NE.0) THEN
c        CALL PRARR1(' CHI : ',CCHI,4*NR*NP,203)
c      ENDIF
      DO 45 I=1,NR
        QS(I) = QS(I) * ALFA / PI
        QPLOT(NR+1-I) = QS(I)
        DQS(I) = DQS(I) * ALFA / PI
        DQPLOT(NR+1-I) = DQS(I)
   45 CONTINUE
      QS0 = QS(1)
      QAXIS = QS(1)
      WRITE(20,*)
      WRITE(20,31) QS(1)
      WRITE(20,32) QS(NR)
   31 FORMAT('  Q ON AXIS = ',f7.4)
   32 FORMAT('  Q AT BOUNDARY = ',f7.4)

      IF (NPR1.NE.0) THEN
        CALL PRARR1(' Q PROFILE : ',QS,NR,203)
      ENDIF
      DQEC = DQS(NR)
C--------------------------- DETERMINE POSITIONS OF EQUIDISTANT CHI'S --
C                            AND CALCULATE MATRIX ELEMENTS -------------
   54 FORMAT('NUMBER OF CHI VALUES : ',I3)
      JS0 = NR - 1
      IF (IAS.EQ.0) THEN
        DO 55 J=1,NCHI
          CHIKN(J) =  PI * REAL(J-1)/REAL(NCHI-1)
          CHI(J)   =  PI * REAL(J-1)/REAL(NCHI-1)
   55   CONTINUE
      ELSE
        DO 56 J=1,NCHI
          CHIKN(J) = 2. * PI * REAL(J-1)/REAL(NCHI)
          CHI(J)   = 2. * PI * REAL(J-1)/REAL(NCHI)
   56   CONTINUE
      ENDIF
C      CALL PRARR1(' CHI : ',CHI,NCHI,201)
C      CALL PRARR2(' CCHI : ',' ',CCHI,4,NR*NP,4,203)
C      CALL PRARR2(' XX : ',' ',XX,4,NR*NP,4,203)
C      CALL PRARR2(' YY : ',' ',YY,4,NR*NP,4,203)
      DO 60 I=1,NR-1
        CALL RADMESH(RADPSI(I),ZPSI,DZPSI,DDZPSI)
        PSIR  = - DPSIKN(I)  /(2.*REAL(NR-1))
C-------------------------- FIRST POINT IS KNOWN -----------------------
        NO = (I-1)*NCHI + 1
        K = 1
        S = -1.
        SCHI((I-1)*NCHI+J) = S
        IJCHI((I-1)*NCHI+J) = K
        N1 =  (I-1)*NP + K
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
        CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >              -1.,S,XCHI(NO),XR,XS,XRS,XRR,XSS)
        CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >              -1.,S,YCHI(NO),YR,YS,YRS,YRR,YSS)
        CALL INTERP(CCHI(1,N1),CCHI(1,N2),CCHI(1,N3),CCHI(1,N4),
     >              -1.,S,DCHI,CHIR,CHIS,CHIRS,CHIRR,CHISS)
C--------------------------------------------- VACUUM DATA ---------
        IF (I.EQ.1) THEN
          VX(1) = XCHI(NO)
          VY(1) = YCHI(NO)
        ENDIF
C-----------------------------------------------------------------------
        EJAC  = (XR*YS-XS*YR)
        EJAC2 = EJAC**2
C----------------------------- DATA FOR VECTOR PLOT TO FILE 21 -----
        PSIX = PSIR * YS / EJAC
        PSIY = - PSIR * XS / EJAC
        CHIX = (CHIR * YS - CHIS * YR) / EJAC
        CHIY = (-CHIR * XS + CHIS * XR)/ EJAC
        WRITE(21,61) SQRT(ZPSI),DCHI,XCHI(NO),YCHI(NO),
     >              PSIX,PSIY,CHIX,CHIY
C----------------------------------------------------------------------
        GRPS2 = PSIR**2 * (XS**2 + YS**2) / EJAC2
        NOG = (NR-1)*NCHI - I*NCHI + 1
        GEM11(NCHI+NOG) = GRPS2 * (CPSURF/RADIUS)**2
        GRPGRC = ( CHIR * PSIR * (XS**2 + YS**2) -
     >             CHIS * PSIR * (XR*XS + YR*YS) ) / EJAC2
        GEM12(NCHI+NOG) = GRPGRC * CPSURF / (RADIUS**2)
        GEM33(NCHI+NOG) =((1.+EPS*XCHI(NO))/(1.+EPS*XAXIS))**2
        XOUT(NCHI+NOG) = XCHI(NO)
        YOUT(NCHI+NOG) = YCHI(NO)
C------------------------------------ CHECK JACOBIAN -----------------
                GRCHI2 = CHIX**2 + CHIY**2
                DUM1 = RBPHI(NR-I+1)**2 /
     >             (QS(NR-I+1)**2 * GEM33(NCHI+NOG))
                DUM2 = GEM11(NCHI+NOG) * GRCHI2 / (RADIUS**2)
                DUM3 = GEM12(NCHI+NOG)
                DUM4 = DUM2 - DUM3*DUM3
		ERRJ = ABS(DUM4-DUM1)
                IF (ERRJ.GT.MAXERR) THEN
                  MAXERR=ERRJ
                  IERR = I
                  JERR = J
                  SERR = SQRT(ZPSI)
                  CERR = DCHI
                ENDIF
CCC                WRITE(20,*) DUM1,DUM4,(dum1-dum4)/dum1
C---------------------------------------------------------------------
        JBASE = 2
        DO 70 K = 1,NP-1
          CHIN =.FALSE.
          DO 80 J = JBASE,NCHI
            ZCHI = CHIKN(J)
            NO = (I-1)*NCHI + J
            IF ((((CCHI(1,(I-1)*NP+K).LE.ZCHI).AND.(CCHI(1,(I-1)*NP+K+1)
     >           .GE.ZCHI)) ).OR.
     >           ((J.EQ.NCHI).AND.(K.EQ.NP-1).AND.(IAS.EQ.0))) THEN
              CHIN = .TRUE.
              NOM = (I-1)*NP + K
              NOP = NOM + 1
              A3 = (CCHI(1,NOM)+CCHI(3,NOM)-CCHI(1,NOP)+CCHI(3,NOP))/4.
              A2 = (- CCHI(3,NOM) + CCHI(3,NOP))/4.
           A1=(-3*CCHI(1,NOM)-CCHI(3,NOM)+3*CCHI(1,NOP)-CCHI(3,NOP))/4.
           A0=( 2*CCHI(1,NOM)+CCHI(3,NOM)+2*CCHI(1,NOP)-CCHI(3,NOP))/4.
     >       - ZCHI
              CALL SOLVP3(A0,A1,A2,A3,S,S2,S3,IFAIL)
	      IF ((IAS.EQ.0).AND.(J.EQ.NCHI)) THEN
	        S = 1.0
		IFAIL=0
	      ENDIF
              IF (IFAIL.EQ.0) THEN
                SCHI((I-1)*NCHI+J) = S
                IJCHI((I-1)*NCHI+J) = K
                N1 = (I-1)*NP + K
                N2 = N1 + 1
                N3 = N2 + NP
                N4 = N1 + NP
                CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                      -1.,S,XCHI(NO),XR,XS,XRS,XRR,XSS)
                CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                      -1.,S,YCHI(NO),YR,YS,YRS,YRR,YSS)
                CALL INTERP(CCHI(1,N1),CCHI(1,N2),CCHI(1,N3),CCHI(1,N4),
     >                      -1.,S,DCHI,CHIR,CHIS,CHIRS,CHIRR,CHISS)
C--------------------------------------------- VACUUM DATA ---------
                IF (I.EQ.1) THEN
                  VX(J) = XCHI(NO)
                  VY(J) = YCHI(NO)
                ENDIF
C-----------------------------------------------------------------------
                EJAC  = (XR*YS-XS*YR)
                EJAC2 = EJAC**2
C-------------------------------- DATA FOR VECTOR PLOT TO FILE 21 -----
                PSIX = PSIR * YS / EJAC
                PSIY = - PSIR * XS / EJAC
                CHIX = (CHIR * YS - CHIS * YR) / EJAC
                CHIY = (-CHIR * XS + CHIS * XR)/ EJAC
                WRITE(21,61) SQRT(ZPSI),DCHI,XCHI(NO),YCHI(NO),
     >                        PSIX,PSIY,CHIX,CHIY
C----------------------------------------------------------------------
                GRPS2 = PSIR**2 * (XS**2 + YS**2) / EJAC2
                NOG = (NR-1)*NCHI - I*NCHI + J
                GEM11(NCHI+NOG) = GRPS2 * (CPSURF/RADIUS)**2
                GRPGRC = ( CHIR * PSIR * (XS**2 + YS**2) -
     >                     CHIS * PSIR * (XR*XS + YR*YS) ) / EJAC2
                GEM12(NCHI+NOG) = GRPGRC * CPSURF / (RADIUS**2)
                GEM33(NCHI+NOG) =((1.+EPS*XCHI(NO))/(1.+EPS*XAXIS))**2
                XOUT(NCHI+NOG) = XCHI(NO)
                YOUT(NCHI+NOG) = YCHI(NO)
C------------------------------------ CHECK JACOBIAN -----------------
                GRCHI2 = CHIX**2 + CHIY**2
                DUM1 = RBPHI(NR-I+1)**2 /
     >             (QS(NR-I+1)**2 * GEM33(NCHI+NOG))
                DUM2 = GEM11(NCHI+NOG) * GRCHI2 / (RADIUS**2)
                DUM3 = GEM12(NCHI+NOG)
                DUM4 = DUM2 - DUM3*DUM3
		ERRJ = ABS(DUM4-DUM1)
                IF (ERRJ.GT.MAXERR) THEN
                  MAXERR=ERRJ
                  IERR = I
                  JERR = J
                  SERR = SQRT(ZPSI)
                  CERR = DCHI
                ENDIF
CCC                WRITE(20,*) DUM1,DUM4,dum1-dum4
C---------------------------------------------------------------------
              ELSE
                WRITE(20,*) 'ERROR IN SOLVP3 I,J,K : ',I,J,K,S,s2,s3
                WRITE(*,*) A0,A1,A2,A3,ZCHI
                WRITE(*,*) CCHI(1,(I-1)*NP+K),CCHI(1,(I-1)*NP+K+1)
              ENDIF
            ELSEIF (CHIN) THEN
              JBASE = J
              GOTO 70
            ENDIF
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
   61 FORMAT(8E16.8)

   62 FORMAT(' MAX. ERROR IN JACOBIAN AFTER MAPPING : ',
     >       1PE10.3,0P2F10.6,2I4)
      WRITE(20,*)
      WRITE(20,*) '***************************************************'
      WRITE(20,62) MAXERR,SERR,CERR,IERR,JERR
      WRITE(20,*) '***************************************************'


C-------------------- WRITE GEOMETRIC QUANTITIES TO TAPE12 -----------
      NMAP = 12
      OPEN(NMAP)
      WRITE(NMAP,8) JS0
      WRITE(NMAP,6) (CS(JS),JS=1,JS0+1)
      WRITE(NMAP,6) (QS(JS),JS=1,JS0+1)
      WRITE(NMAP,7) DQS(1),DQEC
      WRITE(NMAP,6) (DQS(JS),JS=2,JS0+1)
      WRITE(NMAP,6) (CURJ(JS),JS=1,JS0+1)
      WRITE(NMAP,7) DJ0,DJE
      WRITE(NMAP,8) NCHI
      WRITE(NMAP,6) (CHI(JS),JS=1,NCHI)
      WRITE(NMAP,6) (GEM11(JS),JS=NCHI+1,(JS0+1)*NCHI)
      WRITE(NMAP,6) (GEM12(JS),JS=NCHI+1,(JS0+1)*NCHI)
      WRITE(NMAP,7) CPSURF,RADIUS
      WRITE(NMAP,6) (GEM33(JS),JS=NCHI+1,(JS0+1)*NCHI)
      WRITE(NMAP,9) RAXIS
      WRITE(NMAP,6) (P0(JS),JS=1,JS0+1)
      WRITE(NMAP,7) DP0,DPE
      WRITE(NMAP,6) (RBPHI(JS),JS=1,JS0+1)
      WRITE(NMAP,7)  DRBPHI0,DRBPHIE
C----------------------------------------- ADDITIONAL DATA FOR VACUUM --
      WRITE(NMAP,6) (VX(JS),JS=1,NCHI)
      WRITE(NMAP,6) (VY(JS),JS=1,NCHI)
      WRITE(NMAP,9) EPS
      WRITE(NMAP,6) (XOUT(JS),JS=NCHI+1,(JS0+1)*NCHI)
      WRITE(NMAP,6) (YOUT(JS),JS=NCHI+1,(JS0+1)*NCHI)
C----------------------------------------- write profiles to vector file
      WRITE(21,11)  (P0(JS),RBPHI(JS),QS(JS),JS=1,JS0+1)
      WRITE(21,11)  CPSURF
      CLOSE(21)
    6 FORMAT(4E16.8)
    7 FORMAT(2E16.8)
    8 FORMAT(I5)
    9 FORMAT(E16.8)
   11 FORMAT(3E16.8)


      IF (NPL1.NE.0) THEN
      CALL LBLBOT('HELENA EQUILIBRIUM Version 12',34)
      CALL PLOTM(XCHI,YCHI,NR-1,NCHI,IAS)
C---------------------------------------- PLOT PROFILES --------------
      DO 65 I=1,NR
        IF (IAS.EQ.1) THEN
          ILEFT = (I-1)*NP + (NP+1)/2
	ELSE
	  ILEFT = (I-1)*NP + NP
	ENDIF
        IRIGHT= (I-1)*NP + 1
        XL = XX(1,ILEFT)
        XR = XX(1,IRIGHT)
        ZPSI = PSIKN(I)
        IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
          CUPLOT(I) = DGDPSI(ZPSI) + B*XL*(1.+EPS*XL/2.)*DPDPSI(ZPSI)
          CUPLOT(2*NR-I)= DGDPSI(ZPSI)+B*XR*(1.+EPS*XR/2.)*DPDPSI(ZPSI)
        ELSE
          CUPLOT(I) = C*DGDPSI(ZPSI) + B*(1.+EPS*XL)**2 * DPDPSI(ZPSI)
          CUPLOT(2*NR-I)=C*DGDPSI(ZPSI)+B*(1.+EPS*XR)**2 *DPDPSI(ZPSI)
        ENDIF
        CUPLOT(I) = A * CUPLOT(I) / (1.+EPS*XL)
        IF (I.NE.NR) THEN
          CUPLOT(2*NR-I) = A * CUPLOT(2*NR-I) / (1.+EPS*XR)
        ENDIF
   65 CONTINUE
      DO 100 I=1,NR
        IF (IAS.EQ.1) THEN
          ILEFT = (I-1)*NP + (NP+1)/2
	ELSE
	  ILEFT = (I-1)*NP + NP
	ENDIF
        IRIGHT= (I-1)*NP + 1
        XPLOT(I) = XX(1,ILEFT)
        XPLOT(2*NR-I) = XX(1,IRIGHT)
        PPLOT(I) = P0(NR-I+1)/P0(1)
        PPLOT(2*NR-I) = P0(NR-I+1)/P0(1)
        PSIPLOT(I) = PSIKN(I)
        PSIPLOT(2*NR-I) = PSIKN(I)
        QPLOT(2*NR-I) = QPLOT(I)
        DQPLOT(2*NR-I) = DQPLOT(I)
  100 CONTINUE
      IF (NPR1.NE.0) THEN
      WRITE(20,*)
      WRITE(20,*) '***************************************************'
      WRITE(20,*) '* I,     X,          PSI,          P,         Q   *'
      WRITE(20,*) '***************************************************'
      WRITE(20,101) (I,XPLOT(I),PSIPLOT(I),PAX*PPLOT(I),QPLOT(I),
     >               I=1,2*NR-1)
      WRITE(20,*) '***************************************************'
      WRITE(20,*)
      ENDIF
  101 FORMAT(I4,4E14.6)

C      CALL PRARR1('XPLOT : ',XPLOT,2*NR-1,203)
C      CALL PRARR1('PPLOT : ',PPLOT,2*NR-1,203)
C      CALL PRARR1('PSIPLOT : ',PSIPLOT,2*NR-1,203)
C      CALL PRARR1('QPLOT : ',QPLOT,2*NR-1,203)

      WRITE(25,*) '  I,           X,          PSI,        P,       Q '
      WRITE(25,106) (XPLOT(I),PSIPLOT(I),PAXIS*PPLOT(I),QPLOT(I),
     >               I=1,2*NR-1)
  106 FORMAT(4E12.4)
      IF (XPLOT(1).LT.-0.999)       XPLOT(1)     =-0.999
      IF (XPLOT(2*NR-1).GT.0.999)   XPLOT(2*NR-1)= 0.999

      OPEN(UNIT=98, FILE='profilesout.TXT')
       DO I=1,2*NR-1
       WRITE(98,106) XPLOT(I),CUPLOT(I),PAXIS*PPLOT(I),QPLOT(I)
       END DO
      close(UNIT=98)

      CALL LPLOT(3,2,1,XPLOT,PPLOT,2*NR-1,1,'PRESSURE',8,'X',1,'P/P0',4)
      CALL LPLOT(3,3,1,XPLOT,QPLOT,2*NR-1,1,'Q-PROFILE',9,'X',1,'Q',1)
      CALL LBLBOT('HELENA EQUILIBRIUM Version 12',34)
      CALL LPLOT(2,2,1,XPLOT,PSIPLOT,2*NR-1,1,'FLUX',4,'X',1,'PSI',3)
      CALL LPLOT(3,2,1,XPLOT,PPLOT,2*NR-1,1,'PRESSURE',8,'X',1,'P/P0',4)
      CALL LPLOT(2,3,1,XPLOT,CUPLOT,2*NR-1,1,'CURRENT DENSITY',15,'X',1,
     >           'J',1)
      CALL LPLOT(3,3,1,XPLOT,QPLOT,2*NR-1,1,'Q-PROFILE',9,'X',1,'Q',1)
      CALL LBLBOT('HELENA EQUILIBRIUM Version 12',34)

      ENDIF

      RETURN
      END

************************************************************************
*DECK PROFILES
      SUBROUTINE PROFILES(P0,RBPHI,DP,DRBPHI,A)
C-----------------------------------------------------------------------
C SUBROUTINE TO INTEGRATE THE EQUILIBRIUM PROFILES
C USING THE HBT DEFINITIONS
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE NODES
      REAL P0(*),RBPHI(*),DP(*),DRBPHI(*)

C---------------------------- DERIVATIVES TO R = SQRT(PSI) !!!! --------
      DO 10 J=1,NR
        FLUX = PSIKN(NR-J+1)
        IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
          P0(J) =  .5 * A * B * PRES(FLUX)
          DP(J) = - A * B * DPDPSI(FLUX) * SQRT(FLUX)
          DGAM = - 2. * A * DGDPSI(FLUX) * SQRT(FLUX)
          GAM = A * XGAMMA(FLUX)
          RBPHI(J) = P0(J) - EPS * GAM
          DRBPHI(J) = DP(J) - EPS * DGAM
        ELSE
          P0(J) = EPS * A * B * PRES(FLUX)
          DP(J) = -2.*EPS*A*B * DPDPSI(FLUX) * SQRT(FLUX)
          RBPHI(J) = - EPS * A * C * XGAMMA(FLUX)
          DRBPHI(J) = EPS * 2.*A*C * DGDPSI(FLUX) * SQRT(FLUX)
        ENDIF
        RBPHI(J) = SQRT( 1. - 2.*EPS*RBPHI(J) /ALFA**2)
        DRBPHI(J) = - 1./(2.*RBPHI(J)) * 2.*EPS*DRBPHI(J)/ALFA**2
   10 CONTINUE
      RETURN
      END

***********************************************************************
*DECK PRES
      FUNCTION PRES(FLUX)
C-----------------------------------------------------------------------
C PRESSURE PROFILE AS A FUNCTION OF PSI
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF
      NPT=1001
      NINT = INT((NPT-1)*FLUX)+1
      IF (FLUX.GE.1.) NINT=NPT-1
      DPS = 1./REAL(NPT-1)
      SUM = PINT(NINT+1)
      DPSI = NINT*DPS - FLUX
      PSII = (NINT-1)*DPS
      DPRI = DPRES(NINT) + (FLUX-PSII)/DPS*(DPRES(NINT+1)-DPRES(NINT))
      SUM  = SUM + DPSI * (DPRI + DPRES(NINT+1))/2.
      PRES = SUM
      RETURN
      END

************************************************************************
*DECK XGAMMA
      FUNCTION XGAMMA(FLUX)
C-----------------------------------------------------------------------
C SECOND PROFILE AS A FUNCTION OF PSI
C-----------------------------------------------------------------------
      USE PARAM
      USE COMDAT
      USE COMPROF
      NPT=1001
      NINT = INT((NPT-1)*FLUX)+1
      IF (FLUX.GE.1.) NINT=NPT-1
      DPS = 1./REAL(NPT-1)
      SUM = GINT(NINT+1)
      DPSI = NINT*DPS - FLUX
      PSII = (NINT-1)*DPS
      DF2I = DGAM(NINT)+(FLUX-PSII)/DPS*(DGAM(NINT+1)-DGAM(NINT))
      SUM = SUM + DPSI * (DF2I + DGAM(NINT+1))/2.
      XGAMMA = SUM
      RETURN
      END

************************************************************************
*DECK PLOTM
      SUBROUTINE PLOTM(XX,YY,NR,NP,IAS)
C-----------------------------------------------------------------------
C SIMPLE PLOTROUTINE TO PLOT A GRID KNOWN ONLY BY THE POSITIONS, NOT
C THE DERIVATIVES (SEE PLOTGR)
C-----------------------------------------------------------------------
      REAL XX(*),YY(*)
      CHARACTER LABEL*40

      LABEL='STRAIGHT FIELD LINE COORDINATES'
      IF (IAS.EQ.1) THEN
        CALL NFRAME(22,11,1,-1.2,1.2,-2.5,2.5,LABEL,31,'X',1,'Y',1)
      ELSE
        CALL NFRAME(22,11,1,-1.2,1.2,0.,4.0,LABEL,31,'X',1,'Y',1)
      ENDIF
      DO 10 I=1,NR,2
        NO = (I-1)*NP + 1
        CALL LPLOT(2,1,1,XX(NO),YY(NO),-NP,1,' ',1,' ',1,' ',1)
        IF (IAS.EQ.1) THEN
          CALL LPLOT(2,1,1,XX(NO),YY(NO),-NP,NP-1,' ',1,' ',1,' ',1)
        ENDIF
   10 CONTINUE
      DO 20 J=1,NP
        CALL LPLOT(2,1,1,XX(J),YY(J),-NR*NP,NP,' ',1,' ',1,' ',1)
   20 CONTINUE
      RETURN
      END

************************************************************************
*DECK TRIANG
      SUBROUTINE TRIANG(XX,YY,PSI,XAXIS,NR,NP,XELL,XTRIA)
C----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE TRIANGULARITY OF FLUXSURFACES AS A FUNCTION
C OF THE MINOR RADIUS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      REAL XX(4,*),YY(4,*),PSI(*),
     >     RR(8*MAXMNODE+2),RRT(4*MAXMNODE),XTRIA(*),XELL(*)

      WRITE(20,*) '**************************************************'
      WRITE(20,*) '*  ELLIPTICITY AND TRIANGULARITY   (FOURIER CO.) *'
      WRITE(20,*) '**************************************************'
      WRITE(20,*) '* INDEX,  S,    RADIUS,   SHIFT,   ELLIP,   TRIA *'
      WRITE(20,*) '**************************************************'
      DO 10 I=1,NR-1
        DO 20 J=1, NP-1
          N1 = (I-1)*NP + J
          N2 = N1 + 1
          N3 = N2 + NP
          N4 = N1 + NP
          DO 30 NGS=1,4
            R  = -1.
            S  = -1. + 0.5 * REAL(NGS-1)
            CALL INTERP1(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                  R,S,X)
            CALL INTERP1(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                  R,S,Y)
            CALL INTERP1(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                   PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),R,S,PS)
            RRT(4*(J-1)+ NGS) = SQRT((X-XAXIS)**2 + Y*Y)
   30     CONTINUE
   20   CONTINUE
        RRT(4*(NP-1)+1) = SQRT((XX(1,I*NP)-XAXIS)**2)
        DO 35 J=1,4*(NP-1)+1
          RR(J) = RRT(4*(NP-1)-J+2)
   35   CONTINUE
        DO 40 J=1,4*(NP-1)-1
          RR(4*(NP-1)+J+1) = RR(4*(NP-1)-J+1)
   40   CONTINUE
        NPT = 8*(NP-1)
        CALL RFT2(RR,NPT,1)
   41   FORMAT(I3,E16.8)
        XRAD   =    RR(1)/REAL(NPT)
        XSHIFT =   -2.*RR(3)/REAL(NPT) / SQRT(PS)
        XELL(I)  = -4.*RR(5)/REAL(NPT) / SQRT(PS)
        XTRIA(I) =  8.*RR(7)/REAL(NPT) / PS
        WRITE(20,11) I,SQRT(PS),XRAD,XSHIFT,XELL(I),XTRIA(I)
   10 CONTINUE
   11 FORMAT(I3,5E16.8)
c-----------------------------------------------------------------
c alternative to calculate the ellipticity and triangularity by
c finding the point dy/ds=0 on each fluxsurface.
c-----------------------------------------------------------------
      WRITE(20,*) '**************************************************'
      WRITE(20,*) '*  ELLIPTICITY AND TRIANGULARITY (GEOMETRIC CO.) *'
      WRITE(20,*) '**************************************************'
      DO 60 I=1,NR-1
      DO 50 J=1, NP-1
        N1 = (I-1)*NP + J
        N2 = (I-1)*NP + J + 1
        IF  (YY(3,N1)*YY(3,N2).LE.0.) THEN
C------------------------------------- QUAD. EQ FOR S VALUE AT MINIMUM -
        YYM  = YY(1,N1)
        YYMR = YY(3,N1)
        YYP  = YY(1,N2)
        YYPR = YY(3,N2)
        AA =  3. * (YYM + YYMR - YYP + YYPR ) / 4.
        BB =  ( - YYMR + YYPR ) / 2.
        CC =  ( - 3*YYM - YYMR + 3*YYP - YYPR) / 4.
        DET = BB*BB - 4.*AA*CC
        S  = ROOT(AA,BB,CC,DET,1.)
        IF (ABS(S).GT.1.+1.E-8) THEN
          S = ROOT(AA,BB,CC,DET,-1.)
        ENDIF
        CALL CUB1D(XX(1,N1),XX(3,N1),XX(1,N2),XX(3,N2),S,XTOP,DUMMY)
        CALL CUB1D(YY(1,N1),YY(3,N1),YY(1,N2),YY(3,N2),S,YTOP,DUMMY)
        XLEFT = XX(1,(I-1)*NP+1)
        XRIGHT = XX(1,(I-1)*NP+NP)
        XGEO = (XLEFT + XRIGHT)/2.
        XRAD = (XRIGHT-XLEFT)/2.
        XSHIFT = (XAXIS - XGEO)/XRAD
        ELL    = YTOP/XRAD - 1.
        TRI    = -(XTOP-XGEO)/XRAD**2
        PS = PSI(4*(N1-1)+1)
        WRITE(20,51) I,SQRT(PS),XRAD,XSHIFT,ELL,TRI
        ENDIF
   50 CONTINUE
   51 FORMAT(i3,5e16.8)
   60 CONTINUE

      RETURN
      END

************************************************************************
*DECK DIAGNO
      SUBROUTINE DIAGNO(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,IAS,
     >                  XAXIS,YAXIS,IAV,ZVOL,ZVOLP,ZJAR)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE OUTPUT QUANTITIES
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMPROF
      USE COMOUT
      USE COMPRI
      USE COMPLO
      REAL XX(4,*),YY(4,*),PSI(*),ZVOL(*),ZVOLP(*),ZJAR(*)
      REAL XL(NRMMAX),ZPS(NRMMAX),AVC(NRMMAX)
      REAL ZPAR(NRMMAX),ZAR(NRMMAX)

      PI = 2. * ASIN(1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.
cWO calculation of PAREAPED added
      AREA = 0.
      VOLUME = 0.
      BP2VOL = 0.
      PAREA = 0.
	  PAREAPED = 0.
      PAR2  = 0.
      CAREA = 0.
C------------------------------------- NELM ELEMENTS ------------------
      DO 10 I=NR-1,1,-1
      DO 15 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 20 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 30 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                   R,S,X,XR,XS)
            CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                   R,S,Y,YR,YS)
            CALL INTERP2(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                   PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                   R,S,PS,PSR,PSS)
            XJAC =  XR*YS - XS*YR
            AREA = AREA + WR * WS * XJAC
            VOLUME = VOLUME + (1.+EPS*X) * WR * WS * XJAC
            GRPS2 = PSR**2 * (XS**2 + YS**2) / XJAC**2
            BP2 = GRPS2 / (1.+EPS*X)**2
            BP2VOL = BP2VOL + (1.+EPS*X) * BP2 *WR*WS*XJAC
            PRESS = PRES(PS)
            PAREA = PAREA + WR * WS * XJAC * PRESS
			
			IF (PS .LT. 0.96) THEN
             PAREAPED = PAREAPED + WR * WS * XJAC * PRES(0.96)
			ELSE
             PAREAPED = PAREAPED + WR * WS * XJAC * PRES(PS)
			ENDIF
			
            PAR2  = PAR2  + WR * WS * XJAC * PRESS**2
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            CAREA = CAREA + WR * WS * XJAC * ARHS
   30     CONTINUE
   20   CONTINUE
   15 CONTINUE
        IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
          ZJAR(I) = -FACTAS * A * CAREA
          ZPAR(I) = -FACTAS * 0.5 * A * B * PAREA
        ELSE
          ZJAR(I) = -FACTAS * A * CAREA
          ZPAR(I) = -FACTAS * EPS * A * B * PAREA
        ENDIF
        ZAR(I)  = FACTAS * ABS(AREA)
	ZVOL(NR-I+1) = FACTAS * ABS(VOLUME)
   10 CONTINUE
      AREA   = FACTAS * ABS(AREA)
      VOLUME = FACTAS *  ABS(VOLUME)
      RAV = VOLUME / AREA
      IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
        PAREA = FACTAS * 0.5 * A * B * PAREA
		PAREAPED = FACTAS * 0.5 * A * B * PAREAPED
        PAR2  = FACTAS * (0.5*A*B)**2 * PAR2
        CAREA = FACTAS * A * CAREA
      ELSE
        PAREA =  FACTAS * EPS * A * B * PAREA
		PAREAPED =  FACTAS * EPS * A * B * PAREAPED
        PAR2  =  FACTAS * (EPS * A * B)**2 * PAR2
        CAREA =  FACTAS * A * CAREA
      ENDIF

C      WRITE(*,*) 'LUCE TEST : ',BP2VOL/VOLUME,BP2VOL,VOLUME
C      WRITE(*,*) 'LUCE EPS,ALFA : ',EPS,ALFA

      VOLUME = 2 * PI * VOLUME
      XLI = ABS(4*FACTAS*PI * BP2VOL / (RAV * CAREA**2))
      BETA   = 2* (EPS/ALFA**2) * ABS(PAREA) / AREA
      BETASTAR = 2 * (EPS/ALFA**2) * SQRT(ABS(PAR2)/AREA)
      BETAPL = (8*PI/EPS)    * ABS(PAREA) / CAREA**2
	  BETAPED1 = (8*PI/EPS)    * ABS(PAREAPED) / CAREA**2
      CURRENT  = ABS(EPS*CAREA / ALFA)

      DO I=1,NR
        ZVOL(I) = ZVOL(I) * 2. * PI
	ZJAR(I) = ZJAR(I) * EPS / ALFA
      ENDDO

      WRITE(20,*)
      WRITE(20,*) '***************************************'
      WRITE(20,11) XAXIS,YAXIS
      WRITE(20,2) BETAPL
      WRITE(20,3) BETA
      WRITE(20,8) BETASTAR
      WRITE(20,12) 1.256637*BETA/CURRENT
      WRITE(20,4) CURRENT
      WRITE(20,5) AREA
      WRITE(20,6) VOLUME
      WRITE(20,7) XLI
      WRITE(20,9) ALFA
      WRITE(20,13) A,B,C
	  WRITE(20,14) BETAPED1
      WRITE(20,*) '***************************************'
      WRITE(20,*)

      WRITE(TXTOUT(20),2) BETAPL
      WRITE(TXTOUT(21),3) BETA
      WRITE(TXTOUT(22),8) BETASTAR
      WRITE(TXTOUT(23),4) CURRENT
      WRITE(TXTOUT(24),5) AREA
      WRITE(TXTOUT(25),6) VOLUME
      WRITE(TXTOUT(26),7) XLI
      WRITE(TXTOUT(27),9) ALFA
      IF (NPL1.NE.0) THEN
        IX = 500
	IY = 722
	CALL DLCH(IX,IY,'HELENA OUTPUT PARAMETERS : ',27,-2)
	IY=IY-24
        DO JJ=20,27
          CALL DLCH(IX,IY,TXTOUT(JJ),100,-1)
	  IY = IY - 18
        ENDDO
      ENDIF
    2 FORMAT('  POLOIDAL BETA : ',E12.4)
    3 FORMAT('  TOROIDAL BETA : ',E12.4)
    4 FORMAT('  TOTAL CURRENT : ',E12.4)
    5 FORMAT('  TOTAL AREA    : ',E12.4)
    6 FORMAT('  TOTAL VOLUME  : ',E12.4)
    7 FORMAT('  INT. INDUCTANCE : ',E12.4)
    8 FORMAT('  BETA STAR     : ',E12.4)
    9 FORMAT('  POL. FLUX     : ',E12.4)
   11 FORMAT('  MAGNETIC AXIS : ',2F9.5)
   12 FORMAT('  NORM. BETA    : ',F9.5)
   13 FORMAT('  A,B,C         : ',3E12.4)
cWO
   14 FORMAT('  Beta pedestal : ',E12.4)

C------------------------------------- NELM ELEMENTS ------------------
      R = -1.
      DO 40 I=1,NR-1
      XLENGTH = 0.
      AVCUR = 0.
      SUMC3 = 0.
      SUMC5 = 0.
      DO 50 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                 R,S,X,XR,XS)
          CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                 R,S,Y,YR,YS)
          CALL INTERP2(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                 PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                 R,S,PS,PSR,PSS)

          XJAC =  XR*YS - XS*YR
          BIGR = (1. + EPS * X)
	  ZJDCHI = BIGR * XJAC / DABS(PSR)
	  DL = SQRT(XS**2 + YS**2)
	  XLENGTH = XLENGTH +  DL * WS
	  SUMC5 = SUMC5 + WS * ZJDCHI
	  IF (IAV .EQ.0) THEN
	    DLV = DL
	  ELSE
	    DLV = ZJDCHI
	  ENDIF
          SUMC3 = SUMC3 + WS * DLV
          IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
            ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
          ELSE
            ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
          ENDIF
          ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
          IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            AVCUR = AVCUR + WS * ARHS * DLV
   60     CONTINUE
   50   CONTINUE
        ZPS(I) = PS
        XL(I)  = FACTAS * XLENGTH
        AVC(I) = FACTAS * A * EPS * AVCUR/ (ALFA * SUMC3)
	ZVOLP(NR-I+1) = FACTAS * ABS(SUMC5)
   40 CONTINUE

      AVC(NR) = AVC(NR-1) - (AVC(NR-2)-AVC(NR-1))
     >        /(ZPS(NR-2)-ZPS(NR-1))*ZPS(NR-1)

c      XLI = ABS(4*FACTAS*PI * BP2VOL / (RAV * CAREA**2))

C      WRITE(*,*) ' TEST LUCE 2:',XL(1),CAREA

      XLI_LUCE = 2.*FACTAS*PI*ABS(BP2VOL/VOLUME)*XL(1)**2/(CAREA**2)
      XLI_3    = ABS(4.*FACTAS*PI * BP2VOL / (CAREA**2))

      WRITE(20,17) XLI_LUCE
      WRITE(20,18) XLI_3
      WRITE(20,19) XL(1)
   17 FORMAT('  XLI(LUCE) : ',E12.4)
   18 FORMAT('  XLI(3)    : ',E12.4)
   19 FORMAT('  LENGTH    : ',E12.4)


      IF (NPR1.NE.0) THEN
        WRITE(20,*)
        WRITE(20,*) '***********************************************'
     >             //'**********************************'
        WRITE(20,*) '* I   PSI     S      <J>     ERROR   LENGTH    '
     >             //'AREA   VOLUME   CURRENT   BUSSAC   P/P0 *'
        WRITE(20,*) '***********************************************'
     >             //'**********************************'

      DO 70 I=NR-1,1,-1
        IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
          PRS =  0.5 * A * B * PRES(ZPS(I))
        ELSE
          PRS =  EPS * A * B * PRES(ZPS(I))
        ENDIF
        BUSSAC = -(2./EPS)*(PRS - ZPAR(I) / ZAR(I))
     >          / (ZJAR(I)/XL(I))**2
        ERR = AVC(I)/AVC(NR) - CURPHI(ZPS(I))
        WRITE(20,72) I,ZPS(I),SQRT(ZPS(I)),AVC(I)/AVC(NR),
     >               ERR ,XL(I),ZAR(I), ZVOL(NR-I+1), ZJAR(I), BUSSAC,
     >               PRES(ZPS(I))/PRES(ZPS(NR-1))
   70 CONTINUE
      WRITE(20,*) '***********************************************'
     >             //'********'
      WRITE(20,*)
      ENDIF

   72 FORMAT(I5,10E14.6)
      RETURN
      END

************************************************************************
*DECK MOMENTS
      SUBROUTINE MOMENTS(XX,YY,PSI,NR,NP,XIAB,ALFA,EPS,IAS)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE SOME CURRENT MOMENTS
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPLO
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),PSI(4*MAXMNODE)

      PI = 2. * ASIN(-1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      R = -1.
      YM1 = 0.
      YZ1 = 0.
      CUR = 0.
C------------------------------------- INTEGRAL OVER PLASMA BOUNDARY
      DO 10 J=1,NP-1
        NELM = J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 20 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)
          XJAC =  XR*YS - XS*YR
          DL = SQRT(XS**2+YS**2)
          BPL = EPS * PSR / XJAC * SQRT(XS**2+YS**2)/(1.+EPS*X)
          F1  = X * (1. + EPS*X/2.)
          F1Z = Y
          YM1 = YM1 + F1  * ABS(BPL) * DL * WS
          YZ1 = YZ1 + F1Z * ABS(BPL) * DL * WS
          CUR = CUR + BPL * DL * WS
   20   CONTINUE
   10 CONTINUE
      YM1 =  FACTAS * YM1 / (ALFA * XIAB)
      YZ1 =  FACTAS * YZ1/  (ALFA * XIAB)
      WRITE(20,*) ym1,yz1,cur
c-------------------------------- the radius of the current center RJ
      RJ = SQRT(1. + 2.*EPS*YM1) / EPS
      ZJ = YZ1
      YM1a = 0.
      YZ1a = 0.
      YM2 = 0.
      YM3 = 0.
      YM4 = 0.
      YM5 = 0.
      YM6 = 0.
      DO 30 J=1,NP-1
        NELM = J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 40 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)
          XJAC =  XR*YS - XS*YR
          DL = SQRT(XS**2+YS**2)
          BPL = EPS * PSR / XJAC * SQRT(XS**2+YS**2)/(1.+EPS*X)
          RMRJ = (1. + EPS*X - RJ * EPS) / EPS
          F1 = RMRJ * (1. + 0.5*RMRJ/RJ)
          F1Z = Y   * (1. + RMRJ/RJ)
          F2 = F1**2 - F1Z**2
          F3 = F1**3 - 3*F1*F1Z**2 + Y**2/RJ  * F1Z**2
          F4 = F1**4 - 6.*F1**2 * F1Z**2 + F1Z**2 *(3*F1Z**2 - 2*Y**2)
     >       - 0.8 * Y**4 * F1Z**2 / RJ**2
c---------------------------------- function from mathematica mom2.ma
      R1 = (1. + EPS * X)/EPS
      F2=-0.5*R1**2+0.25*R1**4/RJ**2+0.25*RJ**2-1.*R1**2*Y**2/RJ**2
      F3=-0.125*(-1.*R1**6+3.*R1**4*RJ**2-3.*R1**2*RJ**4+1.*RJ**6
     >   +12.*R1**4*Y**2-12.*R1**2*RJ**2*Y**2-8.*R1**2*Y**4)/RJ**3
      F4=0.0625*(1.*R1**8-4.*R1**6*RJ**2+6.*R1**4*RJ**4-4.*R1**2*RJ**6
     >   + 1.*RJ**8 - 24.*R1**6*Y**2 + 48.*R1**4*RJ**2*Y**2
     >   - 24.*R1**2*RJ**4*Y**2
     >   + 48.*R1**4*Y**4-32.*R1**2*RJ**2*Y**4-12.8*R1**2*Y**6)/RJ**4
      F5=-0.03125*(-1.*R1**10 + 5.*R1**8*RJ**2 - 10.*R1**6*RJ**4
     >   + 10.*R1**4*RJ**6-5.*R1**2*RJ**8+1.*RJ**10+40.*R1**8*Y**2
     >   - 120.*R1**6*RJ**2*Y**2 + 120.*R1**4*RJ**4*Y**2
     >   - 40.*R1**2*RJ**6*Y**2 - 160.*R1**6*Y**4
     >   +240.*R1**4*RJ**2*Y**4-80.*R1**2*RJ**4*Y**4+128.*R1**4*Y**6
     >   -64.*R1**2*RJ**2*Y**6 -18.2857*R1**2*Y**8)/RJ**5
      F6=0.015625*(1.*R1**12 - 6.*R1**10*RJ**2 + 15.*R1**8*RJ**4
     >   - 20.*R1**6*RJ**6+15.*R1**4*RJ**8-6.*R1**2*RJ**10+1.*RJ**12
     >   - 60.*R1**10*Y**2 +
     >     240.*R1**8*RJ**2*Y**2 - 360.*R1**6*RJ**4*Y**2 +
     >     240.*R1**4*RJ**6*Y**2-60.*R1**2*RJ**8*Y**2+400.*R1**8*Y**4 -
     >     960.*R1**6*RJ**2*Y**4+720.*R1**4*RJ**4*Y**4 -
     >     160.*R1**2*RJ**6*Y**4-640.*R1**6*Y**6+768.*R1**4*RJ**2*Y**6 -
     >     192.*R1**2*RJ**4*Y**6+274.286*R1**4*Y**8 -
     >     109.714*R1**2*RJ**2*Y**8-24.381*R1**2*Y**10)/RJ**6
c----------------------------------------------------------------------------

          YM1a = YM1a + F1  * ABS(BPL) * DL * WS
          YZ1a = YZ1a + F1Z * ABS(BPL) * DL * WS
          YM2 = YM2 + F2 * ABS(BPL) * DL * WS
          YM3 = YM3 + F3 * ABS(BPL) * DL * WS
          YM4 = YM4 + F4 * ABS(BPL) * DL * WS
          YM5 = YM5 + F5 * ABS(BPL) * DL * WS
          YM6 = YM6 + F6 * ABS(BPL) * DL * WS
   40   CONTINUE
   30 CONTINUE
      YM1a = FACTAS * YM1a / (ALFA * XIAB)
      YZ1a = FACTAS * YZ1a / (ALFA * XIAB)
      YM2  = FACTAS * YM2  / (ALFA * XIAB)
      YM3  = FACTAS * YM3  / (ALFA * XIAB)
      YM4  = FACTAS * YM4  / (ALFA * XIAB)
      YM5  = FACTAS * YM5  / (ALFA * XIAB)
      YM6  = FACTAS * YM6  / (ALFA * XIAB)
      WRITE(20,*) '********************************************'
      WRITE(20,*) '*     CURRENT MOMENTS                      *'
      WRITE(20,*) '********************************************'
      WRITE(20,41) RJ,ZJ
      WRITE(20,42) YM1A, YZ1A
      WRITE(20,43) YM2
      WRITE(20,44) YM3
      WRITE(20,45) YM4
      WRITE(20,46) YM5
      WRITE(20,47) YM6
      WRITE(20,*)
   41 format(' RJ, ZJ : ',2e12.4)
   42 format(' YM1, YZ1 : ',2e12.4)
   43 format(' YM2 : ',e12.4)
   44 format(' YM3 : ',e12.4)
   45 format(' YM4 : ',e12.4)
   46 format(' YM5 : ',e12.4)
   47 format(' YM6 : ',e12.4)
      RETURN
      END

************************************************************************
*DECK BETAPOL
      SUBROUTINE BETAPOL(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                   IGAM,ISOL,BPL,IAS)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE OUTPUT QUANTITIES
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),PSI(4*MAXMNODE)
      REAL XL(NRMMAX),ZPS(NRMMAX)
      REAL ZJAR(NRMMAX),ZPAR(NRMMAX),ZAR(NRMMAX)

      PI = 2. * ASIN(1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      PAREA = 0.
      CAREA = 0.
C------------------------------------- NELM ELEMENTS ------------------
      DO 10 I=NR-1,1,-1
      DO 15 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 20 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 30 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                  R,S,X,XR,XS)
            CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                  R,S,Y,YR,YS)
            CALL INTERP2(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                   PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                  R,S,PS,PSR,PSS)
            XJAC =  XR*YS - XS*YR
            PAREA = PAREA + WR * WS * XJAC * PRES(PS)
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            CAREA = CAREA + WR * WS * XJAC * ARHS
   30     CONTINUE
   20   CONTINUE
   15   CONTINUE
   10 CONTINUE
      IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
        PAREA = 0.5 * FACTAS * A * B * PAREA
        CAREA =  FACTAS * A * CAREA
      ELSE
        PAREA = FACTAS * EPS * A * B * PAREA
        CAREA = FACTAS * A * CAREA
      ENDIF
      BPL = -(8*PI/EPS) * PAREA/ CAREA**2
      RETURN
      END


************************************************************************
*DECK BETAPOL
      SUBROUTINE BETAPOL_PED(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                   IGAM,ISOL,BPL,IAS,PSIPED)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE OUTPUT QUANTITIES
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),PSI(4*MAXMNODE)
      REAL XL(NRMMAX),ZPS(NRMMAX)
      REAL ZJAR(NRMMAX),ZPAR(NRMMAX),ZAR(NRMMAX)

      PI = 2. * ASIN(1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      PAREA = 0.
      CAREA = 0.
C------------------------------------- NELM ELEMENTS ------------------
      DO 10 I=NR-1,1,-1
      DO 15 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 20 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 30 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                  R,S,X,XR,XS)
            CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                  R,S,Y,YR,YS)
            CALL INTERP2(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                   PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                  R,S,PS,PSR,PSS)
            XJAC =  XR*YS - XS*YR

            IF (PS .LT. PSIPED) THEN
             PAREA = PAREA + WR * WS * XJAC * PRES(PSIPED)
			ELSE
             PAREA = PAREA + WR * WS * XJAC * PRES(PS)
			ENDIF
			
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            CAREA = CAREA + WR * WS * XJAC * ARHS
           

   30     CONTINUE
   20   CONTINUE
   15   CONTINUE
   10 CONTINUE
      IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
        PAREA = 0.5 * FACTAS * A * B * PAREA
        CAREA =  FACTAS * A * CAREA
      ELSE
        PAREA = FACTAS * EPS * A * B * PAREA
        CAREA = FACTAS * A * CAREA
      ENDIF
      BPL = -(8*PI/EPS) * PAREA/ CAREA**2
      RETURN
      END



************************************************************************
*DECK CURRENT
      SUBROUTINE CURRENT(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                   IGAM,ISOL,CUR,IAS)
C----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE TOTAL CURRENT
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),PSI(4*MAXMNODE)

      PI = 2. * ASIN(1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      CAREA = 0.
C------------------------------------- NELM ELEMENTS ------------------
      DO 10 I=1,NR-1
      DO 15 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN R -----
        DO 20 NGR=1,4
          R = XGAUSS(NGR)
          WR = WGAUSS(NGR)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 30 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                  R,S,X,XR,XS,XRS,XRR,XSS)
            CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                  R,S,Y,YR,YS,YRS,YRR,YSS)
            CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                  PSI(4*(N4-1)+1),
     >                  R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)
            XJAC =  XR*YS - XS*YR
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = C*DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            CAREA = CAREA + WR * WS * XJAC * ARHS
   30     CONTINUE
   20   CONTINUE
   15 CONTINUE
   10 CONTINUE
      CUR = FACTAS * A * EPS * CAREA
C      WRITE(20,4) CUR
    4 FORMAT(' TOTAL CURRENT (ALFA=1) : ',E12.4)
      RETURN
      END


************************************************************************
*DECK FLXINT
      SUBROUTINE FLXINT(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,
     >                  IGAM,ISOL,IAS,IAV)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE NEW DF2 PROFILE
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      REAL XX(4,MAXMNODE),YY(4,MAXMNODE),PSI(4*MAXMNODE)
      REAL AVC(NRMMAX),XL(NRMMAX),ZPS(NRMMAX)
      REAL DF2TMP(NRMMAX)
      REAL DD1(NRMMAX),DD2(NRMMAX),DD3(NRMMAX),DD4(NRMMAX),ABLTG(3)

      PI = 2. * ASIN(-1.)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
         CUR0 = (1.+B*XAXIS*(1.+EPS*XAXIS/2.))/(1.+EPS*XAXIS)
      ELSE
         CUR0 = (1.+B*(1.+EPS*XAXIS)**2)/(1.+EPS*XAXIS)
      ENDIF
      DF2TMP(1) = 1.
C------------------------------------- NELM ELEMENTS ------------------
      R = -1.
      DO 40 I=1,NR-1
      XLENGTH = 0.
      SUMC1 = 0.
      SUMC2 = 0.
      SUMC3 = 0.
      DO 50 J=1,NP-1
        NELM = (I-1)*(NP-1) + J
        N1 = NODENO(NELM,1)
        N2 = NODENO(NELM,2)
        N3 = NODENO(NELM,3)
        N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)
          XJAC =  XR*YS - XS*YR
          DL = SQRT(XS**2+YS**2)
          BIGR = (1. + EPS * X)
	  ZJDCHI = BIGR * XJAC / DABS(PSR)
          XLENGTH = XLENGTH + DL * WS
          IF (IAV .EQ.0) THEN
	    DLV = DL
	  ELSE
	    DLV = ZJDCHI
	  ENDIF
	  SUMC3 = SUMC3 + WS * DLV
	  SUMC1 = SUMC1 + 1./(1.+EPS*X) * DLV * WS
          IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
             SUMC2 = SUMC2 + X*(1.+EPS*X/2.)/(1.+EPS*X) * DLV * WS
          ELSE
             SUMC2 = SUMC2 + (1.+EPS*X) * DLV * WS
          ENDIF
   60     CONTINUE
   50   CONTINUE
	ZPS(NR-I+1) = PS
        DF2TMP(NR-I+1) = (CUR0*CURPHI(PS)*SUMC3 - B*SUMC2*DPDPSI(PS))
     >       / SUMC1
   40 CONTINUE
c------------------------------------- calc equidistant gamma profile
      DF2TMP(1) = 1.
      ZPS(1) = 0.
      CALL SPLINE(NR,ZPS,DF2TMP,0.,0.,2,DD1,DD2,DD3,DD4)
      DO 70 I=2,NPTS
        IF ((IGAM.EQ.2).OR.(IGAM.EQ.11)) THEN
          SS = REAL(I-1)/REAL(NPTS-1)
          PS = SS*SS
        ELSE
          PS = REAL(I-1)/REAL(NPTS-1)
        ENDIF
        DF2(I) = SPWERT(NR,PS,DD1,DD2,DD3,DD4,ZPS,ABLTG)
   70 CONTINUE
       DF2(1) = 1.
      RETURN
      END

************************************************************************
*DECK FLXINT2
      SUBROUTINE FLXINT2(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,
     >                   CX,CY,IAS,AMPL,SUMDQ)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE NEW DF2 PROFILE
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPIE
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL ZPS(NRMMAX),XPLOT(NRMMAX),QPLIN(NRMMAX)
      REAL DF2OLD(NRMMAX),DF2TMP(NRMMAX),QPLOT(NRMMAX),DELTAQ(NRMMAX)
      REAL DD1(NRMMAX),DD2(NRMMAX),DD3(NRMMAX),DD4(NRMMAX),ABLTG(3)

      SAVE ISAVE

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      CALL PROFILES(DD1,DD2,DD3,DD4,A)
C------------------------------------- NELM ELEMENTS ------------------
      R = -1.
      DO 40 I=1,NR-1
      SUMQ  = 0.
      SUMQR = 0.
      DO 50 J=1,NP-1
        N1 = (I-1)*NP + J
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)
          XJAC =  XR*YS - XS*YR
          BIGR = (1. + EPS * X)
          SUMQ = SUMQ - WS * XJAC / ( BIGR * ABS(PSR))
          ER = XRR*YS + XR*YRS - XRS*YR - XS*YRR
          SUMQR = SUMQR + PSRR * XJAC / ((PSR**2)*BIGR) * WS
          SUMQR = SUMQR - ER / (BIGR*PSR)               * WS
          SUMQR = SUMQR + XJAC*EPS*XR/((BIGR**2)*PSR)   * WS
   60     CONTINUE
   50   CONTINUE
        ZPS(NR-I+1) = PS
        QTMP  = 0.5 * FACTAS * (SUMQ * ALFA / PI) * DD2(NR-I+1)
c        DQTMP = (SUMQR * ALFA /PI) / (0.5/FLOAT(NR-1)) / (2.*SQRT(PS))
c        DQTMP = DD2(NR-I+1)*DQTMP + DD4(NR-I+1)*QTMP
	QPLIN(NR-I+1) = QPRFL(PS)
	DELTAQ(NR-I+1) = QPLIN(NR-I+1) - QTMP
        DF2OLD(NR-I+1) = DGDPSI(PS)
        QPLOT(NR-I+1) = QTMP
        XPLOT(NR-I+1) = SQRT(PS)
   40 CONTINUE
      QTMP0 = DD2(1) * (ALFA/(2.*SQRT(CX*CY)*(1.+EPS*XAXIS)))
      QPLOT(1) = QTMP0
      XPLOT(1) = 0.
      ZPS(1) = 0.
      QPLIN(1) = QPRFL(0.)
      DELTAQ(1) = (QPLIN(1)-QTMP0)
      DF2OLD(1) = 1.

      SUMDQ = 0.
      DF2TMP(1)    = DF2OLD(1)    + 0.25*AMPL*(DELTAQ(1)   -DELTAQ(2))
      ZPS(1) = 0.
      DF2TMP(NR) = DF2OLD(NR)
      DO I=2,NR-1
        DF2TMP(NR-I+1)=DF2OLD(NR-I+1)
     >                + AMPL*(DELTAQ(NR-I)-DELTAQ(NR-I+2))
     >                + AMPL*DELTAQ(NR-I+1)/4.

        DF2TMP(NR-I+1) = DF2TMP(NR-I+1)/DF2TMP(1)
c	WRITE(20,41) XPLOT(NR-I+1),QPLOT(NR-I+1),DELTAQ(NR-I+1),
c     >  DELTAQ(NR-I)-DELTAQ(NR-I+2),DF2TMP(NR-I+1),DF2OLD(NR-I+1)
        SUMDQ = SUMDQ + ABS(DELTAQ(NR-I+1))
      ENDDO
      DF2TMP(1) = 1.
      DF2TMP(NR-1) = DF2OLD(NR-1)+0.25*AMPL*(DELTAQ(NR-2)-DELTAQ(NR-1))

   41 FORMAT(12E12.4)
      SUMDQ = SUMDQ/FLOAT(NR-2)
ccc      WRITE(20,41) ZPS(1),QTMP0,QPRFL(0.),DF2TMP(1),DGDPSI(0.)

c-------------------------------------- still freedom to scale q-profile
ccc      ALFA = ALFA * QPLIN(1)/QPLOT(1)
CCC      WRITE(20,*) ALFA,QPLIN(1),QPLOT(1)

c--------------------------------------- plot q-profile
c      IF (ISAVE.NE.999) THEN
c	ISAVE = 999
c	QMX = FLOAT(INT(MAXVAL(QPLIN)))+1
c        CALL NFRAME(2,1,1,0.,1.,0.,2.0,' ',1,'s',1,'FdF',1)
c        CALL NFRAME(3,1,1,0.,1.,0.,QMX,' ',1,'s',1,'q',1)
c        CALL LPLOT6(2,1,XPLOT,DF2TMP,-NR,'DF2')
c        CALL LINCOL(1)
c        CALL LPLOT6(3,1,XPLOT,QPLIN,-NR,'q-profile')
c        CALL LINCOL(0)
c      ENDIF
c      CALL LPLOT6(2,1,XPLOT,DF2TMP,-NR,'FdF')
c      CALL LPLOT6(3,1,XPLOT,QPLOT,-NR,'q')
c      CALL LINCOL(1)
c      CALL LPLOT6(3,1,XPLOT,QPLIN,-NR,'q-profile')
c      CALL LINCOL(0)
c------------------------------------- calc equidistant gamma profile
      CALL SPLINE(NR,XPLOT,DF2TMP,0.,0.,2,DD1,DD2,DD3,DD4)
      DO 70 I=1,NPTS
        IF ((IGAM.EQ.2).OR.(IGAM.EQ.11)) THEN
          SS = REAL(I-1)/REAL(NPTS-1)
          PS = SS*SS
        ELSE
          PS = REAL(I-1)/REAL(NPTS-1)
	  SS = SQRT(PS)
        ENDIF
        DUMMY = SPWERT(NR,SS,DD1,DD2,DD3,DD4,XPLOT,ABLTG)
	DF2(I) = DUMMY
c        WRITE(20,41) PS,SS,DF2(I),DF2TMP(I)
   70 CONTINUE
      DF21 = DF2(1)
      DO I=1,NPTS
        DF2(I) = DF2(I)/DF21
      ENDDO
      RETURN
      END

************************************************************************
*DECK PROFQ
      SUBROUTINE PROFQ(XAXIS,XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,
     >                   CX,CY,QPRF,ZJZ2,Q95OUT,Q1OUT,IAS)
C----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE THE Q-PROFILE
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPIE
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL ZPS(NRMMAX),QPRF(*),QQ(NRMMAX),ZJZ2(*),AVC(NRMMAX)
      REAL DQ1(NRMMAX),DQ2(NRMMAX),DQ3(NRMMAX),DQ4(NRMMAX),ABLTG(3)
      REAL P0TMP(NRMMAX),DPTMP(NRMMAX),F2TMP(NRMMAX),DFTMP(NRMMAX)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      CALL PROFILES(P0TMP,F2TMP,DPTMP,DFTMP,A)

C------------------------------------- NELM ELEMENTS ------------------
      R = -1.
      DO 40 I=1,NR-1
      SUMQ  = 0.
      AVCUR=0.
      XLENGTH=0.
      DO 50 J=1,NP-1
        N1 = (I-1)*NP + J
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS)
          CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS)
          CALL INTERP2(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS)
          XJAC =  XR*YS - XS*YR
          BIGR = (1. + EPS * X)
          SUMQ = SUMQ - WS * XJAC / ( BIGR * ABS(PSR))

          XLENGTH = XLENGTH + SQRT(XS**2 + YS**2) * WS
          IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
            ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
          ELSE
            ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
          ENDIF
          ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            AVCUR = AVCUR + WS * ARHS * SQRT(XS**2 + YS**2)
   60     CONTINUE
   50   CONTINUE
        ZPS(NR-I+1) = PS
        QQ(NR-I+1) = FACTAS/2. * F2TMP(NR-I+1) * SUMQ * ALFA / PI
        AVC(NR-I+1) = FACTAS * A * EPS * AVCUR/ (ALFA * XLENGTH)
   40 CONTINUE
c------------------------------------- calc equidistant q profile
      QQ(1) = F2TMP(1)* ALFA/(2.*SQRT(CX*CY)*(1.+EPS*XAXIS))
      ZPS(1) = 0.
      CALL SPLINE(NR,ZPS,QQ,0.,0.,2,DQ1,DQ2,DQ3,DQ4)
      DO I=1,NPTS
c        IF ((IGAM.EQ.2).OR.(IGAM.EQ.11)) THEN
          SS = REAL(I-1)/REAL(NPTS-1)
          PS = SS*SS
c        ELSE
c          PS = REAL(I-1)/REAL(NPTS-1)
c        ENDIF
        QPRF(I) = SPWERT(NR,PS,DQ1,DQ2,DQ3,DQ4,ZPS,ABLTG)
      ENDDO
      PS = 0.95
      Q95OUT = SPWERT(NR,PS,DQ1,DQ2,DQ3,DQ4,ZPS,ABLTG)
      Q1OUT  = QQ(NR)
c------------------------------------- calc equidistant current profile
c      AVC(1) = AVC(2) - (AVC(3)-AVC(2))
c     >        /(SQRT(ZPS(3))-SQRT(ZPS(2)))*SQRT(ZPS(2))
      AVC(1) = AVC(2) - (AVC(3)-AVC(2))
     >        /(ZPS(3)-ZPS(2))*ZPS(2)

      CALL SPLINE(NR,ZPS,AVC,0.,0.,2,DQ1,DQ2,DQ3,DQ4)
      DO I=1,NPTS
c        IF ((IGAM.EQ.2).OR.(IGAM.EQ.11)) THEN
          SS = REAL(I-1)/REAL(NPTS-1)
          PS = SS*SS
c        ELSE
c          PS = REAL(I-1)/REAL(NPTS-1)
c        ENDIF
         ZJZ2(I) = SPWERT(NR,PS,DQ1,DQ2,DQ3,DQ4,ZPS,ABLTG)
      ENDDO
      RETURN
      END


************************************************************************
*DECK MERCIER
      SUBROUTINE MERCIER(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,
     >                   DIMERC,DRMERC,HH,QQ,DQ,GEONC,ZJPAR)
C----------------------------------------------------------------------
C SUBROUTINE FOR THE IDEAL AND RESISTIVE MERCIER CRITERION
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPIE
      USE COMPLO
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL P0TMP(NRMMAX),DPTMP(NRMMAX),FTMP(NRMMAX),DFTMP(NRMMAX)
      REAL ZJJ1(NRMMAX),ZJJ2(NRMMAX),ZJJ3(NRMMAX),ZJJ4(NRMMAX)
      REAL ZJJ5(NRMMAX),ZJJ6(NRMMAX),DJJ5(NRMMAX)
      REAL QQ(*),DQ(*),ZPS(NRMMAX),BIGG(NRMMAX)
      REAL DRMERC(*),DIMERC(*),HH(*),GEONC(*),ZJPAR(*)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      R0 = 1.
      B0 = 1.
      RADIUS = EPS * R0
      CPSURF = RADIUS**2 * B0 / ALFA

      CALL PROFILES(P0TMP,FTMP,DPTMP,DFTMP,A)

      PSCALE  = B0**2 * EPS / ALFA**2
      RBSCALE = B0 * R0
      DO I = 1,NR
        P0TMP(I) = P0TMP(I) * PSCALE
        DPTMP(I) = DPTMP(I) * PSCALE
        FTMP(I)  = FTMP(I)  * RBSCALE
        DFTMP(I) = DFTMP(I) * RBSCALE
      ENDDO

      WRITE(20,1)
      WRITE(20,2)
      WRITE(20,1)
      WRITE(20,*) '*I,    SS,     D_I,     D_R,     H,'
     >          //'       Q,     SHEAR, GEONC,    PEEL *'
      WRITE(20,1)
    1 FORMAT(' ',72('*'))
    2 FORMAT(' * IDEAL AND RESISTIVE MERCIER CRITERION ',26(' '),'*')
    3 FORMAT(I3,F8.3,1P2E10.2,0P3F8.3,1PE10.2,1P3E12.4)
c----------------------------------- flux surface integrals
      R = -1.
      DO 40 I=1,NR-1
      SUMQ  = 0.
      SUMJ1 = 0.
      SUMJ2 = 0.
      SUMJ3 = 0.
      SUMJ4 = 0.
      SUMJ5 = 0.
      SUMJ6 = 0.
      SUMJ5R= 0.
      SUMQR = 0.
      SUMB0  = 0.
      SUMBOP = 0.
      SUMOR2 = 0.
      SUMPE  = 0.
      DO 50 J=1,NP-1
        N1 = (I-1)*NP + J
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)

          XJAC =  XR*YS - XS*YR
          XJACR = XRR*YS + XR*YRS - XRS*YR - XS*YRR

          DL = SQRT(XS**2+YS**2)
          BIGR = (1. + EPS * X)
          GRADPS2 = PSR**2 * (XS**2 + YS**2) / XJAC**2
          DSSDR = ABS(PSR) / (2.*SQRT(PS))
c------------------------------------------------ normalisations
          GRADPS2 = GRADPS2 * (CPSURF/RADIUS)**2
	  XJAC    = XJAC * RADIUS**2
	  BIGR    = BIGR
	  DL      = RADIUS * DL
	  PSR     = CPSURF * PSR
          PSRR    = CPSURF * PSRR
	  XJACR   = XJACR  * RADIUS**2
c---------------------------------------------------------------

	  ZJDCHI = BIGR * XJAC / ABS(PSR)

	  SUMJ1 = SUMJ1 - WS * ZJDCHI / (GRADPS2 * BIGR**2)
	  SUMJ2 = SUMJ2 - WS * ZJDCHI / (GRADPS2)
	  SUMJ3 = SUMJ3 - WS * ZJDCHI * (BIGR**2 / GRADPS2)
	  SUMJ4 = SUMJ4 - WS * ZJDCHI / BIGR**2
	  SUMJ5 = SUMJ5 - WS * ZJDCHI
	  SUMJ6 = SUMJ6 - WS * ZJDCHI * GRADPS2 / BIGR**2
          SUMQ  = SUMQ  - WS * XJAC / ( BIGR * ABS(PSR))

          SUMQR = SUMQR + PSRR * XJAC / ((PSR**2)*BIGR) * WS
          SUMQR = SUMQR - XJACR / (BIGR*PSR)            * WS
          SUMQR = SUMQR + XJAC*EPS*XR/((BIGR**2)*PSR)   * WS

	  SUMJ5R = SUMJ5R + XJAC  * BIGR * PSRR / (PSR**2) * WS
	  SUMJ5R = SUMJ5R - XJACR * BIGR / PSR             * WS
	  SUMJ5R = SUMJ5R - XJAC  * EPS *XR / PSR          * WS

c----------------------------------------- D_NC from Hegna's

          B02 = (FTMP(NR-I+1)/BIGR)**2 + GRADPS2/BIGR**2

          SUMB0  = SUMB0  - WS * ZJDCHI * B02
          SUMBOP = SUMBOP - WS * ZJDCHI * B02 / GRADPS2
	  SUMOR2 = SUMOR2 - WS * ZJDCHI / BIGR**2

c----------------------------------------- Peeling mode stability

          SUMPE = SUMPE + WS * ZJDCHI * BIGR / GRADPS2**(1.5)
     >          * ( FTMP(NR-I+1)*DPTMP(NR-I+1) + DFTMP(NR-I+1)*B02 )

c          IF (I .EQ. 1) THEN
c
c            WRITE(*,'(8e16.8)') X,Y,GRADPS2,BIGR,
c     >                 FTMP(NR-I+1)*DPTMP(NR-I+1) + DFTMP(NR-I+1)*B02
c          ENDIF

   60     CONTINUE
   50   CONTINUE
        NI = NR-I+1
        ZPS(NI) = SQRT(PS)
        ZJJ1(NI) = FACTAS * SUMJ1 / (2.*PI)
	ZJJ2(NI) = FACTAS * SUMJ2 / (2.*PI)
        ZJJ3(NI) = FACTAS * SUMJ3 / (2.*PI)
	ZJJ4(NI) = FACTAS * SUMJ4 / (2.*PI)
        ZJJ5(NI) = FACTAS * SUMJ5 / (2.*PI)
	ZJJ6(NI) = FACTAS * SUMJ6 / (2.*PI)

	DJJ5(NI) = FACTAS * SUMJ5R/DSSDR  / (2.*PI)

	QQ(NI) = FACTAS * FTMP(NI) * SUMQ / (2.*PI)
	DQ(NI) = DFTMP(NI)*SUMQ + FTMP(NI)*SUMQR/DSSDR
	DQ(NI) = DQ(NI) * FACTAS / (2.*PI)

	OR2AV = FACTAS * SUMOR2 / (2.*PI) / ZJJ5(NI)
	B02AV = FACTAS * SUMB0  / (2.*PI) / ZJJ5(NI)
	GBAR  = FACTAS * SUMBOP / (2.*PI)

        GEONC(NI) = GBAR / B02AV

        ZJPAR(NI) = (- DPTMP(NI) * FTMP(NI) - DFTMP(NI) * B02AV)
     >            /( 2.*CPSURF*SQRT(PS))

        BIGG(NI) = FACTAS * SUMBOP / (2.*PI)

	DIMERC(NI) = (DPTMP(NI)*FTMP(NI)*ZJJ2(NI)/DQ(NI) - 0.5)**2
     >     +  DPTMP(NI)/DQ(NI)**2 * (DJJ5(NI) - DPTMP(NI)*ZJJ3(NI))
     >     * (FTMP(NI)**2 * ZJJ1(NI) + ZJJ4(NI) )

	HH(NI) = FTMP(NI) * DPTMP(NI)/DQ(NI) * (ZJJ2(NI)
     >    - ZJJ5(NI)*(ZJJ4(NI)+FTMP(NI)**2 * ZJJ1(NI))
     >    / (ZJJ6(NI) + FTMP(NI)**2 * ZJJ4(NI) ) )

        DRMERC(NI) = DIMERC(NI) - (HH(NI) - 0.5)**2

        SHEAR = SQRT(PS) * DQ(NI) / QQ(NI)

        DM = - DIMERC(NI) + 0.25

        PEEL = 0.
        SUMPE = SUMPE * FACTAS / (2.*PI)

        IF (DM .LT. 0.25) THEN
          PEEL  = 1. + 2./DQ(NI) * SUMPE - SQRT(1.- 4.*DM)
        ENDIF

	WRITE(20,3) I,SQRT(PS),-DIMERC(NI),-DRMERC(NI),HH(NI),QQ(NI),
     >              SHEAR,GEONC(NI),DM,PEEL,SUMPE

   40 CONTINUE
      WRITE(20,1)
   41 FORMAT(I3,1P15E11.3)
      IF (NPL1.NE.0) THEN
        CALL LBLTOP('Mercier Criteria',16)
        CALL LPLOT(2,2,1,ZPS(5),DIMERC(5),NR-4,1,'-DI',3,'S',1,' ',1)
        CALL LPLOT(2,3,1,ZPS(5),DRMERC(5),NR-4,1,'-DR',3,'S',1,' ',1)
        CALL LPLOT(3,2,1,ZPS(5),HH(5),NR-4,1,'H',1,'S',1,' ',1)
        CALL LPLOT(3,3,1,ZPS(5),QQ(5),NR-4,1,'Q',1,'S',1,' ',1)
      ENDIF

C----------------------- write data for analysis of Hegna's paper to file
      WRITE(41,*)  ' NR :'
      WRITE(41,*)  NR-2
      WRITE(41,*) ' DI, DR, H, GEONC : '
      DO I=2,NR-1
        WRITE(41,42) I,DIMERC(I),DRMERC(I),HH(I),GEONC(I)
      ENDDO

   42 FORMAT(I4,1P12E12.4)

      RETURN
      END
************************************************************************
*DECK CIRCUL
      SUBROUTINE CIRCUL(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,
     >                  FCIRC,B02AV,B0MAX,RAV)
C----------------------------------------------------------------------
C SUBROUTINE FOR THE FRACTION OF CIRCULATING PARTICLES
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPIE
      USE COMPLO
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL FTMP(NRMMAX),DFTMP(NRMMAX),P0TMP(NRMMAX),DPTMP(NRMMAX)
      REAL B0MAX(*),B02AV(*),FCIRC(*),RAV(*),SUMK(NRMMAX)

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      R0 = 1.
      B0 = 1.
      RADIUS = EPS * R0
      CPSURF = RADIUS**2 * B0 / ALFA

      CALL PROFILES(P0TMP,FTMP,DPTMP,DFTMP,A)

      RBSCALE = B0 * R0
      DO I = 1,NR
        FTMP(I)  = FTMP(I)  * RBSCALE
        DFTMP(I) = DFTMP(I) * RBSCALE
      ENDDO

c---------------------- find Bmax and B^2 average on every surface
      R = -1.
      DO 40 I=1,NR-1

      SUM1 = 0.
      SUM2 = 0.
      SUMR = 0.
      B02MAX = 0.

      DO 50 J=1,NP-1
        N1 = (I-1)*NP + J
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 60 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)

          XJAC =  XR*YS - XS*YR
          BIGR = (1. + EPS * X)
          GRADPS2 = PSR**2 * (XS**2 + YS**2) / XJAC**2
c------------------------------------------------ normalisations
          GRADPS2 = GRADPS2 * (CPSURF/RADIUS)**2
	  XJAC    = XJAC * RADIUS**2
	  BIGR    = BIGR
	  PSR     = CPSURF * PSR
c---------------------------------------------------------------
	  ZJDCHI = BIGR * XJAC / ABS(PSR)
          B02 = (FTMP(NR-I+1)/BIGR)**2 + GRADPS2/BIGR**2
	  SUM1  = SUM1 - WS * ZJDCHI
          SUM2  = SUM2 - WS * ZJDCHI * B02
          SUMR  = SUMR - WS * ZJDCHI * BIGR

          IF (B02.GT.B02MAX) B02MAX = B02
   60     CONTINUE
   50   CONTINUE
        NI = NR-I+1
        SUM1 = FACTAS * SUM1 / (2.*PI)
	SUM2 = FACTAS * SUM2 / (2.*PI)
        SUMR = FACTAS * SUMR / (2.*PI)
	B02AV(NI) = SUM2 / SUM1
        B0MAX(NI) = SQRT(ABS(B02MAX))
	RAV(NI)   = SUMR / SUM1
   40 CONTINUE
   41 FORMAT(I3,1P15E11.3)

c---------------------------------- calculate average term in integral
      NK=101
      R = -1.
      DO 140 I=1,NR-1
        SUM1 = 0.
        DO K=1,NK
          SUMK(K) = 0.
        ENDDO
      DO 150 J=1,NP-1
        N1 = (I-1)*NP + J
        N2 = N1 + 1
        N3 = N2 + NP
        N4 = N1 + NP
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
        DO 160 NGS=1,4
          S  = XGAUSS(NGS)
          WS = WGAUSS(NGS)
          CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                R,S,X,XR,XS,XRS,XRR,XSS)
          CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                R,S,Y,YR,YS,YRS,YRR,YSS)
          CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),PSI(4*(N3-1)+1),
     >                PSI(4*(N4-1)+1),
     >                R,S,PS,PSR,PSS,PSRS,PSRR,PSSS)

          XJAC =  XR*YS - XS*YR
          BIGR = (1. + EPS * X)
	  GRADPS2 = PSR**2 * (XS**2 + YS**2) / XJAC**2
c------------------------------------------------ normalisations
	  GRADPS2 = GRADPS2 * (CPSURF/RADIUS)**2
	  XJAC    = XJAC * RADIUS**2
	  BIGR    = BIGR
	  PSR     = CPSURF * PSR
c---------------------------------------------------------------
	  ZJDCHI = BIGR * XJAC / ABS(PSR)
          B02 = (FTMP(NR-I+1)/BIGR)**2 + GRADPS2/BIGR**2
	  B0 = SQRT(B02)
	  BM = B0MAX(NR-I+1)
	  DO K = 1, NK
	    ZLAM = FLOAT(K-1)/FLOAT(NK-1)
	    SUMK(K) = SUMK(K) - WS * ZJDCHI * SQRT(ABS(1.-ZLAM*B0/BM))
	  ENDDO
          SUM1  = SUM1 - WS * ZJDCHI
  160     CONTINUE
  150   CONTINUE
        NI = NR-I+1

c------------------------------------------ integrate over lambda
        SUM = 0.
        DO K=2,NK-1
	  ZLAM = FLOAT(K-1)/FLOAT(NK-1)
	  SUM = SUM + ZLAM * SUM1 / SUMK(K)
	ENDDO
	FCIRC(NI) = (SUM + 0.5 * SUM1/SUMK(NK)) / FLOAT(NK-1)
	FCIRC(NI) = 0.75 * FCIRC(NI)  * B02AV(NI) / B0MAX(NI)**2
c        WRITE(*,41) I,SQRT(PS),B02AV(NI),B0MAX(NI),FCIRC(NI)
  140 CONTINUE
c      write(*,*) ' DONE CIRCUL'
      RETURN
      END

************************************************************************
*DECK NEO
      SUBROUTINE NEO(FCIRC,TE,TI,ZZNE,Z,F,Q,R,EPS,ZDNE,DTE,DTI,
     >               ZNUE,ZNUI,SIGSPITZ,SIGNEO,ZJBT,HJBT)
C-----------------------------------------------------------------------
C formula for neoclassical quantities from Phys. Plasmas 6 (1999) 2834.
C Te, Ti given in electron Volt, density im m^-3
C-----------------------------------------------------------------------
      FT = 1. - FCIRC

      ZNE  = ZZNE * 1E19
      ZNI  = ZNE  / Z
      ZDNI = ZDNE / Z
      DNI  = ZDNI * 1E19

      RPE = TE / (TI+TE)

      ZLE = 31.3 - LOG(SQRT(ZNE) / TE)
      ZLI = 30.0 - LOG(Z**3 * SQRT(ZNI) / TI**1.5)


      ZNUE = 6.921E-18 * q * R * ZNE * Z    * ZLE /(TE**2 * EPS**1.5)
      ZNUI = 4.90E-18  * q * R * ZNI * Z**4 * ZLI /(TI**2 * EPS**1.5)

      F33TEF = FT / (1. + (0.55-0.1*FT)*SQRT(ZNUE)
     >       + 0.45*(1.-FT)*ZNUE/Z**1.5)

      ZNZ = 0.58 + 0.74/(0.76+Z)

      SIGSPITZ = 1.9012E4 * TE**1.5 / (Z * ZNZ * ZLE)

      X=F33TEF
      SIGNEO = 1. - (1.+0.36/Z)*X + 0.59/Z * X*X - 0.23/Z * X**3
      SIGNEO = SIGNEO * SIGSPITZ

      F31TEF = FT / (1. + (1.-0.1*FT)*SQRT(ZNUE)
     >                  + 0.5*(1.-FT)*ZNUE/Z)

      X = F31TEF
      ZL31 = (1. + 1.4/(Z+1.))*X - 1.9/(Z+1.)*X*X + 0.3/(Z+1)* X**3
     >     + 0.2/(Z+1)*X**4

      F32TEFE = FT/(1. + 0.26*(1.-FT)*SQRT(ZNUE)
     >                 + 0.18*(1.-0.37*FT)*ZNUE/SQRT(Z) )

      F32TEFI = FT/(1. + (1.+0.6*FT)*SQRT(ZNUE)
     >                 + 0.85*(1.-0.37*FT)*ZNUE*(1.+Z))

      X = F32TEFE
      F32EE = (0.05+0.62*Z)/(Z*(1.+0.44*Z)) * (X - X**4)
     >      + 1./(1.+0.22*Z)*(X**2 - X**4 - 1.2*(X**3 - X**4))
     >      + 1.2/(1.+0.5*Z)*X**4

      Y = F32TEFI
      F32EI = - (0.56+1.93*Z)/(Z*(1.+0.44*Z)) * (Y - Y**4)
     >      + 4.95/(1.+2.48*Z) * (Y**2 - Y**4 - 0.55*(Y**3-Y**4))
     >      - 1.2/(1.+0.5*Z)*Y**4

      ZL32 = F32EE + F32EI

      F34TEF = FT/(1.+ (1.-0.1*FT)*SQRT(ZNUE) + 0.5*(1.-0.5*FT)*ZNUE/Z)

      A0 = - 1.17*(1.-FT)/(1.-0.22*FT-0.19*FT*FT)

c----------------------------- sign in front of 0.315 changed !!!

      ANU = ((A0 + 0.25*(1.-FT*FT)*SQRT(ZNUI))/(1.+0.5*SQRT(ZNUI))
     >    + 0.315*ZNUI**2 * FT**6)
     >    /(1. + 0.15*ZNUI**2*FT**6)

      X = F34TEF
      ZL34 = (1. + 1.4/(Z+1.))*X - 1.9/(Z+1.)*X*X + 0.3/(Z+1)* X**3
     >     + 0.2/(Z+1)*X**4

      P = (Z*ZNI*TI+ZNE*TE) * 1.602E-19
      ZJBT = - F * P * (ZL31 * DNI/ZNI + RPE*(ZL31+ZL32)*DTE/TE
     >       + (1.-RPE)*(1.+ZL34/ZL31 * ANU) * ZL31 * DTI/TI)

c------------------------------ Hirshman formula at low collisionality
      XX = (1.-FCIRC)/FCIRC

      DX = 1.414*Z + Z*Z + XX*(0.754 + 2.657*Z + 2*Z*Z)
     >                + XX*XX*(0.348 + 1.243*Z + Z*Z)

      ALFAI = -1.172/(1+0.462*XX)

      IF (DNI.NE.0.) THEN
         ETAI = (dTE/TE) / (dNI/ZNI)
      ELSE
         ETAI = 1.E10
      ENDIF

      A1 = (Z*DNI*(TE+TI) + Z*ZNI*(DTE+DTI)) * 1.6022E-19

      A2E = ETAI/(1.+ETAI) * Z/(Z+1.) * A1
      A2I = A2E

      HL31 = F * XX * (0.754 + 2.210*Z + Z*Z
     >                 + XX * (0.348 + 1.243*Z + Z*Z) ) / DX

      HL32 = - F * XX * (0.884 + 2.074*Z) / DX

      HJBT = -HL31 * ( A1 + ALFAI/Z * A2I) - HL32 * A2E

      RETURN
      END

      SUBROUTINE EQDSK(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,IAS,CUR,
     >                 R0,B0,QPROF)
!-----------------------------------------------------------------------
! WRITE EQDSK FILE, BUT WITHOUT THE FLUXSURFACES, WHICH ARE SET TO 0.
!-----------------------------------------------------------------------
      USE PARAM
      USE COMPIE
      REAL XX(4,*),YY(4,*),PSI(*),QPROF(*)
      REAL P0TMP(NRMMAX),DPTMP(NRMMAX),FTMP(NRMMAX),DFTMP(NRMMAX)

      DO I=1,NR
        FLUX = FLOAT(I-1)/FLOAT(NR-1)
        IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
          P0TMP(I) =  .5 * A * B * PRES(FLUX)
          DPTMP(I) = -.5 * A * B * DPDPSI(FLUX)
          DGAM = - A * DGDPSI(FLUX)
          GAM  =   A * XGAMMA(FLUX)
          FTMP(I) = P0TMP(I)  - EPS * GAM
          DFTMP(I) = DPTMP(I) - EPS * DGAM
        ELSE
          P0TMP(I) = EPS * A * B * PRES(FLUX)
          DPTMP(I) = -EPS* A * B * DPDPSI(FLUX)
          FTMP(I) = - EPS * A * C * XGAMMA(FLUX)
          DFTMP(I) =  EPS * A * C * DGDPSI(FLUX)
        ENDIF
        FTMP(I)  = SQRT( 1. - 2.*EPS*FTMP(I) /ALFA**2)
        DFTMP(I) = - 1./(2.*FTMP(I)) * 2.*EPS*DFTMP(I)/ALFA**2
      ENDDO

      RADIUS = EPS * R0
      CPSURF = RADIUS**2 * B0 / ALFA

      ZMU0 = 4E-7*PI

      PSCALE  = B0**2 * EPS / (ZMU0 * ALFA**2)
      RBSCALE = B0 * R0

      DO I = 1,NR
        P0TMP(I) = P0TMP(I) * PSCALE
        DPTMP(I) = DPTMP(I) * PSCALE / CPSURF
        FTMP(I)  = FTMP(I)  * RBSCALE
        DFTMP(I) = FTMP(I)*(DFTMP(I) * RBSCALE / CPSURF)
      ENDDO

      NP2 = NP
      IF (IAS.EQ.1) NP2 = (NP-1)/2 + 1
      REAST = XX(1,1)  * EPS*R0
      RWEST = XX(1,NP2)* EPS*R0
      RGRIDL = REAST-RWEST

      ZMID   = YY(1,NR*NP)*EPS*R0
      RMAXIS = (1.+XX(1,NR*NP)*EPS)*R0
      ZMAXIS = ZMID
      XIP = CUR * EPS*R0*B0/ZMU0

      OPEN(11,FILE='EQDSK')
   11 FORMAT(5E17.9)
      WRITE(11,*) ' HELENA PRODUCED EQDSK FILE ',NR,NP
      WRITE(11,11) 0.,0.,R0,RGRIDL,ZMID
      WRITE(11,11) RMAXIS,ZMAXIS,0.,0.,B0
      WRITE(11,11) XIP,0.,0.,RMAXIS,0.
      WRITE(11,11) ZMAXIS,0.,0.,0.,0.

      WRITE(11,11) (FTMP(I),I=1,NR)
      WRITE(11,11) (P0TMP(I),I=1,NR)
      WRITE(11,11) (DFTMP(I),I=1,NR)
      WRITE(11,11) (DPTMP(I),I=1,NR)
      WRITE(11,11) (0.,I=1,NR*NP)
      WRITE(11,11) (QPROF(I),I=1,NR)

      IF (IAS.EQ.1) THEN
      WRITE(11,*) NP-1,1
      WRITE(11,11) (R0*(1.+EPS*XX(1,J)),EPS*R0*YY(1,J),J=1,NP-1)
      ELSE
      WRITE(11,*) 2*(NP-1),1
      WRITE(11,11) (R0*(1.+EPS*XX(1,J)),EPS*R0*YY(1,J),J=1,NP)
      WRITE(11,11) (R0*(1.+EPS*XX(1,NP-J)),-EPS*R0*YY(1,NP-J),J=1,NP-2)
      ENDIF
      CLOSE(11)
!---------------------------------------------- END OF EQDSK WRITING
      RETURN
      END


************************************************************************
*DECK WORLD
      SUBROUTINE WORLD(XX,YY,PSI,NR,NP,A,B,C,EPS,ALFA,IGAM,ISOL,IAS,
     >                 R0,B0,ZN0,ZEFF,RPE,ETAEI,
     >                 FCIRC,QPROF,DQPROF,RAV,GEONC,
     >                 ZJPAR,ZJAR,DIMERC,DRMERC,HH)
C----------------------------------------------------------------------
C SUBROUTINE TO TRANSLATE BACK TO REAL WORLD VALUES USING THE MAJOR
C RADIUS (m) , VACUUM MAGNETIC FIELD (T) AT THE GEOMETRIC CENTRE,
C AND THE CENTRAL DENSITY (10^19 m-3), RPE = Te/(Te+Ti)
C CALCULATES BOOTSTRAP CURRENT PROFILE.
C----------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE GAUSINT
      USE COMOUT
      USE COMPROF
      USE COMPIE
      USE COMPLO
      REAL XX(4,*),YY(4,*),PSI(*),FCIRC(*),QPROF(*),RAV(*),GEONC(*)
      REAL DIMERC(*),DRMERC(*),HH(*),DQPROF(*),ZJPAR(*),ZJAR(*)
      REAL P0TMP(NRMMAX),DPTMP(NRMMAX),FTMP(NRMMAX),DFTMP(NRMMAX)
      REAL TE(NRMMAX),TI(NRMMAX),ZNI(NRMMAX),ZJBT(NRMMAX),SR(NRMMAX)
      REAL ZNUE(NRMMAX),ZNUI(NRMMAX),SIGNEO(NRMMAX),SIGSPITZ(NRMMAX)
      REAL HJBT(NRMMAX),AVC(NRMMAX),SP(NRMMAX),STAB(NRMMAX)
      REAL DNC(NRMMAX),DQDPSI(NRMMAX),QP(5), XMID(NRMMAX)
      REAL ZJBTCUT(NRMMAX), ZJPARCUT(NRMMAX)
      INTEGER IRIGHT

      FACTAS = 2.
      IF (IAS.EQ.1) FACTAS=1.

      RADIUS = EPS * R0
      CPSURF = RADIUS**2 * B0 / ALFA

      CALL PROFILES(P0TMP,FTMP,DPTMP,DFTMP,A)

      ZMU0 = 4E-7*PI
      ZE   = 1.6022E-19

      PSCALE  = B0**2 * EPS / (ZMU0 * ALFA**2)
      RBSCALE = B0 * R0

      DO I = 1,NR
        P0TMP(I) = P0TMP(I) * PSCALE
        DPTMP(I) = DPTMP(I) * PSCALE
        FTMP(I)  = FTMP(I)  * RBSCALE
        DFTMP(I) = DFTMP(I) * RBSCALE
      ENDDO

      WRITE(20,*)
      WRITE(20,*) '*************************************************'
      WRITE(20,*) ' REAL WORLD OUTPUT :'
      WRITE(20,*) '*************************************************'
      WRITE(20,1)  R0
      WRITE(20,2)  B0
      WRITE(20,3)  ZN0
      WRITE(20,4)  ZEFF
      WRITE(20,6)  RPE
      WRITE(20,*) '*************************************************'
    1 FORMAT(' MAJOR RADIUS    : ',F9.4,' [m]')
    2 FORMAT(' MAGNETIC FIELD  : ',F9.4,' [T]')
    3 FORMAT(' CENTRAL DENSITY : ',F9.4,' 10^19 [m^-3]')
    4 FORMAT(' ZEFF            : ',F9.4)
    6 FORMAT(' TE/(TE+TI)      : ',F9.4)

      WRITE(20,*) '**************************************************'
      WRITE(20,*) 'R[m], S,P [Pa], Ne [10^19m^-3], Te [eV],  Ti [eV],'
      WRITE(20,*) '**************************************************'

      DO I=2,NR-1
        PS = PSI(4*(NR-I)*NP + 1)
        SS = SQRT(PS)

        IRIGHT= (NR-I)*NP + 1
        XMID(I) = R0 * ( 1.d0 + EPS * XX(1,IRIGHT))

c------------------------------------ to be elaborated
c        ZNI(I) =  ZN0 * (1. - PS)
c	DNI    = -ZN0

ccc        ETAEI = 1.
        PWR = 1./(ETAEI+1.)

        ZNI(I) = ZN0 * (P0TMP(I)/P0TMP(1)) ** PWR
        DNI = ZN0 * PWR * ABS((P0TMP(I)/P0TMP(1))) ** (PWR-1.)
     >      * DPTMP(I)/P0TMP(1) / (2.*SS)

ccc	ZNI(I) =  ZN0
ccc	DNI = 0.
        TE(I)  = RPE      * P0TMP(I) / (ZE * ZNI(I)*1E19)
	TI(I)  = (1.-RPE) * P0TMP(I) / (ZE * ZNI(I)*1E19)

	DTE = RPE/(ZE*1E19) * ( DPTMP(I)/(2.*SS * ZNI(I))
     >                        - P0TMP(I)*DNI/ZNI(I)**2)
	DTI = (1.-RPE)/(ZE*1E19) * ( DPTMP(I)/(2.*SS * ZNI(I))
     >                             - P0TMP(I)*DNI/ZNI(I)**2)

	DTE = DTE / CPSURF
	DTI = DTI / CPSURF
	DNI = DNI / CPSURF

	RR = R0 * RAV(I)

	CALL NEO(FCIRC(I),TE(I),TI(I),ZNI(I),ZEFF,
     >           FTMP(I),QPROF(I),RR,
     >           EPS*SS,DNI,DTE,DTI,
     >           ZZNUE,ZZNUI,ZSIGSP,ZSIGNEO,ZZJBT,ZHJBT)

        ZJBT(I) = ZZJBT
        HJBT(I) = ZHJBT
        SIGNEO(I) = ZSIGNEO
        SIGSPITZ(I) = ZSIGSP
        ZNUE(I) = ZZNUE
        ZNUI(I) = ZZNUI
        SR(I) = SS

        ZJPAR(I) = ZJPAR(I) / (ZMU0 * R0 / B0**2)



	WRITE(20,7) XMID(I), PS, SQRT(PS),P0TMP(I),ZNI(I),TE(I),
     >      TI(I),ZJBT(I),HJBT(I),HJBT(I)*ZMU0*R0/B0**2,ZJPAR(I),ZJAR(I)
      ENDDO
c------------------------------- flux averaged current profile
      R = -1.
      DO 40 I=2,NR-1
        XLENGTH = 0.
        AVCUR = 0.
        DO 50 J=1,NP-1
          NELM = (I-1)*(NP-1) + J
          N1 = NODENO(NELM,1)
          N2 = NODENO(NELM,2)
          N3 = NODENO(NELM,3)
          N4 = NODENO(NELM,4)
C------------------------------------- 4 POINT GAUSSIAN INT. IN S -----
          DO 60 NGS=1,4
            S  = XGAUSS(NGS)
            WS = WGAUSS(NGS)
            CALL INTERP2(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                   R,S,X,XR,XS)
            CALL INTERP2(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                   R,S,Y,YR,YS)
            CALL INTERP1(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                   PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                   R,S,PS)
            XJAC =  XR*YS - XS*YR
            XLENGTH = XLENGTH + SQRT(XS**2 + YS**2) * WS
            IF ((IGAM.GE.1).AND.(IGAM.LE.4)) THEN
              ARHS = DGDPSI(PS) + B*X*(1.+EPS*X/2.)*DPDPSI(PS)
            ELSE
              ARHS = C*DGDPSI(PS) + B*(1.+EPS*X)**2 * DPDPSI(PS)
            ENDIF
            ARHS =  ARHS / (1.+EPS*X)
C---------------------------------------------- SOLOVIEV RHS ----------
            IF (ISOL.EQ.1) ARHS = (1.+ B*X*(1.+EPS*X/2.))/(1.+EPS*X)
            AVCUR = AVCUR + WS * ARHS * SQRT(XS**2 + YS**2)
   60     CONTINUE
   50   CONTINUE
        XLENGTH = FACTAS * XLENGTH
        AVC(NR-I+1) = FACTAS * A * EPS * AVCUR/ (ALFA * XLENGTH)
	AVC(NR-I+1) = AVC(NR-I+1) * B0 / (ZMU0 * R0)
   40 CONTINUE
      WRITE(20,*)
      WRITE(20,*) '************************************************'
     >            //'***********************'
      WRITE(20,*) '*   S,       Q,      Fcirc,     NU_e,     NU_i, '
     >            //'  SIG(Spitz),SIG(neo) *'
      WRITE(20,*) '************************************************'
     >            //'***********************'
      DO I=2,NR-1
        WRITE(20,7) SR(I),QPROF(I),FCIRC(I),ZNUE(I),ZNUI(I),
     >              SIGSPITZ(I),SIGNEO(I),AVC(NR-I+1)
      ENDDO

      DO I=2,NR-1
      ZJBTCUT(I) = ZJBT(I)*0.5d0*(1.d0+TANH((SR(I)-0.9d0)/0.05d0))
      ZJPARCUT(I) = ZJPAR(I)*0.5d0*(1.d0-TANH((SR(I)-0.9d0)/0.02d0))
      END DO

      OPEN(UNIT=15, FILE='profilesreal.TXT')
      DO I=2,NR-1
      WRITE(15,7) SR(I), ZJPAR(I), P0TMP(I), ZJBT(I)
      END DO
      close(UNIT=15)

      OPEN(UNIT=15, FILE='currents.TXT')
      DO I=2,NR-1
      WRITE(15,7) SR(I), ZJPARCUT(I), ZJBTCUT(I),
     >            ZJPARCUT(I)+ZJBTCUT(I)
      END DO
      close(UNIT=15)

      OPEN(UNIT=15, FILE='curr_psi.TXT')
      DO I=2,NR-1
      WRITE(15,7) SR(I), ZJPARCUT(I)+ZJBTCUT(I)
      END DO
      close(UNIT=15)

    7 FORMAT(1P12E14.6)



      IF (NPL1.NE.0) THEN
      CALL LBLTOP('Real World Data',15)
      CALL LPLOT(4,4,1,SR,P0TMP,NR-1,1,'Pressure',8,'S',1,'P [Pa]',6)
      CALL LPLOT(5,4,1,SR(2),TE(2),NR-2,1,'Temperatures',12,'S',1,
     >           'T [keV]',7)
      CALL LINCOL(1)
      CALL LPLOT(5,4,1,SR(2),TI(2),-NR+2,1,' ',1,'S',1,' ',1)
      CALL LINCOL(0)
      CALL LPLOT(6,4,1,SR(2),ZNI(2),NR-2,1,'Density',7,'S',1,
     >           'n [m^-3]',8)
      CALL LPLOT(4,5,2,SR(2),ZNUE(2),NR-2,1,'Collisionality',14,
     >           'S',1,' ',8)
      CALL LINCOl(1)
      CALL LPLOT(4,5,2,SR(2),ZNUI(2),-NR+2,1,'Collisionality',14,
     >           'S',1,' ',8)
      CALL LINCOL(0)
      CALL LPLOT(5,5,2,SR(2),SIGSPITZ(2),NR-2,1,'Conductivity',12,
     >           'S',1,' ',8)
      CALL LINCOl(1)
      CALL LPLOT(5,5,2,SR(2),SIGNEO(2),-NR+2,1,'Conductivity',12,
     >           'S',1,' ',8)
      CALL LINCOL(0)
      CALL LPLOT(4,6,1,SR(2),ZJBT(2),NR-2,1,'Bootstrap current',19,
     >           'S',1,'J [A/m^2]',9)
      CALL LINCOL(1)
      CALL LPLOT(4,6,1,SR(2),HJBT(2),-NR+2,1,' ',1,'S',1,' ',1)
      CALL LINCOL(0)
c------------------------------------- fujitsu does not support maxval()
      ZJMAX = -1.e20
      DO I=2,NR-1
        ZJMAX = MAX(ZJMAX,ZJPAR(I))
      ENDDO
      CALL NFRAME(5,6,1,0.,1.,0.,ZJMAX,
     >            'Current Profile',12,'S',1,'J.B [TA/m^2]',11)
      CALL LPLOT(5,6,1,SR(2),ZJPAR(2),-NR+2,1,'Current profile',12,
     >           'S',1,'J [A/m^2]',9)
      CALL LINCOL(1)
      CALL LPLOT(5,6,1,SR(2),ZJBT(2),-NR+2,1,' ',1,'S',1,' ',1)
      CALL LINCOL(0)
      ENDIF

C----------------------- Stability Condition Hegna (Phys. Plas. 1999)
      WRITE(20,*)
      WRITE(20,*) '*****************************************'
      WRITE(20,*) '*  NeoClassical Stability Condition     *'
      WRITE(20,*) '*****************************************'
      IP = 0
      DO I=2,NR-1
C----------------------- Undo scaling, i.e. back to dimensionless
        P0TMP(I) = P0TMP(I) * ZMU0/B0**2
        DPTMP(I) = DPTMP(I) * ZMU0/B0**2
        FTMP(I)  = FTMP(I)  / RBSCALE
        DFTMP(I) = DFTMP(I) / RBSCALE

        ZJBT(I)  = ZJBT(I) * ZMU0 * R0 / B0**2
        DQDPSI(I) = DQPROF(I) / (2.*SR(I)* CPSURF/(R0**2 * B0) )
	IF (DQDPSI(I) .NE. 0. ) THEN
            DNC(I) =  ZJBT(I) /DQDPSI(I) * GEONC(I)
	ELSE
	    DNC(I)  =  0.
	ENDIF
	IF (DIMERC(I).GT.0.) THEN
	  IP = IP +1
	  ALFAS  = 0.5 + SQRT(DIMERC(I))
          ALFAL  = 0.5 - SQRT(DIMERC(I))
          EF = -DRMERC(I) - HH(I)**2
	  STAB(IP)   = DNC(I) + (-DRMERC(I))/(ALFAS - HH(I))
	  SP(IP) = SR(I)
	  WRITE(20,101) SR(I),QPROF(I),ZJBT(I),DIMERC(I),
     >                  DNC(I),STAB(IP)
	ELSE
	  WRITE(20,102) SR(I),QPROF(I),ZJBT(I),DIMERC(I),DNC(I)
	ENDIF
      ENDDO
  102 FORMAT(' S =',F7.3,' q =',F7.3,' J.B=',F7.3,' D_I =',1PE10.3,
     >       ' D_nc =',E10.3,' Mercier Unstable')
  101 FORMAT(' S =',F7.3,' q =',F7.3,' J.B=',F7.3,' D_I =',1PE10.3,
     >       ' D_nc =',E10.3,' D_nc+D_R/(a-H) =',E10.3)
C
C----------------------- write data for analysis of Hegna's paper to file
C
      WRITE(41,*)  ' R0, B0, CPSURF, EPS :'
      WRITE(41,41) R0,B0,CPSURF,EPS
      WRITE(41,*)  ' I, S, P0, DP0, F, DF, Q, DQDPSI, J.B, D_nc : '
      DO I=2,NR-1
        WRITE(41,42) I,SR(I),P0TMP(I),DPTMP(I),FTMP(I),DFTMP(I),         WRITE(ID,392)  I,DPRI,I,DF2I,I,ZJZI,I,QPRF(I)
     >               QPROF(I),DQDPSI(I),ZJBT(I),DNC(I)
      ENDDO
   41 FORMAT(1P12E12.4)
   42 FORMAT(I4,1P12E12.4)
   43 FORMAT(12F8.4)
      NQP = 3
      QP(1) = 1.5
      QP(2) = 2.0
      QP(3) = 3.0
      WRITE(20,*)
      WRITE(20,*) '     SP,     Q,     SHEAR,      D_I,      D_R,'
     >         //'    E+F,     H,      D_NC,       DSP,      ZJBT'
      DO I=1,NQP
        DO J=2,NR-1
           IF ((QP(I).GT.QPROF(J)).AND.(QP(I).LT.QPROF(J+1))) THEN
             DLT = (QP(I)-QPROF(J))/(QPROF(J+1)-QPROF(J))
             SSP = SR(J)     + DLT*(SR(J+1)-SR(J))
             DQP = DQPROF(J) + DLT*(DQPROF(J+1)-DQPROF(J))
             DNP = DNC(J)    + DLT*(DNC(J+1)-DNC(J))
             DIP = DIMERC(J) + DLT*(DIMERC(J+1)-DIMERC(J))
             DRP = DRMERC(J) + DLT*(DRMERC(J+1)-DRMERC(J))
             HP  = HH(J)     + DLT*(HH(J+1)-HH(J))
             ZJP = ZJBT(J)   + DLT*(ZJBT(J+1)-ZJBT(J))
             EFP = -DRP - HP**2
             EF2 = -DIP - HP + 0.25
             SHRP = DQP/QP(I)  *SSP
             DSP = 0.
             IF (DIP.GE.0.) DSP = 0.5-SQRT(DIP) - HP

             WRITE(20,43) SSP,QP(I),SHRP,-DIP,-DRP,EFP,EF2,
     >                    HP,DNP,DSP,ZJP
           ENDIF
         ENDDO
       ENDDO


c      CALL LPLOT6(1,1,QPROF(20:NR-1),GEONC(20:NR-1),NR-20,'GEONC')
c      CALL LPLOT6(1,1,QPROF(20:NR-1),ZJBT(20:NR-1),NR-20,'ZJBT')
c      CALL LPLOT6(1,1,QPROF(20:NR-1),DQDPSI(20:NR-1),NR-20,'DQDPSI')
c      CALL LPLOT6(1,1,QPROF(30:NR-1),DNC(30:NR-1),NR-30,'DNC')

      CALL LPLOT(6,6,1,SR(5),DNC(5),NR-6,1,'D_nc, D_R/(a_s-H)',17,
     >           'S',1,' ',1)
      CALL LINCOL(1)
      CALL LPLOT(6,6,1,SP,STAB,-IP,1,'D_nc, D_R/(a_s-H)',17,
     >           'S',1,' ',1)
      CALL LINCOL(0)
      RETURN
      END
************************************************************************
*DECK ERRSOL
      SUBROUTINE ERRSOL(XX,YY,PSI,NR,NP)
C-----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE THE ERROR IN PSI AND GRAD(PSI) ON THE NODES
C INPUT IS A RESTORED PSI VECTOR
C-----------------------------------------------------------------------
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL ABSPSI(10,11),ABSPSX(10,11),ABSPSY(10,11),ABPSXY(10,11)
      SMALL = 1E-8
      DO 5 K=1,10
        DO 6 L=1,10
          ABSPSI(K,L) = 0.
          ABSPSX(K,L) = 0.
          ABSPSY(K,L) = 0.
          ABPSXY(K,L) = 0.
    6   CONTINUE
    5 CONTINUE
      ERRMPS = 0.
      ERRMPX = 0.
      ERRMPY = 0.
      ERRPXY = 0.
      ERRPS = 0.
      ERRPSX = 0.
      ERRPSY = 0.
      DO 60 L=1,10
        S = -1. + 2.*(L-1)/10.
      DO 10 I=1,NR-1
        DO 20 J=1,NP-1
          DO 30 K=1,10
            R = -1. +  2.*(K-1)/10.
            NODE = (I-1)*NP + J
            N1 = NODE
            N2 = NODE + 1
            N3 = NODE + 1 + NP
            N4 = NODE + NP
            CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >                  R,S,X,XR,XS,XRS,XRR,XSS)
            CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >                  R,S,Y,YR,YS,YRS,YRR,YSS)
            CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >                  PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >                  R,S,EQPSI,PSIR,PSIS,PSRS,PSRR,PSSS)
            CALL SOLOV(X,Y,SOLPSI,SPSIX,SPSIY,SPSIXY)
            ABSPSI(K,L) = ABSPSI(K,L) + ABS(SOLPSI - EQPSI)
            EJAC = XR*YS - XS*YR
            PSIX =  ( YS * PSIR - YR * PSIS) / EJAC
            PSIY =  (-XS * PSIR + XR * PSIS) / EJAC
            ABSPSX(K,L) = ABSPSX(K,L) + ABS(PSIX-SPSIX)
            ABSPSY(K,L) = ABSPSY(K,L) + ABS(PSIY-SPSIY)
            IF (ABS(SOLPSI-EQPSI).GT.ERRMPS) THEN
              ERRMPS = ABS(SOLPSI-EQPSI)
            ENDIF
            IF (ABS(PSIX-SPSIX).GT.ERRMPX) THEN
              ERRMPX = ABS(PSIX-SPSIX)
            ENDIF
            IF (ABS(PSIY-SPSIY).GT.ERRMPY) THEN
              ERRMPY = ABS(PSIY-SPSIY)
              RMAX = R
              SMAX = S
              IMAX = I
              JMAX = J
            ENDIF
            RY = - XS / EJAC
            RX =   YS / EJAC
            SY =   XR / EJAC
            SX = - YR / EJAC
            ER = XRR*YS + XR*YRS - XRS*YR - XS*YRR
            ES = XRS*YS + XR*YSS - XSS*YR - XS*YRS
            RXY = - (ER*RY + ES*SY) * YS / EJAC**2
     >            + (YRS*RY + YSS*SY) / EJAC
            SXY =   (ER*RY + ES*SY) * YR / EJAC**2
     >            - (YRR*RY + YRS*SY) / EJAC
            PSIXY = PSRR*RX*RY + PSRS*(RX*SY+SX*RY) + PSIR * RXY
     >            + PSSS*SX*SY + PSIS*SXY
            IF (ABS(PSIXY-SPSIXY).GT.ERRPXY) THEN
              ERRPXY = ABS(PSIXY-SPSIXY)
            ENDIF
            ABPSXY(K,L) = ABPSXY(K,L) + ABS(PSIXY-SPSIXY)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      DO 40 K=1,10
        R = -1. + 2.*(K-1)/10.
        ABSPSI(K,L) = ABSPSI(K,L) /((NR-1)*(NP-1))
        ABSPSX(K,L) = ABSPSX(K,L) /((NR-1)*(NP-1))
        ABSPSY(K,L) = ABSPSY(K,L) /((NR-1)*(NP-1))
        ERRPS = ERRPS + ABSPSI(K,L)
        ERRPSX = ERRPSX + ABSPSX(K,L)
        ERRPSY = ERRPSY + ABSPSY(K,L)
   40 CONTINUE
   60 CONTINUE
      WRITE(20,*) 'MAX ERR IN PSI :',ERRMPS
      WRITE(20,*) 'MAX ERR IN PSIX:',ERRMPX
      WRITE(20,*) 'MAX ERR IN PSIY:',ERRMPY
      WRITE(20,*) ' AT I,J : ',I,J
      WRITE(20,*) '    R,S : ',RMAX,SMAX
      WRITE(20,*) 'MAX ERR IN PSIXY :',ERRPXY
      WRITE(20,*) 'AVERAGE ERROR IN PSI : ',ERRPS/100.
      WRITE(20,*) 'AVERAGE ERROR IN PSIX : ',ERRPSX/100.
      WRITE(20,*) 'AVERAGE ERROR IN PSIY : ',ERRPSY/100.
      RETURN
      END

************************************************************************
*DECK SOLOV
      SUBROUTINE SOLOV(X,Y,PSI,PSIX,PSIY,PSIXY)
C-----------------------------------------------------------------------
C SOLOVIEV EQUILIBRIUM
C-----------------------------------------------------------------------
      USE COMDAT

      Y = -Y
      PAR4 = QUAD
      PAR3 = TRIA
      PAR2 = EPS
      PAR1 = ELLIP
      PSI = (1.0 - 0.25*PAR2**2) * (1.+PAR2*X)**2 * (Y/PAR1)**2
     >    + (X - 0.5*PAR2*(1.0 - X*X))**2
     >    + (4.- PAR2**2) * (X - 0.5*PAR2*(1.0-X**2)) *PAR4*Y/(2*PAR2)

      PSIX = 2*(X - 0.5*PAR2*(1-X*X)) * (1.+PAR2*X) +
     >       (2*PAR2*(1.-0.25*PAR2**2) * (1.+PAR2*X) ) * (Y/PAR1)**2
     >     +(4.- PAR2**2) * (1.+ 0.5*PAR2*(2*X))      *PAR4*Y/(2*PAR2)

      PSIY = 2*Y / PAR1**2 * (1.-0.25*PAR2**2)*(1.+ PAR2*X)**2
     >     + (4.- PAR2**2) * (X - 0.5*PAR2*(1.0-X**2))*PAR4/(2*PAR2)

      PSIY = - PSIY

      PSIXY= 2*Y / PAR1**2 * (1.-0.25*PAR2**2)*(1.+ PAR2*X)*2*PAR2
     >     + (4.- PAR2**2) * (1.+ 0.5*PAR2*(2*X))     *PAR4/(2*PAR2)

      PSIXY = - PSIXY

      RETURN
      END

************************************************************************
*DECK INFO
C***********************************************************************
C***********************************************************************
C**                                                                   **
C**  BELOW FOLLOW THE SOURCES OF THE HGOLIB ROUTINES AS USED IN       **
C**  HELENA :                                                         **
C**                     - GRID2NV      - ZERO                         **
C**                     - RFT2         - RTRAN2                       **
C**                     - RFT          - FFT2                         **
C**                     - PRARR1       - FSUM2                        **
C**                     - PRARR2                                      **
C**                                                                   **
C***********************************************************************
C***********************************************************************

************************************************************************
*DECK PRARR1
      SUBROUTINE PRARR1(NAME,ARRAY,ISIZE,LL)
C
C     ******************************************************************
C     * PRINTS 1D ARRAY(ISIZE) WITH THE TITLE NAME.                    *
C     * IF ISIZE.LT.0, INDEX I STARTS COUNTING FROM 0 (RATHER THAN 1). *
C     * THE ARGUMENT LL HAS THE DOUBLE FUNCTION OF COMMUNICATING THE   *
C     * OUTPUT UNIT IOUT=LL/10 AND THE FORMAT SWITCH L=LL-IOUT*10.     *
C     * IOUT=0 (L=LL): FILE IS "OUTPUT".                               *
C     * L=1: E-FORMAT, L=2: F-FORMAT (WIDTH OF 132 CHARACTERS),        *
C     * L=3: E-FORMAT, L=4: F-FORMAT (WIDTH OF 80 CHARACTERS).         *
C     ******************************************************************
C
      CHARACTER*(*) NAME
      DIMENSION ARRAY(*)
C
      IOUT=LL/10
      L=LL-IOUT*10
      IF(L.EQ.0) THEN
         RETURN
      ELSEIF(L.EQ.1) THEN
         KSTEP=8
         ASSIGN 11 TO IFM
         IF(ISIZE.GE.100) ASSIGN 111 TO IFM
      ELSEIF(L.EQ.2) THEN
         KSTEP=8
         ASSIGN 12 TO IFM
      ELSEIF(L.EQ.3) THEN
         KSTEP=4
         ASSIGN 13 TO IFM
         IF(ISIZE.GE.100) ASSIGN 113 TO IFM
      ELSEIF(L.EQ.4) THEN
         KSTEP=4
         ASSIGN 14 TO IFM
      ENDIF
C
      I0=0
      IS=ISIZE
      IF(ISIZE.LT.0) THEN
         I0=1
         IS=-ISIZE
      ENDIF
C
      IF(IOUT.EQ.0) THEN
         WRITE(   *,10) NAME
         DO 20 K=1,IS,KSTEP
         KPLUS=MIN(K+KSTEP-1,IS)
   20    WRITE(   *,IFM) (ARRAY(I),I-I0,I=K,KPLUS)
      ELSE
         WRITE(IOUT,10) NAME
         DO 200 K=1,IS,KSTEP
         KPLUS=MIN(K+KSTEP-1,IS)
  200    WRITE(IOUT,IFM) (ARRAY(I),I-I0,I=K,KPLUS)
      ENDIF
      RETURN
C
   10 FORMAT(/1X,A/)
   11 FORMAT(1X,8(1PE12.4,'(',0P,I2,')'))
  111 FORMAT(1X,8(1PE11.4,'(',0P,I3,')'))
   12 FORMAT(1X,8(F11.5,'(',I3,')'))
   13 FORMAT(1X,5(1PE12.4,'(',0P,I2,')'))
  113 FORMAT(1X,5(1PE11.4,'(',0P,I3,')'))
   14 FORMAT(1X,5(F11.5,'(',I3,')'))
      END
C
************************************************************************
*DECK PRARR2
      SUBROUTINE PRARR2(NAME,INDS,ARRAY,ISIZE,JSIZE,IMAX,LL)
C
C     ******************************************************************
C     * PRINTS 2D ARRAY(ISIZE,JSIZE) WITH TITLE AND INDICES INDICATED  *
C     * BY THE CHARACTER VARIABLES NAME AND INDS.                      *
C     * ISIZE SHOULD NOT EXCEED THE DIMENSION IMAX DECLARED IN THE     *
C     * CALLING PROGRAM.                                               *
C     * IF ISIZE.LT.0, INDEX I STARTS COUNTING FROM 0 (RATHER THAN 1). *
C     * IF JSIZE.LT.0, INDEX J STARTS COUNTING FROM 0 (RATHER THAN 1). *
C     * THE ARGUMENT LL HAS THE DOUBLE FUNCTION OF COMMUNICATING THE   *
C     * OUTPUT UNIT IOUT=LL/10 AND THE FORMAT SWITCH L=LL-IOUT*10.     *
C     * IOUT=0 (L=LL): FILE IS "OUTPUT".                               *
C     * L=1: E-FORMAT, L=2: F-FORMAT (WIDTH OF 132 CHARACTERS),        *
C     * L=3: E-FORMAT, L=4: F-FORMAT (WIDTH OF 80 CHARACTERS).         *
C     ******************************************************************
C
      CHARACTER*(*) NAME
      CHARACTER*(*) INDS
      CHARACTER*10 INDS1
      DIMENSION ARRAY(IMAX,*)
C
      INDS1=INDS
      IOUT=LL/10
      L=LL-IOUT*10
      IF(L.EQ.0) THEN
         RETURN
      ELSEIF(L.EQ.1) THEN
         KSTEP=8
         ASSIGN 11 TO IFM1
         ASSIGN 21 TO IFM2
      ELSEIF(L.EQ.2) THEN
         KSTEP=12
         ASSIGN 12 TO IFM1
         ASSIGN 22 TO IFM2
      ELSEIF(L.EQ.3) THEN
         KSTEP=5
         ASSIGN 13 TO IFM1
         ASSIGN 23 TO IFM2
      ELSEIF(L.EQ.4) THEN
         KSTEP=8
         ASSIGN 14 TO IFM1
         ASSIGN 24 TO IFM2
      ENDIF
C
      I0=0
      IS=ISIZE
      IF(ISIZE.LT.0) THEN
         I0=1
         IS=-ISIZE
      ENDIF
      J0=0
      JS=JSIZE
      IF(JSIZE.LT.0) THEN
         J0=1
         JS=-JSIZE
      ENDIF
C
      IF(IOUT.EQ.0) THEN
         WRITE(   *,10) NAME
         DO 20 K=1,JS,KSTEP
         KPLUS=MIN(K+KSTEP-1,JS)
         WRITE(   *,IFM1) INDS1,(J-J0,J=K,KPLUS)
         WRITE(   *,*)
         DO 20 I=1,IS
   20    WRITE(   *,IFM2) I-I0,(ARRAY(I,J),J=K,KPLUS)
      ELSE
         WRITE(IOUT,10) NAME
         DO 200 K=1,JS,KSTEP
         KPLUS=MIN(K+KSTEP-1,JS)
         WRITE(IOUT,IFM1) INDS1,(J-J0,J=K,KPLUS)
         WRITE(IOUT,*)
         DO 200 I=1,IS
  200    WRITE(IOUT,IFM2) I-I0,(ARRAY(I,J),J=K,KPLUS)
      ENDIF
      RETURN
C
   10 FORMAT(/1X,A)
   11 FORMAT(/2X,A10,I7,7I13)
   12 FORMAT(/2X,A10,I3,11I9)
   13 FORMAT(/2X,A10,I7,4I13)
   14 FORMAT(/2X,A10,I3,7I9)
   21 FORMAT(1X,I3,2X,1P,8E13.4)
   22 FORMAT(1X,I3,2X,12F9.5)
   23 FORMAT(1X,I3,2X,1P,5E13.4)
   24 FORMAT(1X,I3,2X,8F9.5)
      END

************************************************************************
*DECK ZERO
      SUBROUTINE ZERO(X1,Y1,X2,Y2,FUNC,ERR,X,Y,IZERO,LL)
C
C     ******************************************************************
C     * THE ZERO Y=FUNC(X)=0 ON THE INTERVAL (X1,X2) IS FOUND.         *
C     * THE FUNCTION FUNC(X) SHOULD BE PROVIDED BY AN EXTERNAL.        *
C     * UPON RETURN IZERO IS THE NUMBER OF ITERATIONS WHICH WERE       *
C     * REQUIRED TO OBTAIN ABS(Y).LE.ERR.                              *
C     * DIAGNOSTIC INFORMATION IS PRINTED IF L.NE.0.                   *
C     * THE ARGUMENT LL HAS THE DOUBLE FUNCTION OF COMMUNICATING THE   *
C     * OUTPUT UNIT IOUT=LL/10 AND THE PRINT SWITCH L=LL-IOUT*10.      *
C     * IOUT=0 (L=LL): FILE IS "OUTPUT".                               *
C     * MODIFIED BY JAN REM TO IMPROVE CONVERGENCE 25/08/84.           *
C     ******************************************************************
C
      EXTERNAL FUNC
C
      IOUT=LL/10
      L=LL-IOUT*10
      NIZERO=50
      IF((Y1.LT.0..AND.Y2.LT.0.).OR.(Y1.GT.0..AND.Y2.GT.0.)) THEN
         IF(IOUT.EQ.0) WRITE(   *,5)
         IF(IOUT.NE.0) WRITE(IOUT,5)
         RETURN
      ENDIF
      IF(L.NE.0) THEN
         IF(IOUT.EQ.0) WRITE(   *,6)
         IF(IOUT.NE.0) WRITE(IOUT,6)
      ENDIF
      X1S=X1
      Y1S=Y1
      X2S=X2
      Y2S=Y2
      IF(X1.GT.X2) THEN
         X1=X2S
         Y1=Y2S
         X2=X1S
         Y2=Y1S
      ENDIF
      SIG=1.
      IF(Y1.GE.0.) THEN
         SIG=-1.
         Y1=-Y1
         Y2=-Y2
      ENDIF
C
C     ***BEGIN LOOP ON IZERO***
      IZERO=0
   10 X0=X1-(X2-X1)*Y1/(Y2-Y1)
      Y0=SIG*FUNC(X0)
      IZERO=IZERO+1
      IF(L.NE.0) THEN
         IF(IOUT.EQ.0) WRITE(   *,11) IZERO,X1,Y1,X2,Y2,X0,Y0
         IF(IOUT.NE.0) WRITE(IOUT,11) IZERO,X1,Y1,X2,Y2,X0,Y0
      ENDIF
      IF(ABS(Y0).LE.ERR) GOTO 20
      IF(ABS(Y0).GE.0.2*MIN(-Y1,Y2)) THEN
         A=((Y2-Y0)/(X2-X0)-(Y0-Y1)/(X0-X1))/(X2-X1)
         B=(Y2-Y1)/(X2-X1)-A*(X2+X1)
         C=Y0-A*X0*X0-B*X0
         XN=(-B+SQRT(B*B-4.*A*C))/(2.*A)
         YN=SIG*FUNC(XN)
         IZERO=IZERO+1
         IF(L.NE.0) THEN
            IF(IOUT.EQ.0) WRITE(   *,12) IZERO,X1,Y1,X2,Y2,XN,YN
            IF(IOUT.NE.0) WRITE(IOUT,12) IZERO,X1,Y1,X2,Y2,XN,YN
         ENDIF
         IF(ABS(YN).LE.ERR) GOTO 30
         IF(YN.LT.0.) THEN
            X1=XN
            Y1=YN
            IF(Y0.GT.0.) THEN
               X2=X0
               Y2=Y0
            ENDIF
         ELSE
            X2=XN
            Y2=YN
            IF(Y0.LT.0.) THEN
               X1=X0
               Y1=Y0
            ENDIF
         ENDIF
      ELSEIF(Y0.GT.ERR) THEN
         X2=X0
         Y2=Y0
      ELSE
         X1=X0
         Y1=Y0
      ENDIF
      IF(IZERO.LT.NIZERO) GOTO 10
C     ***END LOOP ON IZERO***
      IF(IOUT.EQ.0) WRITE(   *,13) NIZERO
      IF(IOUT.NE.0) WRITE(IOUT,13) NIZERO
C
   20 X=X0
      Y=SIG*Y0
      X1=X1S
      Y1=Y1S
      X2=X2S
      Y2=Y2S
      RETURN
C
   30 X=XN
      Y=SIG*YN
      X1=X1S
      Y1=Y1S
      X2=X2S
      Y2=Y2S
      RETURN
C
C     * FORMATS.
    5 FORMAT(/1X,'***SUBROUTINE ZERO: Y1 AND Y2 VIOLATE REQUIREMENTS')
    6 FORMAT(/30X,'==== SUBROUTINE ZERO ====')
   11 FORMAT(1X,'IZERO=',I3,3X,'X1=',1PE12.4,3X,'Y1=',1PE12.4,
     A       3X,'X2=',1PE12.4,3X,'Y2=',1PE12.4,3X,'X0=',1PE12.4,
     B       3X,'Y0=',1PE12.4)
   12 FORMAT(1X,'IZERO=',I3,3X,'X1=',1PE12.4,3X,'Y1=',1PE12.4,
     A       3X,'X2=',1PE12.4,3X,'Y2=',1PE12.4,3X,'XN=',1PE12.4,
     B       3X,'YN=',1PE12.4)
   13 FORMAT(/1X,'***SUBROUTINE ZERO: NO CONVERGENCE FOR X IN ',I3,
     A       ' STEPS')
      END
C
************************************************************************
*DECK RFT2
      SUBROUTINE RFT2(DATA,NR,KR)
C
C     ******************************************************************
C     * REAL FOURIER TRANSFORM.                                        *
C     * INPUT:  NR REAL COEFFICIENTS                                   *
C     *             DATA(1),DATA(1+KR),....,DATA(1+(NR-1)*KR).         *
C     * OUTPUT: NR/2+1 COMPLEX COEFFICIENTS                            *
C     *            (DATA(1),      DATA(1+KR))                          *
C     *            (DATA(1+2*KR), DATA(1+3*KR))                        *
C     *             .............................                      *
C     *            (DATA(1+NR*KR),DATA(1+(NR+1)*KR).                   *
C     * THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST *
C     * (NR+1)*KR+1 ELEMENTS. (I.E., NR+2 IF INCREMENT KR=1).          *
C     * LASL ROUTINE MAY 75, CALLING FFT2 AND RTRAN2.                  *
C     ******************************************************************
C
      DIMENSION DATA(*)
      CALL FFT2(DATA(1),DATA(KR+1),NR/2,-(KR+KR))
      CALL RTRAN2(DATA,NR,KR,1)
      RETURN
      END
C
C
************************************************************************
*DECK RTRAN2
      SUBROUTINE RTRAN2(DATA,NR,KR,KTRAN)
C
C     ******************************************************************
C     * INTERFACE BETWEEN RFT2, RFI2, AND FFT2.                        *
C     * THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST *
C     * (NR+1)*KR+1 ELEMENTS.                                          *
C     * LASL ROUTINE MAY 75, CALLED FROM RFT2 AND RFI2.                *
C     ******************************************************************
C
      DIMENSION DATA(*)
      KS=2*KR
      N=NR/2
      NMAX=N*KS+2
      KMAX=NMAX/2
      THETA=1.5707963267949/N
      DC=2.*SIN(THETA)**2
      DS=SIN(2.*THETA)
      WS=0.
      IF(KTRAN.LE.0) THEN
         WC=-1.0
         DS=-DS
      ELSE
         WC=1.0
         DATA(NMAX-1)=DATA(1)
         DATA(NMAX-1+KR)=DATA(KR+1)
      ENDIF
      DO 10 K=1,KMAX,KS
         NK=NMAX-K
         SUMR=.5*(DATA(K)+DATA(NK))
         DIFR=.5*(DATA(K)-DATA(NK))
         SUMI=.5*(DATA(K+KR)+DATA(NK+KR))
         DIFI=.5*(DATA(K+KR)-DATA(NK+KR))
         TR=WC*SUMI-WS*DIFR
         TI=WS*SUMI+WC*DIFR
         DATA(K)=SUMR+TR
         DATA(K+KR)=DIFI-TI
         DATA(NK)=SUMR-TR
         DATA(NK+KR)=-DIFI-TI
         WCA=WC-DC*WC-DS*WS
         WS=WS+DS*WC-DC*WS
         WC=WCA
   10 CONTINUE
      RETURN
      END
C
************************************************************************
*DECK FFT2
      SUBROUTINE FFT2 (DATAR,DATAI,N,INC)
C
C     ******************************************************************
C     * FFT2 FORTRAN VERSION CLAIR NIELSON MAY 75.                     *
C     ******************************************************************
C
      DIMENSION DATAR(*), DATAI(*)
      KTRAN=ISIGN(-1,INC)
      KS=IABS(INC)
      IP0=KS
      IP3=IP0*N
      IREV=1
      DO 20 I=1,IP3,IP0
         IF(I.LT.IREV) THEN
            TEMPR=DATAR(I)
            TEMPI=DATAI(I)
            DATAR(I)=DATAR(IREV)
            DATAI(I)=DATAI(IREV)
            DATAR(IREV)=TEMPR
            DATAI(IREV)=TEMPI
         ENDIF
         IBIT=IP3/2
   10    IF(IREV.GT.IBIT) THEN
            IREV=IREV-IBIT
            IBIT=IBIT/2
            IF(IBIT.GE.IP0) GOTO 10
         ENDIF
   20    IREV=IREV+IBIT
      IP1=IP0
      THETA=REAL(KTRAN)*3.1415926535898
   30 IF(IP1.GE.IP3) GOTO 60
      IP2=IP1+IP1
      SINTH=SIN(.5*THETA)
      WSTPR=-2.*SINTH*SINTH
      WSTPI=SIN(THETA)
      WR=1.
      WI=0.
      DO 50 I1=1,IP1,IP0
         DO 40 I3=I1,IP3,IP2
            J0=I3
            J1=J0+IP1
            TEMPR=WR*DATAR(J1)-WI*DATAI(J1)
            TEMPI=WR*DATAI(J1)+WI*DATAR(J1)
            DATAR(J1)=DATAR(J0)-TEMPR
            DATAI(J1)=DATAI(J0)-TEMPI
            DATAR(J0)=DATAR(J0)+TEMPR
   40       DATAI(J0)=DATAI(J0)+TEMPI
         TEMPR=WR
         WR=WR*WSTPR-WI*WSTPI+WR
   50    WI=WI*WSTPR+TEMPR*WSTPI+WI
      IP1=IP2
      THETA=.5*THETA
      GOTO 30
   60 RETURN
      END

************************************************************************
*DECK GRID2NV
      SUBROUTINE GRID2NV(TIN,TOUT,JPTS,ACC,IGRD,LL)
C------------------------------------------------------------------------
C  THE FUNCTION TIN(TOUT), GIVEN ON THE GRID TOUT=2*PI*(J-1)/JPTS,
C  IS INVERTED TO GIVE TOUT(TIN) ON THE GRID TIN=2*PI*(I-1)/JPTS.
C  THIS IS DONE BY DETERMINING THE ZEROS OF THE FUNCTION
C     Y(T)=T+SUM(GF(M)*SIN(M*T))-2*PI*(I-1)/JPTS,
C  WHERE GF(M) ARE THE FOURIER COEFFICIENTS OF G(T)=TIN(T)-T.
C  DIAGNOSTIC INFORMATION IS PRINTED IF L.NE.0.
C  THE ARGUMENT LL HAS THE DOUBLE FUNCTION OF COMMUNICATING THE
C  OUTPUT UNIT IOUT=LL/10 AND THE PRINT SWITCH L=LL-IOUT*10.
C  IOUT=0 (L=LL): FILE IS "OUTPUT".
C-----------------------------------------------------------------------
      PARAMETER (JMAX=1024,NINV=100)
      DIMENSION TIN(*),TOUT(*),T(JMAX+1),G(JMAX+1),GFCOS(JMAX/2-1),
     A          GFSIN(JMAX/2-1)
      EQUIVALENCE(T(1),G(1))
C
      PI=3.1415926535898
      MHARM=JPTS/2-1
C
      DO 9 JJ=2,JPTS
        IF (TIN(JJ-1).GT.TIN(JJ))  TIN(JJ)=TIN(JJ)+2*PI
    9 CONTINUE
      IOUT=LL/10
      L=LL-IOUT*10
      IF(L.NE.0) THEN
         IF(IOUT.EQ.0) WRITE(   *,3)
         IF(IOUT.NE.0) WRITE(IOUT,3)
      ENDIF
      DO 10 I=1,JPTS
   10 G(I)=TIN(I)-2.*PI*(I-1.)/JPTS
      CALL RFT(G,GFNUL,GFCOS,GFSIN,JPTS,MHARM)
      IF(L.NE.0) THEN
         CALL PRARR1('TIN(J) : ',TIN,JPTS,IOUT*10+L)
         CALL PRARR1('G(I):',T,JPTS,IOUT*10+L)
         IF(IOUT.EQ.0) WRITE(20,56) GFNUL
         IF(IOUT.NE.0) WRITE(IOUT,56) GFNUL
         CALL PRARR1('GFCOS(M):',GFCOS,MHARM,IOUT*10+L)
         CALL PRARR1('GFSIN(M):',GFSIN,MHARM,IOUT*10+L)
      ENDIF
      DO 20 I=1,JPTS+1
        T(I)=2.*PI*(I-1.)/JPTS
   20 CONTINUE
      J1=1
      IGRDNV=1
      IFIRST=0
      ICIRC= - (INT(TIN(1)/(2*PI)+10000) -  9999)
      IF (ABS(TIN(1)).LT.1E-12) ICIRC=0
      DO 80 I=1,JPTS
        J=J1
        T1=T(J) + ICIRC*2*PI
        CALL FSUM2(SUM1,T1,GFNUL,GFCOS,GFSIN,MHARM)
        Y1=T1+SUM1-T(I)
   30   CONTINUE
          T0=T1
          Y0=Y1
          IF (ABS(Y0).LE.ACC) THEN
            TOUT(I)=T0
            GOTO 80
          ENDIF
          IF (J.NE.JPTS+1) GOTO 31
            IF (IFIRST.EQ.0)  THEN
              J=1
              ICIRC=ICIRC+1
              IFIRST=1
            ELSE
              WRITE(IOUT,70)
              GOTO 90
            ENDIF
   31     J=J+1
          T1=T(J) + ICIRC*2*PI
          CALL FSUM2(SUM1,T1,GFNUL,GFCOS,GFSIN,MHARM)
          Y1=T1+SUM1-T(I)
          IF(SIGN(1.,Y0).EQ.SIGN(1.,Y1)) GOTO 30
        J1=J-1
        DO 40 N=1,NINV
          T2=T0-(T1-T0)*Y0/(Y1-Y0)
          CALL FSUM2(SUM2,T2,GFNUL,GFCOS,GFSIN,MHARM)
          Y2=T2+SUM2-T(I)
          IF(L.NE.0) THEN
             IF(IOUT.EQ.0) WRITE(   *,25) N,T0,T1,Y0,Y1,T2,Y2,J
             IF(IOUT.NE.0) WRITE(IOUT,25) N,T0,T1,Y0,Y1,T2,Y2,J
          ENDIF
          IF(ABS(Y2).LE.ACC) GOTO 50
          IF(SIGN(1.,Y2).EQ.SIGN(1.,Y1)) GOTO 45
          T0=T2
          Y0=Y2
          GOTO 40
   45     T1=T2
          Y1=Y2
   40   CONTINUE
   50   TOUT(I)=T2
        IF(L.NE.0) THEN
          IF(IOUT.EQ.0) WRITE(   *,55) I,N
          IF(IOUT.NE.0) WRITE(IOUT,55) I,N
        ENDIF
        IF(N.GT.IGRDNV) IGRDNV=N
   80 CONTINUE
   90 RETURN
C
C     * FORMATS.
    3 FORMAT(///1X,'SUBROUTINE GRIDINV')
   25 FORMAT(1X,'N=',I3,' T0=',F10.5,' T1=',F10.5,
     A       ' Y0=',F10.5,' Y1=',F10.5,' T2=',F10.5,' Y2=',F10.5,
     B       ' J=',I3)
   55 FORMAT(1X,'I=',I3,5X,'N=',I3)
   56 FORMAT(/1X,'GFNUL = ',1PE12.4)
   70 FORMAT(/1X,'***SUBROUTINE GRIDINV: NO ZERO FOUND ')
      END

************************************************************************
*DECK FSUM2
      SUBROUTINE FSUM2(F,T,FFNUL,FFCOS,FFSIN,MHARM)
C-----------------------------------------------------------------------
C FOURIER SYNTHESIS OF GENERAL  FUNCTION F(T) AT SINGLE POINT T.
C-----------------------------------------------------------------------
      DIMENSION FFCOS(*),FFSIN(*)
      CO=COS(T)
      SI=SIN(T)
      C=1.
      S=0.
      SUM=.5*FFNUL
      DO 10 M=1,MHARM
        CA=C*CO-S*SI
        S=S*CO+C*SI
        C=CA
        SUM=SUM+FFCOS(M)*C + FFSIN(M)*S
   10 CONTINUE
      F=SUM
      RETURN
      END

************************************************************************
*DECK RFT
      SUBROUTINE RFT(F,FFNUL,FFCOS,FFSIN,JPTS,MHARM)
C-----------------------------------------------------------------------
C  CALCULATES FOURIER COSINE AND SINE COEFFICIENTS FFCOS AND
C  FFSIN OF THE ARRAY FF CORRESPONDING TO THE  FUNCTION
C  F(T)=.5*FFNUL+SUM(FFCOS(M)*COS(M*T)+FFSIN(M)*SIN(M*T))
C  WHERE MHARM.LE.JPTS/2-1, FFNUL=FF(0) AND T=2*PI*(J-1)/JPTS.
C  THE INPUT ARRAY F(J) IS NOT DESTROYED BY CALLING RFTCOS.
C  TYPICAL USE IS FOR MHARM MUCH SMALLER THAN JPTS/2-1, SO THAT
C  RFT2 CANNOT BE USED DIRECTLY.
C-----------------------------------------------------------------------
      PARAMETER (JMAX=1024)
      DIMENSION F(*),FFCOS(*),FFSIN(*),FSTORE(JMAX+2)
      DO 10 J=1,JPTS
   10 FSTORE(J)=F(J)
      CALL RFT2(FSTORE,JPTS,1)
      FAC=2./JPTS
      FFNUL=FSTORE(1)*FAC
      DO 20 M=1,MHARM
        FFCOS(M)=FSTORE(2*M+1)*FAC
        FFSIN(M) = - FSTORE(2*M+2)*FAC
   20 CONTINUE
      RETURN
      END

************************************************************************
*DECK RFI2
      SUBROUTINE RFI2(DATA,NR,KR)
C-----------------------------------------------------------------------
C  INVERSE OF RFT2.
C  WHEN USING RFI2 IT IS NECESSARY TO HAVE VANISHING IMAGINARY
C  PARTS OF THE FIRST AND LAST ELEMENT OF THE INPUT VECTOR:
C  DATA(1+KR)=DATA(1+(NR+1)*KR)=0.
C  THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST
C    (NR+1)*KR+1 ELEMENTS.
C  LASL ROUTINE MAY 75, CALLING RTRAN2 AND FFT2.
C -----------------------------------------------------------------------
C
      DIMENSION DATA(*)
      CALL RTRAN2(DATA,NR,KR,-1)
      MR=NR*KR
      FNI=2./NR
      DO 10 I=1,MR,KR
   10 DATA(I)=FNI*DATA(I)
      CALL FFT2(DATA(1),DATA(KR+1),NR/2,(KR+KR))
      RETURN
      END

**************************************************************************
*DECK MNEWTAX
      SUBROUTINE mnewtax(ntrial,x,n,tolx,tolf,errx,errf)
C-------------------------------------------------------------------------
C ROUTINE TO SOLVE TWO NONLINEAR EQUATIONS USING NEWTONS METHOD FROM
C NUMERICAL RECIPES.
C LU DECOMPOSITION REPLACED BY EXPLICIT SOLUTION OF 2X2 MATRIX.
C-------------------------------------------------------------------------
      USE PARAM
      USE CORNERS
      USE FAXIS
      PARAMETER (NEQ=2)
      REAL X(NEQ),FVEC(NEQ),FJAC(NEQ,NEQ)
      INTEGER n,ntrial
      REAL tolf,tolx
      INTEGER i,k,indx(NEQ)
      REAL d,errf,errx,p(NEQ)

      do 14  k=1,ntrial
c----------------------------- usrfun iserted here -----------------------
        R = X(1)
        S = X(2)
        N1 = NODENO(NAXIS,1)
        N2 = NODENO(NAXIS,2)
        N3 = NODENO(NAXIS,3)
        N4 = NODENO(NAXIS,4)
        CALL INTERP(PSI(4*(N1-1)+1),PSI(4*(N2-1)+1),
     >            PSI(4*(N3-1)+1),PSI(4*(N4-1)+1),
     >            R,S,ZPSI,ZPSIR,ZPSIS,ZPSIRS,ZPSIRR,ZPSISS)
        FVEC(1) = ZPSIR
        FVEC(2) = ZPSIS
        FJAC(1,1) = ZPSIRR
        FJAC(1,2) = ZPSIRS
        FJAC(2,1) = ZPSIRS
        FJAC(2,2) = ZPSISS
c-------------------------------------------------------------------------
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        temp = p(1)
        dis = fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)
        p(1) = (fjac(2,2)*p(1)-fjac(1,2)*p(2))/dis
        p(2) = (fjac(1,1)*p(2)-fjac(2,1)*temp)/dis
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          x(i)=x(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *N*1V45_Lt+V'.


************************************************************************
*DECK SPLINE
      SUBROUTINE SPLINE(N,X,Y,ALFA,BETA,TYP,A,B,C,D)
C-----------------------------------------------------------------------
C     INPUT:
C
C     N     ANZAHL DER KNOTEN
C     X     ARRAY DER X-WERTE
C     Y     ARRAY DER Y-WERTE
C     ALFA  RANDBEDINGUNG IN X(1)
C     BETA        "       IN X(N)
C     TYP   =  0  NOT-A-KNOT SPLINE
C              1  ALFA, BETA 1. ABLEITUNGEN VORGEGEBEN
C              2    "    "   2.     "           "
C              3    "    "   3.     "           "
C
C     BEMERKUNG: MIT TYP = 2 UND ALFA = BETA = 0 ERHAELT MAN
C           EINEN NATUERLICHEN SPLINE
C
C     OUTPUT:
C
C     A, B, C, D     ARRAYS DER SPLINEKOEFFIZIENTEN
C       S = A(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2+ D(I)*(X-X(I))**3
C
C     BEI ANWENDUNGSFEHLERN WIRD DAS PROGRAMM MIT ENTSPRECHENDER
C     FEHLERMELDUNG ABGEBROCHEN
C-----------------------------------------------------------------------
C
      INTEGER  N, TYP
      REAL     X(N), Y(N), ALFA, BETA, A(N), B(N), C(N), D(N)
      INTEGER  I, IERR
      REAL     H(1001)

C
      IF((TYP.LT.0).OR.(TYP.GT.3)) THEN
         WRITE(*,*) ' ERROR IN ROUTINE SPLINE: WRONG TYPE'
         STOP
      ENDIF
C
      IF (N.LT.3) THEN
         WRITE(*,*) ' ERROR IN ROUTINE SPLINE: N < 3 '
         STOP
      ENDIF
C
C
C     BERECHNE DIFFERENZ AUFEINENDERFOLGENDER X-WERTE UND
C     UNTERSUCHE MONOTONIE
C
      DO I = 1, N-1
        H(I) = X(I+1)- X(I)
        IF ( H(I).LE. 0.0) THEN
         WRITE(*,*) ' NON-MONOTONIC COORDINATE IN SPLINE: X(I-1)>=X(I)'
         STOP
       ENDIF
      ENDDO
C
C     AUFSTELLEN DES GLEICHUNGSSYSTEMS
C
      DO I = 1, N-2
         A(I) = 3.0 * ((Y(I+2)-Y(I+1)) / H(I+1) - (Y(I+1)-Y(I)) / H(I))
         B(I) = H(I)
         C(I) = H(I+1)
         D(I) = 2.0 * (H(I) + H(I+1))
      ENDDO
C
C     BERUECKSICHTIGEN DER RANDBEDINGUNGEN
C
C     NOT-A-KNOT
C
      IF (TYP.EQ.0) THEN
         A(1)   = A(1) * H(2) / (H(1) + H(2))
         A(N-2) = A(N-2) * H(N-2) / (H(N-1) + H(N-2))
         D(1)   = D(1) - H(1)
         D(N-2) = D(N-2) - H(N-1)
         C(1)   = C(1) - H(1)
         B(N-2) = B(N-2) - H(N-1)
      ENDIF
C
C     1. ABLEITUNG VORGEGEBEN
C
      IF (TYP.EQ.1) THEN
         A(1)   = A(1) - 1.5 * ((Y(2)-Y(1)) / H(1) - ALFA)
         A(N-2) = A(N-2) - 1.5 * (BETA - (Y(N)-Y(N-1)) / H(N-1))
         D(1)   = D(1) - 0.5 * H(1)
         D(N-2) = D(N-2) - 0.5 * H(N-1)
      ENDIF
C
C     2. ABLEITUNG VORGEGEBEN
C
      IF (TYP.EQ.2) THEN
         A(1)   = A(1) - 0.5 * ALFA * H(1)
         A(N-2) = A(N-2) - 0.5 * BETA * H(N-1)
      ENDIF
C
C     3. ABLEITUNG VORGEGEBEN
C
      IF (TYP.EQ.3 ) THEN
         A(1)   = A(1) + 0.5 * ALFA * H(1) * H(1)
         A(N-2) = A(N-2) - 0.5 * BETA * H(N-1)* H(N-1)
         D(1)   = D(1) + H(1)
         D(N-2) = D(N-2) + H(N-1)
      ENDIF
C
C     BERECHNUNG DER KOEFFIZIENTEN
C
      CALL SGTSL(N-2,B,D,C,A,IERR)

      IF(IERR.NE.0) THEN
         WRITE(*,21)
         STOP
      ENDIF
C
C     UEBERSCHREIBEN DES LOESUNGSVEKTORS
C
      CALL DCOPY(N-2,A,1,C(2),1)
C
C     IN ABHAENGIGKEIT VON DEN RANDBEDINGUNGEN WIRD DER 1. UND
C     DER LETZTE WERT VON C KORRIGIERT
C
      IF (TYP.EQ.0) THEN
         C(1) = C(2) + H(1) * (C(2)-C(3)) / H(2)
         C(N) = C(N-1) + H(N-1) * (C(N-1)-C(N-2)) / H(N-2)
      ENDIF
C
      IF (TYP.EQ.1) THEN
         C(1) = 1.5*((Y(2)-Y(1)) / H(1) - ALFA) / H(1) - 0.5 * C(2)
         C(N) = -1.5*((Y(N)-Y(N-1)) / H(N-1)-BETA) / H(N-1)-0.5*C(N-1)
      ENDIF
C
      IF (TYP.EQ.2) THEN
         C(1) = 0.5 * ALFA
         C(N) = 0.5 * BETA
      ENDIF
C
      IF (TYP.EQ.3) THEN
         C(1) = C(2) - 0.5 * ALFA * H(1)
         C(N) = C(N-1) + 0.5 * BETA * H(N-1)
      ENDIF
C
      CALL DCOPY(N,Y,1,A,1)
C
      DO I = 1, N-1
         B(I) = (A(I+1)-A(I)) / H(I) - H(I) * (C(I+1)+2.0 * C(I)) / 3.0
         D(I) = (C(I+1)-C(I)) / (3.0 * H(I))
      ENDDO
C
      B(N) = (3.0 * D(N-1) * H(N-1) + 2.0 * C(N-1)) * H(N-1) + B(N-1)
C
      RETURN
C
   21 FORMAT(1X,'ERROR IN SGTSL: MATRIX SINGULAR')
      END

************************************************************************
*DECK SPWERT
      REAL FUNCTION SPWERT(N,XWERT,A,B,C,D,X,ABLTG)
C-----------------------------------------------------------------------
C     INPUT:
C
C     N           ANZAHL DER KNOTENPUNKTE
C     XWERT       STELLE AN DER FUNKTIONSWERTE BERECHNET WERDEN
C     A, B, C, D  ARRAYS DER SPLINEKOEFFIZIENTEN (AUS SPLINE)
C     X           ARRAY DER KNOTENPUNKTE
C
C     OUTPUT:
C
C     SPWERT   FUNKTIONSWERT AN DER STELLE XWERT
C     ABLTG(I) WERT DER I-TEN ABLEITUNG BEI XWERT
C-----------------------------------------------------------------------
C
      INTEGER  N
      REAL     XWERT, A(N), B(N), C(N), D(N), X(N), ABLTG(3)
      INTEGER  I, K, M
C
C     SUCHE PASSENDES INTERVALL (BINAERE SUCHE)
C
      I = 1
      K = N
C
   10 M = (I+K) / 2
C
      IF(M.NE.I) THEN
         IF(XWERT.LT.X(M)) THEN
            K = M
         ELSE
            I = M
         ENDIF
         GOTO 10
      ENDIF
C
      XX = XWERT - X(I)
C
      ABLTG(1) = (3.0 * D(I) * XX + 2.0 * C(I)) * XX + B(I)
      ABLTG(2) = 6.0 * D(I) * XX + 2.0 * C(I)
      ABLTG(3) = 6.0 * D(I)
C
      SPWERT = ((D(I)*XX + C(I))*XX + B(I))*XX + A(I)
C
      RETURN
      END

************************************************************************
*DECK HELBAL
      SUBROUTINE HELBAL(ZVOL,ZVOLP,XAXIS)
C-----------------------------------------------------------------------
C PROGRAM TO DETERMINE THE BALLOONING STABILITY OF HELENA EQUILIBRIA
C         - READS THE MAPPING FILE AS USED BY CASTOR
C         - CALCULATES STABILITY INDEX USING SUYDAM METHOD
C         - FOR SYMMETRIC AND ASYMMETRIC PLASMA BOUNDARY SHAPES
C
C         VERSION : 1                 DATE : 28-09-95
C-----------------------------------------------------------------------
C STATUS :
C
C 28/9/95   - tested for symmetric soloviev equilibrium with E=1.41
C             eps=0.382, compares well with Turnbull results.
C           - also tested for sym. soloviev asymmetric helena
C
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMPIO
      USE COMDAT
      USE COMMAP
      USE COMSPL
      USE COMB02
      USE COMPQ
      USE COMNAM
      REAL ZVOL(*),ZVOLP(*),XAXIS
      CHARACTER*25 BAL,INIBAL

      PI = 2.*ASIN(1.)
      CALL INIT(IAS)
C---------------------- READ HELENA MAPPING FILE -----------------------
      CALL IODSK
C---------------------- CALCULATE B-FIELD ON GRID POINTS ---------------
      CALL BFIELD(IAS)
C---------------------- CALCULATE P AND Q COEFF. ON ALL GRIDPOINTS -----
      CALL PQ(IAS)
C---------------------- CALCULATE STABILITY INDEX ----------------------
      WRITE(20,*)
      WRITE(20,*) '****************************************************'
     >            //'***************************'
      WRITE(20,*) '* I, FLUX,  RHO,   Q,    SHEAR,   SHEAR1, ALPHA,'//
     >              '  ALPHA1,  FMARG,  BALLOONING *'
      WRITE(20,*) '****************************************************'
     >            //'***************************'
      DO 10 IPSI=2,NPSI
         FACT = 1.
         CALL SUYDAM(IPSI,0.,TBB,TBF,NCPQ,1.,BAL)
C---------------------- FIND DISTANCE FROM STABILITY BOUNDARY ----------
         INIBAL = BAL
         IF (BAL.EQ.' STABLE') THEN
            DELF = 2.
         ELSE
            DELF = 0.5
         ENDIF
         FACT = DELF * FACT
         DO 20 NIT=1,10
           CALL SUYDAM(IPSI,0.,TBB,TBF,NCPQ,FACT,BAL)
           IF (BAL.EQ.INIBAL) THEN
              FACT = DELF * FACT
           ELSE
              GOTO 30
           ENDIF
   20   CONTINUE
C--------------------- UPPER AND LOWER LIMII ESTABLISHED --------------
   30   CONTINUE
        IF (INIBAL.EQ.' STABLE') THEN
          FUNST = FACT
          FSTAB = FACT / 2.
        ELSE
          FSTAB = FACT
          FUNST = FACT * 2.
        ENDIF
c------------------------------------- BISECTION TO FIND FACTOR
        DO 40 NIT = 1, 10
           FACT = (FUNST + FSTAB)/2.
           CALL SUYDAM(IPSI,0.,TBB,TBF,NCPQ,FACT,BAL)
           IF (BAL.EQ.' STABLE') THEN
             FSTAB = FACT
           ELSE
             FUNST = FACT
           ENDIF
   40   CONTINUE
        SHEAR = CS(IPSI)/QS(IPSI) * DQS(IPSI)
        ALPHA = - 2.*QS(IPSI)**2 * P2(IPSI)/EPS

        SHEAR1 = 2.*ZVOL(IPSI)/QS(IPSI) * DQS(IPSI)
     >         / (2.*CS(IPSI)*ZVOLP(IPSI))
c-------------------------------------------- Lao's definition
c        ALPHA1 = -P2(IPSI)/(2.*CS(IPSI)) * ZVOLP(IPSI) / CPSURF**2
c     >         * SQRT(ZVOL(IPSI)/(4.*PI*(1.+EPS*XAXIS))) * EPS**3
c
c--------------------------------------------- use rho as radius
        RHO = SQRT(ZVOL(IPSI)/ZVOL(NR))
        DRHODS =  CS(IPSI)/ (RHO*ZVOL(NR)) * ZVOLP(IPSI)
        ALPHA2 = ALPHA / DRHODS

c----------------------------------------- Lao corrected? Needs check
c        ALPHA11 = -P2(IPSI)/(2.*CS(IPSI)) * ZVOLP(IPSI) / CPSURF**2
c     >          / (PI*PI*EPS) * RHO * ZVOL(NR) *EPS**4
c----------------------------------------------------------------------

	FMARG = (FSTAB+FUNST)/2.
C----------------------------------- temporary fix for negative shear and IAS=1
        IF ((SHEAR<0.).AND.(IAS.EQ.1)) THEN
	  FMARG = 0.
	  INIBAL =  ' STABLE'
        ENDIF
        WRITE(20,11) IPSI,CS(IPSI)**2,RHO,QS(IPSI),SHEAR,SHEAR1,
     >               ALPHA,ALPHA2,FMARG, INIBAL
   10 CONTINUE
      WRITE(20,*) '****************************************************'
     >            //'*****************'
   11 FORMAT(' ',I3,' ',F8.5,' ',F8.5,' ',F8.5,' ',F8.5,' ',F8.5,' ',
     >       F8.5,' ',F8.5,' ',F8.3,'  ',A25)

      END

************************************************************************
*DECK INIT
      SUBROUTINE INIT(IAS)
C-----------------------------------------------------------------------
C SUBROUTINE TO INITIALIZE THE INPUT VARIABLES, ALSO WRITES THE HEADER
C OF THE OUTPUT FILE
C-----------------------------------------------------------------------
      USE COMNAM

      TBB = -100.
      IF (IAS.EQ.0) TBB = 0.
      TBF =  100.
      WRITE(20,*) '*******************************'
      WRITE(20,*) '*   BALLOONING STABILITY      *'
      WRITE(20,*) '*   PROGRAM HELBAL VERSION 1  *'
      WRITE(20,*) '*******************************'
      WRITE(20,*)  TBB,TBF
      RETURN
      END
************************************************************************
*DECK PQ
      SUBROUTINE PQ(IAS)
C-----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE P AND Q COEFFICENTS IN THE BALLOONING
C EQUATION IN ALL GRID POINTS.
C THE TERMS DEPENDING ON THE EXTENDED BALLOONING ANGLE ARE SEPERATED.
C THE BALLOONING EQUATION (POGUSTE AND YURCHENKO) READS:
C       DT((P0 + T * P1 + T^2 P2) DT(F)) - (Q0 + T * Q1) F = 0
C
C NOTE : THE FINAL ARRAYS CP AND CQ HAVE A DIFFERENT INDEX FROM THE
C GEMIJ ARRAYS. THE J INDEX RUNS FROM 0 TO 2PI (INCLUDING 2PI), THE
C NUMBER OF POINTS IN POLOIDAL DIRECTION IS NCHI+1 FOR IAS=1 AND
C 2*NCHI-1 FOR IAS=0.
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMPIO
      USE COMMAP
      USE COMSPL
      USE COMB02
      USE COMPQ
      USE COMNAM

      NCPQ = 2*NCHI-1
      IF (IAS.EQ.1) NCPQ = NCHI+1
      DO 10 I=2,NPSI
        SPS2 = 2.*CPSURF*CS(I)
        ZQ = QS(I)
        ZT = RBPHI(I)
        ZDP = P2(I)
        ZDQ = DQS(I)
        DO 20 J=1,NCHI
          IJ = (I-1)*NCHI + J
          IF (IAS.EQ.0) THEN
            IJ1 = (I-1)*(2*NCHI-1) + J
          ELSE
            IJ1 = (I-1)*(NCHI+1) + J
          ENDIF
C----------------------------------- LOWER INDEX GEOMETRIC COEFF. -------
c          G33 = GEM33(IJ)
c          G11 = SPS2**2 *(1. + ZQ**2/ZT**2
c     >                       * GEM33(IJ) * GEM12(IJ)**2 ) / GEM11(IJ)
c          G12 = - SPS2 * ZQ**2 / ZT**2 * GEM12(IJ) * GEM33(IJ)
c          G22 = ZQ**2 / ZT**2 * GEM11(IJ) * GEM33(IJ)
c          ZJ  = SPS2 * ZQ * GEM33(IJ) / ZT
C-----------------------------------------------------------------------
          BETA =  - GEM12(IJ)/GEM11(IJ)
          CP0(IJ1) = 1. / (GEM33(IJ)*GEM11(IJ))
     >             + ZQ**2 * GEM11(IJ) *BETA**2 /(GEM33(IJ)*B02(IJ))
          CP1(IJ1) = 2.*BETA * GEM11(IJ) * (ZDQ*ZQ)
     >                               /(SPS2*GEM33(IJ)*B02(IJ))
          CP2(IJ1) = GEM11(IJ) * ZDQ**2
     >                               /(SPS2**2 *GEM33(IJ)*B02(IJ))
          CQ0(IJ1) = -ZDP* ZQ**2 * GEM33(IJ) / (SPS2*ZT*B02(IJ))**2
     >             *( (2.*ZDP+DSB02(IJ))*B02(IJ) +
     >                SPS2*BETA*GEM11(IJ)/GEM33(IJ) * DTB02(IJ))
          CQ1(IJ1) = ZDP*ZDQ*ZQ /(SPS2 * B02(IJ))**2 * DTB02(IJ)
C-----------------------------------------------------------------------
          IF (IAS.EQ.0) THEN
             IJ2 = I*(2*NCHI-1) - J + 1
             CP0(IJ2) = CP0(IJ1)
             CP1(IJ2) =-CP1(IJ1)
             CP2(IJ2) = CP2(IJ1)
             CQ0(IJ2) = CQ0(IJ1)
             CQ1(IJ2) =-CQ1(IJ1)
          ELSEIF ((IAS.EQ.1).AND.(J.EQ.1)) THEN
             IJ2 = I*(NCHI+1)
             CP0(IJ2) = CP0(IJ1)
             CP1(IJ2) = CP1(IJ1)
             CP2(IJ2) = CP2(IJ1)
             CQ0(IJ2) = CQ0(IJ1)
             CQ1(IJ2) = CQ1(IJ1)
          ENDIF
   20   CONTINUE
   10 CONTINUE
      RETURN
      END

************************************************************************
*DECK BFIELD
      SUBROUTINE BFIELD(IAS)
C-----------------------------------------------------------------------
C SUBROUTINE TO CALCULATE :
C      - THE TOTAL MAGNETIC FIELD SQUARED
C      - THE RADIAL AND POLOIDAL DERIVATIVE OF THE TOTAL FIELD SQUARED
C    ON THE GRID POINTS OF THE STRAIGHT FIELD LINE COORDINATE SYSTEM.
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMPIO
      USE COMMAP
      USE COMSPL
      USE COMB02
      USE COMNAM
      REAL SP(NPSIMAX),S1(NPSIMAX),S2(NPSIMAX),S3(NPSIMAX),S4(NPSIMAX)
C---------------------------------------- B0 ON GRID POINTS -----------
      DO 10 I=1,NPSI
        DO 20 J=1,NCHI
           IJ = (I-1)*NCHI + J
           B02(IJ) = (GEM11(IJ) + RBPHI(I)**2)/GEM33(IJ)
   20   CONTINUE
   10 CONTINUE
C---------------------------------------- D(B02)/D(THETA) -------------
      DO 30 I=1,NPSI
        IJSTART = (I-1)*NCHI + 1
        CALL DERIV(B02(IJSTART),DTB02(IJSTART),NCHI,IAS)
   30 CONTINUE
C---------------------------------------- D(B02)/DS --------------------
      DO 40 J=1,NCHI
        DO 50 I=1,NPSI
          IJ = (I-1)*NCHI + J
          SP(I) = B02(IJ)
   50   CONTINUE
        CALL SPLINE(NPSI,CS,SP,0.0,0.0,2,S1,S2,S3,S4)
        DO 60 I=1,NPSI
          IJ = (I-1)*NCHI + J
          DSB02(IJ) = S2(I)
   60   CONTINUE
   40 CONTINUE
      RETURN
      END
************************************************************************
*DECK IODSK
      SUBROUTINE IODSK
C-----------------------------------------------------------------------
C    - READS HELENA MAPPING FILE
C    - SCALES QUANTITIES WITH Q ON AXIS
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMPIO
      USE COMMAP
      USE COMSPL
      USE COMNAM
      REAL C1(NPSIMAX),dummy(3)

      NPSI = JS0 + 1
      NG = NPSI*NCHI

      DO 30 JC=1,NCHI
        GEM11(JC) = 0.
        GEM33(JC) = RAXIS**2
   30 CONTINUE
C
      DO 40 JC=1,NCHI
         CALL DCOPY(NPSI-1,GEM12(NCHI+JC),NCHI,C1,1)
         CALL SPLINE(NPSI-1,CS(2),C1,0.0,0.0,2,Q1,Q2,Q3,Q4)
         GEM12(JC) = SPWERT(NPSI-1,0.0,Q1,Q2,Q3,Q4,CS(2),DUMMY)
   40 CONTINUE

C------------------------------------------------------------------
C     SCALE QUANTITIES WITH VALUE OF Q ON AXIS (TOTAL CURRENT)
C------------------------------------------------------------------
      SCALEQ = QS(1)/QAXIS
C
      CPSURF = CPSURF*SCALEQ
      CALL DSCAL(NPSI,SCALEQ**2,P0,1)
      CALL DSCAL(NPSI*NCHI,-SCALEQ,GEM12,1)
      CALL DSCAL(NPSI*NCHI,SCALEQ**2,GEM11,1)
C
      RBPHI02 = RBPHI(1)**2
C
      DP0     = DP0*SCALEQ**2
      DPE     = DPE*SCALEQ**2
      DRBPHI0 = DRBPHI0*SCALEQ**2
      DRBPHIE = DRBPHIE*SCALEQ/SQRT(1.+RBPHI02/RBPHI(NPSI)**2*
     >          (1./SCALEQ**2-1.))
C
      DO 50 J=1,NPSI
         WURZEL   = SQRT(1.+RBPHI02/RBPHI(J)**2*(1./SCALEQ**2-1.))
         QS(J)    = WURZEL * QS(J)
         RBPHI(J) = WURZEL * SCALEQ * RBPHI(J)
   50 CONTINUE
C
      WRITE(20,*)  NPSI,NCHI
      WRITE(20,51) SCALEQ
      WRITE(20,52) CPSURF
      WRITE(20,53) (QS(JJ),JJ=1,JS0+1)
      WRITE(20,54) (P0(JJ),JJ=1,NPSI)
      WRITE(20,55) (RBPHI(JJ),JJ=1,NPSI)
C
C     SPLINES
C
      H12 = CS(NPSI-3) - CS(NPSI-2)
      H13 = CS(NPSI-3) - CS(NPSI-1)
      H14 = CS(NPSI-3) - CS(NPSI)
      H23 = CS(NPSI-2) - CS(NPSI-1)
      H24 = CS(NPSI-2) - CS(NPSI)
      H34 = CS(NPSI-1) - CS(NPSI)

      DQ1 = QS(NPSI-3) * H24 * H34 / (H12 * H13 * H14)
     >    - QS(NPSI-2) * H14 * H34 / (H12 * H23 * H24)
     >    + QS(NPSI-1) * H14 * H24 / (H13 * H23 * H34)
     >    - QS(NPSI) * (H24*H34 + H14*H34 + H14*H24) / (H14 * H24 * H34)

      CALL SPLINE(NPSI,CS,QS,DQ0,DQ1,1,Q1,Q2,Q3,Q4)
      CALL SPLINE(NPSI,CS,P0,DP0,DPE,1,P1,P2,P3,P4)
      CALL SPLINE(NPSI,CS,RBPHI,DRBPHI0,DRBPHIE,1,RBP1,RBP2,RBP3,RBP4)
C
      CALL DCOPY(NPSI,Q2,1,DQS,1)
      RETURN
C
   51 FORMAT(/' AFTER SCALE: SCALEQ=',1P,E12.4,0P)
   52 FORMAT(/' CPSURF = ',1P,E12.4,0P)
   53 FORMAT(/' QS'/(1X,1P,5E16.8,0P))
   54 FORMAT(/' P0'/(1X,1P,5E16.8,0P))
   55 FORMAT(/' RBPHI'/(1X,1P,5E16.8,0P))
      END

************************************************************************
*DECK KGS
      SUBROUTINE KGS(T,T0,CP,CQ,IPSI)
C-----------------------------------------------------------------------
C SUBROUTINE TO EVALUATE THE P AND Q COEFFICIENT, USED IN THE SUYDAM
C ROUTINE.
C NOTE : THE ROUTINE ASUMES THAT ONLY VALUES ON GRID POINTS ARE
C REQUESTED. NO INTERPOLATION IS DONE!
C-----------------------------------------------------------------------
      USE COMMAX
      USE COMPQ
      TWOPI = 4.*ASIN(1.)
      DT = TWOPI/REAL(NCPQ-1)
      TM = MOD(MOD(T,TWOPI)+TWOPI, TWOPI)
C-------------------------- NCPQ=2*NCHI-1 (IAS=0) OR NCPQ=NCHI+1 (IAS=1)
      JM = INT((NCPQ-1) * TM/TWOPI) + 1
      IJM = (IPSI-1)*NCPQ + JM
      FRAC = (TM - REAL(JM-1)*DT)/DT
      CP0D = CP0(IJM) + (CP0(IJM+1)-CP0(IJM))*FRAC
      CP1D = CP1(IJM) + (CP1(IJM+1)-CP1(IJM))*FRAC
      CP2D = CP2(IJM) + (CP2(IJM+1)-CP2(IJM))*FRAC
      CQ0D = CQ0(IJM) + (CQ0(IJM+1)-CQ0(IJM))*FRAC
      CQ1D = CQ1(IJM) + (CQ1(IJM+1)-CQ1(IJM))*FRAC
      CP = CP0D + CP1D*(T-T0) + CP2D*(T-T0)**2
      CQ = CQ0D + CQ1D*(T-T0)
      RETURN
      END

************************************************************************
*DECK SUYDAM
      SUBROUTINE SUYDAM(IPSI,T0,TBB,TBF,NCPQ,FACT,BAL)
C-----------------------------------------------------------------------
C  SUBROUTINE TO VERIFY THE SIGN OF THE BALLOONING ENERGY
C  FUNCTIONAL ACCORDING TO THE SUYDAM FINITE DIFFERENCE SCHEME.
C-----------------------------------------------------------------------
      REAL KP1,KP2
      CHARACTER*25 BAL

      PI  = 2*ASIN(1.)
      BAL = ' STABLE'
      DBT = 4.*ASIN(1.)/REAL(NCPQ-1)

      N  = INT((TBF - TBB) / DBT)
      ALP= 1.
      T2 = TBB - DBT/2.

      KP2 = 0.
      GP2 = 0.
C
      DO 10 I=1,N
        T1  = T2
        KP1 = KP2
        GP1 = GP2
        T2  = T1+DBT
        CALL KGS(T2,T0,KP2,GP2TEMP,IPSI)
        GP2 = GP2TEMP * FACT
        A11=(KP2+KP1)/DBT+0.25*(GP1+GP2)*DBT
        A01=-KP1/DBT+0.25*GP1*DBT
        IF (I.EQ.1) A01 = 0.
        ALP=A11-(A01*A01)/ALP
        IF(ALP.LE.0.) THEN
          TS = (T1+T2)/2.
          WRITE(BAL,11) TS
          RETURN
        ENDIF
   10 CONTINUE
   11 FORMAT(' UNSTABLE AT T = ',F8.3)
   20 RETURN
      END

**********************************************************************
*DECK DERIV
      SUBROUTINE DERIV(ARRIN,DARR,NCHI,IAS)
C---------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE THETA DERIVATE USING FFT.
C ARRIN : THE INPUT ARRAY
C DARR  : THE RESULTING THETA DERIVATIVE
C NCHI  : THE NUMBER OF POLOIDAL POINTS
C IAS   : 0 FOR UP/DOWN SYMMETRIC EQUILIBRIA, 1 FOR ASYMMETRIC EQUIL.
C---------------------------------------------------------------------
      USE COMMAX
      REAL       ARRIN(*),DARR(*),FF(2*NCHIMAX+2),DF(2*NCHIMAX+2)
      INTEGER    INDEX(2*NCHIMAX)

      PI = 2.* ASIN(1.)

      DO 10 J=1,NCHI
         INDEX(J) = J
   10 CONTINUE

      IF (IAS.EQ.0) THEN
        DO 20 J=NCHI+1,2*NCHI-2
          INDEX(J) = 2*NCHI-J
   20   CONTINUE
      ENDIF

      IF (IAS.EQ.0) THEN
         N=2*(NCHI-1)
      ELSE
         N=NCHI
      ENDIF
      DO 40 J=1,N
        FF(J) = ARRIN(INDEX(J))
   40 CONTINUE

      CALL RFT2(FF,N,1)

      DO 50 J = 1,N/2
        DF(2*J-1) = - FLOAT(J-1) * FF(2*J)
        DF(2*J)   =   FLOAT(J-1) * FF(2*J-1)
   50 CONTINUE

      DF(2) = 0.
      DF(N+2) = 0.

      CALL RFI2(DF,N,1)

      DO 60 J=1,NCHI
        DARR(J) = DF(J)
   60 CONTINUE
      END

************************************************************************
*DECK DBLELM
      SUBROUTINE DBLELM(XX,YY,PSI,NODENO,NR,NP,INDEX)
C-----------------------------------------------------------------------
C SUBROUTINE TO ADJUST EXISTING GRID TO HALF THE RADIAL SIZE OF
C RADIAL ELEMENT INDEX, I.E. BETWEEN INDEX I=INDEX,INDEX+1
C-----------------------------------------------------------------------
      USE PARAM
      REAL XX(4,*),YY(4,*),PSI(*)
      REAL XN(4,NPMAX),YN(4,NPMAX),PSN(4,NPMAX)
      INTEGER NODENO(MAXMNODE,4)

      RR =  0.
      SS = -1.
      DO J=1,NP-1
        NN = (INDEX-1)*(NP-1) + J
        N1 = NODENO(NN,1)
        N2 = NODENO(NN,2)
        N3 = NODENO(NN,3)
        N4 = NODENO(NN,4)
        CALL INTERP(XX(1,N1),XX(1,N2),XX(1,N3),XX(1,N4),
     >	            RR,SS,X,XR,XS,XRS,XRR,XSS)
        CALL INTERP(YY(1,N1),YY(1,N2),YY(1,N3),YY(1,N4),
     >	            RR,SS,Y,YR,YS,YRS,YRR,YSS)
        CALL INTERP(PSI(4*N1-3),PSI(4*N2-3),PSI(4*N3-1),
     >              PSI(4*N4-3),RR,SS,PS,PSR,PSS,PSRS,PSRR,PSSS)
        XN(1,J) = X
	XN(2,J) = XR / 1.
	XN(3,J) = XS
	XN(4,J) = XRS / 1.
        YN(1,J) = Y
	YN(2,J) = YR  / 1.
	YN(3,J) = YS
	YN(4,J) = YRS / 1.
        PSN(1,J) = PS
	PSN(2,J) = PSR / 1.
	PSN(3,J) = PSS
	PSN(4,J) = PSRS / 1.
      ENDDO
      DO K=1,4
        XN(K,NP) = XN(K,1)
	YN(K,NP) = YN(K,1)
	PSN(K,NP) = PSN(K,1)
      ENDDO
      DO I=NR,INDEX+1,-1
        DO J=1,NP
	  NODE = (I-1)*NP+J
	  DO K=1,4
	    XX(K,NODE+NP) = XX(K,NODE)
	    YY(K,NODE+NP) = YY(K,NODE)
	    PSI(4*(NODE+NP-1)+K) = PSI(4*(NODE-1)+K)
	  ENDDO
	ENDDO
      ENDDO
      DO J=1,NP
        NODE = INDEX*NP + J
	DO K=1,4
	  XX(K,NODE) = XN(K,J)
	  YY(K,NODE) = YN(K,J)
	  PSI(4*(NODE-1)+K) = PSN(K,J)
        ENDDO
      ENDDO
      NR = NR + 1
      DO J=1, NP
        NODE = INDEX*NP + J
        XX(2,NODE) = XX(2,NODE) / 1.5
        YY(2,NODE) = YY(2,NODE) / 1.5
        PSI(4*(NODE-1)+2) = PSI(4*(NODE-1)+2) / 1.5
        XX(4,NODE) = XX(4,NODE) / 1.5
        YY(4,NODE) = YY(4,NODE) / 1.5
        PSI(4*(NODE-1)+4) = PSI(4*(NODE-1)+4) / 1.5
      ENDDO
      DO J=1, NP
        NODE = (INDEX+2)*NP + J
        XX(2,NODE) = XX(2,NODE) / 1.5
        YY(2,NODE) = YY(2,NODE) / 1.5
        PSI(4*(NODE-1)+2) = PSI(4*(NODE-1)+2) / 1.5
        XX(4,NODE) = XX(4,NODE) / 1.5
        YY(4,NODE) = YY(4,NODE) / 1.5
        PSI(4*(NODE-1)+4) = PSI(4*(NODE-1)+4) / 1.5
      ENDDO
      NO = 0
      DO N=1,NR-1
        DO M=1,NP-1
          NO = NO + 1
          NODENO(NO,1) = (N-1)*NP + M
          NODENO(NO,2) = NODENO(NO,1) + 1
          NODENO(NO,3) = NODENO(NO,2) + NP
          NODENO(NO,4) = NODENO(NO,1) + NP
        ENDDO
      ENDDO


      RETURN
  111 FORMAT(6E12.4)
      END
************************************************************************
*DECK GSCHECK
      SUBROUTINE GSCHECK
C-----------------------------------------------------------------------
C    - CHECKS THE ERROR IN THE GRAD-SHAFRANOV EQUATION USING THE
C      STRAIGHT FIELD LINE SYSTEM
C-----------------------------------------------------------------------
      USE PARAM
      USE COMMAX
      USE COMPIO
      USE COMMAP
      USE COMSPL
      USE COMNAM
      REAL GP(NPSIMAX),DGT(NCHIMAX),dummy(3)
      REAL GP1(NPSIMAX), GP2(NPSIMAX), GP3(NPSIMAX),  GP4(NPSIMAX)
      REAL ERROR(NPSIMAX-1,NCHIMAX),ZC(10),MAXERR
      REAL ERROR2(NPSIMAX-1,NCHIMAX),FF(NCHI),FP(NCHI)

      NPSI = JS0 + 1
C
      DO JC=1,NCHI
         CALL DCOPY(NPSI-1,GEM12(NCHI+JC),NCHI,C1,1)
         CALL SPLINE(NPSI-1,CS(2),C1,0.0,0.0,2,Q1,Q2,Q3,Q4)
         GEM12(JC) = SPWERT(NPSI-1,0.0,Q1,Q2,Q3,Q4,CS(2),DUMMY)
      ENDDO

      CALL SPLINE(NPSI,CS,P0,DP0,DPE,1,P1,P2,P3,P4)
      CALL SPLINE(NPSI,CS,RBPHI,DRBPHI0,DRBPHIE,1,RBP1,RBP2,RBP3,RBP4)

      MAXERR = -1.e20
      DO J=1,NCHI
        DO I=1,NPSI
          IJ1 = (I-1)*NCHI + J
	  GP(I) = GEM11(IJ1) * QS(I) / RBPHI(I)
	ENDDO
        DGP0 = (GP(2)-GP(1))/(CS(2)-CS(1))
	DGP1 = (GP(NPSI)-GP(NPSI-1))/(CS(NPSI)-CS(NPSI-1))
        CALL SPLINE(NPSI,CS,GP,DGP0,DGP1,2,GP1,GP2,GP3,GP4)
        CALL SPLINE(NPSI,CS,GP,0.,0.,2,GP1,GP2,GP3,GP4)

        DO I=2,NPSI-1
	  SPS2 = 2.*CS(I) * CPSURF
	  IJ1 = (I-1)*NCHI + 1
	  IJT = (I-1)*NCHI + J
          CALL DERIV(GEM12(IJ1),DGT,NCHI,IAS)
          DGP = (GEM11(IJT+NCHI)*QS(I+1)/RBPHI(I+1) -
     >          GEM11(IJT-NCHI)*QS(I-1)/RBPHI(I-1))/
     >           (CS(I+1)-CS(I-1))
c          DGT2 = (GEM12(IJT+1)-GEM12(IJT-1))/(CHI(J+1)-CHI(J-1))
	  ERR = ABS(RBPHI(I)/QS(I)*DGP + DGT(J) * SPS2
     >          + (RBPHI(I)*RBP2(I) + GEM33(IJT)*P2(I)))
          IF (ERR.GT.MAXERR) THEN
	    MAXERR = ERR
	    IM = I
	    JM = J
	  ENDIF
C          WRITE(20,21) I,J,CS(I),CHI(J),ERR,ALOG(ERR)
	  ERROR(I-1,J) = ALOG(ERR)
	  ERROR2(I-1,J) = ERR
        ENDDO
      ENDDO
      WRITE(20,*)
      WRITE(20,23) MAXERR,IM,JM
      WRITE(20,*)
   21 FORMAT(2I3,2F8.3,3e12.4)
   22 FORMAT(2I3,F8.3,3e12.4)
   23 FORMAT(' MAX ERROR IN GS AFTER MAPPING : ',E12.4,2I5)

c      CALL CPLOTX(21,1,1,CS(2),CHI,NPSI-1,NCHI,1,1,ERROR,NPSIMAX-1,
c     >            ZC,-10,'ERROR IN GS',11,'X',1,'Y',1,1.,4,0)
c      CALL P3PLOT(21,1,CS(2),CHI,NPSI-1,NCHI,ERROR,NPSIMAX-1,
c     >            10.,40.,'ERROR IN GS',11)

c      WRITE(20,*)
c      WRITE(20,*) ' ERROR IN GS : HARMONICS'
c      WRITE(20,*)
c      DO J=1,NCHI
c        FP(J) = FLOAT(J-1)
c      ENDDO
c      DO I=2,NPSI
c        DO J=1,NCHI
c	  FF(J) = ERROR2(I-1,J)
c	ENDDO
c        CALL RFT2(FF,NCHI,1)
c	DO M=1,NCHI/2
c          FF(M) = 2. * FF(M) / REAL(NCHI)
c	  WRITE(20,22) I,M,CS(I),FF(M)
c        ENDDO
c     ENDDO

      RETURN
C
      END

*DECK SGTSL
*** FROM NETLIB, TUE AUG 28 08:28:34 EDT 1990 ***
C
      SUBROUTINE SGTSL(N,C,D,E,B,INFO)
      INTEGER N,INFO
      REAL C(1),D(1),E(1),B(1)
C
C     SGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
C     SIDE WILL FIND THE SOLUTION.
C
C     ON ENTRY
C
C        N       INTEGER
C                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C        C       REAL(N)
C                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
C                ON OUTPUT C IS DESTROYED.
C
C        D       REAL(N)
C                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C                ON OUTPUT D IS DESTROYED.
C
C        E       REAL(N)
C                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
C                ON OUTPUT E IS DESTROYED.
C
C        B       REAL(N)
C                IS THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       IS THE SOLUTION VECTOR.
C
C        INFO    INTEGER
C                = 0 NORMAL VALUE.
C                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
C                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
C                    THIS IS DETECTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     NO EXTERNALS
C     FORTRAN ABS
C
C     INTERNAL VARIABLES
C
      INTEGER K,KB,KP1,NM1,NM2
      REAL T
C     BEGIN BLOCK PERMITTING ...EXITS TO 100
C
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = 0.0E0
            E(N) = 0.0E0
C
            DO 30 K = 1, NM1
               KP1 = K + 1
C
C              FIND THE LARGEST OF THE TWO ROWS
C
               IF (ABS(C(KP1)) .LT. ABS(C(K))) GO TO 10
C
C                 INTERCHANGE ROW
C
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
C
C              ZERO ELEMENTS
C
               IF (C(K) .NE. 0.0E0) GO TO 20
                  INFO = K
C     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = 0.0E0
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (C(N) .NE. 0.0E0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
C
C           BACK SOLVE
C
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END

!      SUBROUTINE CPUTIME(TIME1)
!-----------------------------------------------------------------------
! DUMMY ROUTINE TO REPLACE THE STANDARD CALL CPU_TIME FOR MACHINES
! WHERE THIS IS NOT AVAILABLE. X05BAF() IS A NAG ROUTINE.
!-----------------------------------------------------------------------
!      REAL TIME1
!      TIME1 = X05BAF()
!      RETURN
!      END
