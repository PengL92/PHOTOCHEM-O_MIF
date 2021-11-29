       program atm_chem
 
C     This program was created by James Kasting and modified mainly by
c     Alex Pavlov (AP) and Kara Kerelove (KK). The user-friendly version 
c     of this code was created by Antigona Segura (AS) (2005).
c     James Kasting (JK), Peng Liu (PL) and Jingjun Liu (JL) developed
c     the version from Stanton et al (2008). JL deleted the chloride 
c     and made O2/CO2 as long-lived species. JK, PL, and JL added the
c     oxygen isotope module to calculate O-MIF signals in all O-atoms
c     species following the logic in Liu et al (2019).
c     Some modifications are identified with the author's initials.
c
c     The code is mostly written in f77 but is compiled in f90 and it 
c     contains some f90 features.
c     
c     The UV fluxes for stars that are not the Sun were provided by
c     Martin Cohen from new IUE (+model) data. 
c
c     Check the notes along the main program before use it. Look for 
c     the word 'NOTE:'
C
C     This program has been updated to the JPL '92 rate constant
C     recommendations. Most of rate constants have also been updated
C     according to the latest studies collected in NIST database.
C
C         THIS PROGRAM IS DESIGNED SPECIFICALLY FOR AN O2 LEVEL
C     FROM 1.E-3 TO 1 PAL.
C         THIS PROGRAM IS A ONE-DIMENSIONAL MODEL OF THE PRIMORDIAL
C     ATMOSPHERE.  THE MIXING RATIOS OF THE LONG-LIVED SPECIES
C     ARE CALCULATED FROM THE EQUATION
C
C     DF/DT  =  (1/N)*D/DZ(KN*DF/DZ) + P/N - LF
C
C     WHERE
C     F = MIXING RATIO (USOL)
C     K = EDDY DIFFUSION COEFFICIENT (EDD)
C     N = TOTAL NUMBER DENSITY (DEN)
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTONS
C     METHOD, E.G. IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.
C
C          THE LIST OF SUBROUTINES IS AS FOLLOWS:
C     (1) GRID   -  SETS UP THE ALTITUDE GRID
C     (2) RATES  -  DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE
C     (3) PHOTO  -  COMPUTES PHOTOLYSIS RATES (CALLS TWOSTR)
C     (4) DENSTY -  COMPUTES ATMOSPHERIC DENSITIES FROM HYDROSTATIC
C                   EQUILIBRIUM AND INITIALIZES ABSORBER PROFILES
C     (5) DIFCO  -  COMPUTES DK = K*N BETWEEN GRID POINTS
C     (6) OUTPUTP -  PRINTS OUT RESULTS
C     (7) DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C                  POINTS BY CALLING CHEMPL
C     (8) CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C                  FOR ONE SPECIES AT ALL GRID POINTS
C     (9) LTNING -  COMPUTES LIGHTNING PRODUCTION RATES FOR O2 AND
C                   N2 BASED ON CHAMEIDES' RESULTS
C    (10) SCALE  -  CAN BE USED TO SCALE THE EQUATIONS BY ARBITRARY
C                   FACTORS TO IMPROVE CONDITIONING OF MATRIX
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (1) TBDY   -  COMPUTES 3-BODY REACTION RATES
C     (2) E1     - EXPONENTIAL INTEGRAL OF ORDER ONE
C

C ***** REACTION LIST *****
C     1)  H2O + O(1D) = 2OH
C     2)  H2 + O(1D) = OH + H
C     3)  H2 + O = OH + H
C     4)  H2 + OH = H2O + H
C     5)  H + O3 = OH + O2
C     6)  H + O2 + M = HO2 + M
C     7)  H + HO2 = H2 + O2
C     8)  H + HO2 = H2O + O
C     9)  H + HO2 = OH + OH
C    10)  OH + O = H + O2
C    11)  OH + HO2 = H2O + O2
C    12)  OH + O3 = HO2 + O2
C    13)  HO2 + O = OH + O2
C    14)  HO2 + O3 = OH + 2O2
C    15)  HO2 + HO2 = H2O2 + O2
C    16)  H2O2 + OH = HO2 + H2O
C    17)  O + O + M = O2 + M
C    18)  O + O2 + M = O3 + M
C    19)  O + O3 = 2O2
C    20)  OH + OH = H2O + O
C    21)  O(1D) + N2 = O(3P) + N2
C    22)  O(1D) + O2 = O(3P) + O2
C    23)  O2 + HV = O(3P) + O(1D)
C    24)  O2 + HV = O(3P) + O(3P)
C    25)  H2O + HV = H + OH
C    26)  O3 + HV = O2 + O(1D)
C    27)  O3 + HV = O2 + O(3P)
C    28)  H2O2 + HV = 2OH
C    29)  CO2 + HV = CO + O(3P)
C    30)  CO + OH = CO2 + H
C    31)  CO + O + M = CO2 + M
C    32)  H + CO + M = HCO + M
C    33)  H + HCO = H2 + CO
C    34)  HCO + HCO = H2CO + CO
C    35)  OH + HCO = H2O + CO
C    36)  O + HCO = H + CO2
C    37)  O + HCO = OH + CO
C    38)  H2CO + HV = H2 + CO
C    39)  H2CO + HV = HCO + H
C    40)  HCO + HV = H + CO
C    41)  H2CO + H = H2 + HCO
C    42)  CO2 + HV = CO + O(1D)
C    43)  H + H + M = H2 + M
C    44)  HCO + O2 = HO2 + CO
C    45)  H2CO + OH = H2O + HCO
C    46)  H + OH + M = H2O + M
C    47)  OH + OH + M = H2O2 + M
C    48)  H2CO + O = HCO + OH
C    49)  H2O2 + O = OH + HO2
C    50)  HO2 + HV = OH + O
C    51)  CH4 + HV  =  1CH2 + H2
C    52)  CH3OOH + HV  =  H3CO + OH
C    53)  N2O + HV  =  N2 + O
C    54)  HNO2 + HV  = NO + OH
C    55)  HNO3 + HV  = NO2 + OH
C    56)  NO + HV  =  N + O
C    57)  NO2 + HV  =  NO + O
C    58)  CH4 + OH  =  CH3 + H2O
C    59)  CH4 + O(1D)  =  CH3 + OH
C    60)  CH4 + O(1D)  =  H2CO + H2
C    61)  1CH2 + CH4  =  2 CH3
C    62)  1CH2 + O2  =  H2CO + O
C    63)  1CH2 + N2  =  3CH2 + N2
C    64)  3CH2 + H2  =  CH3 + H
C    65)  3CH2 + CH4  =  2 CH3
C    66)  3CH2 + O2  =  H2CO + O
C    67)  CH3 + O2 + M  =  CH3O2 + M
C    68)  CH3 + OH  =  H2CO + H2
C    69)  CH3 + O  =  H2CO + H
C    70)  CH3 + O3  =  H2CO + HO2
C    71)  CH3O2 + HO2  =  CH3OOH + O2
C    72)  CH3O2 + CH3O2  =  2 H3CO + O2
C    73)  CH3O2 + NO  =  H3CO + NO2
C    74)  H3CO + O2  =  H2CO + HO2
C    75)  H3CO + O  =  H2CO + OH
C    76)  H3CO + OH  =  H2CO + H2O
C    77)  N2O + O(1D)  =  NO + NO
C    78)  N2O + O(1D)  =  N2 + O2
C    79)  N + O2  =  NO + O
C    80)  N + O3  =  NO + O2
C    81)  N + OH  =  NO + H
C    82)  N + NO  =  N2 + O
C    83)  NO + O3  =  NO2 + O2
C    84)  NO + O + M  =  NO2 + M
C    85)  NO + HO2  =  NO2 + OH
C    86)  NO + OH + M  =  HNO2 + M
C    87)  NO2 + O  =  NO + O2
C    88)  NO2 + OH + M  =  HNO3 + M
C    89)  NO2 + H  =  NO + OH
C    90)  HNO3 + OH  =  H2O + NO3
C    91)  HO2 + NO2 + M  =  HO2NO2 + M
C    92)  HO2NO2 + OH  =  NO2 + H2O + O2
C    93)  HO2NO2 + O  =  NO2 + OH + O2
C    94)  HO2NO2 + M  =  HO2 + NO2 + M
C    95)  HO2NO2 + HV  =  HO2 + NO2
C    96)  CH3OOH + OH  =  CH3O2 + H2O
C    97)  CH3O2 + OH  = H3CO + HO2
C    98)  O3 + NO2  =  O2 + NO3
C    99)  NO2 + NO3  =  NO + NO2 + O2
C   100)  O + NO3  =  O2 + NO2
C   101)  NO + NO3  =  NO2 + NO2
C   102)  OH + NO3  =  HO2 + NO2
C   103)  HO2 + NO3  =  HNO3 + O2
C   104)  NO2 + O + M  =  NO3 + M
C   105)  NO3 + hv  =  NO2 + O
C   106)  NO3 + NO2 + M  =  N2O5 + M
C   107)  N2O5 + hv  =  NO2 + NO3
C   108)  N2O5 + M  =  NO2 + NO3 + M
C   109)  N2O5 + H2O  =  2 HNO3
C************** Added sulfur chemistry ********
C   110)  SO   + HV   =     S    +     O
C   111)  SO2  + HV   =     SO   +     O
C   112)  H2S  + HV   =     HS   +     H
C   113)  SO   + O2   =     O    +     SO2
C   114)  SO   + HO2  =     SO2  +     OH
C   115)  SO   + O    =     SO2
C   116)  SO   + OH   =     SO2  +     H
C   117)  SO2  + OH   =     HSO3
C   118)  SO2  + O    =     SO3
C   119)  SO3  + H2O  =     H2SO4
C   120)  HSO3 + O2   =     HO2  +     SO3
C   121)  HSO3 + OH   =     H2O  +     SO3
C   122)  HSO3 + H    =     H2   +     SO3
C   123)  HSO3 + O    =     OH   +     SO3
C   124)  H2S  + OH   =     H2O  +     HS
C   125)  H2S  + H    =     H2   +     HS
C   126)  H2S  + O    =     OH   +     HS
C   127)  HS   + O    =     H    +     SO
C   128)  HS   + O2   =     OH   +     SO
C   129)  HS   + HO2  =     H2S  +     O2
C   130)  HS   + HS   =     H2S  +     S
C   131)  HS   + HCO  =     H2S  +     CO
C   132)  HS   + H    =     H2   +     S
C   133)  HS   + S    =     H    +     S2
C   134)  S    + O2   =     SO   +     O
C   135)  S    + OH   =     SO   +     H
c-as This reaction was deleted and subtituted by 182a).Note: 182a is from the old order
C   136)  S    + HCO  =     HS   +     CO
c   136a) SO2  + HV   =     S    +     O2
C   137)  S    + HO2  =     HS   +     O2
C   138)  S    + HO2  =     SO   +     OH
C   139)  HS   + H2CO =     H2S  +     HCO
C   140)  SO2  +     HV  =      SO21
C   141)  SO2  +     HV  =      SO23
C   142)  H2SO4 +     HV =       SO2   +   OH  +    OH
C   143)  SO3   +    HV  =      SO2  +     O
C   144)  SO21  +    M   =      SO23 +     M
C   145)  SO21  +    M   =      SO2  +     M
C   146)  SO21  +    HV  =      SO23 +     HV
C   147)  SO21  +    HV  =      SO2  +     HV
C   148)  SO21  +    O2  =      SO3  +     O
C   149)  SO21  +    SO2 =      SO3  +     SO
C   150)  SO23  +    M   =      SO2  +     M
C   151)  SO23  +    HV  =      SO2  +     HV
C   152)  SO23  +    SO2 =      SO3  +     SO
C   153)  SO    +    NO2 =      SO2  +     NO
C   154)  SO    +    O3  =      SO2  +     O2
C   155)  SO2   +    HO2 =      SO3  +     OH
C   156)  HS    +    O3  =      HSO  +     O2
C   157)  HS    +    NO2 =      HSO  +     NO
C   158)  S     +    O3  =      SO   +     O2
C   159)  SO    +    SO  =      SO2  +     S
C   160)  SO3   +    SO  =      SO2  +     SO2
C   161)  S     +    CO2 =      SO   +     CO
C   162)  SO    +    HO2 =      HSO  +     O2
C   163)  SO    +    HCO =      HSO  +     CO
C   164)  H     +    SO  =      HSO
C   165)  HSO   +    HV  =      HS   +     O
C   166)  HSO   +    OH  =      H2O  +     SO
C   167)  HSO   +    H   =      HS   +     OH
C   168)  HSO   +    H   =      H2   +     SO
C   169)  HSO   +    HS  =      H2S  +     SO
C   170)  HSO   +    O   =      OH   +     SO
C   171)  HSO   +    S   =      HS   +     SO
c-as Reactions added for N2O 
c   172)  N2 + O1D = N2O
c   173)  N2O + H = NO + NO + OH
c   174)  N2O + NO = NO2 + N2
C***********************************************************
C
C***********************************************************
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.  THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS IN FIVE 10-DIGIT COLUMNS
C     STARTING IN COLUMN 11, I.E.
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NSP1 = NSP + 1 (INCLUDES HV)
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C               IS SOLVED
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C        PHOTOLYSIS REACTIONS ARE IDENTIFIED BY THE SYMBOL HV (NOT
C     COUNTED IN EVALUATING NSP).  THREE-BODY REACTIONS ARE WRITTEN
C     IN TWO-BODY FORM, SO THE DENSITY FACTOR MUST BE INCLUDED IN
C     THE RATE CONSTANT.
C        THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICE
C
C     ISPEC(NSP2) = VECTOR CONTAINING THE HOLLERITH NAMES OF THE
C                  CHEMICAL SPECIES.  THE LAST ENTRY MUST BE HV.
C     JCHEM(5,NR) = MATRIX OF CHEMICAL REACTIONS.  THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     ILOSS(2,NSP,NMAX) = MATRIX OF LOSS PROCESSES.  ILOSS(1,I,L)
C                         HOLDS REACTION NUMBER J, ILOSS(2,I,L) HOLDS
C                         REACTANT NUMBER.
C     IPROD(NSP,NMAX) = MATRIX OF PRODUCTION PROCESSES.  IPROD(I,L)
C                       HOLDS REACTION NUMBER J.
C     NUML(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF ILOSS
C     NUMP(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF IPROD
C
c      PARAMETER(NZ=64, NQ=34, NQT=NQ+1)
c      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
c      PARAMETER(NR=220, NF=34)
c      PARAMETER(NSP=55, NSP1=NSP+1, NSP2=NSP+2, NMAX=70)

       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNEQ_LDA.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
!C-PL AS I ADDED SGFLUX(NQ),SMFLUX(NQ),VDEP(NQ) TO comBBLOK.inc, delete their statement as follow
      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
     2  ,USAVE(NQ,NZ),R(NZ),U(NQ)!SGFLUX(NQ),SMFLUX(NQ),VDEP(NQ),
      DIMENSION UOLD(28,100),TOLD(100),EDDOLD(100),DENOLD(100),
     2  SO4AEROLD(100),AERSOLOLD(100,3),WFALLOLD(100,3),RPAROLD(100,3)
      DIMENSION DPU(NZ,3),DPL(NZ,3)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)                                    
      DIMENSION TSAV(NZ),TDUM(NZ),TP1(NZ) 
      DIMENSION water(NZ),FLOW(NQT),fluxsave(108),sfxsave(10)

      CHARACTER :: STARR*3,DIRDATA*4, AA*11,DIRIO*2
C-PL ADDED NEW VARIABLES FOR ISOTOPE CALCULATION
C-JL ADD HQO AND HQONO2 
      DIMENSION DELTAH2CO17(NZ),DELTAO17(NZ),DELTAH2O17(NZ),
     2 DELTAOH17(NZ),DELTAHO2A17(NZ),DELTAH2O217(NZ),DELTAO3A17(NZ),
     3 DELTAO3B17(NZ),DELTACO17(NZ),DELTACH3OOHA17(NZ),
     4 DELTACH3OOHB17(NZ),DELTACH3O2A17(NZ),DELTACH3O2B17(NZ),
     5 DELTAN2O17(NZ),DELTANO17(NZ),DELTANO217(NZ),DELTAHNO217(NZ),
     6 DELTAHNO317(NZ),DELTAHO2NO2A17(NZ),DELTAHO2NO2B17(NZ),
     7 DELTANO317(NZ),DELTAN2O5A17(NZ),DELTAN2O5B17(NZ),DELTAO217(NZ),
     8 DELTASO17(NZ),DELTASO217(NZ),DELTAH2SO417(NZ),DELTAHSO17(NZ),
     9 DELTACO217(NZ),DELTAO317(NZ),DELTACH3OOH17(NZ),
     1 DELTACH3O217(NZ),DELTAHO2NO217(NZ),DELTAN2O517(NZ), 
     2 DELTAHO2B17(NZ),DELTAHO2NO2C17(NZ),DELTAHO217(NZ)

C-JL ADD 18O ISOTOPE 
      DIMENSION DELTAH2CO18(NZ),DELTAO18(NZ),DELTAH2O18(NZ),
     2 DELTAOH18(NZ),DELTAHO2A18(NZ),DELTAH2O218(NZ),DELTAO3A18(NZ),
     3 DELTAO3B18(NZ),DELTACO18(NZ),DELTACH3OOHA18(NZ),
     4 DELTACH3OOHB18(NZ),DELTACH3O2A18(NZ),DELTACH3O2B18(NZ),
     5 DELTAN2O18(NZ),DELTANO18(NZ),DELTANO218(NZ),DELTAHNO218(NZ),
     6 DELTAHNO318(NZ),DELTAHO2NO2A18(NZ),DELTAHO2NO2B18(NZ),
     7 DELTANO318(NZ),DELTAN2O5A18(NZ),DELTAN2O5B18(NZ),DELTAO218(NZ),
     8 DELTASO18(NZ),DELTASO218(NZ),DELTAH2SO418(NZ),DELTAHSO18(NZ),
     9 DELTACO218(NZ),DELTAO318(NZ),DELTACH3OOH18(NZ),
     1 DELTACH3O218(NZ),DELTAHO2NO218(NZ),DELTAN2O518(NZ), 
     2 DELTAHO2B18(NZ),DELTAHO2NO2C18(NZ),DELTAHO218(NZ)

      DIMENSION C2DELTAH2CO(NZ),C2DELTAO(NZ),C2DELTAH2O(NZ),
     2 C2DELTAOH(NZ),C2DELTAHO2(NZ),C2DELTAH2O2(NZ),C2DELTAO3A(NZ),
     3 C2DELTAO3B(NZ),C2DELTACO(NZ),C2DELTACH3OOHA(NZ),
     4 C2DELTACH3OOHB(NZ),C2DELTACH3O2A(NZ),C2DELTACH3O2B(NZ),
     5 C2DELTAN2O(NZ),C2DELTANO(NZ),C2DELTANO2(NZ),C2DELTAHNO2(NZ),
     6 C2DELTAHNO3(NZ),C2DELTAHO2NO2A(NZ),C2DELTAHO2NO2B(NZ),
     7 C2DELTANO3(NZ),C2DELTAN2O5A(NZ),C2DELTAN2O5B(NZ),C2DELTAO2(NZ),
     8 C2DELTASO(NZ),C2DELTASO2(NZ),C2DELTAH2SO4(NZ),
     9 C2DELTAHSO(NZ),C2DELTACO2(NZ),C2DELTAO3(NZ),C2DELTACH3OOH(NZ),
     1 C2DELTACH3O2(NZ),C2DELTAHO2NO2(NZ),C2DELTAN2O5(NZ)


      CHARACTER*5 label1
      CHARACTER*40 CHEMI(5,NRI)
      DIMENSION DELTA18O(NQI)
      DIMENSION DELTA17O(NQI),CAP17(NQI),USOLI(NQI,NZ),ZKM(NZ)
      REAL DELTAO3,DELTACH3OOH,DELTACH3O2,DELTAHO2NO2,DELTAN2O5,
     2 CAPO317, CAPCH3OOH17,CAPC3O217,CAPHO2NO217,CAPN2O517
     3 ,LIGHTNO2I

      COMMON/ISOTOP/ARI(NRI,NZ),ILOSSI(2,NSPI,NMAXI),
     2  JCHEMI(5,NRI),NUMLI(NSPI),NUMPI(NSPI),TPI(NQI),TLI(NQI),
     3  YPI(NQI,NZ),YLI(NQI,NZ),SRI(NQI),TLOSSI(NQI),PHIDEPI(NQI),
     4 ISPECI(NSPI2),DIZ(NSPI2,NZ),IPRODI(NSPI,NMAXI),PSO4AER
     5 ,LBOUNDI(NQI),FLOWI(NQI),FUPI(NQI),CONI(NQI)
      COMMON/RATESIS/RAINGCI(NQI,NZ),DDI(NQI,NZ),DKI(NQI,NZ),
     2  DUI(NQI,NZ),DLI(NQI,NZ),DHUI(NQI,NZ),DHLI(NQI,NZ),HIZ(NQI,NZ),
     3  VDEPI(NQI),RATI(NRI)

      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,LOQ,
     4  LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2

      INTEGER OTYPE 
      REAL GPPXOYI,GPPCDEI,SYMRATIO17(NZ),SYMRATIO18(NZ) 
     
c Name of the star
      INCLUDE 'INCLUDECHEM/comSTR.inc'
      INCLUDE 'INCLUDECHEM/comFLUXPHOTO.inc'
c DIRP contains DIRCOUP, this variable is defined as a character 
c in comDIRP.inc
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comEBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comGBLOK.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comPRESS1.inc'
      INCLUDE 'INCLUDECHEM/comQBLOK.inc'
      INCLUDE 'INCLUDECHEM/comRBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'

c in may 17 2019, JL move CO2 from inert to long lived 
        DATA LH2CO,LO,LH2O,LOH,LHO2,LH2O2,LO3,LH,LH2,LCH4,LCO,
     2  LCH3OOH,LCH3O2,LN2O,LNO,LNO2,LHNO2,LHNO3,LHO2NO2,LNO3,LN2O5,
     3  LO2,LH2S,LHS,LSO,LSO2,LH2SO4,LHSO,LCO2,LSO4AER,LCH21,LCH23,
     4  LO1D,LCH3,LH3CO,LHCO,LN,LS,LSO21,LSO23,LHSO3,LSO3,LS2,LN2/
     5  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     6  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     7  44/

C-AP
C
C   NO PREDISSOCIATION COEFFICIENTS (ALLEN AND FREDERICK, 1982)
C
      DATA ANO/-1.790868E+1, -1.924701E-1, -7.217717E-2, 5.648282E-2,
     2  4.569175E-2, 8.353572E-3, 3*0.,
     3  -1.654245E+1, 5.836899E-1, 3.449436E-1, 1.700653E-1,
     4  -3.324717E-2, -4.952424E-2, 1.579306E-2, 1.835462E-2,
     5  3.368125E-3/
C
      DATA BNO/7.836832E+3, -1.549880E+3, 1.148342E+2, -3.777754E+0,
     2  4.655696E-2, 1.297581E+4, -2.582981E+3, 1.927709E+2,
     3  -6.393008E+0, 7.949835E-2/
C
      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/
      DATA RNO2/60*0., .79, .83, .66, .15, 4*0./
C
C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/
C
C ***** SOLUBILITY (GIORGI AND CHAMEIDES) *****
      DATA H/1.3E+04, 1.0E-99, 1.0E-99, 1.0E+05, 3.3E+04,
c            h2co      o         h2o      OH       HO2
     2       2.0E+05, 1.03e-2, 1.0E-99, 1.0E-99, 1.0E-99,
c             h2o2     o3
     3       1.0E-03, 2.0E+05, 3.3E+04, 2.5E-02, 1.9E-03,
     4       7.0E-03, 7.0E+11, 7.0E+11, 7.0E+11, 7.0E-03,
     5       7.0E+11, 3.2E-4,   0.14,    1.E+5,   1.9E-3,   
c                     o2       
     6       1.E+4, 7E+11, 9E+3,3.72E-2/
c                                co2 
       
C-AP  I have changed the Henry constant similar to the Archean code
C-AP  only for the sulfur species  from Archean. Note that Henry(SO2) 
C-AP  in the table is higher because we use effective Henry constant
C-AP  from Archean.

c  Temperature from the US Standard Atmosphere 1976. Used when the 
c  code is not coupled to the climate model.
C  JK   Data are estimated above 64 km. I'm adjusting the temperature
C       near the tropopause downward in order to get better 
C       statospheric H2O
      DATA T/288.15, 281.65, 275.15, 268.66, 262.17,
     &       255.68, 249.19, 242.70, 236.21, 229.73,
     &       223.25, 216.77, 216.65, 216.65, 216.65,
     &       216.65, 216.65, 216.65, 216.65, 216.65,
     &       216.65, 217.58, 218.57, 219.57, 220.56, 
     &       221.55, 222.54, 223.54, 224.53, 225.52,
     &       226.51, 227.50, 228.49, 230.97, 233.74, 
     &       236.51, 239.28, 242.05, 244.82, 247.58,
     &       250.35, 253.14, 255.88, 258.64, 261.40,
     &       264.16, 266.96, 269.68, 270.65, 270.65,
     &       270.65, 270.65, 269.03, 266.27, 263.52, 
     &       260.77, 258.02, 255.27, 252.52, 249.77, 
     &       247.02, 244.27, 241.53, 230.78, 230.50,
     &       230.50, 227.00, 225.10, 222.00, 219.60,
     &       216.00, 214.30, 212.00, 210.30, 208.00,
     &       206.40, 204.00, 202.50, 200.00, 198.60,
     &       196.00, 194.70, 192.00, 190.80, 189.00,
     &       186.90, 186.90, 186.90, 186.90, 186.90,
     &       188.00, 190.00, 191.00, 192.50, 194.00,
     &       195.00, 196.50, 198.00, 200.00, 202.00/
C

      DATA EDD/64*0., 5.174E+05, 5.674E+05, 6.224E+05, 6.826E+05,
     2  7.487E+05, 8.212E+05, 9.007E+05, 9.879E+05, 28*1.E6/
 
c NOTE: Boundary conditions are defined here
C ***** UPPER BOUNDARY FLUXES *****
      DATA SMFLUX/NQ*0./
C
C ***** EFFUSION VELOCITIES *****
      DATA VEFF/NQ*0./
C
C ***** UPPER BOUNDARY CONDITIONS *****
c      DATA MBOUND/0, 1, 8*0, 1, 24*0/
C JK 6/2019 Ignore downward fluxes of CO and O
c in may 6, JL removed CL               14
      DATA MBOUND/0, 1, 8*0, 1, 3*0, 0,6*0, 0 ,6*0, 0  ,0/
C                    O       co     NO      O2     CO2  aer
C   0 = CONSTANT EFFUSION VELOCITY (VEFF)
C   1 = CONSTANT FLUX (SMFLUX)
C
C ***** LOWER BOUNDARY FLUXES *****
c NOTE: SGFLUX is redefined later in the main code when the O2 
c mixing ratio is less than 0.21 or when the star is not the Sun
      DATA SGFLUX/NQ*0./
C
C ***** DEPOSITION VELOCITIES (NQ)*****
c      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1., 3*0., 0.2, 1.,
c     2  0., 2*0.2, 6*0.5, 0., 0.5, 1., 3*0.5, 0.02, 1, 3.E-4, 
c     3  3*1./

c in May 6 2018, JL rip off all the CL species 
      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1., 3*0., 0.2, 1.,
c               h2co o  h2o  oh  ho2 h2o2 o3  h        ch3ooh ch3o2
     2  0., 2*0.2, 5*0.5, 1.e-4,0.02, 1., 3.E-4, 3*1.,0./
c      n2o                o2    h2s   hs   so         co2

C ***** LOWER BOUNDARY CONDITIONS (NQT)*****
c NOTE: Mixing ratios for present Earth for H2,CO,CH4, N2O and CH3Cl
c are defined after the model parameters.
c fixed mixing ratios for present run
c                             H2 CH4 CO     N2O 
      DATA LBOUND/2*0, 1, 5*0, 1, 1, 1, 2*0, 1,0,4*0, 2*0, 1, 
     2  3*0, 0, 2*0, 1, 0/
c            so2
c                                                
C-PL FOR THE LOW O2 RUNS
c       DATA LBOUND/2*0, 1, 5*0, 2, 2, 2, 2*0, 2, 5*0, 2*0, 1,6*0,1,0/
c                                             n2o         o2    co2

C
C   0 = CONSTANT DEPOSITION VELOCITY (VDEP)
C   1 = CONSTANT MIXING RATIO
C   2 = CONSTANT UPWARD FLUX (SGFLUX)
C
C Modify Henry constants for sensitivity run comparison to Archean atmosphere
c      H(LHO2) = 9.E3
c      H(LH2O2) = 6.2E5
c      H(LH2CO) =  4.25E4

C============== FILE SECTION =============================

      DIRDATA = 'DATA'
      DIRIO = 'IO'

c **** INPUT FILES
c
c This file contains the chemical reactions used in the code
      OPEN(unit=61, file= DIRDATA//'/primo3s.chm')
c-as Next file contains solar flux and atmospheric data
c-as it MUST be read for all the cases  
      OPEN(unit=62, file= DIRDATA//'/photos.pdat', status='old') 
      OPEN(unit=63,file= DIRDATA//'/h2so4.pdat',status='old')
      OPEN(unit=64,file= DIRDATA//'/eddy.pdat',status='old')

c Files with the input parameters
c NOTE: Check this two before runnig the program
       OPEN(unit=65,file= DIRIO//'/input_atmchem.dat')
       OPEN(unit=66,file= DIRIO//'/planet.dat')  

c-as Next file was formerly named atm_chem_coefs.dat.
c-as I has the same format as atm_composition.out   
      OPEN(unit=67, file= DIRIO//'/atm_composition.dat') 

c Next file is generated by the climate model contains altitude,
c temperature and water. Formerly called photo_input.dat.
c Only used when ICOUPLE=1
      OPEN(unit=71,file= DIRIO//'/fromClima2Photo.dat') 
      
c-as  Unit 72 is an input and output file that is shared with the climate code
c-as  when ICOUPLE= 1. It is OPEN later in the program to WRITE on it.
c Only used when ICOUPLE=1    
      OPEN(unit=72,file= DIRIO//'/mixing_ratios.dat')

      OPEN(unit=73, file= DIRDATA//'/faruvs.pdat')
      
c  NOTE: IMPORTANT files to read the far UV of all the stars (including 
c  the Sun) and the UV fluxes for stars others than the Sun.
c        74     fluxesKGF_photo.pdat
c        75     M star flux (name it as you like)
c        76     far UV flux (name depends on the star)

C **** OUTPUT FILES
c Main output
      OPEN(unit=90,file= DIRIO//'/outchem.dat')

c  This is an output file containing altitude, H2O, O3.
c  To be used as input of the climate model, formerly called Pass2SurfMP.dat
      OPEN(unit=84,file= DIRIO//'/fromPhoto2Clima.dat')

c This file contains total hydrogen mixing ratio vs. altitude (cm)
C These OUTPUT files are opened along the program 
      OPEN(unit=86,file= DIRIO//'/total_hydrogen.out')
      OPEN(unit=87,file= DIRIO//'/diffusion_coefs.out')

c In may 9 2019, Jl add the reaction rate table to the code to make life easier 
      OPEN(UNIT=15,file=DIRIO//'/int.rates.out.dat')
      OPEN(UNIT=16,file=DIRIO//'/int.iso_rates.out.dat')     
c in may 16 2019, JL add a clean printout for h2so4 to see possible disarray of input 
      open(UNIT=158, file=DIRIO//'/h2so4.out.dat') 

C These OUTPUT files are opened along the program 
c	UNIT   	NAME  
c        19     mixing_ratios.dat (main program)

c=================================================================

C********* SET MODEL PARAMETERS *****

C     ZY = SOLAR ZENITH ANGLE (IN DEGREES)
C     AGL = DIURNAL AVERAGING FACTOR FOR PHOTORATES
C     ISEASON = TELLS WHETHER P AND T VARY WITH TIME (THEY DON'T FOR
C               ISEASON < 3)
C     IZYO2 = TELLS WHETHER SOLAR ZENITH ANGLE VARIES WITH TIME (0 SAYS
C             IT DOESN'T; 1 SAYS IT DOES)
C     IO2 = 0 FOR ALLEN AND FREDERICK O2 SCHUMANN-RUNGE COEFFICIENTS
C         = 1 FOR EXPONENTIAL SUM FITS (FOR LOW-O2 ATMOSPHERES)
C     INO = 0 FOR ALLEN AND FREDERICK NO PREDISSOCIATION COEFFICIENTS
C         = 1 FOR MODIFIED CIESLIK AND NICOLET FORMULATION
C     EPSJ = AMOUNT BY WHICH TO PERTURB THINGS FOR JACOBIAN CALCULATION
C     ZTROP = TROPOPAUSE HEIGHT (ABOVE WHICH H2O BEHAVES AS A NONCONDENS
C             ABLE GAS
C     STARR - Character variable to choose a star, it can be:
c             Sun, F2V, K2V,dMV 
c             DO NOT FORGET quotation marks
c     ICOUPLE - 1 = Coupled to the climate model              
c               0 = Not coupled
C     FCO2 = CO2 mixing ratio when ICOUPLE = 0
C     FO2 = O2 mixing ratio when ICOUPLE = 0
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED
      read(65,555)
      read(65,*)AA,STARR
      read(65,*)AA,FLUXFAC
      read(65,*)AA,INIT
      read(65,*)AA,TSTOP
      read(65,*)AA,DT
      read(65,*)AA,NSTEPS
      read(65,*)AA,ZY
      read(65,*)AA,AGL
      read(65,*)AA,ISEASON
      read(65,*)AA,IZYO2
      read(65,*)AA,IO2
      read(65,*)AA,INO
      read(65,*)AA,EPSJ
      read(65,*)AA,ZTROP
      read(65,*)AA,FCO2
      read(65,*)AA,FO2   

      read(65,*)AA,GPPOXY
      read(65,*)AA,GPPCDE
      read(65,*)AA,COMBIN
      read(65,*)AA,RECOMB
  555  format(3/)
      close(65)
      ICOUPLE = 0

*** Read the flux from a star
      call readstar(FLUXFAC)

***** READ THE PLANET PARAMETER DATAFILE *****
      READ(66,502) G,FSCALE,ALB,DELZ,ZTROP,JTROP
 502  FORMAT(F5.1/,F4.2/,F5.3/,E5.1/,E5.1/,I2)
      close(66)
C
***** READ THE ATMOSPHERIC COMPOSITION
       READ(67,400) Usol,TDUM,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 400  FORMAT(1P8E12.5)
      close(67)

c      do J = 1,NZ
c      do I = 1,28 
c      USOL(I,J) = Uold(I, J)
c      end do 
c      USOL(29,J) = 3.55E-4 
c      end do 

C-PL CHANGE O2 RATIO HERE TO MAKE SURE FO2 GET THE RIGHT NUMBER 05/24/2019
       DO I=1,NZ
       USOL(LO2,I)=FO2 * USOL(LO2,I)/USOL(LO2,1)
       USOL(LCO2,I)=FCO2 * USOL(LCO2,I)/USOL(LCO2,1)
       O2(I) = USOL(LO2,I)
       CO2(I) = USOL(LCO2,I)
       END DO  

      if (ICOUPLE.eq.1) then
	DO J=1, NZ
c Read the temperature and water profiles from the climate code
	  READ(71,*) Z(J),T(J),water(J)
	END DO
      close(71)
      endif
 351  FORMAT(I3,1PE10.3, 1PE12.3, 1PE12.3)

C-KK  Surface mixing ratios to share with the climate code
c  FAR is not needed on this code but it must be transfered to 
c  the climate model.
        READ(72,*) FAR                  !Argon
	READ(72,*) FCH4			!Methane
	READ(72,*) JTROP		!Tropopause layer
        READ(72,*) O3OLD                !Former O3 column depth
      close(72)

C   Initial constant mixing ratios used for present Earth
c       if (STARR=="Sun".and.FO2.gt.0.20) then
c           READ(67,400) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR

       if(INIT.eq.1) goto 77
c      Surface mixing ratios
         if(LBOUND(9).eq.1)FH2 = 5.5E-7
         if(LBOUND(11).eq.1)FCO = 9.0E-8
         if(LBOUND(14).eq.1)FN2O = 3.0E-7
c         if(LBOUND(23).eq.1)FCH3CL = 5.0E-10

       DO i = 1, NZ
        if(LBOUND(9).eq.1) then
           USOL(LH2, i) = USOL(LH2,i)*(FH2/USOL(LH2,1))
        endif
       if(LBOUND(10).eq.1) then
           USOL(LCH4,i) = USOL(LCH4,i)*(FCH4/USOL(LCH4,1))
        endif
        if(LBOUND(11).eq.1)then
            USOL(LCO,i) = USOL(LCO,i)*(FCO/USOL(LCO,1))
        endif
        if(LBOUND(14).eq.1) then
            USOL(LN2O,i) = USOL(LN2O,i)*(FN2O/USOL(LN2O,1))
        endif
c        if(LBOUND(23).eq.1)then
c            USOL(LCH3CL,i) = USOL(LCH3CL,i)*(FCH3CL/USOL(LCH3CL,1))
c        endif
        END DO
   77   continue       
c        endif
       
C-KK	Constant flux for long-lived gases. For details see the model 
c       description on: Segura et al. (2003) Astrobiology, 3(4), 689-708.
C-KK	Only used at O2 levels lower than 1 PAL and for planets around 
C-KK    other stars
C-KK	9 - H2   10 - CH4   11 - CO  14 - N2O   23 - CH3Cl
         if(LBOUND(9).eq.2)SGFLUX(9) = -8.61E9
         if(LBOUND(10).eq.2)SGFLUX(10) = 2.02E11
c        if(LBOUND(10).eq.2)SGFLUX(10) = 1.5E11  !value for AD Leo
         if(LBOUND(11).eq.2)SGFLUX(11) = 3.46E11
         if(LBOUND(14).eq.2)SGFLUX(14) = 1.12E9
c JL SO2
         if(LBOUND(26).eq.2)SGFLUX(26) = 6.7E9
c in may 6, 2019, JL removed Cl 
c         if(LBOUND(23).eq.2)SGFLUX(23) = 3.99E8


       write(90,*)"SOLAR ZENITH ANGLE (deg) = ", ZY
	 write(90,*)"FIXED SPECIES MIXING RATIOS:"
	 write(90,*)"   O2 = ",FO2," CO2 = ",FCO2," CH4 = ",FCH4
C	 print*,"   H2 = ",FH2," CO = ",FCO," N2O = ",FN2O

c
c in May 6, 2019, JL remake species list without cl 
C   Long-lived species
       ISPEC(1) = 4HH2CO
       ISPEC(2) = 1HO
       ISPEC(3) = 3HH2O
       ISPEC(4) = 2HOH
       ISPEC(5) = 3HHO2
       ISPEC(6) = 4HH2O2
       ISPEC(7) = 2HO3
       ISPEC(8) = 1HH
       ISPEC(9) = 2HH2
       ISPEC(10) = 3HCH4
       ISPEC(11) = 2HCO
       ISPEC(12) = 6HCH3OOH
       ISPEC(13) = 5HCH3O2
       ISPEC(14) = 3HN2O
       ISPEC(15) = 2HNO
       ISPEC(16) = 3HNO2
       ISPEC(17) = 4HHNO2
       ISPEC(18) = 4HHNO3
       ISPEC(19) = 6HHO2NO2
       ISPEC(20) = 3HNO3
       ISPEC(21) = 4HN2O5
       ISPEC(22) = 2HO2      
       ISPEC(23) = 3HH2S
       ISPEC(24) = 2HHS
       ISPEC(25) = 2HSO
       ISPEC(26) = 3HSO2
       ISPEC(27) = 5HH2SO4
       ISPEC(28) = 3HHSO
       ISPEC(29) = 3HCO2
C
C   TRIDIAGONAL SOLVER
       ISPEC(30) = 6HSO4AER
C
C   Short-lived species
       ISPEC(31) = 4HCH21
       ISPEC(32) = 4HCH23
       ISPEC(33) = 3HO1D
       ISPEC(34) = 3HCH3
       ISPEC(35) = 4HH3CO
       ISPEC(36) = 3HHCO
       ISPEC(37) = 1HN
       ISPEC(38) = 1HS
       ISPEC(39) = 4HSO21
       ISPEC(40) = 4HSO23
       ISPEC(41) = 4HHSO3
       ISPEC(42) = 3HSO3
C
C   Inert species
       ISPEC(43) = 2HS2
       ISPEC(44) = 2HN2
       ISPEC(45) = 2HHV
       ISPEC(46) = 1HM
      do i=1,NSP
       NUML(i) = 0
       NUMP(i) = 0
      enddo


C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ (61,200) JCHEM
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      close(61)
 
C NOTE: This command will print the list of chemical reactions in the 
c       output file 
c      WRITE(90,201)(J,(JCHEM(M,J),M=1,5),J=1,NR)
 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
      KJAC = LDA*NEQ
c      WRITE(90,202) NQ,NZ,KJAC
 202  FORMAT(//1X,'NQ=',I2,5X,'NZ=',I3,5X,'KJAC=',I7)
C     
    
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
c      Print *,:'Reading in the reactions'
      DO 5 J=1,NR
c      print *,'J =',J
      DO 5 M=1,5
      IF(JCHEM(M,J).EQ.' ' ) GO TO 5
C-AP Change NSP1 to NSP2 because we added M to the species list
      DO 6 I=1,NSP2
      IF(JCHEM(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      GO TO 25
   5  CONTINUE


C 
C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
c      print *,'Exchanging names for numbers'
      DO 7 M=1,2
      N = 3-M
      DO 7 J=1,NR
c      print *,'J =',J

      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
      NUML(I) = NUML(I) + 1
      IF(NUML(I).GT.NMAX) GO TO 20
      K = NUML(I)
      ILOSS(1,I,K) = J
      ILOSS(2,I,K) = JCHEM(N,J)
   7  CONTINUE
   
C
      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE
 

c      print *, 'I =', I
c      print *, 'NUMP =', NUMP
c      print *, 'J=',J
c       print *, 'K=',K
C-AP
C
c** Read the eddy diffusion profile
	do J=1,NZ
	 READ(64,245) EDD(J),Z(J)
	end do 
 245  FORMAT(2(1PE11.3, 1X))
      close(64)
      LTIMES = 0       !Counter for the photorate subroutine
      ISULF = 0
      VOLFLX = 3.E9    !Volcanic flux of SO2 (I think! JK)
c      print *, 'Z before call',Z

      CALL GRID
c      print *, 'Z after call',Z
      DZ = Z(2) - Z(1)

      CALL DENSTY(O2,CO2,P0,P)

      CALL RATES

      CALL DIFCO(O2,CO2)


c  If the star is other than the Sun or there is a flare
c  the fluxes are saved here
      do i=1,108
       fluxsave(i)=FLUX(i)
      enddo     

c  Read the data for the photochemistry and the UV flux of the Sun
      CALL READPHOTO
      if(STARR.ne.'Sun') then
        do i=1,108
         FLUX(i)= fluxsave(i)
        enddo       
      endif
   
c** Water calculation

      JTROP = ZTROP/DZ + 0.01

      CALL PSATRAT(H2O)
c      print *, 'hello world' 
       DO 23 J=1,JTROP 
  23    USOL(LH2O,J) = H2O(J)
 
       
       if(ICOUPLE.eq.1) then
C-KK	This is added to make sure that tropospheric water is being
C-KK	handled consistently. The #s are imported from the climate model.
        DO J = 1, NZ
         USOL(LH2O,J) = water(J)
        END DO
       endif       !end of water calculation


      CALL LTNING(FO2,FCO2,P0)
c      print *, 'After call to LTNING'
      CALL AERTAB
      NZ1 = NZ - 1
      HA = 1.38E-16*T(NZ)/(1.67E-24*28.*980.)


      DTINV = 1./DT
      TIME = 0.
C
C ***** PRINT OUT INITIAL DATA *****
c in may 10 2019, JL common out call for test 
      CALL OUTPUTP(0,NSTEPS,0.,FLOW)	
C    
C ***** SET JACOBIAN DIMENSIONING PARAMETERS *****
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   PRINT OUT RESULTS EVERY NPR TIME STEPS
      NPR = NSTEPS
      PRN = NPR
C
C   DO PHOTORATES EVERY MP TIME STEPS
      NPHOT = 0
      MP = 3
      PM = MP
      NN = 0 
 
C
C ***** START THE TIME-STEPPING LOOP *****
C	STARSHIP
      DO 1 N=1,NSTEPS
      TIME = TIME + DT
      NN = NN + 1
      MS = (N-1)/MP
      SM = (N-1)/PM
      IF(NN.EQ.NSTEPS) SM = MS
      IF(SM-MS.GT.0.01) GO TO 18
      IF(N.GT.1 .AND. TIME.LT.1.E4) GO TO 18
c      print *, 'JTROP =', JTROP 
c      print *, 'DEN =', DEN
C
C   STORE ABSORBERS USED TO BLOCK OUT SOLAR UV RADIATION
      DO 35 I=1,NZ
      H2O(I) = ABS(USOL(LH2O,I))
      O3(I) = ABS(USOL(LO3,I))
      O2(I) = USOL(LO2,I) !FO2
      CO2(I) = USOL(LCO2,I) !FCO2
C-AP
      FSO2(I) = ABS(USOL(LSO2,I))
      H2S(I) = ABS(USOL(LH2S,I))
C-AP
      CH4(I) = ABS(USOL(LCH4,I))


c      print *, 'H2O(I) =', H2O 
c      print *, 'O3(I) = ', O3 
c      print *, 'FSO2(I) = ', FSO2
c      print *, 'H2S(I) = ', H2S 
c      print *, 'CH4(I) = ', CH4
  35  CONTINUE
       
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1    
!      CALL DENSTY(O2,CO2,P0)   		
!      CALL DIFCO(O2,CO2)
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)
      CALL AERCON(H2O)
c      print *, ' po3d(1) = ', PO3D(1) 

C
C ****** TIME-DEPENDENT BOUNDARY CONDITIONS ********
C
C  UPPER BOUNDARY
C  Escape of hydrogen: VEFF(H) = (Bi/Ha)/Nt, Nt=total number density
      BOVERH = DI(LH,NZ)*DEN(NZ)/HA
      VEFF(LH) = BOVERH/DEN(NZ)
      BOVERH2 = DI(LH2,NZ)*DEN(NZ)/HA
      VEFF(LH2) = BOVERH2/DEN(NZ)
      BOVERCH4 = DI(LCH4,NZ)*DEN(NZ)/HA
      VEFF(LCH4) = BOVERCH4/DEN(NZ)
      BOVERH2O = DI(LH2O,NZ)*DEN(NZ)/HA
      VEFF(LH2O) = BOVERH2O/DEN(NZ)
C

 
      IF(NN.EQ.NSTEPS) write(90, 63) BOVERH,VEFF(LH),BOVERH2,VEFF(LH2)
     & ,BOVERCH4,VEFF(LCH4),BOVERH2O,VEFF(LH2O),VEFF(LO2),VEFF(LCO2)
  63  FORMAT(/'Information on upper boundary conditions'/'BOVERH=',
     2  1PE10.3,' VEFF(LH)= ',E10.3,' BOVERH2=',
     3   E10.3,' VEFF(LH2)= ',E10.3,/'BOVERCH4= ',E10.3,' VEFF(LCH4)= '
     4  ,E10.3,' BOVERH2O- ',E10.3,' VEFF(LH2O)= ',E10.3,/'VEFF(LO2)= '
     5  ,E10.3,3x,'VEFF(LCO2)= ',E10.3)
C
C PL 07/2019 Calculate SMFLUX(LCO) 
      VO2 = (PO2(NZ) + PO2D(NZ)) * HA
      VCO2 = (PCO2(NZ) + PCO2D(NZ)) * HA
      VEFF(LO2) = VO2
      VEFF(LCO2) = VCO2
      SMFLUX(LCO) = - VCO2*CO2(NZ)*DEN(NZ)   
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ) - 2.*VO2*O2(NZ)*DEN(NZ)
     2  + 2.*VEFF(LCH4)*USOL(LCH4,NZ)*DEN(NZ)

C
      NMP = NSTEPS - MP
      IF (NN.gt.1 .and. nn.LT.NSTEPS) GO TO 18

      write(90,97)
  97  FORMAT(//1X,'PHOTOLYSIS RATES')
      write(90,98) 
  98  FORMAT(/5X,'Z',7X,'PO2',6X,'PO2D',5X,'PCO2',5X,'PCO2D',4X,
     2  'PH2O',5X,'PO3',6X,'PO3D',5X,'PH2O2',4X,'PHCO',5X,'PH2',
     3  6X,'PHO2')
      write(90,99)(Z(I),PO2(I),PO2D(I),PCO2(I),PCO2D(I),PH2O(I),
     2  PO3(I),PO3D(I),PH2O2(I),PHCO(I),PH2(I),PHO2(I),I=1,NZ,
     3  3)

  99  FORMAT(2X,1P12E9.2)
      write(90,198)
 198  FORMAT(/5X,'Z',6X,'PCH4',5X,'PCH3OOH',2X,'PN2O',5X,'PHNO3',4X,
     2  'PNO',6X,'PNO2',5X,'PHNO4',4X,'PCCL3F',3X,'PCCL2F2',2X,
     3  'PCCL4',4X,'PCH3CL')
      write(90,99) (Z(I),PCH4(I),PMOOH(I),PN2O(I),PHNO3(I),PNO(I),
     2  PNO2(I),PHNO4(I),PCCL3F(I),PCCL2F2(I),PCCL4(I),PCH3CL(I),
     3  I=1,NZ,3)
      write(90,197)
 197  FORMAT(/5X,'Z',6X,'PMCCL3',3X,'PCL2',5X,'PHOCL',4X,'PNOCL',4X,
     2  'PCLONO',3X,'PCLONO2',2X,'PCLO2',4X,'PHCL')
      write(90,199) (Z(I),PMCCL3(I),PCL2(I),PHOCL(I),PNOCL(I),PCLONO(I),
     2  PCLONO2(I),PCLO2(I),PHCL(I),I=1,NZ,3)
 199  FORMAT(2X,1P9E9.2)
      write(90,298)
 298  format(/5x,'Z',6x,'PNO3',5x,'PN2O5',4x,'PCL2O2',3x,'PSO2',5x,
     2  'PH2S',5x,'PSO21',4x,'PSO23',4x,'PH2SO4')
      write(90,299) (z(i),pno3(i),pn2o5(i),pcl2o2(i),pso2(i),ph2s(i),
     2  pso21(i),pso23(i),ph2so4(i),i=1,nz,3)
 299  format(2x,1p9e9.2)
  18  CONTINUE

      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL SEDMNT(FSULF,IDO)
c      print *,'After SEDMNT'
      DO J=1,NZ
        AERSOL(J,1) = SO4AER(J)*DEN(J)/CONVER(J,1)
      ENDDO                                                               
C

C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,LDA
      DO 17 K=1,NEQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NEQ
  19  RHS(K) = 0.

C
C     (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
c      print *
c      print *,'Before call to DOCHEM'

c      print *, 'USOl =', USOL
      CALL DOCHEM(FVAL,IDO)
c      print *, 'DEN = ', DEN
c      print *, 'DEN = ', DEN 
      DO 9 I=1,NQ
      DO 9 J=1,NZ
      K = I + (J-1)*NQ
      RHS(K) = FVAL(I,J)
   9  USAVE(I,J) = USOL(I,J)
C
      DO 3 I=1,NQ
      DO 11 J=1,NZ
      R(J) = EPSJ * ABS(USOL(I,J))
  11  USOL(I,J) = USAVE(I,J) + R(J)
      CALL DOCHEM(FV,0)
C
      DO 12 M=1,NQ
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQ
  12  DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)
C
      DO 10 J=1,NZ
  10  USOL(I,J) = USAVE(I,J)
   3  CONTINUE
C
C   COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO 13 I = 1,NQ
      DO 14 J=2,NZ1
      K = I + (J-1)*NQ
      RHS(K) = RHS(K) - DD(I,J)*USOL(I,J) 
     2  + (DU(I,J) + DHU(I,J))*USOL(I,J+1) 
     3  + (DL(I,J) - DHL(I,J))*USOL(I,J-1)      
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DD(I,J)
      DJAC(KU,K+NQ) = - DU(I,J) - DHU(I,J)
  14  DJAC(KL,K-NQ) = - DL(I,J) + DHL(I,J)
  13  CONTINUE
C
C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQ
      U(K) = USOL(K,1)
      LB = LBOUND(K)
C
C   CONSTANT DEPOSITION VELOCITY
      IF(LB.NE.0) GO TO 16
      RHS(K) = RHS(K) + (DU(K,1) + DHU(K,1))*USOL(K,2) - DU(K,1)*U(K)
     2  - (VDEP(K)/DZ - HI(K,1)/(2.*DZ))*U(K)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(K,1) + VDEP(K)/DZ
     2  - HI(K,1)/(2.*DZ)
      DJAC(KU,K+NQ) = - DU(K,1) - DHU(K,1)
      GO TO 15
C
C   CONSTANT MIXING RATIO
  16  IF(LB.NE.1) GO TO 31
      RHS(K) = 0.
      DO 36 M=1,NQ
      MM = KD + K - M
  36  DJAC(MM,M) = 0.
      DJAC(KU,K+NQ) = 0.
      DJAC(KD,K) = DTINV + DU(K,1)
      GO TO 15
C
C   CONSTANT UPWARD FLUX
C-PL fixed sign error in last term of DJAC from SH
  31  CONTINUE
      RHS(K) = RHS(K) + (DU(K,1) + DHU(K,1))*USOL(K,2) - DU(K,1)*U(K)
     2   + HI(K,1)*U(K)/(2.*DZ) + SGFLUX(K)/DEN(1)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(K,1) - HI(K,1)/(2.*DZ)
      DJAC(KU,K+NQ) = - DU(K,1) - DHU(K,1)
  15  CONTINUE
C
C ***** UPPER BOUNDARY CONDITIONS *****
      DO 30 I=1,NQ
      U(I) = USOL(I,NZ)
      K = I + NZ1*NQ
      MB = MBOUND(I)
C
C   CONSTANT EFFUSION VELOCITY
      IF(MB.NE.0) GO TO 29
      RHS(K) = RHS(K) + (DL(I,NZ) - DHL(I,NZ))*USOL(I,NZ1) 
     2  - DL(I,NZ)*U(I) - (VEFF(I)/DZ + HI(I,NZ)/(2.*DZ))*U(I)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(I,NZ) + VEFF(I)/DZ
     2  + HI(I,NZ)/(2.*DZ)
      DJAC(KL,K-NQ) = - DL(I,NZ) + DHL(I,NZ)
      GO TO 30
C
C   CONSTANT DOWNWARD FLUX
  29  CONTINUE
      RHS(K) = RHS(K) + (DL(I,NZ) - DHL(I,NZ))*USOL(I,NZ1)
     2  - DL(I,NZ)*U(I) - HI(I,NZ)*U(I)/(2.*DZ) 
     3  - SMFLUX(I)/DEN(NZ)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(I,NZ) + HI(I,NZ)/(2.*DZ) 
      DJAC(KL,K-NQ) = - DL(I,NZ) + DHL(I,NZ)
  30  CONTINUE
C
      DO 33 J=1,NZ
      IF(Z(J).GT.ZTROP) GO TO 34
      K = 3 + (J-1)*NQ
      RHS(K) = 0.
      DO 32 M=1,NQ
      MM = M - 3 + KD
  32  DJAC(MM,K) = 0.
      DJAC(KD,K) = DTINV
      DJAC(KU,K+NQ) = 0.
      IF(J.EQ.1) GO TO 33
      DJAC(KL,K-NQ) = 0.
  33  CONTINUE
  34  CONTINUE
C
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INDEX)
      IF(INDEX.NE.0.) write(90,103)N,INDEX
 103  FORMAT(/1X,'N =',I3,5X,'INDEX =',I3)
      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)
C

C   COMPUTE NEW CONCENTRATIONS
      EMAX = 0.
      DO 26 I=1,NQ
      DO 26 J=1,NZ
      K = I + (J-1)*NQ
C-KK	For all runs, inc. standard O2
      IF((I.EQ.LH2S).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHS).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHSO).AND.(Z(J).GT.2.E6)) GOTO 26
      IF((I.EQ.LCH3CL).AND.(Z(J).GT.2.8E6)) GOTO 26
      IF((I.EQ.LH).AND.(Z(J).LT.3.E6)) GOTO 26
      IF (I.EQ.LSO) GOTO 26
      IF (I.EQ.LCL2O2) GOTO 26
      IF((I.EQ.LN2O5).AND.(Z(J).GT.1.E6)) GOTO 26
      IF((I.EQ.LHNO3).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LHO2NO2).AND.(Z(J).GT.5.E6)) GOTO 26
      IF((I.EQ.LHNO2).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LNO3).AND.(Z(J).GT.7.E6)) GOTO 26
      IF((I.EQ.LCLONO2).AND.(Z(J).GT.5.E6)) GOTO 26
      IF((I.EQ.LCH3OOH).AND.(Z(J).GT.7.E6)) GOTO 26

C
      REL(I,J) = RHS(K)/USOL(I,J)
      EREL = ABS(REL(I,J))
      EMAX = AMAX1(EMAX,EREL)
      IF(EREL.LT.EMAX) GO TO 26
      IS = I
      JS = J
      UMAX = USOL(I,J)
      RMAX = RHS(K)
  26  USOL(I,J) = USOL(I,J) + RHS(K)

C-PL  DO NOT LET [H2O] CHANGE AGTER CHEMICAL REACTION
C
      DO 4 J=1,NZ
      IF(Z(J).LT.ZTROP) USOL(3,J) = H2O(J)
   4  CONTINUE

C-AP Adding tridiagonal solver
C       TRIDIAGONAL INVERSION *****
      L=1
      I = NQ + L
      IF(I.EQ.LSO4AER) MZ = 50
C-AP      IF(I.EQ.LS8TESTAER) MZ = 40
C-AP      IF(I.EQ.LHCAER) MZ = NZ
      MZ1 = MZ - 1
      MZP1 = MZ + 1
C
C   COMPUTE ADVECTION TERMS FOR PARTICLES
      DPU(1,L) = WFALL(2,L)*DEN(2)/DEN(1)/(2.*DZ)
      DPL(NZ,L) = WFALL(NZ1,L)*DEN(NZ1)/DEN(NZ)/(2.*DZ)
      DO 38 J=2,NZ1
      DPU(J,L) = WFALL(J+1,L)*DEN(J+1)/DEN(J)/(2.*DZ)
  38  DPL(J,L) = WFALL(J-1,L)*DEN(J-1)/DEN(J)/(2.*DZ)
C                                                                              
C
C   TA = LOWER DIAGONAL, TB = DIAGONAL, TC = UPPER DIAGONAL, TY =
C   RIGHT-HAND SIDE
      DO 70 J=1,NZ
      TA(J) = 0.
      TB(J) = 0.
      TC(J) = 0.
  70  TY(J) = 0.
C
      DO 44 J=1,MZ
      TB(J) = YL(I,J)
  44  TY(J) = YP(I,J)/DEN(J)
C
      DO 45 J=2,MZ1
      TA(J) = - DL(I,J) + DPL(J,L)
      TB(J) = TB(J) + DD(I,J)
  45  TC(J) = - DU(I,J) - DPU(J,L)
C                                                                              
C   BOUNDARY CONDITIONS
      TA(MZ) = - DL(I,MZ) + DPL(MZ,L)
      TB(MZ) = TB(MZ) + DL(I,MZ) + 0.5*WFALL(MZ,L)/DZ
      TB(1) = TB(1) + DU(I,1) + (.01 - 0.5*WFALL(1,L))/DZ
      TC(1) = - DU(I,1) - DPU(1,L)
C
      CALL SGTSL(MZ,TA,TB,TC,TY,NFLAG)
C-AP      STOP
      IF (NFLAG.NE.0) write(90,401) N,NFLAG,I
 401  FORMAT(//1X,'TRIDIAGONAL SOLVER FAILED AT N =',I3,2X,
     2  'NFLAG =',I2,2X,'SPECIES #',I2)
C
      IF(I.EQ.LSO4AER) THEN
        DO 59 J=1,MZ
   59     SO4AER(J) = TY(J)
      ENDIF
C
C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION
      IF(I.EQ.LSO4AER .AND. MZ.NE.NZ) THEN
        DO 61 J=MZP1,NZ
          SO4AER(J) = SO4AER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
   61     SO4AER(J) = AMAX1(SO4AER(J),1E-100)
      ENDIF

      if(NN.eq.NSTEPS) GOTO 2010
C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
      IF(EMAX.GT.0.20)  DT = 0.7*DTSAVE
      IF(EMAX.GT.0.15)  DT = 0.9*DTSAVE
      IF(EMAX.LT.0.10)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.05)  DT = 1.3*DTSAVE
      IF(EMAX.LT.0.03)  DT = 1.5*DTSAVE
      IF(EMAX.LT.0.01)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.003) DT = 5.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 10.*DTSAVE

c Adjusting time step to assure that TIME<=TSTOP
      TIME1 = TIME + DT
      if(TIME1.ge.TSTOP) then
        DT = ABS(TIME-TSTOP)
        NN = NSTEPS -1
      endif
      DTINV = 1./DT
C
2010  ISP = ISPEC(IS)
      ZMAX = Z(JS)
      IF(SM-MS.GT.0.01) GO TO 317
      write(90,100)N,EMAX,ISP,ZMAX,
     & UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
C-AP
C   COMPUTE ATMOSPHERIC OXIDATION STATE

      DO 42 I=1,NQ
      SR(I) = 0.
      DO 43 J=1,JTROP
  43  SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ
      PHIDEP(I) = VDEP(I)*USOL(I,1)*DEN(1)
  42  TLOSS(I) = SR(I) + PHIDEP(I)


      SR(LSO4AER) = 0.
      DO 48 J=1,JTROP
      SR(LSO4AER) = SR(LSO4AER) + RAINGC(LH2SO4,J)*SO4AER(J)*DEN(J)*DZ
  48  CONTINUE
      PHIDEP(LSO4AER) = (WFALL(1,1) + .01) * SO4AER(1) * DEN(1)
      TLOSS(LSO4AER) = SR(LSO4AER) + PHIDEP(LSO4AER)
C
C   COMPUTE SULFUR BUDGET AND READJUST SO2 (H2S) OUTGASSING RATE IF SO
C   DESIRED (PROGRAM IS SET UP FOR PURE SO2 OUTGASSING)
      SLOSS = TLOSS(LH2S) + TLOSS(LHS) + TLOSS(LSO) +
     2  TLOSS(LSO2) + TLOSS(LH2SO4) + TLOSS(LHSO) + 
     3  TLOSS(LSO4AER)
      SLOSSP = SLOSS - TLOSS(LSO2)
      IF (ISULF.EQ.0 .OR. TIME.LT.1.E6) GO TO 316
      IF (LBOUND(LSO2).EQ.2) SGFLUX(LSO2) = SGFLUX(LSO2) * VOLFLX/SLOSS
      IF (LBOUND(LH2S).EQ.2) SGFLUX(LH2S) = SGFLUX(LH2S) * VOLFLX/SLOSS
 316  CONTINUE
C
 317  CONTINUE
C
C   RETRY TIME STEP IF EMAX EXCEEDS 30 PERCENT
      IF(EMAX.LT.0.3) GO TO 28
      write(90,*)'RETRY TIME STEP BECAUSE EMAX=',EMAX
      DT = 0.5*DTSAVE
      TIME = TIME - DTSAVE
      if(TIME.lt.0) TIME = 0.
      DO 27 I=1,NQ
      DO 27 J=1,NZ
  27  USOL(I,J) = USAVE(I,J)
  28  CONTINUE
C
      NS = N/NPR
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF(SN-NS.GE.1.E-4) GO TO 37 !C-PL prevent print redundant info. for first step 
C
      write(90,*)
      write(90,*)'BEFOREOUT'
      write(90,*) 'N=',N,'NN=',NN,'SN=',SN,'NS=',NS


c in may 10 2019, common out call for test 
      CALL OUTPUTP(NN,NSTEPS,TIME,FLOW)
  37  CONTINUE
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
      IF(TIME.GE.TSTOP) NN = NSTEPS - 1
      if(TIME.lt.TSTOP.and.N.eq.NSTEPS) then
       write(*,'(A,I4,A)')'Time not reached after',NSTEPS,
     &  ' of the photochemical model. The run is stopping now.'
       STOP
      endif
          
   1  CONTINUE      
     
C ***** END THE TIME-STEPPING LOOP *****

  22  CONTINUE
      if(ICOUPLe.eq.1) NSTEPS =N

      OPEN(unit=81,file= DIRIO//'/atm_composition.out')
      WRITE(81,399) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 399  FORMAT(1P8E12.5)
      close(81) 
       

c Transfer results to the climate model (ICOUPLE=1)
      if (ICOUPLE.eq.1) then
       DO 255 I=1,NZ
c Transfer O3 and H2O
        WRITE(84,254) Z(I),PRESS(I),O3(I),H2O(I) 
 254    FORMAT(1PE9.3,3(E10.2))
 255   CONTINUE
       close(84)
      endif
C-KK  Need to find coldtrap by locating water mixing ratio minimum. 

	Jcold = 0

	DO J = 1, NZ
 	   IF (JCOLD .EQ. 0) THEN
           IF (T(J) .LT. T(J+1)) JCOLD = J
	   END IF
	END DO

	OPEN(unit=19,file= DIRIO//'/mixing_ratios.out')
c     Transfer surface mixing ratios
        WRITE(19,*) FAR
	WRITE(19,*) USOL(LCH4,1) 
	WRITE(19,*) FCO2
	WRITE(19,*) O2(1)
	WRITE(19,*) Jcold
        WRITE(19,*) O3COL
       close(19)

      GO TO 21
  20  write(90,300)I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
      GO TO 21
  25  write(90,301)IERR
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)
C
  21  CONTINUE
c      close(82)


C-PL Added the oxygen isotope calculation 
C-SH Added code from Alex Pavlov's version 

      PRINT*, ''
      PRINT *, '----------START O ISOTOPE CALCULATION----------'
      PRINT*, ''
      PRINT 631, GPPOXY,GPPCDE
 631  FORMAT(1X,' Gross Primary Productivity FOR O2 CO2= '
     2 ,1PE8.2,2X,1PE8.2,1X,'cm-2s-1')
      IF (COMBIN.EQ.2) THEN
      PRINT*, ' TURN OFF Combinatorial Effect'
      ELSE 
      PRINT*, ' TURN ON Combinatorial Effect'
      END IF

      IF (RECOMB.EQ.1) THEN
      PRINT*, ' TURN OFF Symmetry Effect'
      ELSE 
      PRINT*, ' TURN ON Symmetry Effect'
      END IF

      USOLI = 0.

      
c    JL RE-INITIATE PARAMETERS FROM SUBROUTINE OXYGEN
C
C-------START 17O-------------------------
      OTYPE = 2

      CALL OXYGEN(USOLI,FLOW,OTYPE,DENCOR) !USOLI IS THE OUT&IN_PUT VAR

      GPPOXYI = GPPOXY*DENCOR
      GPPCDEI = GPPCDE*DENCOR 
C-JL DETERMINE THE ASYMMETRIC TO SYMMETRIC OZONE RATIO 
      DO J =1,NZ 
      SYMRATIO17(J) = USOLI(LO2Q,J)/USOLI(LOQO,J) 
      END DO 

      DO J=1,NZ !C-PL DELTA HERE IS FOR DIFFERENT SPECIES IN DIFFERENT LAYER
      DELTAH2CO17(J) = (USOLI(LH2CQ,J)/DENCOR-USOL(LH2CO,J))
     2 /USOL(LH2CO,J)*1000
      DELTAO17(J) = (USOLI(LQ,J)/DENCOR-USOL(LO,J))/USOL(LO,J)*1000
      DELTAH2O17(J) = (USOLI(LH2Q,J)/DENCOR-USOL(LH2O,J))
     2 /USOL(LH2O,J)*1000
      DELTAOH17(J) = (USOLI(LQH,J)/DENCOR-USOL(LOH,J))/USOL(LOH,J)*1000
C-JL     DELTAHO217(J) = (USOLI(LHOQ,J)/2.-USOL(LHO2,J))/USOL(LHO2,J)*1000 
      DELTAHO2A17(J) = (USOLI(LHOQ,J)/DENCOR-USOL(LHO2,J))
     2 /USOL(LHO2,J)*1000
      DELTAHO2B17(J) = (USOLI(LHQO,J)/DENCOR-USOL(LHO2,J))
     2 /USOL(LHO2,J)*1000
      DELTAHO217(J) = (DELTAHO2A17(J)+DELTAHO2B17(J))/2.

      DELTAH2O217(J) = (USOLI(LH2OQ,J)/DENCOR/2.-USOL(LH2O2,J))
     2 /USOL(LH2O2,J)*1000
      DELTAO3A17(J) = (USOLI(LO2Q,J)/DENCOR/2.-USOL(LO3,J))
     2 /USOL(LO3,J)*1000
      DELTAO3B17(J) = (USOLI(LOQO,J)/DENCOR-USOL(LO3,J))
     2 /USOL(LO3,J)*1000
      DELTAO317(J) =((USOLI(LO2Q,J)+USOLI(LOQO,J))/DENCOR/3.
     2 -USOL(LO3,J))/USOL(LO3,J)*1000! 2./ 2./3.*DELTAO3A17(J)+1./3.*DELTAO3B17(J)    !TOTAL DELTA O3
      DELTACO17(J) = (USOLI(LCQ,J)/DENCOR-USOL(LCO,J))/USOL(LCO,J)*1000
      DELTACH3OOHA17(J) = (USOLI(LCH3OQH,J)/DENCOR-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOHB17(J) = (USOLI(LCH3QOH,J)/DENCOR-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOH17(J)=(DELTACH3OOHA17(J)+DELTACH3OOHB17(J))/2. ! TOTAL DELTA CH3OOH
      DELTACH3O2A17(J) = (USOLI(LCH3OQ,J)/DENCOR-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O2B17(J) = (USOLI(LCH3QO,J)/DENCOR-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O217(J)  = (DELTACH3O2A17(J)+DELTACH3O2B17(J))/2. !TOTAL CH3O2
      DELTAN2O17(J) = (USOLI(LN2Q,J)/DENCOR-USOL(LN2O,J))
     2 /USOL(LN2O,J)*1000
      DELTANO17(J) = (USOLI(LNQ,J)/DENCOR-USOL(LNO,J))/USOL(LNO,J)*1000
      DELTANO217(J) = (USOLI(LNOQ,J)/DENCOR/2.-USOL(LNO2,J))
     2 /USOL(LNO2,J)*1000
      DELTAHNO217(J) = (USOLI(LHNOQ,J)/DENCOR/2.-USOL(LHNO2,J))
     2 /USOL(LHNO2,J)*1000
      DELTAHNO317(J) = (USOLI(LHNO2Q,J)/DENCOR/3.-USOL(LHNO3,J))
     2 /USOL(LHNO3,J)*1000

C     DELTAHO2NO2A17(J) = (USOLI(LHOQNO2,J)/2.-USOL(LHO2NO2,J))
C     2 /USOL(LHO2NO2,J)*1000
C JL 
      DELTAHO2NO2A17(J) = (USOLI(LHOQNO2,J)/DENCOR-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO2B17(J) = (USOLI(LHO2NOQ,J)/DENCOR/2.-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO2C17(J) = (USOLI(LHQONO2,J)/DENCOR-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
C JL REDEINE FROM 1/3 TO 1/4.
      DELTAHO2NO217(J) = 1./4.*DELTAHO2NO2A17(J)+1./2.*DELTAHO2NO2B17(J)
     2 + 1./4.*DELTAHO2NO2C17(J) !TOTAL HO2NO2
      DELTANO317(J) = (USOLI(LNO2Q,J)/DENCOR/3.-USOL(LNO3,J))
     2 /USOL(LNO3,J)*1000
      DELTAN2O5A17(J) = (USOLI(LN2O4Q,J)/DENCOR/4.-USOL(LN2O5,J))
     2 /USOL(LN2O5,J)*1000
      DELTAN2O5B17(J) = (USOLI(LN2QO4,J)/DENCOR-USOL(LN2O5,J)) 
     2 /USOL(LN2O5,J)*1000
      DELTAN2O517(J)  = 4./5.*DELTAN2O5A17(J) + 1./5.*DELTAN2O5B17(J) !TOTAL N2O5
      DELTAO217(J) = (USOLI(LOQ,J)/DENCOR/2.-USOL(LO2,J))
     2 /USOL(LO2,J)*1000
      DELTASO17(J) = (USOLI(LSQ,J)/DENCOR-USOL(LSO,J))/USOL(LSO,J)*1000
      DELTASO217(J) = (USOLI(LSOQ,J)/DENCOR/2.-USOL(LSO2,J))
     2 /USOL(LSO2,J)*1000
      DELTAH2SO417(J) = (USOLI(LH2SO3Q,J)/DENCOR/4.-USOL(LH2SO4,J))
     2 /USOL(LH2SO4,J)*1000
      DELTAHSO17(J) = (USOLI(LHSQ,J)/DENCOR-USOL(LHSO,J))
     2 /USOL(LHSO,J)*1000
      DELTACO217(J) = (USOLI(LCOQ,J)/DENCOR/2.-USOL(LCO2,J))
     2 /USOL(LCO2,J)*1000

      ENDDO 

      
      DO 420 I=1,NQI
      SRI(I) = 0.
      DO 430 J=1,JTROP
  430 SRI(I) = SRI(I) + RAINGCI(I,J)*USOLI(I,J)*DEN(J)*1.e5
      PHIDEPI(I) = VDEPI(I)*USOLI(I,1)*DEN(1)
  420 TLOSSI(I) = SRI(I) + PHIDEPI(I)

      DO I=1,NQI
      IF (LBOUNDI(I).EQ.0) THEN
      TLOSSI(I) = SRI(I) + PHIDEPI(I)
      ELSE IF (LBOUNDI(I).NE.0.AND.FLOWI(I).LT.0) THEN
      TLOSSI(I) = SRI(I) - FLOWI(I) 
      ELSE     
      TLOSSI(I) = SRI(I)
      END IF
      END DO

      DO J=1,NZ
      SO4LOSS = SO4LOSS + CONSO4(J)*DEN(J)*1.e5
      ENDDO
      
      PRINT*,'FOR 17O' 
      print*, 'RAINOUT RATE, PHIDEP, AND LOWER B.C.'
      print*, '   FOLLOWED BY TP, TL, FUP, FLOW, CON'
      print *

      print 1800, (ISPECI(I),I=1,10)
      print 1810, 'SR ',(SRI(I),I=1,10)
      print 1810, 'Phi',(PHIDEPI(I),I=1,10)
      print 1911, 'LBC',(LBOUNDI(I),I=1,10)
      print 1810, 'TPI ',(TPI(I),I=1,10)
      print 1810, 'TLI ',(TLI(I),I=1,10)
      print 1810, 'FUP ',(FUPI(I),I=1,10)
      print 1810, 'FLOW',(FLOWI(I),I=1,10)
      print 1810, 'CON ',(CONI(I),I=1,10)
      print*, ''

      print 1800, (ISPECI(I),I=11,20)
      print 1810, 'SR ',(SRI(I),I=11,20)
      print 1810, 'Phi',(PHIDEPI(I),I=11,20)
      print 1911, 'LBC',(LBOUNDI(I),I=11,20)
      print 1810, 'TPI ',(TPI(I),I=11,20)
      print 1810, 'TLI ',(TLI(I),I=11,20)
      print 1810, 'FUP ',(FUPI(I),I=11,20)
      print 1810, 'FLOW',(FLOWI(I),I=11,20)
      print 1810, 'CON ',(CONI(I),I=11,20)
      print*, ''

      print 1800, (ISPECI(I),I=21,30)
      print 1810, 'SR ',(SRI(I),I=21,30)
      print 1810, 'Phi',(PHIDEPI(I),I=21,30)
      print 1911, 'LBC',(LBOUNDI(I),I=21,30)
      print 1810, 'TPI ',(TPI(I),I=21,30)
      print 1810, 'TLI ',(TLI(I),I=21,30)
      print 1810, 'FUP ',(FUPI(I),I=21,30)
      print 1810, 'FLOW',(FLOWI(I),I=21,30)
      print 1810, 'CON ',(CONI(I),I=21,30)
      print*, ''

      print 1800, (ISPECI(I),I=31,NQI)
      print 1810, 'SR ',(SRI(I),I=31,NQI)
      print 1810, 'Phi',(PHIDEPI(I),I=31,NQI)
      print 1911, 'LBC',(LBOUNDI(I),I=31,NQI)
      print 1810, 'TPI ',(TPI(I),I=31,NQI)
      print 1810, 'TLI ',(TLI(I),I=31,NQI)
      print 1810, 'FUP ',(FUPI(I),I=31,NQI)
      print 1810, 'FLOW',(FLOWI(I),I=31,NQI)
      print 1810, 'CON ',(CONI(I),I=31,NQI)
      print*, ''
 1800 format(7x,10(2x,A8,2x))!JL(4x,10(A8))
 1810 format(A5,1P10E12.4)
 1911 FORMAT(A4,7x,10(I1,11X))

      print*, 'delta 17O values versus altitude'
      PRINT 102, 'Z(km)',(ISPECI(I),I=1,10)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
      PRINT 104,J-0.5,
     2 DELTAH2CO17(J),DELTAO17(J),DELTAH2O17(J),
     3 DELTAOH17(J),DELTAHO217(J),DELTAH2O217(J),DELTAO3A17(J),
     4 DELTAO3B17(J),DELTACO17(J),DELTACH3OOHA17(J)
      ENDDO
      print*, ''
      PRINT 102, 'Z(km)',(ISPECI(I),I=11,20)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 104,J-0.5,
     1 DELTACH3OOHB17(J),DELTACH3O2A17(J),DELTACH3O2B17(J),
     2 DELTAN2O17(J),DELTANO17(J),DELTANO217(J),DELTAHNO217(J),
     3 DELTAHNO317(J),DELTAHO2NO2A17(J),DELTAHO2NO2B17(J)
      ENDDO
      print*, ''
      PRINT 1102, 'Z(km)',(ISPECI(I),I=21,30)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 1104,J-0.5,
     1 DELTANO317(J),DELTAN2O5A17(J),DELTAN2O5B17(J),DELTAO217(J),
     2 DELTASO17(J),DELTASO217(J),DELTAH2SO417(J),DELTAHSO17(J),
     3 DELTACO217(J), DELTAHO2B17(J)
      ENDDO
C JL MODIFY THE PRINT OUT 
      PRINT *, ''
      PRINT 1103, 'Z(km)',ISPECI(31),'SYM-17'
      DO J = 1,NZ 
      PRINT 104, J-0.5, DELTAHO2NO2C17(J),SYMRATIO17(J) 
      END DO 

      print*, ''
      PRINT*, 'Isotope with more than one geometric configurations'
      PRINT*,'DELTA'
      PRINT 102, 'Z(km)','O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I',
     2 'HO2I'
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 104,J-0.5,
     1 DELTAO317(J), DELTACH3OOH17(J) ,DELTACH3O217(J),
     2 DELTAHO2NO217(J),DELTAN2O517(J),DELTAHO217(J)
      ENDDO
  102 format(2x,A8,10(A8,2x))
 1103 FORMAT(2X,A8,2(A8,2X))
 1102 format(2x,A8,2x,11(A8,2x))
  104 format(2x,F6.1,1x,1P10E10.2)
 1104 format(2x,F6.1,1x,1P11E10.2)         
c 1031 format(A8,A5,A3,1PE10.2)
c 1034 format(I2,1x,1P12E10.2)
c 1035 format(3x,12(A8,2x))
 
C-PL-JL COMPUTE CONSERVATION OF 17O 
      OXYDEPI = 0.
      OXYRANI = 0. 
      OXYUPI  = 0.
      OXYLOSI = 0.
      OXYDEPI = - (FLOWI(LH2CQ) + FLOWI(LQ) + FLOWI(LQH) + FLOWI(LHOQ)
     2  + FLOWI(LH2OQ) + FLOWI(LO2Q) + FLOWI(LOQO) + FLOWI(LCH3OQH)
     3  + FLOWI(LCH3QOH) + FLOWI(LCH3OQ) + FLOWI(LCH3QO) + FLOWI(LNQ) 
     4  + FLOWI(LNOQ) + FLOWI(LHNOQ)+ FLOWI(LHNO2Q) + FLOWI(LHOQNO2) 
     5  + FLOWI(LHO2NOQ) + FLOWI(LNO2Q) + FLOWI(LN2O4Q)
     6  + FLOWI(LN2QO4) + FLOWI(LSQ) +FLOWI(LSOQ) 
     7  + FLOWI(LH2SO3Q) + FLOWI(LHSQ)+FLOWI(LHQO)+FLOW(LHQONO2)
     8  + FLOWI(LOQ) + FLOWI(LCOQ) + FLOWI(LCQ) + FLOWI(LN2Q)) !

      OXYRANI = SRI(LH2CQ) + SRI(LQ) + SRI(LQH) + SRI(LHOQ)
     2  + SRI(LH2OQ) + SRI(LO2Q) + SRI(LOQO) + SRI(LCH3OQH)
     3  + SRI(LCH3QOH) + SRI(LCH3OQ) + SRI(LCH3QO) + SRI(LNQ) 
     4  + SRI(LNOQ) + SRI(LHNOQ)+ SRI(LHNO2Q) + SRI(LHOQNO2) 
     5  + SRI(LHO2NOQ) + SRI(LNO2Q) + SRI(LN2O4Q)
     6  + SRI(LN2QO4) + SRI(LSQ) +SRI(LSOQ) 
     7  + SRI(LH2SO3Q) + SRI(LHSQ) 
     8  + SRI(LCQ) + SRI(LN2Q) + SRI(LOQ) + SRI(LCOQ)+SRI(LHQO)
     9  + SRI(LHQONO2)!+ PSO4AER

       OXYUPI = FUPI(LH2CQ) + FUPI(LQ) + FUPI(LQH) + FUPI(LHOQ)
     2  + FUPI(LH2OQ) + FUPI(LO2Q) + FUPI(LOQO) + FUPI(LCH3OQH)
     3  + FUPI(LCH3QOH) + FUPI(LCH3OQ) + FUPI(LCH3QO) + FUPI(LNQ) 
     4  + FUPI(LNOQ) + FUPI(LHNOQ)+ FUPI(LHNO2Q) + FUPI(LHOQNO2) 
     5  + FUPI(LHO2NOQ) + FUPI(LNO2Q) + FUPI(LN2O4Q)
     6  + FUPI(LN2QO4) + FUPI(LSQ) +FUPI(LSOQ) 
     7  + FUPI(LH2SO3Q) + FUPI(LHSQ) 
     8  + FUPI(LCQ) + FUPI(LN2Q) + FUPI(LOQ) + FUPI(LCOQ)
     9  + FUPI(LHQO) + FUPI(LHQONO2)

      OXYLOSI = OXYDEPI + OXYRANI + OXYUPI + TPI(LH2Q)-TLI(LH2Q)
 
      OXYPROI = 0.
C      OXYPROI = FLOWI(LCQ) + FLOWI(LN2Q)!SGFLUX(LCO) + SGFLUX(LN2O)
C_JL move the FLOW terms to the OXYDEPI
      OXYPROI =  2.*GPPOXYI+2.*GPPCDEI    

      CONOXYI = OXYLOSI - OXYPROI
      print 177, OXYLOSI,OXYPROI,CONOXYI,
     2 CONOXYI/Min(OXYLOSI,OXYPROI) *100. 
 177  FORMAT(/1X,'CONSERVATION OF OXYGENI:',/5X,'OXYLOSI =',1PE10.3,
     2  2X,'OXYPROI =',E10.3,2X,'CONOXYI =',E10.3,
     3  3X,'ERRI = ',E10.3,' %' )
      print 178, OXYDEPI,OXYRANI,OXYUPI,TPI(LH2Q)-TLI(LH2Q)
 178  FORMAT(/5X,'OXYDEPI =',1PE10.3,2X,'OXYRANI =',E10.3,
     2 2X,'OXYUPI = ',E10.3,2X,'TPI-TLI(H2Q) =',2E10.3)

      PRINT 189, FLOWI(LCQ),FLOWI(LN2Q),2.*GPPOXYI+FLOWI(LOQ)
     2 ,2*GPPOXYI,2.*GPPCDEI+FLOWI(LCOQ),2*GPPCDEI
 189  FORMAT(5X,'FLOW(LCQ) =',1PE10.3,2X,'FLOW(LN2Q) =',E10.3,
     2 /5X,'GPP-FLOW(LOQ)=',E10.3,2X,'GPPOQ=',E10.3,
     2 2X,'GPP-FLOW(LCOQ)=',E10.3,2X,'GPPCOQ=',E10.3)

      print 187, FLOW(LO2)*2.,FLOWI(LOQ)
 187  format(/5X,'2*FLOW(O2) =',1PE12.5,2X,'FLOW(OQ) =',1PE12.5)
      print 188, FLOW(LO2)/SL(LO2,1), FLOW(LCO2)/SL(LCO2,1)
 188  format(5X,'Vdep O2/CO2 from Main code =',1P2E18.8)

      print 289, TPH2OI,TLH2OI,TPH2OI-TLH2OI,TPI(LH2Q)-TPH2OI,
     2 TLI(LH2Q)-TLH2OI,TPI(LH2Q)-TPH2OI-(TLI(LH2Q)-TLH2OI),
     3 TPI(LH2Q),TLI(LH2Q),TPI(LH2Q)-TLI(LH2Q)
 289  FORMAT(/5X,' H2Q     TPI       TLI       TPI-TLI ',
     2       /5X, 'TROPOS',1P3E10.3,
     2       /5X, 'STRATO',1P3E10.3,
     3       /5X, 'TOTAL ',1P3E10.3)


      print*, ''
      print*, '--------out of isotope subroutine for 17O-------------'
      print*, ''
    
C-PL WE ONLY CAL. 17O FRACOXYGEN=0.5305(Crockford et al.,2018)      
      fracoxy =0.528! 0.5305


C-PL PRINT OUT INTEGRATED REACTION RATE FOR EVERY SPECIES 06/12/2019
C

      USOLI = 0.
C-JL------------------ DOING 18O HERE---------------- 
C      go to 1997
      OTYPE = 3
      CALL OXYGEN(USOLI,FLOW,OTYPE,DENCOR)

      GPPOXYI = GPPOXY*DENCOR
      GPPCDEI = GPPCDE*DENCOR   
C-JL DETERMINE THE ASYMMETRIC TO SYMMETRIC OZONE RATIO 
      DO J =1,NZ
      SYMRATIO18(J) = USOLI(LO2Q,J)/USOLI(LOQO,J)
      END DO


      DO J=1,NZ !C-PL DELTA HERE IS FOR DIFFERENT SPECIES IN DIFFERENT LAYER
      DELTAH2CO18(J) = (USOLI(LH2CQ,J)/DENCOR-USOL(LH2CO,J))
     2 /USOL(LH2CO,J)*1000
      DELTAO18(J) = (USOLI(LQ,J)/DENCOR-USOL(LO,J))/USOL(LO,J)*1000
      DELTAH2O18(J) = (USOLI(LH2Q,J)/DENCOR-USOL(LH2O,J))
     2 /USOL(LH2O,J)*1000
      DELTAOH18(J) = (USOLI(LQH,J)/DENCOR-USOL(LOH,J))/USOL(LOH,J)*1000
C-JL     DELTAHO217(J) = (USOLI(LHOQ,J)/2.-USOL(LHO2,J))/USOL(LHO2,J)*1000 
      DELTAHO2A18(J) = (USOLI(LHOQ,J)/DENCOR-USOL(LHO2,J))
     2 /USOL(LHO2,J)*1000
      DELTAHO2B18(J) = (USOLI(LHQO,J)/DENCOR-USOL(LHO2,J))
     2 /USOL(LHO2,J)*1000
      DELTAHO218(J) = (DELTAHO2A18(J)+DELTAHO2B18(J))/2.

      DELTAH2O218(J) = (USOLI(LH2OQ,J)/DENCOR/2.-USOL(LH2O2,J))
     2 /USOL(LH2O2,J)*1000
      DELTAO3A18(J) = (USOLI(LO2Q,J)/DENCOR/2.-USOL(LO3,J))
     2 /USOL(LO3,J)*1000
      DELTAO3B18(J) = (USOLI(LOQO,J)/DENCOR-USOL(LO3,J))
     2 /USOL(LO3,J)*1000
      DELTAO318(J)  =((USOLI(LO2Q,J)+USOLI(LOQO,J))/DENCOR/3.
     2 -USOL(LO3,J))/USOL(LO3,J)*1000! 2./3.*DELTAO3A18(J)+1./3.*DELTAO3B18(J)    !TOTAL DELTA O3
      DELTACO18(J) = (USOLI(LCQ,J)/DENCOR-USOL(LCO,J))/USOL(LCO,J)*1000
      DELTACH3OOHA18(J) = (USOLI(LCH3OQH,J)/DENCOR-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOHB18(J) = (USOLI(LCH3QOH,J)/DENCOR-USOL(LCH3OOH,J))
     2 /USOL(LCH3OOH,J)*1000
      DELTACH3OOH18(J)=(DELTACH3OOHA17(J)+DELTACH3OOHB17(J))/2. ! TOTAL DELTA CH3OOH
      DELTACH3O2A18(J) = (USOLI(LCH3OQ,J)/DENCOR-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O2B18(J) = (USOLI(LCH3QO,J)/DENCOR-USOL(LCH3O2,J))
     2 /USOL(LCH3O2,J)*1000
      DELTACH3O218(J)  = (DELTACH3O2A18(J)+DELTACH3O2B18(J))/2. !TOTAL CH3O2
      DELTAN2O18(J) = (USOLI(LN2Q,J)/DENCOR-USOL(LN2O,J))
     2 /USOL(LN2O,J)*1000
      DELTANO18(J) = (USOLI(LNQ,J)/DENCOR-USOL(LNO,J))/USOL(LNO,J)*1000
      DELTANO218(J) = (USOLI(LNOQ,J)/DENCOR/2.-USOL(LNO2,J))
     2 /USOL(LNO2,J)*1000
      DELTAHNO218(J) = (USOLI(LHNOQ,J)/DENCOR/2.-USOL(LHNO2,J))
     2 /USOL(LHNO2,J)*1000
      DELTAHNO318(J) = (USOLI(LHNO2Q,J)/DENCOR/3.-USOL(LHNO3,J))
     2 /USOL(LHNO3,J)*1000

C     DELTAHO2NO2A17(J) = (USOLI(LHOQNO2,J)/2.-USOL(LHO2NO2,J))
C     2 /USOL(LHO2NO2,J)*1000
C JL 
      DELTAHO2NO2A18(J) = (USOLI(LHOQNO2,J)/DENCOR-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO2B18(J) = (USOLI(LHO2NOQ,J)/DENCOR/2.-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
      DELTAHO2NO2C18(J) = (USOLI(LHQONO2,J)/DENCOR-USOL(LHO2NO2,J))
     2 /USOL(LHO2NO2,J)*1000
C JL REDEINE FROM 1/3 TO 1/4.
      DELTAHO2NO218(J) = 1./4.*DELTAHO2NO2A18(J)+1./2.*DELTAHO2NO2B18(J)
     2 + 1./4.*DELTAHO2NO2C18(J) !TOTAL HO2NO2
      DELTANO318(J) = (USOLI(LNO2Q,J)/DENCOR/3.-USOL(LNO3,J))
     2 /USOL(LNO3,J)*1000
      DELTAN2O5A18(J) = (USOLI(LN2O4Q,J)/DENCOR/4.-USOL(LN2O5,J))
     2 /USOL(LN2O5,J)*1000
      DELTAN2O5B18(J) = (USOLI(LN2QO4,J)/DENCOR-USOL(LN2O5,J)) 
     2 /USOL(LN2O5,J)*1000
      DELTAN2O518(J)  = 4./5.*DELTAN2O5A18(J) + 1./5.*DELTAN2O5B18(J) !TOTAL N2O5
      DELTAO218(J) = (USOLI(LOQ,J)/DENCOR/2.-USOL(LO2,J))
     2 /USOL(LO2,J)*1000
      DELTASO18(J) = (USOLI(LSQ,J)/DENCOR-USOL(LSO,J))/USOL(LSO,J)*1000
      DELTASO218(J) = (USOLI(LSOQ,J)/DENCOR/2.-USOL(LSO2,J))
     2 /USOL(LSO2,J)*1000
      DELTAH2SO418(J) = (USOLI(LH2SO3Q,J)/DENCOR/4.-USOL(LH2SO4,J))
     2 /USOL(LH2SO4,J)*1000
      DELTAHSO18(J) = (USOLI(LHSQ,J)/DENCOR-USOL(LHSO,J))
     2 /USOL(LHSO,J)*1000
      DELTACO218(J) = (USOLI(LCOQ,J)/DENCOR/2.-USOL(LCO2,J))
     2 /USOL(LCO2,J)*1000
      ENDDO 

      
      DO 421 I=1,NQI
      SRI(I) = 0.
      DO 431 J=1,JTROP
  431 SRI(I) = SRI(I) + RAINGCI(I,J)*USOLI(I,J)*DEN(J)*1.e5
      PHIDEPI(I) = VDEPI(I)*USOLI(I,1)*DEN(1)
  421 TLOSSI(I) = SRI(I) + PHIDEPI(I)

      DO I=1,NQI
      IF (LBOUNDI(I).EQ.0) THEN
      TLOSSI(I) = SRI(I) + PHIDEPI(I)
      ELSE IF (LBOUNDI(I).NE.0.AND.FLOWI(I).LT.0) THEN
      TLOSSI(I) = SRI(I) - FLOWI(I) 
      ELSE     
      TLOSSI(I) = SRI(I)
      END IF
      END DO

      DO J=1,NZ
      SO4LOSS = SO4LOSS + CONSO4(J)*DEN(J)*1.e5
      ENDDO
      
      PRINT*,'-------------FOR 18O------------' 
      print*, 'RAINOUT RATE, PHIDEP, AND LOWER B.C.'
      print*, '   FOLLOWED BY TP, TL, FUP, FLOW, CON'
      print *

      print 2800, (ISPECI(I),I=1,10)
      print 2810, 'SR ',(SRI(I),I=1,10)
      print 2810, 'Phi',(PHIDEPI(I),I=1,10)
      print 2911, 'LBC',(LBOUNDI(I),I=1,10)
      print 2810, 'TPI ',(TPI(I),I=1,10)
      print 2810, 'TLI ',(TLI(I),I=1,10)
      print 2810, 'FUP ',(FUPI(I),I=1,10)
      print 2810, 'FLOW',(FLOWI(I),I=1,10)
      print 2810, 'CON ',(CONI(I),I=1,10)
      print*, ''

      print 2800, (ISPECI(I),I=11,20)
      print 2810, 'SR ',(SRI(I),I=11,20)
      print 2810, 'Phi',(PHIDEPI(I),I=11,20)
      print 2911, 'LBC',(LBOUNDI(I),I=11,20)
      print 2810, 'TPI ',(TPI(I),I=11,20)
      print 2810, 'TLI ',(TLI(I),I=11,20)
      print 2810, 'FUP ',(FUPI(I),I=11,20)
      print 2810, 'FLOW',(FLOWI(I),I=11,20)
      print 2810, 'CON ',(CONI(I),I=11,20)
      print*, ''

      print 2800, (ISPECI(I),I=21,30)
      print 2810, 'SR ',(SRI(I),I=21,30)
      print 2810, 'Phi',(PHIDEPI(I),I=21,30)
      print 2911, 'LBC',(LBOUNDI(I),I=21,30)
      print 2810, 'TPI ',(TPI(I),I=21,30)
      print 2810, 'TLI ',(TLI(I),I=21,30)
      print 2810, 'FUP ',(FUPI(I),I=21,30)
      print 2810, 'FLOW',(FLOWI(I),I=21,30)
      print 2810, 'CON ',(CONI(I),I=21,30)
      print*, ''

      print 2800, (ISPECI(I),I=31,NQI)
      print 2810, 'SR ',(SRI(I),I=31,NQI)
      print 2810, 'Phi',(PHIDEPI(I),I=31,NQI)
      print 2911, 'LBC',(LBOUNDI(I),I=31,NQI)
      print 2810, 'TPI ',(TPI(I),I=31,NQI)
      print 2810, 'TLI ',(TLI(I),I=31,NQI)
      print 2810, 'FUP ',(FUPI(I),I=31,NQI)
      print 2810, 'FLOW',(FLOWI(I),I=31,NQI)
      print 2810, 'CON ',(CONI(I),I=31,NQI)
      print*, ''
 2800 format(7x,10(2x,A8,2x))!JL(4x,10(A8))
 2810 format(A5,1P10E12.4)
 2911 FORMAT(A4,7x,10(I1,11X))

      print*, 'delta 18O values versus altitude'
      PRINT 302, 'Z(km)',(ISPECI(I),I=1,10)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
      PRINT 304,J-0.5,
     2 DELTAH2CO18(J),DELTAO18(J),DELTAH2O18(J),
     3 DELTAOH18(J),DELTAHO218(J),DELTAH2O218(J),DELTAO3A18(J),
     4 DELTAO3B18(J),DELTACO18(J),DELTACH3OOHA18(J)
      ENDDO
      print*, ''
      PRINT 302, 'Z(km)',(ISPECI(I),I=11,20)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 304,J-0.5,
     1 DELTACH3OOHB18(J),DELTACH3O2A18(J),DELTACH3O2B18(J),
     2 DELTAN2O18(J),DELTANO18(J),DELTANO218(J),DELTAHNO218(J),
     3 DELTAHNO318(J),DELTAHO2NO2A18(J),DELTAHO2NO2B18(J)
      ENDDO
      print*, ''
      PRINT 2102, 'Z(km)',(ISPECI(I),I=21,30)
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 2104,J-0.5,
     1 DELTANO318(J),DELTAN2O5A18(J),DELTAN2O5B18(J),DELTAO218(J),
     2 DELTASO18(J),DELTASO218(J),DELTAH2SO418(J),DELTAHSO18(J),
     3 DELTACO218(J), DELTAHO2B18(J)
      ENDDO
C JL MODIFY THE PRINT OUT 
      PRINT *, ''
      PRINT 2103, 'Z(km)',ISPECI(31), 'SYM-18'
      DO J = 1,NZ 
      PRINT 304, J-0.5, DELTAHO2NO2C18(J),SYMRATIO18(J) 
      END DO 

      print*, ''
      print*,'-------18O-------------------'
      PRINT*, 'Isotope with more than one geometric configurations'
      PRINT*,'DELTA'
      PRINT 302, 'Z(km)','O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I',
     2 'HO2I'
      DO J=1,NZ
      !if(J.le.10.or.modulo(J,4).eq.1) 
       PRINT 304,J-0.5,
     1 DELTAO318(J), DELTACH3OOH18(J) ,DELTACH3O218(J),
     2 DELTAHO2NO218(J),DELTAN2O518(J),DELTAHO218(J)
      ENDDO
  302 format(2x,A8,10(A8,2x))
 2103 FORMAT(2X,A8,2(A8,2X))
 2102 format(2x,A8,2x,11(A8,2x))
  304 format(2x,F6.1,1x,1P10E10.2)
 2104 format(2x,F6.1,1x,1P11E10.2)         
c 1031 format(A8,A5,A3,1PE10.2)
c 1034 format(I2,1x,1P12E10.2)
c 1035 format(3x,12(A8,2x))
 
C-JL COMPUTE CONSERVATION OF 18O
      OXYDEPI = 0.
      OXYRANI = 0. 
      OXYUPI  = 0.
      OXYLOSI = 0.
      OXYDEPI = - (FLOWI(LH2CQ) + FLOWI(LQ) + FLOWI(LQH) + FLOWI(LHOQ)
     2  + FLOWI(LH2OQ) + FLOWI(LO2Q) + FLOWI(LOQO) + FLOWI(LCH3OQH)
     3  + FLOWI(LCH3QOH) + FLOWI(LCH3OQ) + FLOWI(LCH3QO) + FLOWI(LNQ) 
     4  + FLOWI(LNOQ) + FLOWI(LHNOQ)+ FLOWI(LHNO2Q) + FLOWI(LHOQNO2) 
     5  + FLOWI(LHO2NOQ) + FLOWI(LNO2Q) + FLOWI(LN2O4Q)
     6  + FLOWI(LN2QO4) + FLOWI(LSQ) +FLOWI(LSOQ) 
     7  + FLOWI(LH2SO3Q) + FLOWI(LHSQ)+FLOWI(LHQO)+FLOW(LHQONO2)
     8  + FLOWI(LOQ) + FLOWI(LCOQ) + FLOWI(LCQ) + FLOWI(LN2Q)) !

      OXYRANI = SRI(LH2CQ) + SRI(LQ) + SRI(LQH) + SRI(LHOQ)
     2  + SRI(LH2OQ) + SRI(LO2Q) + SRI(LOQO) + SRI(LCH3OQH)
     3  + SRI(LCH3QOH) + SRI(LCH3OQ) + SRI(LCH3QO) + SRI(LNQ) 
     4  + SRI(LNOQ) + SRI(LHNOQ)+ SRI(LHNO2Q) + SRI(LHOQNO2) 
     5  + SRI(LHO2NOQ) + SRI(LNO2Q) + SRI(LN2O4Q)
     6  + SRI(LN2QO4) + SRI(LSQ) +SRI(LSOQ) 
     7  + SRI(LH2SO3Q) + SRI(LHSQ) 
     8  + SRI(LCQ) + SRI(LN2Q) + SRI(LOQ) + SRI(LCOQ)+SRI(LHQO)
     9  + SRI(LHQONO2)!+ PSO4AER

       OXYUPI = FUPI(LH2CQ) + FUPI(LQ) + FUPI(LQH) + FUPI(LHOQ)
     2  + FUPI(LH2OQ) + FUPI(LO2Q) + FUPI(LOQO) + FUPI(LCH3OQH)
     3  + FUPI(LCH3QOH) + FUPI(LCH3OQ) + FUPI(LCH3QO) + FUPI(LNQ) 
     4  + FUPI(LNOQ) + FUPI(LHNOQ)+ FUPI(LHNO2Q) + FUPI(LHOQNO2) 
     5  + FUPI(LHO2NOQ) + FUPI(LNO2Q) + FUPI(LN2O4Q)
     6  + FUPI(LN2QO4) + FUPI(LSQ) +FUPI(LSOQ) 
     7  + FUPI(LH2SO3Q) + FUPI(LHSQ) 
     8  + FUPI(LCQ) + FUPI(LN2Q) + FUPI(LOQ) + FUPI(LCOQ)
     9  + FUPI(LHQO) + FUPI(LHQONO2)

      OXYLOSI = OXYDEPI + OXYRANI + OXYUPI + TPI(LH2Q)-TLI(LH2Q)
 
      OXYPROI = 0.
C      OXYPROI = FLOWI(LCQ) + FLOWI(LN2Q)!SGFLUX(LCO) + SGFLUX(LN2O)

C_JL move the FLOW terms to the OXYDEPI
      OXYPROI =  2.*GPPOXYI+2.*GPPCDEI

      CONOXYI = OXYLOSI - OXYPROI
      print *
      print *,'-----18O budget-----------------'
      print 277, OXYLOSI,OXYPROI,CONOXYI,
     2 CONOXYI/MIN(OXYLOSI,OXYPROI) *100. 
 277  FORMAT(/1X,'CONSERVATION OF OXYGENI:',/5X,'OXYLOSI =',1PE10.3,
     2  2X,'OXYPROI =',E10.3,2X,'CONOXYI =',E10.3,
     3  3X,'ERRI = ',E10.3,' %' )
      print 278, OXYDEPI,OXYRANI,OXYUPI,TPI(LH2Q)-TLI(LH2Q)
 278  FORMAT(/5X,'OXYDEPI =',1PE10.3,2X,'OXYRANI =',E10.3,
     2 2X,'OXYUPI = ',E10.3,2X,'TPI-TLI(H2Q) =',2E10.3)

      PRINT 489, FLOWI(LCQ),FLOWI(LN2Q),2.*GPPOXYI+FLOWI(LOQ)
     2 ,2*GPPOXYI,2.*GPPCDEI+FLOWI(LCOQ),2*GPPCDEI
 489  FORMAT(5X,'FLOW(LCQ) =',1PE10.3,2X,'FLOW(LN2Q) =',E10.3,
     2 /5X,'GPP-FLOW(LOQ)=',E10.3,2X,'GPPOQ=',E10.3,
     2 2X,'GPP-FLOW(LCOQ)=',E10.3,2X,'GPPCOQ=',E10.3)

      print 287, FLOW(LO2)*2.,FLOWI(LOQ)
 287  format(/5X,'2*FLOW(O2) =',1PE12.5,2X,'FLOW(OQ) =',1PE12.5)
      print 288, FLOW(LO2)/SL(LO2,1), FLOW(LCO2)/SL(LCO2,1)
 288  format(5X,'Vdep O2/CO2 from Main code =',1P2E18.8)

      print 389, TPH2OI,TLH2OI,TPH2OI-TLH2OI,TPI(LH2Q)-TPH2OI,
     2 TLI(LH2Q)-TLH2OI,TPI(LH2Q)-TPH2OI-(TLI(LH2Q)-TLH2OI),
     3 TPI(LH2Q),TLI(LH2Q),TPI(LH2Q)-TLI(LH2Q)
 389  FORMAT(/5X,' H2Q     TPI       TLI       TPI-TLI ',
     2       /5X, 'TROPOS',1P3E10.3,
     2       /5X, 'STRATO',1P3E10.3,
     3       /5X, 'TOTAL ',1P3E10.3)

      DO J=1,NZ
      C2DELTAH2CO(J)=1000.*LOG(DELTAH2CO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAH2CO18(J)/1000.+1.)

      C2DELTAO(J)=1000.*LOG(DELTAO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAO18(J)/1000.+1.)

      C2DELTAH2O(J)=1000.*LOG(DELTAH2O17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAH2O18(J)/1000.+1.)

      C2DELTAOH(J)=1000.*LOG(DELTAOH17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAOH18(J)/1000.+1.)

      C2DELTAHO2(J)=1000.*LOG(DELTAHO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHO218(J)/1000.+1.)

      C2DELTAH2O2(J)=1000.*LOG(DELTAH2O217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAH2O218(J)/1000.+1.)

      C2DELTAO3A(J)=1000.*LOG(DELTAO3A17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAO3A18(J)/1000.+1.)

      C2DELTAO3B(J)=1000.*LOG(DELTAO3B17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAO3B18(J)/1000.+1.)

      C2DELTAO3(J)=1000.*LOG(DELTAO317(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAO318(J)/1000.+1.) 
  
      C2DELTACO(J)=1000.*LOG(DELTACO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACO18(J)/1000.+1.)

      C2DELTACH3OOHA(J)=1000.*LOG(DELTACH3OOHA17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3OOHA18(J)/1000.+1.)

      C2DELTACH3OOHB(J)=1000.*LOG(DELTACH3OOHB17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3OOHB18(J)/1000.+1.)

      C2DELTACH3OOH(J)=1000.*LOG(DELTACH3OOH17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3OOH18(J)/1000.+1.)

      C2DELTACH3O2A(J)=1000.*LOG(DELTACH3O2A17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3O2A18(J)/1000.+1.)

      C2DELTACH3O2B(J)=1000.*LOG(DELTACH3O2B17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3O2B18(J)/1000.+1.)

      C2DELTACH3O2(J)=1000.*LOG(DELTACH3O217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACH3O218(J)/1000.+1.)

      C2DELTAN2O(J)=1000.*LOG(DELTAN2O17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAN2O18(J)/1000.+1.)

      C2DELTANO(J)=1000.*LOG(DELTANO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTANO18(J)/1000.+1.)

      C2DELTANO2(J)=1000.*LOG(DELTANO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTANO218(J)/1000.+1.)

      C2DELTAHNO2(J)=1000.*LOG(DELTAHNO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHNO218(J)/1000.+1.)

      C2DELTAHNO3(J)=1000.*LOG(DELTAHNO317(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHNO318(J)/1000.+1.)

      C2DELTAHO2NO2A(J)=1000.*LOG(DELTAHO2NO2A17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHO2NO2A18(J)/1000.+1.)

      C2DELTAHO2NO2B(J)=1000.*LOG(DELTAHO2NO2B17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHO2NO2B18(J)/1000.+1.)

      C2DELTAHO2NO2(J)=1000.*LOG(DELTAHO2NO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHO2NO218(J)/1000.+1.)

      C2DELTANO3(J)=1000.*LOG(DELTANO317(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTANO318(J)/1000.+1.)

      C2DELTAN2O5A(J)=1000.*LOG(DELTAN2O5A17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAN2O5A18(J)/1000.+1.)

      C2DELTAN2O5B(J)=1000.*LOG(DELTAN2O5B17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAN2O5B18(J)/1000.+1.)

      C2DELTAN2O5(J)=1000.*LOG(DELTAN2O517(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAN2O518(J)/1000.+1.)

      C2DELTAO2(J)=1000.*LOG(DELTAO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAO218(J)/1000.+1.)

      C2DELTASO(J)=1000.*LOG(DELTASO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTASO18(J)/1000.+1.)

      C2DELTASO2(J)=1000.*LOG(DELTASO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTASO218(J)/1000.+1.)

      C2DELTAH2SO4(J)=1000.*LOG(DELTAH2SO417(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAH2SO418(J)/1000.+1.)

      C2DELTAHSO(J)=1000.*LOG(DELTAHSO17(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTAHSO18(J)/1000.+1.)

      C2DELTACO2(J)=1000.*LOG(DELTACO217(J)/1000.+1.)
     2 -fracoxy*1000.*LOG(DELTACO218(J)/1000.+1.)
      ENDDO

      PRINT *,'cap-DELTA 17O profiles'
      print *
      PRINT 1804,'Z(km)',(ISPECI(I),I=1,10)
 1804 format(16x,11(A10))
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      PRINT 889, ZKM(J),
     2 C2DELTAH2CO(J),C2DELTAO(J),C2DELTAH2O(J),C2DELTAOH(J),
     3 C2DELTAHO2(J), C2DELTAH2O2(J),C2DELTAO3A(J),C2DELTAO3B(J),
     4 C2DELTACO(J),C2DELTACH3OOHA(J)
 889  FORMAT(18X,F6.1,1P10E10.2)
      ENDDO

      PRINT*, ''
      PRINT 1804,'Z(km)',(ISPECI(I),I=11,20)
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      PRINT 889, ZKM(J),
     2 C2DELTACH3OOHB(J),C2DELTACH3O2A(J),C2DELTACH3O2B(J),
     3 C2DELTAN2O(J), C2DELTANO(J),C2DELTANO2(J),C2DELTAHNO2(J),
     4 C2DELTAHNO3(J),C2DELTAHO2NO2A(J),C2DELTAHO2NO2B(J)
      ENDDO

      PRINT*, ''
      PRINT 1904,'Z(km)',(ISPECI(I),I=21,29)
 1904 format(16x,12(A10))
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      PRINT 889, ZKM(J),
     2 C2DELTANO3(J),C2DELTAN2O5A(J),C2DELTAN2O5B(J),C2DELTAO2(J),
     3 C2DELTASO(J), C2DELTASO2(J),C2DELTAH2SO4(J),C2DELTAHSO(J),
     4 C2DELTACO2(J)
      ENDDO

! PRINT THE SPECIES HAVE TWO ISOTOPIC SPECIES
      PRINT*, ''
      print 1806, 'Z(km)','O3I','CH3OOHI','CH3O2I','HO2NO2I','N2O5I'
     2 ,'HO2I' 
 1806 format(17x,A7,1x,A7,5X,A7,2X,A7,4X,A7,2X,A7,2X,A7)
      DO J=1,NZ
      ZKM = Z(J)/1.E5
      PRINT 889, ZKM(J),
     2 C2DELTAO3(J),C2DELTACH3OOH(J),C2DELTACH3O2(J),
     3 C2DELTAHO2NO2(J),C2DELTAN2O5(J),C2DELTAHO2(J)
      ENDDO

 1036 format(1P5E10.2,1x,1p4e10.2,1x,0p,4(f6.2,1x))
 1066 format(1P6E10.2,1x,1p5e10.2,1x,0p,5(f6.2,1x))
 1076 format(F6.1,1P4E10.2)
 1086 format(F6.1,1P9E10.2)
 1096 format(1P4E10.2,1x,0p,4(f6.1,1x),4(f6.2,1x),4(f6.1,1x))
      print*, ''
      
      PRINT*, "Success!"
 1039 FORMAT('[',F4.1,' , ',1PE10.2,' , ',1PE10.2,' , ',1PE10.2,']')

C 1997 continue 

      OPEN(UNIT=90,FILE='OXYGEN_MIF/Reac_List.Oxy')
      READ(90,205) CHEMI
 205  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      DO 702 I=1,NSPI
         ISP = ISPECI(I)
         WRITE(16,703) ISP,TPI(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TPI = ',1PE9.2)
       DO 704 N=1,NRI 
          IF(JCHEMI(3,N).EQ.I .OR. JCHEMI(4,N).EQ.I .OR. 
     2       JCHEMI(5,N).EQ.I)THEN
         IF(RATI(N).NE.0.) WRITE(16,705) N,(CHEMI(J,N),J=1,5),RATI(N)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
      WRITE(16,706) ISP,TLI(I)
 706  FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TLI = ',1PE9.2)
       DO 707 N=1,NRI 
          IF(JCHEMI(1,N).EQ.I .OR. JCHEMI(2,N).EQ.I)THEN
      IF(RATI(N).NE.0.) WRITE(16,705) N,(CHEMI(J,N),J=1,5),RATI(N)
          ENDIF
 707   CONTINUE
 702  CONTINUE


      STOP
      END program atm_chem

c---------------------------------------------------------------


      SUBROUTINE OUTPUTP(N,NSTEPS,TIME,FLOW)
       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'
      INCLUDE 'INCLUDECHEM/comSATBLK.inc'
      INCLUDE 'INCLUDECHEM/comRRATS.inc'
      INCLUDE 'INCLUDECHEM/comPRESS1.inc'
c in may 10 2019, JL include Rblcok for JCHEM
      INCLUDE 'INCLUDECHEM/comRBLOK.inc'
c in may 10 2019, JL include chem; other two are place holder
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)  
      DIMENSION FUP(NQT),FLOW(NQT),CON(NQT),FLUXCH(NQT,NZ)
     2  ,ZF(NZ)
c in may 16 2019, JL include the following to print h2so4.pdat
      DIMENSION TTAB(51),PH2O(51,34),PH2SO4(51,34)
      REAL LIGHTNO2,LIGHTNO2I
C
c in may 10 2019, to include primo3s to this subroutine
      OPEN(unit=61, file= 'DATA'//'/primo3s.chm', status='old')
c in may 16 2019, JL add a clean printout for h2so4 to see possible disarray of input 
c      open(UNIT=1958, file='DIRIO'//'/h2so4.out.dat') 

      JS=N

      ISKIP = 4
      JSKIP = ISKIP
      IF(N.EQ.NSTEPS) ISKIP = 2
      TIMEY = TIME/3600./24./365.
      write(90,100)TIME,TIMEY
 100  FORMAT(/1X,'TIME =',E11.4,5X,'TIMEY =',1pe13.4,1X,'YEARS')
      write(90,101)NPHOT
 101  FORMAT(/1X,'NPHOT =',I3)
C
      write(90,105)
 105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES'/)
      IROW = 12
      LR = NQ/IROW + 1
      RL = FLOAT(NQ)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
 110  FORMAT(/5X,'Z',8X,13(A8,1X))
      DO 20 I=1,3
  20  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
      DO 21 I=4,NZ,ISKIP
  21  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
 120  FORMAT(1X,1P13E9.2)

      write(90,114) (ISPEC(K),K=K1,K2)
 114  FORMAT(5X,'Z',8X,13(A8,1X))
      IF (N.EQ.0) GO TO 8
      write(90,140)
 140  FORMAT(1X,'TP, TL')
      write(90,145)(TP(K),K=K1,K2)
      write(90,145)(TL(K),K=K1,K2)
 145  FORMAT(10X,1P12E9.2)
   8  CONTINUE
C
C-AP
      write(90,106)
 106  FORMAT(//1X,'MIXING RATIO OF AEROSOL'/)
      write(90,185)
 185  FORMAT(5X,'Z',6X,'SO4AER')
      DO 18 J=1,3
  18  write(90,182) Z(J),SO4AER(J)
      DO 19 J=4,NZ,ISKIP
  19  write(90,182) Z(J),SO4AER(J)
 182  FORMAT(1X,1P2E9.2)
C-AP
      IF (N.EQ.0) RETURN
C
      write(90,183) TP(LSO4AER)
      write(90,184) TL(LSO4AER)
 183  FORMAT(/2X,'TP',6X,1P1E9.2)
 184  FORMAT(2X,'TL',6X,1P1E9.2)
C
      O3_ATMCM = O3COL/2.687e19
      write(90,150) O3COL
 150  FORMAT(//1X,'OZONE COLUMN DEPTH = ',1PE11.4)
      write(90,1150) O3_ATMCM
 1150 FORMAT(1X,'OZONE COLUMN DEPTH = ',F6.3,' ATM CM')
C-AP
      write(90,152) H2SCOL,SO2COL
 152  FORMAT(/1X,'SULFUR COLUMN DEPTHS:  H2S =',1PE10.3,2X,'SO2 =',
     2  E10.3,2X)
C-AP
	DO i = 1, NZ
	 IF (USOL(3,i) .LT. USOL(3,i+1)) THEN
		JCOLD = i
		GOTO 352
	 END IF
	END DO
 352  CONTINUE
      write(90,151) JCOLD, USOL(3,JCOLD)
 151  FORMAT(/1X,I3,' FH2O AT COLD TRAP =',1PE10.3)
      IF(N.LT.NSTEPS) RETURN
C
C ***** PRINT ON LAST ITERATION ONLY *****
      DO 1 I=1,NZ
   1  ZF(I) = Z(I) + 0.5*DZ
C
      DO 3 K=1,NQ
      DO 2 I=1,NZ
   2  SL(K,I) = USOL(K,I)*DEN(I)
!C-PL Added molecular diffusion terms below because they were 
!     missing (7/25/19)
      DO 4 I=1,NZ1
   4  FLUXCH(K,I) = - DK(K,I)*(USOL(K,I+1) - USOL(K,I))/DZ- 0.5*
     2  (HI(K,I)*DEN(I)*USOL(K,I) + HI(K,I+1)*DEN(I+1)*USOL(K,I+1))
   3  CONTINUE

C
      K = LSO4AER
      J = 1
      DO 5 I=1,NZ1
      FLUXCH(K,I) = -DK(K,I)*(SO4AER(I+1) - SO4AER(I))/DZ
     2  - 0.5*(WFALL(I,J)*DEN(I)*SO4AER(I) + WFALL(I+1,J)*DEN(I+1)
     3         *SO4AER(I+1))
   5  CONTINUE
C
!      print *,'Calculating lower boundary flux, FLOW, for O2'
      DO 15 K=1,NQT
      FLOW(K) = FLUXCH(K,1) - (YP(K,1) - YL(K,1)*SL(K,1))*DZ
      FUP(K) = FLUXCH(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*SL(K,NZ))*DZ
      CON(K) = TP(K) - TL(K) + FLOW(K) - FUP(K)
 502  format('K=',I2,' FLOW(K)=',1PE10.3,' FUP(K)=',E10.3)
  15  CONTINUE

C-AP
      FLOW(3) = FLUXCH(3,11)
      CON(3) = TP(3) - TL(3) + FLOW(3) - FUP(3)
      DO 6 I=1,10
   6  FLUXCH(3,I) = 0.
C
      write(90,125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)
      DO 9 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 22 I=1,NZ,ISKIP
       write(90,120) Z(I),(SL(K,I),K=K1,K2)
  22  CONTINUE
   9  CONTINUE
 355  FORMAT(1PE10.3,1PE12.3)
c     close(83)

      ISKIP = 4
      write(90,155)
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED SPECIES'/)
      ZFL = 0.
      ZFT = ZF(NZ)
      DO 10 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQT
C-AP      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,120) ZFL,(FLOW(K),K=K1,K2)
      DO 23 I=1,NZ,ISKIP
  23  write(90,120) ZF(I),(FLUXCH(K,I),K=K1,K2)
      write(90,120) ZFT,(FUP(K),K=K1,K2)
  10  CONTINUE
C-PL RECALCULATE TLOSS   07/2019
      DO I=1,NQ
      IF (LBOUND(I).EQ.0) THEN
      TLOSS(I) = SR(I) + PHIDEP(I)
      ELSE IF (LBOUND(I).NE.0.AND.FLOW(I).LT.0) THEN
      TLOSS(I) = SR(I) -FLOW(I) 
      ELSE     
      TLOSS(I) = SR(I)
      END IF
      END DO

      write(90,175)
 175  FORMAT(/1X,'RAINOUT RATE, PHIDEP, TLOSS AND LOWER B.C.'/)
      write(90,176)
 176  FORMAT(1X,'FOLLOWED BY TP, TL, FUP, FLOW, CON'/)
      DO 13 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQT
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,145) (SR(K),K=K1,K2)
      write(90,145)(PHIDEP(K),K=K1,K2)
      write(90,145)(TLOSS(K),K=K1,K2) 
      write(90,146) (LBOUND(K),K=K1,K2)
 146  FORMAT(14X,12(I1,8X))
      write(90,145)
      write(90,145) (TP(K),K=K1,K2)
      write(90,145) (TL(K),K=K1,K2)
      write(90,145) (FUP(K),K=K1,K2)
      write(90,145) (FLOW(K),K=K1,K2)
      write(90,145) (CON(K),K=K1,K2)     
  13  CONTINUE


C-PL COMPUTE CONSERVATION OF OXYGEN  07/2019 
      OXYDEP = 0.
      OXYDEP = - (FLOW(LH2CO) + FLOW(LO) + FLOW(LOH) + FLOW(LHO2)*2.
     2  + FLOW(LH2O2)*2. + FLOW(LO3)*3. + FLOW(LCH3OOH)*2.
     3  + FLOW(LCH3O2)*2. + FLOW(LNO) + FLOW(LNO2)*2. + FLOW(LHNO2)*2.
     4  + FLOW(LHNO3)*3. + FLOW(LHO2NO2)*4. + FLOW(LNO3)*3. 
     5  + FLOW(LN2O5)*5. + FLOW(LSO) + FLOW(LSO2)*2. 
     6  + FLOW(LH2SO4) * 4. + FLOW(LHSO)  )
      OXYRAN = SR(LH2CO) + SR(LO) + SR(LOH) + SR(LHO2) * 2.
     2  + SR(LH2O2)*2. + SR(LO3)*3. + SR(LCH3OOH)*2.
     3  + SR(LCH3O2)*2. + SR(LNO) + SR(LNO2)*2. + SR(LHNO2)*2.
     4  + SR(LHNO3)*3. + SR(LHO2NO2)*4. + SR(LNO3)*3. 
     5  + SR(LN2O5)*5. + SR(LSO) + SR(LSO2)*2. 
     6  + SR(LH2SO4) * 4. + SR(LHSO) + TP(LSO4AER) * 4.
     7  + SR(LCO) + SR(LN2O) + SR(LO2)*2. + SR(LCO2) *2.
       OXYUP = FUP(LH2CO) + FUP(LO) + FUP(LOH) + FUP(LHO2) * 2.
     2  + FUP(LH2O2)*2. + FUP(LO3)*3. + FUP(LCH3OOH)*2.
     3  + FUP(LCH3O2)*2. + FUP(LNO) + FUP(LNO2)*2. + FUP(LHNO2)*2.
     4  + FUP(LHNO3)*3. + FUP(LHO2NO2)*4. + FUP(LNO3)*3. 
     5  + FUP(LN2O5)*5. + FUP(LSO) + FUP(LSO2)*2. 
     6  + FUP(LH2SO4) * 4. + FUP(LHSO) 
     7  + FUP(LCO) + FUP(LN2O) + FUP(LO2)*2. + FUP(LCO2) *2.
      OXYLOS = OXYDEP + OXYRAN + OXYUP + TP(LH2O)-TL(LH2O)
      OXYPRO = 0.
      OXYPRO = FLOW(LCO) + FLOW(LN2O)!SGFLUX(LCO) + SGFLUX(LN2O)     
      OXYPRO = OXYPRO + 2.*GPPOXY + FLOW(LO2)*2.
     2                + 2.*GPPCDE + FLOW(LCO2)*2.

      CONOXY = OXYLOS - OXYPRO
      write(90,177)OXYLOS,OXYPRO,CONOXY,CONOXY/MIN(OXYLOS,OXYPRO)*100.
 177  FORMAT(/1X,'CONSERVATION OF OXYGEN:',/5X,'OXYLOS =',1PE10.3,
     2  2X,'OXYPRO =',E10.3,2X,'CONOXY =',E10.3,
     3  3X,'ERR = ',E10.3,'%' )
      write(90,178) OXYDEP,OXYRAN,OXYUP,TP(LH2O)-TL(LH2O)
 178  FORMAT(/5X,'OXYDEP =',1PE10.3,2X,'OXYRAN =',E10.3,
     2 2X,'OXYUP =',E10.3,2X,'TP-TL(H2O) =',2E10.3)

      write(90,189) FLOW(LCO),FLOW(LN2O),2.*GPPOXY+FLOW(LO2)*2.
     2 ,GPPOXY,2.*GPPCDE+FLOW(LCO2)*2,GPPCDE
 189  FORMAT(5X,'FLOW(LCO) =',1PE10.3,2X,'FLOW(LN2O) =',E10.3,
     2 /5X,'GPP-FLOW(LO2)=',E10.3,2X,'GPPOXY =',E10.3,
     3 2X,'GPP-FLOW(LCO2)=',E10.3,2X,'GPPCDE =',E10.3)

      write(90,288) TPH2O,TLH2O,TPH2O-TLH2O,
     2 TP(LH2O)-TPH2O,TL(LH2O)-TLH2O,TP(LH2O)-TPH2O-(TL(LH2O)-TLH2O),
     3 TP(LH2O),TL(LH2O),TP(LH2O)-TL(LH2O)
 288  FORMAT(/5X,' H2O     TP        TL        TP-TL ',
     2       /5X, 'TROPOS',1P3E10.3,
     2       /5X, 'STRATO',1P3E10.3,
     3       /5X, 'TOTAL ',1P3E10.3)

      write(90,285) ABS(FLOW(LO2)),SL(LO2,1)
 285  FORMAT(5x,'FLOW(O2)=',F20.3,/5X,'SL(O2)=',F25.3)
      write(90,286) ABS(FLOW(LCO2)),SL(LCO2,1)
 286  FORMAT(5x,'FLOW(CO2)=',F20.3,/5X,'SL(CO2)=',F25.3)

      write(90,388)
 388  FORMAT(/1X,'RAINOUT RATE, FLOW, FUP,TLOSS')

      write(90,110)  ISPEC(LH2CO),ISPEC(LO),ISPEC(LOH),ISPEC(LHO2),
     2 ISPEC(LH2O2),ISPEC(LO3),ISPEC(LCO),ISPEC(LCH3OOH),ISPEC(LCH3O2),
     3 ISPEC(LN2O)
      write(90,149) 'RAIN',SR(LH2CO),SR(LO),SR(LOH),SR(LHO2)*2.,
     2 SR(LH2O2)*2,SR(LO3)*3.,SR(LCO),SR(LCH3OOH)*2.,SR(LCH3O2)*2.,
     3 SR(LN2O)
      write(90,149) 'FLOW',FLOW(LH2CO),FLOW(LO),FLOW(LOH),
     2 FLOW(LHO2)*2.,FLOW(LH2O2)*2.,FLOW(LO3)*3.,FLOW(LCO),
     3 FLOW(LCH3OOH)*2.,FLOW(LCH3O2)*2.,FLOW(LN2O)
      write(90,149) ' FUP',FUP(LH2CO),FUP(LO),FUP(LOH),
     2 FUP(LHO2)*2.,FUP(LH2O2)*2.,FUP(LO3)*3.,FUP(LCO),
     3 FUP(LCH3OOH)*2.,FUP(LCH3O2)*2.,FUP(LN2O)
      write(90,149) 'TLOS',TLOSS(LH2CO),TLOSS(LO),TLOSS(LOH),
     2 TLOSS(LHO2)*2.,TLOSS(LH2O2)*2.,TLOSS(LO3)*3.,TLOSS(LCO),
     3 TLOSS(LCH3OOH)*2.,TLOSS(LCH3O2)*2.,TLOSS(LN2O)
      write(90,*)

      write(90,110)  ISPEC(LNO),ISPEC(LNO2),ISPEC(LHNO2),ISPEC(LHNO3),
     2 ISPEC(LHO2NO2),ISPEC(LNO3),ISPEC(LN2O5),ISPEC(LO2),ISPEC(LSO),
     3 ISPEC(LSO2)
      write(90,149) 'RAIN',SR(LNO),SR(LNO2)*2.,SR(LHNO2)*2.,
     2 SR(LHNO3)*3.,SR(LHO2NO2)*4.,SR(LNO3)*3.,SR(LN2O5)*5.,
     3 SR(LO2)*2.,SR(LSO),SR(LSO2)*2.
      write(90,149) 'FLOW',FLOW(LNO),FLOW(LNO2)*2.,FLOW(LHNO2)*2.,
     2 FLOW(LHNO3)*3.,FLOW(LHO2NO2)*4.,FLOW(LNO3)*3.,FLOW(LN2O5)*5.
     3 ,FLOW(LO2)*2.,FLOW(LSO),FLOW(LSO2)*2.
      write(90,149) ' FUP',FUP(LNO),FUP(LNO2)*2.,FUP(LHNO2)*2.,
     2 FUP(LHNO3)*3.,FUP(LHO2NO2)*4.,FUP(LNO3)*3.,FUP(LN2O5)*5.,
     3 FUP(LO2)*2.,FUP(LSO),FUP(LSO2)*2.
      write(90,149) 'TLOS',TLOSS(LNO),TLOSS(LNO2)*2.,TLOSS(LHNO2)*2.,
     2 TLOSS(LHNO3)*3.,TLOSS(LHO2NO2)*4.,TLOSS(LNO3)*3.,
     3 TLOSS(LN2O5)*5.,TLOSS(LO2)*2.,TLOSS(LSO),TLOSS(LSO2)*2.
      write(90,*)

      write(90,110)  ISPEC(LH2SO4),ISPEC(LHSO),
     2 ISPEC(LCO2),ISPEC(LSO4AER)
      write(90,149) 'RAIN',SR(LH2SO4)*4.,SR(LHSO),
     2 SR(LCO2)*2.,SR(LSO4AER)*4
      write(90,149) 'FLOW',FLOW(LH2SO4)*4.,FLOW(LHSO),
     2 FLOW(LCO2)*2.,FLOW(LSO4AER)*4.
      write(90,149) ' FUP',FUP(LH2SO4)*4.,
     2 FUP(LHSO),FUP(LCO2)*2.,FUP(LSO4AER)*4.
      write(90,149) 'TLOS',TLOSS(LH2SO4)*4.,TLOSS(LHSO),
     2 TLOSS(LCO2)*2.,TLOSS(LSO4AER)*4.

 149  FORMAT(A5,5X,1P12E9.2)


C-AP
      write(90,179)
  179 FORMAT(/1X,'INTEGRATED REACTION RATES'/)
C      write(90,180) RAT
C  180 FORMAT(1X,1P10E10.3)
C
      WRITE(90,181)
      IROW = 10
      LR = NR/IROW + 1
      RL = FLOAT(NR)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 17 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) THEN
        K2 = NR
        WRITE(90,186) K1,(RAT(K),K=K1,K2),K2
  186   FORMAT(I3,2X,1P10E10.3,2X,I3)
        GO TO 17
      ENDIF
      WRITE(90,180) K1,(RAT(K),K=K1,K2),K2
  180 FORMAT(I3,2X,1P10E10.3,2X,I3)
   17 CONTINUE
      WRITE(90,181)
  181 FORMAT(9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,
     2    '8',9X,'9',8X,'10')
C
      write(90,160)
 160  FORMAT(/1X,'ATMOSPHERIC PARAMETERS AND PH EQ SPECIES')
      NPE = NSP - NQ
      NQ1 = NQ + 1
      LR = NPE/IROW + 1
      RL = FLOAT(NPE)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 12 L=1,LR
      K1 = NQ1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NSP
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 24 I=1,NZ,ISKIP
  24  write(90,120) Z(I),(SL(K,I),K=K1,K2)
  12  CONTINUE
C
C-JL ADD PRINTOUT FOR PRESSURE 
      write(90,190)
 190  FORMAT(/1X,'ATMOSPHERIC PARAMETERS')
      write(90,195)
 195  FORMAT(/4X,'Z',9X,'P'9X,'T',9X,'EDD',7X,'DEN')
      write(90,200)(Z(I),P(I),T(I),EDD(I),DEN(I),I=1,NZ,ISKIP)
 200  FORMAT(1X,1P5E10.3)
C
      write(90,230)
 230  FORMAT(/1X,'SULFATE AEROSOL PARAMETERS')
      write(90,235)
 235  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'FSULF',4X,
     2  'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'H2SO4S',4X,'H2SO4',5X,
     3  'CONSO4',4X,'CONVER')
      write(90,240) (Z(I),AERSOL(I,1),RPAR(I,1),
     & WFALL(I,1),FSULF(I),
     2  TAUSED(I,1),TAUEDD(I),TAUC(I,1),H2SO4S(I),USOL(LH2SO4,I),
     3  CONSO4(I),CONVER(I,1),I=1,NZ,ISKIP)
 240  FORMAT(1X,1P12E10.3)


C
C   Calculate total hydrogen mixing ratio and write to a file
      write(86,122)
 122  format(5x,'zkm',6x,'totH2O',4x,'totH',6x,'totH2',4x,'totCH4',
     2  4x,'tothyd')
      DO I=1,NZ
      zkm = Z(I)/1.E5
      totH2O = USOL(LH2O,I)
      totH = USOL(LH,I)*0.5
      totH2 = USOL(LH2,I)
      totCH4 = USOL(LCH4,I)*2.
      tothyd = totH2O + totH + totH2 + totCH4
      if(I.eq.JTROP) tothytrop = tothyd
      if(I.eq.90) tothyhom = tothyd
      write(86,121) zkm,totH2O,totH,totH2,totCH4,tothyd
      END DO
      close(86)  
 121  FORMAT(1x,1p6e10.3)
C
C Calculate hydrogen escape (in units of H2) based on tothyd at the top
C of the model. Do this both approximately and accurately.
C Approximate
      H2escape1 = 2.5e13*tothyhom
C Accurate
      H2escape2 = FUP(LH2O) + 0.5*FUP(LH) + FUP(LH2) + 2.*FUP(LCH4)
C
      write(90,300)
 300  format(/'Hydrogen escape rate in units of H2 molecules')
      write(90,301) H2escape1
 301  format(5x,'Approximate escape rate =',1pe10.3)
      write(90,302) H2escape2
 302  format(5x,'Escape rate calculated from FUP =',1pe10.3)
      write(90,303) tothyhom
 303  format(5x,'Total hydrogen mixing ratio at the homopause = ',
     2  1pe10.3)
      write(90,304) tothytrop
 304  format(5x,'Total hydrogen mixing ratio at the tropopause =',
     2  1pe10.3)

c in may 9 2019, JL add reaction rate table to the code to make life easier 

      if (N.EQ.NSTEPS) THEN
      REWIND 61
      READ(61,2000)CHEM
2000  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      DO 702 I=1,NSP
         ISP = ISPEC(I)
         WRITE(15,703) ISP,TP(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TP = ',1PE9.2)
       DO 704 Nj=1,NR 
          IF(JCHEM(3,Nj).EQ.I .OR. JCHEM(4,Nj).EQ.I .OR. 
     2       JCHEM(5,Nj).EQ.I)THEN
           IF(RAT(Nj).NE.0.) WRITE(15,705) Nj,(CHEM(J,Nj),J=1,5),RAT(Nj)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
      WRITE(15,706) ISP,TL(I)
 706     FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TL = ',1PE9.2)
       DO 707 Nj=1,NR 
          IF(JCHEM(1,Nj).EQ.I .OR. JCHEM(2,Nj).EQ.I)THEN
          IF(RAT(Nj).NE.0.) WRITE(15,705) Nj,(CHEM(J,Nj),J=1,5),RAT(Nj)
          ENDIF
 707   CONTINUE
 702  CONTINUE
      close(61)
      END IF

      RETURN
      END
