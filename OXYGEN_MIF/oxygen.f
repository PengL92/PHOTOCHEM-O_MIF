C NEED TO CHANGE OUTPUT FORMAT AT THE BOTTOM OF THIS SUBROUTINE C-PL
C******************************************************
      SUBROUTINE OXYGEN(USOLI,FLOW,OTYPE,DENCOR)
      
!      INCLUDE 'parameters_add_CO2.txt'

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/parNEQ_LDA.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
      INCLUDE '../INCLUDECHEM/parNF.inc'

      INCLUDE '../INCLUDECHEM/comABLOK.inc' !NEW:JTROP DEL:PMOD
      INCLUDE '../INCLUDECHEM/comBBLOK.inc' 
      INCLUDE '../INCLUDECHEM/comDBLOK.inc' !NEW:DEQ
      INCLUDE '../INCLUDECHEM/comFBLOK1.inc'
      INCLUDE '../INCLUDECHEM/comGBLOK.inc' !NEW:H(NQ)
      INCLUDE '../INCLUDECHEM/comRBLOK.inc' !A CONVERT TO AR
      INCLUDE '../INCLUDECHEM/comSULBLK.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'
      INCLUDE '../INCLUDECHEM/isotope_name.inc'
      INTEGER OTYPE 
      PARAMETER(NZ1=NZ-1, NQP = 40)!JL:NSLS = # TO SHORT LIVED SPECIES
c      PARAMETER(NEQI=NQI*NZ, LDAI=3*NQI+1, NMAXI=70)
c      PARAMETER(NRI=353, NSPI=78, NSPI1=NSPI+1, NSPI2=NSPI+2)
C     NQI - number of the long lived isotope species
C     NSLS - number of the long lived and short lived species
C     NRI - number of reactions with isotopes
C     NSPI - number of chemical species
C     NSPI1 = NSPI + 1 (INCLUDES HV)
C     NSPI2 = NSPI + 2 (INCLUDES "M") 
C-APN
      DIMENSION FVALI(NQI,NZ),FVI(NQI,NZ),DJACI(LDAI,NEQI),
     2  RHSI(NEQI),IPVT(NEQI),      FLOW_D(10000),USAVEOQ(NZ),
     3  USAVEI(NQI,NZ),RI(NZ),UI(NQI),RELI(NQI,NZ),FLOW(NQT)
!    , DELTAS4(NZ),  DDD
!     4  DELTAS3(NZ),DELTAS8(NZ)
C-APN
      DIMENSION USOLI(NQI,NZ),
     2 FLUXI(NQI,NZ),FLUXCHI(NQI,NZ),ZF(NZ)
      COMMON/RATESIS/RAINGCI(NQI,NZ),DDI(NQI,NZ),DKI(NQI,NZ),
     2  DUI(NQI,NZ),DLI(NQI,NZ),DHUI(NQI,NZ),DHLI(NQI,NZ),HIZ(NQI,NZ),
     3  VDEPI(NQI),RATI(NRI)
      COMMON/CBLOK/O3(NZ),H2O(NZ),O2(NZ),CO2(NZ),DZ,ZTROP
C
      DATA LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,
     4  LOQ,LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2/ 
     1  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     2  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     3  44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
     4  64,65,66,67,68,69,70,71,72,73,74,75,76,77,78/

      COMMON/NIBLOK/LH2CQ,LQ,LH2Q,LQH,LHOQ,LH2OQ,LO2Q,LOQO,LCQ,LCH3OQH,
     2  LCH3QOH,LCH3OQ,LCH3QO,LN2Q,LNQ,LNOQ,LHNOQ,LHNO2Q,
     3  LHOQNO2,LHO2NOQ,LNO2Q,LN2O4Q,LN2QO4,LOQ,
     4  LSQ,LSOQ,LH2SO3Q,LHSQ,LCOQ,LHQO,LHQONO2,LQ1D, 
     5  LH3CQ,LHCQ,LSO1Q,LSO3Q,LHSO2Q,LSO2Q,LIH2CO,LIO,LIH2O,LIOH,
     6  LIHO2,LIH2O2,LIO3,LIH,LIH2,LICH4,LICO,LICH3OOH,LICH3O2,LIN2O,
     7  LINO,LINO2,LIHNO3,LIHO2NO2,LINO3,LIN2O5,LIO2,LIH2S,LIHS,LISO,
     8  LISO2,LIHSO,LICO2,LIO1D,LICH21,LICH23,LICH3,LIH3CO,LIHCO,LIN,
     9  LIS,LISO21,LISO23,LIHSO3,LISO3,LIN2
      DIMENSION SLI(NSPI,NZ),SGFLUXI(NQI),
     2  SMFLUXI(NQI),VEFFI(NQI),MBOUNDI(NQI),
     3  DIZSAVE(NSPI2,NZ),PRODSIS3(NZ)
      integer NK,CT,qq1,qq2,qq3
      REAL SO4LOSS,DENCOR
      REAL S_cond_sink
C 
      OPEN(UNIT=90,FILE='OXYGEN_MIF/Reac_List.Oxy')
C     LONG-LIVED SPECIES
       ISPECI(1) = 4HH2CQ  !H2CO
       ISPECI(2) = 1HQ     !O
       ISPECI(3) = 3HH2Q   !H2O
       ISPECI(4) = 2HQH    !OH
       ISPECI(5) = 3HHOQ   !HO2A
       ISPECI(6) = 4HH2OQ  !H2O2
       ISPECI(7) = 3HO2Q   !O3A
       ISPECI(8) = 3HOQO   !O3B
       ISPECI(9) = 2HCQ    !CO
       ISPECI(10) = 6HCH3OQH !CH3OOHA NOTE6
       ISPECI(11) = 6HCH3QOH !CH3OOHB
       ISPECI(12) = 5HCH3OQ  !CH3O2A
       ISPECI(13) = 5HCH3QO  !CH3O2B  NOTE6
       ISPECI(14) = 3HN2Q    !N2O
       ISPECI(15) = 2HNQ     !NO
       ISPECI(16) = 3HNOQ    !NO2 
       ISPECI(17) = 4HHNOQ   !HNO2
       ISPECI(18) = 5HHNO2Q  !HNO3
       ISPECI(19) = 6HHOQNO2 !HO2NO2A
       ISPECI(20) = 6HHO2NOQ !HO2NO2B
       ISPECI(21) = 4HNO2Q   !NO3
       ISPECI(22) = 5HN2O4Q  !N2O5A NOTE12
       ISPECI(23) = 5HN2QO4  !N2O5B NOTE: HERE Q IS IN BETWEEN TWO N ATOMS
       ISPECI(24) = 2HOQ     !O2 
       ISPECI(25) = 2HSQ     !SO
       ISPECI(26) = 3HSOQ    !SO2
       ISPECI(27) = 6HH2SO3Q !H2SO4
       ISPECI(28) = 3HHSQ    !HSO 
       ISPECI(29) = 3HCOQ    !CO2
C  MARCH 1 2020, JL ADD TWO NEW SPECIES 
       ISPECI(30) = 3HHQO    ! HO2B
       ISPECI(31) = 6HHQONO2 ! HO2NO2C

C     SHORT-LIVED SPECIES     
       ISPECI(32) = 3HQ1D   !O1D  
       ISPECI(33) = 4HH3CQ  !H3CO
       ISPECI(34) = 3HHCQ   !HCO
       ISPECI(35) = 4HSO1Q  !SO21
       ISPECI(36) = 4HSO3Q  !SO23
       ISPECI(37) = 5HHSO2Q !HSO3
       ISPECI(38) = 4HSO2Q  !SO3 

C     INERT SPECIES
       ISPECI(39) = 4HH2CO
       ISPECI(40) = 1HO
       ISPECI(41) = 3HH2O
       ISPECI(42) = 2HOH
       ISPECI(43) = 3HHO2
       ISPECI(44) = 4HH2O2
       ISPECI(45) = 2HO3
       ISPECI(46) = 1HH
       ISPECI(47) = 2HH2
       ISPECI(48) = 3HCH4
       ISPECI(49) = 2HCO
       ISPECI(50) = 6HCH3OOH
       ISPECI(51) = 5HCH3O2
       ISPECI(52) = 3HN2O
       ISPECI(53) = 2HNO
       ISPECI(54) = 3HNO2
       ISPECI(55) = 4HHNO3
       ISPECI(56) = 6HHO2NO2
       ISPECI(57) = 3HNO3
       ISPECI(58) = 4HN2O5
       ISPECI(59) = 2HO2      
       ISPECI(60) = 3HH2S
       ISPECI(61) = 2HHS
       ISPECI(62) = 2HSO
       ISPECI(63) = 3HSO2  !SKIP H2SO4
       ISPECI(64) = 3HHSO  !
       ISPECI(65) = 3HCO2
       ISPECI(66) = 3HO1D  
       ISPECI(67) = 4HCH21
       ISPECI(68) = 4HCH23
       ISPECI(69) = 3HCH3
       ISPECI(70) = 4HH3CO
       ISPECI(71) = 3HHCO
       ISPECI(72) = 1HN
       ISPECI(73) = 1HS
       ISPECI(74) = 4HSO21
       ISPECI(75) = 4HSO23
       ISPECI(76) = 4HHSO3
       ISPECI(77) = 3HSO3 !SHIP S2
       ISPECI(78) = 2HN2
       ISPECI(79) = 2HHV
       ISPECI(80) = 1HM
C
      DO J=1,NZ
       DO I=1,NSPI
        SLI(I,J) = 0. ! numbers for all species
       ENDDO
      ENDDO

      DO J=1,NZ
       DO I=1,NSPI2
        DIZ(I,J) = 0.
        DIZSAVE(I,J) = 0. !initial numbers for HV and M
       ENDDO
      ENDDO
c********************************
C-JL INITIATIZE 17O AND 18O FROM THEIR NATURAL ABUNDANCE 
C FROM FEILBER ET AL. 2013 CPL
      IF(OTYPE.EQ.1) THEN ! ONE ISOTOPE 
      DENCOR = 1.0
      ELSE IF(OTYPE.EQ.2) THEN ! 17O
      DENCOR =0.00076915! 4.000E-4
      ELSE IF(OTYPE.EQ.3) THEN !18O
      DENCOR =0.00410467! 2.000E-3
      
      END IF 

C********** INITIALIZATION OF ISOTOPE SPECIES*****************

      DO J=1,NZ
        SLI(LH2CQ,J) = SL(LH2CO,J)*DENCOR !number density
        SLI(LQ,J) = SL(LO,J)*DENCOR           
        SLI(LH2Q,J) = SL(LH2O,J)*DENCOR           
        SLI(LQH,J) = SL(LOH,J)*DENCOR           
        SLI(LHOQ,J) = SL(LHO2,J)*DENCOR!.JL ADD HQO SO REMOVE *2  
        SLI(LH2OQ,J) = SL(LH2O2,J)*DENCOR*2.      
        SLI(LO2Q,J) = SL(LO3,J)*DENCOR*2. ! QO2=O2Q 
        SLI(LOQO,J) = SL(LO3,J)*1.*DENCOR 
        SLI(LCQ,J) = SL(LCO,J)*DENCOR      
        SLI(LCH3OQH,J) = SL(LCH3OOH,J)*1.*DENCOR       
        SLI(LCH3QOH,J) = SL(LCH3OOH,J)*1.*DENCOR     
        SLI(LCH3OQ,J) = SL(LCH3O2,J)*1.*DENCOR          
        SLI(LCH3QO,J) = SL(LCH3O2,J)*1.*DENCOR
        SLI(LN2Q,J) = SL(LN2O,J)*DENCOR
        SLI(LNQ,J) = SL(LNO,J)*DENCOR
        SLI(LNOQ,J) = SL(LNO2,J)*DENCOR*2. 
        SLI(LHNOQ,J) = SL(LHNO2,J)*DENCOR*2. 
        SLI(LHNO2Q,J) = SL(LHNO3,J)*DENCOR*3. 
        SLI(LHOQNO2,J) = SL(LHO2NO2,J)*DENCOR !JL ADD HQONO2 SO REMOVE *2.
        SLI(LHO2NOQ,J) = SL(LHO2NO2,J)*DENCOR*2. !2true
        SLI(LNO2Q,J) = SL(LNO3,J)*DENCOR*3.   !3 true
        SLI(LN2O4Q,J) = SL(LN2O5,J)*DENCOR*4. !4 true
        SLI(LN2QO4,J) = SL(LN2O5,J)*1.*DENCOR
        SLI(LOQ,J) = SL(LO2,J)*DENCOR*2.!*(1.-.01) 
        SLI(LSQ,J) = SL(LSO,J)*DENCOR
        SLI(LSOQ,J) = SL(LSO2,J)*DENCOR*2. !2true 
        SLI(LH2SO3Q,J) = SL(LH2SO4,J)*DENCOR*4. !4true 
        SLI(LHSQ,J) = SL(LHSO,J)*DENCOR
        SLI(LCOQ,J) = SL(LCO2,J)*DENCOR*2.!*(1.+.08) !2 true
C MARCH 1, 2020 JL
        SLI(LHQO,J) = SL(LHO2,J)*DENCOR
        SLI(LHQONO2,J) =SL(LHO2NO2,J)*DENCOR
C
C JL  HERE IS THE INERT SPECIES 
C**************************************************          
        SLI(LIH2CO,J) = SL(LH2CO,J)           
        SLI(LIO,J) = SL(LO,J)
        SLI(LIH2O,J) = SL(LH2O,J)           
        SLI(LIOH,J) = SL(LOH,J)
        SLI(LIHO2,J) = SL(LHO2,J)           
        SLI(LIH2O2,J) = SL(LH2O2,J)
        SLI(LIO3,J) = SL(LO3,J)         
        SLI(LIH,J) = SL(LH,J)
        SLI(LIH2,J) = SL(LH2,J)          
        SLI(LICH4,J) = SL(LCH4,J)
        SLI(LICO,J) = SL(LCO,J)        
        SLI(LICH3OOH,J) = SL(LCH3OOH,J)
        SLI(LICH3O2,J) = SL(LCH3O2,J)           
        SLI(LIN2O,J) = SL(LN2O,J)
        SLI(LINO,J) = SL(LNO,J)           
        SLI(LINO2,J) = SL(LNO2,J)
        SLI(LIHNO3,J) = SL(LHNO3,J)           
        SLI(LIHO2NO2,J) = SL(LHO2NO2,J)
        SLI(LINO3,J) = SL(LNO3,J)           
        SLI(LIN2O5,J) = SL(LN2O5,J)
        SLI(LIO2,J) = SL(LO2,J)      
        SLI(LIH2S,J) = SL(LH2S,J)
        SLI(LIHS,J) = SL(LHS,J)          
        SLI(LISO,J) = SL(LSO,J)
        SLI(LISO2,J) = SL(LSO2,J) 
        SLI(LIHSO,J) = SL(LHSO,J)
        SLI(LICO2,J) = SL(LCO2,J) 
        SLI(LIO1D,J) = SL(LO1D,J)
        SLI(LICH21,J)= SL(LCH21,J)
        SLI(LICH23,J)= SL(LCH23,J)          
        SLI(LICH3,J) = SL(LCH3,J)
        SLI(LIH3CO,J) = SL(LH3CO,J)          
        SLI(LIHCO,J) = SL(LHCO,J)
        SLI(LIN,J) = SL(LN,J)           
        SLI(LIS,J) = SL(LS,J)
        SLI(LISO21,J) = SL(LSO21,J)           
        SLI(LISO23,J) = SL(LSO23,J)
        SLI(LIHSO3,J) = SL(LHSO3,J)           
        SLI(LISO3,J) = SL(LSO3,J)
        SLI(LIN2,J) = SL(LN2,J)           
      ENDDO

C Change OQ at the surface
!       SLI(LOQ,1) = SLI(LOQ,1)*1.01
      !print 1811, SLI(I,1),I=37,75)
 1811 format(6(1P5E9.2,/))
! 1811 FORMAT(1X,E9.2,1X,I2,1X,E9.2,1X,I2,1X,E9.2)

C******Initialization of the boundary conditions *******
C*********Lower Boundary*******************************
C*****Ground Fluxes********************************************
C
         SGFLUXI(LCQ) = FLOW(LCO)*DENCOR     
         SGFLUXI(LN2Q) = FLOW(LN2O)*DENCOR

C******Deposition velocities ****************************
        VDEPI(LH2CQ) = VDEP(LH2CO)
        VDEPI(LQ) = VDEP(LO)            
        VDEPI(LH2Q) = VDEP(LH2O)            
        VDEPI(LQH) = VDEP(LOH)            
        VDEPI(LHOQ) = VDEP(LHO2)    
        VDEPI(LH2OQ) = VDEP(LH2O2)           
        VDEPI(LO2Q) = VDEP(LO3) 
        VDEPI(LOQO) = VDEP(LO3) 
        VDEPI(LCQ) = VDEP(LCO)       
        VDEPI(LCH3OQH) = VDEP(LCH3OOH)        
        VDEPI(LCH3QOH) = VDEP(LCH3OOH)    
        VDEPI(LCH3OQ) = VDEP(LCH3O2)            
        VDEPI(LCH3QO) = VDEP(LCH3O2) 
        VDEPI(LN2Q) = VDEP(LN2O) 
        VDEPI(LNQ) = VDEP(LNO) 
        VDEPI(LNOQ) = VDEP(LNO2) 
        VDEPI(LHNOQ) = VDEP(LHNO2) 
        VDEPI(LHNO2Q) = VDEP(LHNO3) 
        VDEPI(LHOQNO2) = VDEP(LHO2NO2) 
        VDEPI(LHO2NOQ) = VDEP(LHO2NO2) 
        VDEPI(LNO2Q) = VDEP(LNO3) 
        VDEPI(LN2O4Q) = VDEP(LN2O5) 
        VDEPI(LN2QO4) = VDEP(LN2O5) 
        VDEPI(LOQ) = ABS (FLOW(LO2)/SL(LO2,1))!*7.3 !VDEP(LO2) ! (GPP/N(o2))1.1e14/5.1E18
        VDEPI(LSQ) = VDEP(LSO) 
        VDEPI(LSOQ) = VDEP(LSO2) 
        VDEPI(LH2SO3Q) = VDEP(LH2SO4) 
        VDEPI(LHSQ) = VDEP(LHSO) 
        VDEPI(LCOQ) = ABS (FLOW(LCO2)/SL(LCO2,1)) !VDEP(LCO2)
c        print *,'sl(lo2,1) =',sl(lo2,1)
c        print *,'sl(lco2,1) =',sl(lco2,1)
C MARCH 1,2020 JL 
        VDEPI(LHQO) = VDEP(LHO2) 
        VDEPI(LHQONO2) = VDEP(LHO2NO2)

C******UPPER BOUNDARY***********************************
C********UPWARD FLUX************************************
        SMFLUXI(LH2CQ) = SMFLUX(LH2CO)*DENCOR 
        SMFLUXI(LQ) =  SMFLUX(LO)*DENCOR !4.*FUP(LO2)           
        SMFLUXI(LH2Q) = SMFLUX(LH2O)*DENCOR           
        SMFLUXI(LQH) = SMFLUX(LOH)*DENCOR           
        SMFLUXI(LHOQ) = SMFLUX(LHO2)*DENCOR   
        SMFLUXI(LH2OQ) = SMFLUX(LH2O2)*DENCOR          
        SMFLUXI(LO2Q) = SMFLUX(LO3)*DENCOR
        SMFLUXI(LOQO) = SMFLUX(LO3)*DENCOR
        SMFLUXI(LCQ) = SMFLUX(LCO)*DENCOR      
        SMFLUXI(LCH3OQH) = SMFLUX(LCH3OOH)*DENCOR       
        SMFLUXI(LCH3QOH) = SMFLUX(LCH3OOH)*DENCOR   
        SMFLUXI(LCH3OQ) = SMFLUX(LCH3O2)*DENCOR           
        SMFLUXI(LCH3QO) = SMFLUX(LCH3O2)*DENCOR
        SMFLUXI(LN2Q) = SMFLUX(LN2O)*DENCOR
        SMFLUXI(LNQ) = SMFLUX(LNO)*DENCOR
        SMFLUXI(LNOQ) = SMFLUX(LNO2)*DENCOR
        SMFLUXI(LHNOQ) = SMFLUX(LHNO2)*DENCOR
        SMFLUXI(LHNO2Q) = SMFLUX(LHNO3)*DENCOR
        SMFLUXI(LHOQNO2) = SMFLUX(LHO2NO2)*DENCOR
        SMFLUXI(LHO2NOQ) = SMFLUX(LHO2NO2)*DENCOR
        SMFLUXI(LNO2Q) = SMFLUX(LNO3)*DENCOR
        SMFLUXI(LN2O4Q) = SMFLUX(LN2O5)*DENCOR
        SMFLUXI(LN2QO4) = SMFLUX(LN2O5)*DENCOR
        SMFLUXI(LOQ) = SMFLUX(LO2)*DENCOR
        SMFLUXI(LSQ) = SMFLUX(LSO)*DENCOR
        SMFLUXI(LSOQ) = SMFLUX(LSO2)*DENCOR
        SMFLUXI(LH2SO3Q) = SMFLUX(LH2SO4)*DENCOR
        SMFLUXI(LHSQ) = SMFLUX(LHSO)*DENCOR
        SMFLUXI(LCOQ) = SMFLUX(LCO2)*DENCOR
C MARCH 1 2020 JL 
        SMFLUXI(LHQO) = SMFLUX(LHO2) *DENCOR
        SMFLUXI(LHQONO2) = SMFLUX(LHO2NO2) *DENCOR

C********Effusion velocity******************************
        VEFFI(LH2CQ) = VEFF(LH2CO) 
        VEFFI(LQ) = VEFF(LO)           
        VEFFI(LH2Q) = VEFF(LH2O)           
        VEFFI(LQH) = VEFF(LOH)           
        VEFFI(LHOQ) = VEFF(LHO2)   
        VEFFI(LH2OQ) = VEFF(LH2O2)          
        VEFFI(LO2Q) = VEFF(LO3)
        VEFFI(LOQO) = VEFF(LO3)
        VEFFI(LCQ) = VEFF(LCO)      
        VEFFI(LCH3OQH) = VEFF(LCH3OOH)       
        VEFFI(LCH3QOH) = VEFF(LCH3OOH)   
        VEFFI(LCH3OQ) = VEFF(LCH3O2)           
        VEFFI(LCH3QO) = VEFF(LCH3O2)
        VEFFI(LN2Q) = VEFF(LN2O)
        VEFFI(LNQ) = VEFF(LNO)
        VEFFI(LNOQ) = VEFF(LNO2)
        VEFFI(LHNOQ) = VEFF(LHNO2)
        VEFFI(LHNO2Q) = VEFF(LHNO3)
        VEFFI(LHOQNO2) = VEFF(LHO2NO2)
        VEFFI(LHO2NOQ) = VEFF(LHO2NO2)
        VEFFI(LNO2Q) = VEFF(LNO3)
        VEFFI(LN2O4Q) = VEFF(LN2O5)
        VEFFI(LN2QO4) = VEFF(LN2O5)
        VEFFI(LOQ) = VEFF(LO2)
        VEFFI(LSQ) = VEFF(LSO)
        VEFFI(LSOQ) = VEFF(LSO2)
        VEFFI(LH2SO3Q) = VEFF(LH2SO4)
        VEFFI(LHSQ) = VEFF(LHSO)
        VEFFI(LCOQ) = VEFF(LCO2)
        
C MARCH 1,2020  JL 
        VEFFI(LHQO) = VEFF(LHO2) 
        VEFFI(LHQONO2) = VEFF(LHO2NO2) 
 
C********* Specify what type of the boundary conditions to use**
C*******Lower boundary*******************************
        LBOUNDI(LH2CQ) = LBOUND(LH2CO) 
        LBOUNDI(LQ) = LBOUND(LO)           
        LBOUNDI(LH2Q) = LBOUND(LH2O)           
        LBOUNDI(LQH) = LBOUND(LOH)           
        LBOUNDI(LHOQ) = LBOUND(LHO2)   
        LBOUNDI(LH2OQ) = LBOUND(LH2O2)          
        LBOUNDI(LO2Q) = LBOUND(LO3)
        LBOUNDI(LOQO) = LBOUND(LO3)
        LBOUNDI(LCQ) =  2!LBOUND(LCO) !2      
        LBOUNDI(LCH3OQH) = LBOUND(LCH3OOH)       
        LBOUNDI(LCH3QOH) = LBOUND(LCH3OOH)   
        LBOUNDI(LCH3OQ) = LBOUND(LCH3O2)           
        LBOUNDI(LCH3QO) = LBOUND(LCH3O2)
        LBOUNDI(LN2Q) = 2!LBOUND(LN2O) !2
        LBOUNDI(LNQ) = LBOUND(LNO)
        LBOUNDI(LNOQ) = LBOUND(LNO2)
        LBOUNDI(LHNOQ) = LBOUND(LHNO2)
        LBOUNDI(LHNO2Q) = LBOUND(LHNO3)
        LBOUNDI(LHOQNO2) = LBOUND(LHO2NO2)
        LBOUNDI(LHO2NOQ) = LBOUND(LHO2NO2)
        LBOUNDI(LNO2Q) = LBOUND(LNO3)
        LBOUNDI(LN2O4Q) = LBOUND(LN2O5)
        LBOUNDI(LN2QO4) = LBOUND(LN2O5)
        LBOUNDI(LOQ) = 0!LBOUND(LO2) !0
        LBOUNDI(LSQ) = LBOUND(LSO)
        LBOUNDI(LSOQ) = LBOUND(LSO2)
        LBOUNDI(LH2SO3Q) = LBOUND(LH2SO4)
        LBOUNDI(LHSQ) = LBOUND(LHSO)
        LBOUNDI(LCOQ) = 0!LBOUND(LCO2)!0

C MARCH 1, 2020 JL Changed by Jim (6/23)--Jingjun had 1's
        LBOUNDI(LHQO) = LBOUND(LHO2) 
        LBOUNDI(LHQONO2) = LBOUND(LHO2NO2) 

C********Upper boundary******************************
        MBOUNDI(LH2CQ) = MBOUND(LH2CO) 
        MBOUNDI(LQ) = MBOUND(LO)           
        MBOUNDI(LH2Q) = MBOUND(LH2O)           
        MBOUNDI(LQH) = MBOUND(LOH)           
        MBOUNDI(LHOQ) = MBOUND(LHO2)   
        MBOUNDI(LH2OQ) = MBOUND(LH2O2)          
        MBOUNDI(LO2Q) = MBOUND(LO3)
        MBOUNDI(LOQO) = MBOUND(LO3)
        MBOUNDI(LCQ) = MBOUND(LCO)      
        MBOUNDI(LCH3OQH) = MBOUND(LCH3OOH)       
        MBOUNDI(LCH3QOH) = MBOUND(LCH3OOH)   
        MBOUNDI(LCH3OQ) = MBOUND(LCH3O2)           
        MBOUNDI(LCH3QO) = MBOUND(LCH3O2)
        MBOUNDI(LN2Q) = MBOUND(LN2O)
        MBOUNDI(LNQ) = MBOUND(LNO)
        MBOUNDI(LNOQ) = MBOUND(LNO2)
        MBOUNDI(LHNOQ) = MBOUND(LHNO2)
        MBOUNDI(LHNO2Q) = MBOUND(LHNO3)
        MBOUNDI(LHOQNO2) = MBOUND(LHO2NO2)
        MBOUNDI(LHO2NOQ) = MBOUND(LHO2NO2)
        MBOUNDI(LNO2Q) = MBOUND(LNO3)
        MBOUNDI(LN2O4Q) = MBOUND(LN2O5)
        MBOUNDI(LN2QO4) = MBOUND(LN2O5)
        MBOUNDI(LOQ) = MBOUND(LO2)
        MBOUNDI(LSQ) = MBOUND(LSO)
        MBOUNDI(LSOQ) = MBOUND(LSO2)
        MBOUNDI(LH2SO3Q) = MBOUND(LH2SO4)
        MBOUNDI(LHSQ) = MBOUND(LHSO)
        MBOUNDI(LCOQ) = MBOUND(LCO2)

C   MARCH 1,2020 JL Changed by Jim (6/23)--Jingjun had 1's
        MBOUNDI(LHQO) = MBOUND(LHO2) 
        MBOUNDI(LHQONO2) = MBOUND(LHO2NO2) 

C   0 = CONSTANT DEPOSITION VELOCITY (VDEP)
C   1 = CONSTANT MIXING RATIO
C   2 = CONSTANT UPWARD FLUX (SGFLUX)

C**** ZERO everything
      DO I=1,NSPI
       NUMPI(I) = 0
       NUMLI(I) = 0
      ENDDO

C ***** READ ISOTOPE CHEMISTRY DATA FILE *****
      READ(90,200)JCHEMI
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
c      PRINT 201,(J,(JCHEMI(M,J),M=1,5),J=1,NRI)
 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
      KJACI = LDAI*NEQI
c      PRINT 202,NQI,NZ,KJACI
 202  FORMAT(//1X,'NQI=',I2,5X,'NZ=',I3,5X,'KJACI=',I7)

C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEMI *****
      DO 5 J=1,NRI
      DO 5 M=1,5
      IF(JCHEMI(M,J).EQ.1H ) GO TO 5
      DO 6 I=1,NSPI2
      IF(JCHEMI(M,J).NE.ISPECI(I)) GO TO 6
      JCHEMI(M,J) = I
      GO TO 5
   6  CONTINUE
      IERRI = J
      GO TO 25
   5  CONTINUE

C Jim debug
c      PRINT 2201,(J,(JCHEMI(M,J),M=1,5),J=1,NRI)
 2201 FORMAT(1X,I3,1H),5X,I8,4H +  ,I8,7H  =    ,I8,4H +  ,I8,4X,I8)

C
C ***** FILL UP ISOTOPE CHEMICAL PRODUCTION AND LOSS MATRICES *****
      DO 7 M=1,2 ! FOR REACTANTS
      N = 3-M
      DO 7 J=1,NRI
      I = JCHEMI(M,J)
      IF(I.LT.1.OR.I.GT.NSPI) GO TO 7
      NUMLI(I) = NUMLI(I) + 1
        IF(NUMLI(I).GT.NMAXI) GO TO 20
      K = NUMLI(I)
      ILOSSI(1,I,K) = J          !1:STORE THE REACTIONS NUMBER
      ILOSSI(2,I,K) = JCHEMI(N,J)!2:STORE THE OTHER REACTANT 
   7  CONTINUE

C Jim debug
c      print *
c      print 3300
 3300 format('ILOSSI_1')
c      DO I=1,NSPI
c      PRINT 3301,I,(ILOSSI(1,I,K),K=1,35)
 3301 FORMAT('I=',I3,1X,35I4)
c      enddo
c      PRINT *
c     print 3302
 3302 format('ILOSSI_2')
c      do I=1,NSPI
c      PRINT 3301,I,(ILOSSI(2,I,K),K=1,35)
c      ENDDO

      DO 8 M=3,5 !FOR PRODUCTS
      DO 8 J=1,NRI
      I = JCHEMI(M,J)
      IF(I.LT.1.OR.I.GT.NSPI) GO TO 8
      NUMPI(I) = NUMPI(I) + 1
      IF(NUMPI(I).GT.NMAXI) GO TO 20
      K = NUMPI(I)
      IPRODI(I,K) = J
   8  CONTINUE
C
C Jim debug
c      print *
c      print 3303
 3303 format('IPRODI')
c      DO I=1,NSPI
c      PRINT 3301,I,(IPRODI(I,K),K=1,35)
c      enddo

       DO I=1,NSPI
        DO J=1,NZ
         DIZ(I,J) = SLI(I,J)
         DIZSAVE(I,J) = SLI(I,J)
        ENDDO
       ENDDO
       DO J=1,NZ
        DIZ(NSPI2-1,J) = 1.   ! HV
        DIZ(NSPI2,J) = DEN(J) ! M
        DIZSAVE(NSPI2-1,J) = 1.
        DIZSAVE(NSPI2,J) = DEN(J) 
       ENDDO

C Initialize isotopic long-lived species from number densities
      DO I=1,NQI
       DO J=1,NZ
        USOLI(I,J) = DIZ(I,J)/DEN(J) 
       ENDDO
      ENDDO


C*********************************
C
C-JL AVOID MUTIPLE PRINTOUT
      IF(OTYPE .EQ.3) GO TO 2020
C-PL PRINT OUT THE INITIAL MIXING RATIO FOR ALL THE SPECIES
      print*,''
      print 199,'INITIAL MIXING RATIOS OF ISOTOPIC SPECIES'
      PRINT*, ''
      print 206,'Alt(km)',(ISPECI(I),I=1,10) !197
      DO J=1,3
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=1,10) !198
      END DO      
      DO J=4,NZ,2
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=1,10) 
      ENDDO
      print 206,'Alt(km)',(ISPECI(I),I=1,10)

      print*,''
      print 206,'Alt(km)',(ISPECI(I),I=11,20)
      DO J=1,3
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=11,20)
      END DO   
      DO J=4,NZ,2
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=11,20)
      ENDDO
      print 206,'Alt(km)',(ISPECI(I),I=11,20)
      print*,''

      print 206,'Alt(km)',(ISPECI(I),I=21,29)
      DO J=1,3
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=21,29)
      END DO   
      DO J=4,NZ,2
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=21,29)
      ENDDO
      print 206,'Alt(km)',(ISPECI(I),I=21,29)

C JL DON'T FORGET TO PRINT OUT NEW LONG LIVED SPECIES 
      print*,''
      print 206,'Alt(km)',(ISPECI(I),I=30,31)
      DO J=1,3
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=30,31)
      END DO   
      DO J=4,NZ,2
      print 203,j-0.5,(DIZ(I,J)/DEN(J),I=30,31)
      ENDDO
      print 206,'Alt(km)',(ISPECI(I),I=30,31)

  203 format(F6.1,2x,1P10E9.2)
  206 format(A7,1x,10(1x,A8))!,2X
  199 format(A42)
      print *

C**********************************
 2020 CONTINUE 
      DZ = Z(2) - Z(1) 




      CALL RATESI(OTYPE)


C*************************************

      EPSJ = 1.E-7  ! Fractional amount by which to perturb things for the Jacobian

C   SET JACOBIAN PARAMETERS
      KD = 2*NQI + 1
      KU = KD - NQI
      KL = KD + NQI
C
      NN = 0
      NSTEPS = 50 !pengpeng

      qq1 = 0  
      qq2 = 0
      qq3 = 0
      qq11 = 0.
      qq22 = 0.
      qq33 = 0.
      CT  = 0  

C ***** START ITERATIVE LOOP *****
      DO 1 N=1,NSTEPS
      NN = NN + 1

C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,LDAI
      DO 17 K=1,NEQI
  17  DJACI(J,K) = 0.
      DO 19 K=1,NEQI
  19  RHSI(K) = 0.
      DO I=1,NQI
      DO J=1,NZ
       USOLI(I,J) = ABS(USOLI(I,J))
      ENDDO
      ENDDO
C
C     (DJACI IS EQUAL TO  - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS


      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL DOCHEMI(FVALI,IDO,USOLI,OTYPE,DENCOR)
C      print *, 'LO2 =',USOLI(LOQ,1)



      DO 9 I=1,NQI
      DO 9 J=1,NZ
      K = I + (J-1)*NQI
      RHSI(K) = FVALI(I,J)
   9  USAVEI(I,J) = USOLI(I,J)
C
      DO 3 I=1,NQI
      DO 11 J=1,NZ
      RI(J) = EPSJ * ABS(USOLI(I,J))
  11  USOLI(I,J) = USAVEI(I,J) + RI(J)


      CALL DOCHEMI(FVI,0,USOLI,OTYPE,DENCOR)     

      DO 12 M=1,NQI
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQI
  12  DJACI(MM,K) = (FVALI(M,J) - FVI(M,J))/RI(J)
C
      DO 10 J=1,NZ
  10  USOLI(I,J) = USAVEI(I,J)
   3  CONTINUE

C   COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO 13 I = 1,NQI
      DO 14 J=2,NZ1
      K = I + (J-1)*NQI
      RHSI(K) = RHSI(K) - DDI(I,J)*USOLI(I,J)
     2  + (DUI(I,J) + DHUI(I,J))*USOLI(I,J+1)
     3  + (DLI(I,J) - DHLI(I,J))*USOLI(I,J-1)
      DJACI(KD,K) = DJACI(KD,K) + DDI(I,J)
      DJACI(KU,K+NQI) = - DUI(I,J) - DHUI(I,J)
  14  DJACI(KL,K-NQI) = - DLI(I,J) + DHLI(I,J)
  13  CONTINUE

C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQI
      UI(K) = USOLI(K,1)
      LB = LBOUNDI(K)
C
C   CONSTANT DEPOSITION VELOCITY
      IF(LBOUNDI(K).NE.0) GO TO 16
      RHSI(K)=RHSI(K)+(DUI(K,1)+DHUI(K,1))*USOLI(K,2) 
     2  - DUI(K,1)*UI(K)
     3  - (VDEPI(K)/DZ - HIZ(K,1)/(2.*DZ))*UI(K)
      DJACI(KD,K) = DJACI(KD,K) + DUI(K,1) + VDEPI(K)/DZ
     2  - HIZ(K,1)/(2.*DZ)
      DJACI(KU,K+NQI) = - DUI(K,1) - DHUI(K,1)
      GO TO 15     
C
C   CONSTANT MIXING RATI3O
  16  IF(LBOUNDI(K).NE.1) GO TO 31
      RHSI(K) = 0.
      DO 36 M=1,NQI
      MM = KD + K - M
  36  DJACI(MM,M) = 0.
      DJACI(KU,K+NQI) = 0.
      DJACI(KD,K) = DUI(K,1)
      GO TO 15
C
C   CONSTANT UPWARD FLUX
  31  CONTINUE
      RHSI(K)=RHSI(K)+(DUI(K,1)+DHUI(K,1))*USOLI(K,2)-DUI(K,1)*UI(K)
     2   + HIZ(K,1)*UI(K)/(2.*DZ) + SGFLUXI(K)/DEN(1)/DZ
      DJACI(KD,K) = DJACI(KD,K) + DUI(K,1) + HIZ(K,1)/(2.*DZ)
      DJACI(KU,K+NQI) = - DUI(K,1) - DHUI(K,1)
  15  CONTINUE

C ***** UPPER BOUNDARY CONDITIONS *****
      DO 30 I=1,NQI
      UI(I) = USOLI(I,NZ)
      K = I + NZ1*NQI
      MB = MBOUNDI(I)
C
C   CONSTANT EFFUSION VELOCITY
      IF(MBOUNDI(I).NE.0) GO TO 29
      RHSI(K)=RHSI(K)+(DLI(I,NZ)-DHLI(I,NZ))*USOLI(I,NZ1)
     2  - DLI(I,NZ)*UI(I) - (VEFFI(I)/DZ + HIZ(I,NZ)/(2.*DZ))*UI(I)
      DJACI(KD,K) = DJACI(KD,K) + DLI(I,NZ) + VEFFI(I)/DZ
     2  + HIZ(I,NZ)/(2.*DZ)
      DJACI(KL,K-NQI) = - DLI(I,NZ) + DHLI(I,NZ)
      GO TO 30
C
C   CONSTANT DOWNWARD FLUX
  29  CONTINUE
      RHSI(K) = RHSI(K) + (DLI(I,NZ) - DHLI(I,NZ))*USOLI(I,NZ1)
     2  - DLI(I,NZ)*UI(I) - HIZ(I,NZ)*UI(I)/(2.*DZ)
     3  - SMFLUXI(I)/DEN(NZ)/DZ
      DJACI(KD,K)=DJACI(KD,K) + DLI(I,NZ) + HIZ(I,NZ)/(2.*DZ)
      DJACI(KL,K-NQI) = - DLI(I,NZ) + DHLI(I,NZ)
  30  CONTINUE

C-PL TO FIX MIXING RATIO OF H2O 06/18
      DO 33 J=1,NZ
      IF(Z(J).GT.ZTROP) GO TO 34
C      K = 3 + (J-1)*NQ

C-MARCH 6,2020, JL find a bug (Correct! JK)
      K = 3 + (J-1)*NQI
      RHSI(K) = 0.
  33  CONTINUE
  34  CONTINUE



C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****

      CALL SGBFA(DJACI,LDAI,NEQI,NQI,NQI,IPVT,INDEX)
      IF(INDEX.NE.0) PRINT 103,N,INDEX
 103  FORMAT(/1X,'N =',I3,5X,'INDEX =',I9)
      CALL SGBSL(DJACI,LDAI,NEQI,NQI,NQI,IPVT,RHSI,0)


C
C   COMPUTE NEW CONCENTRATIONS 
      EMAX = 0.
      DO 26 I=1,NQI
      DO 26 J=1,NZ
      K = I + (J-1)*NQI
C-AP
      RELI(I,J) = RHSI(K)/USOLI(I,J)
      EREL = ABS(RELI(I,J))
      EMAX = AMAX1(EMAX,EREL)
      IF(EREL.LT.EMAX) GO TO 26
      IS = I
      JS = J
      UMAX = USOLI(I,J)
      RMAX = RHSI(K)
  26  USOLI(I,J) = USOLI(I,J) + RHSI(K)

      ISP = ISPECI(IS)
      ZMAX = Z(JS)
      PRINT 100, N,EMAX,ISP,ZMAX,UMAX,RMAX
 100  FORMAT(1X,'N =',I3,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X)
      !PRINT 504, DEN(JS),SLI(IS,JS)
 504  FORMAT(9X,'DEN = ',1PE9.2,' DEN(I) = ',1PE9.2)

          
 !CHANGING TIME DOHERE
      GO TO 299
      IF (NN.EQ.1) THEN
       DO J=1,NZ
       USAVEOQ(J) = USOLI(LOQ,J)
       END DO
      END IF 

      IF (NN.GT.10.AND.MOD(NN,3).EQ.0.and.qsss.ne.100) THEN
      CT = CT + 1
C-JL ADD DENCOR HERE
      FLOW_D(CT)= (FLOWI(LOQ) - 2.*FLOW(LO2))!*DENCOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PRINT *, 'ERR',FLOW_D(CT)/(2.*FLOW(LO2))
        IF (ABS(FLOW_D(CT)/(2.*FLOW(LO2))).GT.1.E-5) THEN

         IF (FLOW_D(1).GT.0) THEN
           IF (FLOW_D(CT).GT.0..AND.qq22.LT.1.) THEN

           qq1 = qq1 + 1 
           qq11 = qq11 +1. 
           DO J=1,NZ                  
           USOLI(LOQ,J)= USAVEOQ(J) * (1.- 0.01*qq11)
           END DO

           ELSE IF (FLOW_D(qq1+1).LT.0..and.FLOW_D(qq1).GT.0.
     2     .AND.qs.ne.100  .AND.qq3.lt.1.) THEN
               IF (FLOW_D(CT).LT.FLOW_D(qq1).OR.qq2.eq.0) THEN 
           qq2  = qq2 + 1 
           qq22 = qq22+1.           
           DO J=1,NZ
           USOLI(LOQ,J)= USAVEOQ(J)*(1.- 0.01*(qq11-1.) + 0.001*qq22)
           END DO
           DMIN = MIN(FLOW_D(CT),FLOW_D(CT-1))
               ELSE
               qs = 100
               END IF

           ELSE IF (qs.eq.100..AND.qss.ne.100) THEN
                IF (FLOW_D(CT).LT.FLOW_D(CT-1).OR.qq3.eq.0) THEN !FLOW_D(qq2-1)
           qq3  = qq3  + 1
           qq33 = qq33 + 1.
           DO J=1,NZ
           USOLI(LOQ,J)= USAVEOQ(J) * 
     2             (1.- 0.01*(qq11-1.) + 0.001*(qq22-1.) + 0.0001*qq33)
           END DO
               ELSE
               qss = 100
        PRINT *, 'YOU GOT THE BEST TUNED USOLI(OQ)'              
               END IF
           END IF
         ELSE
          DO I=1,NZ
         USOLI(LOQ,J)=USOLI(LOQ,J)
          END DO
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        IF (FLOW_D(1).LT.0) THEN

!           USOLI(LOQ,J)= USOLI(LOQ,J) * (1.- 0.1)         
!        ELSE 
!           USOLI(LOQ,J)= USOLI(LOQ,J) * (1.+ 0.01)
        ELSE
        PRINT *, 'YOU GOT THE RIGHT USOLI(OQ)'
        qsss = 100
        END IF
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (NN.GT.10.AND.MOD(NN,3).EQ.0) THEN
!      PRINT 287,'CT,q1/2/3', CT,qq1,qq2,qq3
      end if
  287 FORMAT(a10,5X,4I2,3X,1PE10.4)
  299  CONTINUE


!      USOLI(LOQ,1)  = USOLI(LOQ,1) + 2.4*0.21E-5
!      USOLI(LCOQ,1) = USOLI(LCOQ,1)+ 8.*4.E-4*1.E-4

!      USOLI(LOQ,1)  = USOLI(LOQ,2) !floating
!      USOLI(LCOQ,1) = USOLI(LCOQ,2)
      OQSUR  = (USOLI(LOQ,1)/2.-USOL(LO2,1))/USOL(LO2,1)*1000
      COQSUR = (USOLI(LCOQ,1)/2.-USOL(LCO2,1))/USOL(LCO2,1)*1000
!      PRINT 285,USOLI(LOQ,1),USOLI(LCOQ,1)
!      PRINT 286, FLOWI(LOQ),2.*FLOW(LO2)
  285 FORMAT(20X,'OQ,COQ=',1P4E12.5)
  286 FORMAT(20X,'SUROQ,SURCOQ=',1P4E12.5)
C
C Ask aboutn sulfur budget in the main code
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
   1  CONTINUE
C ***** END THE TIME-STEPPING LOOP *****
  22  CONTINUE

C************OUTPUT****************************
! 185  FORMAT(6X,'H2SI',7X,'HSI',9X,'SIO',8X,'SIO2',
!     & 8X,'H2SIO4',8X,'HSIO')
! 186  FORMAT(6X,'H2S',8X,'HS',10X,'SO',9X,'SO2',
!     & 9X,'H2SO4',9X,'HSO')
! 187  FORMAT(6X,'SIO21',7X,'SIO23',8X,'HSIO3',6X,'SIO3',7X,
!     & 'SIS2',7X,'SIS3')
! 188  FORMAT(6X,'SO21',8X,'SO23',9X,'HSO3',7X,'SO3',8X,
!     & 'S3',9X,'S4')

 185  FORMAT(6X,'H2CQ',10X,'Q',9X,'H2Q',10X,'QH',
     & 9X,'HOQ',9X,'H2OQ',9X,'O2Q',9X,'OQO',9X,'CQ')
 186  FORMAT(6X,'H2CO',10X,'O',9X,'H2O',10X,'OH',
     & 9X,'HO2',9X,'H2O2',9X,'O3',11X,'O3',9X,'CO')

 187  FORMAT(6X,'CH3OQH',6X,'CH3QOH',6X,'CH3OQ',7X,'CH3QO',7X,
     & 'N2Q',10X,'NQ',10X,'NOQ',9X,'HNOQ',7X,'HNO2Q')
 188  FORMAT(6X,'CH3OOH',6X,'CH3OOH',6X,'CH3O2',7X,'CH3O2',7X,
     & 'N2O',10X,'NO',10X,'NO2',9X,'HNO2',7X,'HNO3')

 189  FORMAT(6X,'HOQNO2',6X,'HO2NOQ',7X,'NO2Q',7X,'N2O4Q',7X,'N2QO4',
     & 8X,'OQ',11X,'SQ',9X,'SOQ',7X,'H2SO3Q',8X,'HSQ',9X,'COQ')  
 192  FORMAT(6X,'HO2NO2',6X,'HO2NO2',7X,'NO3', 8X,'N2O5', 8X,'N2O5',
     & 9X,'O2',11X,'SO',9X,'SO2',8X,'H2SO4'8X,'HSO',9X,'CO2')

 193  FORMAT(6X,'Q1D',9X,'H3CQ',10X,'HCQ',7X,'SQ21',8X,'SQ23',
     & 8X,'HSO2Q',9X,'SO2Q')  
 194  FORMAT(6X,'O1D',9X,'H3CO',10X,'HCO', 7X,'SO21', 8X,'SO23',
     & 9X,'HSO3',9X,'SO3')
  
 182  FORMAT(3X,1P9E12.5)
 183  FORMAT(3X,1P11E12.5)
 184  FORMAT(3X,1P2E12.5)
  66  PRINT 106
C-AP      PRINT *, 'ITERATION =', IT 
! 106  FORMAT(//1X,'MIXING RATIO OF ISOTOPES scaled by S*/Stot')
 106  FORMAT(//1X,'MIXING RATIO OF ISOTOPES AND NORMAL SPECIES') 
! 107  FORMAT(//1X,'NUMBER DENSITY OF ISOTOPES scaled by S*/Stot')
 107  FORMAT(//1X,'NUMBER DENSITY OF ISOTOPES AND NORMAL SPECIES')
C-PL WHEN STEP=1, ONLY UPDATE SHORT-LIVED SPECIES') ! SEE DOCHEMI
      PRINT 185
      PRINT 182, USOLI(LH2CQ,1),USOLI(LQ,1),
     2 USOLI(LH2Q,1),USOLI(LQH,1),USOLI(LHOQ,1),
     3 USOLI(LH2OQ,1),USOLI(LO2Q,1),USOLI(LOQO,1),
     4 USOLI(LCQ,1)
      PRINT 182, USOLI(LH2CQ,10),USOLI(LQ,10),
     2 USOLI(LH2Q,10),USOLI(LQH,10),USOLI(LHOQ,10),
     3 USOLI(LH2OQ,10),USOLI(LO2Q,10),USOLI(LOQO,10),
     4 USOLI(LCQ,10)
      PRINT 186
      PRINT 182, USOL(LH2CO,1)*DENCOR,USOL(LO,1)*DENCOR,
     2 USOL(LH2O,1)*DENCOR,USOL(LOH,1)*DENCOR,USOL(LHO2,1)*DENCOR,
     3 USOL(LH2O2,1)*DENCOR,USOL(LO3,1)*DENCOR,USOL(LO3,1)*DENCOR,
     4 USOL(LCO,1)*DENCOR
      PRINT 182, USOL(LH2CO,10)*DENCOR,USOL(LO,10)*DENCOR,
     2 USOL(LH2O,10)*DENCOR,USOL(LOH,10)*DENCOR,USOL(LHO2,10)*DENCOR,
     3 USOL(LH2O2,10)*DENCOR,USOL(LO3,10)*DENCOR,USOL(LO3,10)*DENCOR,
     4 USOL(LCO,10)*DENCOR
C-PL PRINT OUT THE MIXING RATIO FOR ALL ISOTOPIC SPECIES
C    WE WANT THE REAL NUMBER FOR ISOTOPIC SPECIES, SO YOU DO NOT NEED
C    *2 OR 3
      print*, ''
      PRINT 187
      PRINT 182, USOLI(LCH3OQH,1),USOLI(LCH3QOH,1),
     2 USOLI(LCH3OQ,1),USOLI(LCH3QO,1),USOLI(LN2Q,1),
     3 USOLI(LNQ,1),USOLI(LNOQ,1),USOLI(LHNOQ,1),  
     4 USOLI(LHNO2Q,1)                                                                                   
      PRINT 182, USOLI(LCH3OQH,10),USOLI(LCH3QOH,10),
     2 USOLI(LCH3OQ,10),USOLI(LCH3QO,10),USOLI(LN2Q,10),
     3 USOLI(LNQ,10),USOLI(LNOQ,10),USOLI(LHNOQ,10),
     4 USOLI(LHNO2Q,10)
      PRINT 188
      PRINT 182, USOL(LCH3OOH,1)*DENCOR,USOL(LCH3OOH,1)*DENCOR,
     2 USOL(LCH3O2,1)*DENCOR,USOL(LCH3O2,1)*DENCOR,USOL(LN2O,1)*DENCOR,
     3 USOL(LNO,1)*DENCOR,USOL(LNO2,1)*DENCOR,USOL(LHNO2,1)*DENCOR,
     4 USOL(LHNO3,1)*DENCOR
      PRINT 182, USOL(LCH3OOH,10)*DENCOR,USOL(LCH3OOH,10)*DENCOR,
     2 USOL(LCH3O2,10)*DENCOR,USOL(LCH3O2,10)*DENCOR,
     3 USOL(LN2O,10)*DENCOR,USOL(LNO,10)*DENCOR,USOL(LNO2,10)*DENCOR,
     4 USOL(LHNO2,10)*DENCOR,USOL(LHNO3,10)*DENCOR
      print*, ''
      PRINT 189
      PRINT 183, USOLI(LHOQNO2,1),USOLI(LHO2NOQ,1),
     2 USOLI(LNO2Q,1),USOLI(LN2O4Q,1),USOLI(LN2QO4,1),
     3 USOLI(LOQ,1),USOLI(LSQ,1),USOLI(LSOQ,1),
     4 USOLI(LH2SO3Q,1),USOLI(LHSQ,1),USOLI(LCOQ,1)
C     5 USOLI(LHQO,1),USOLI(LHQONO2,1)
      PRINT 183, USOLI(LHOQNO2,10),USOLI(LHO2NOQ,10),
     2 USOLI(LNO2Q,10),USOLI(LN2O4Q,10),
     3 USOLI(LN2QO4,10),USOLI(LOQ,10),USOLI(LSQ,10),
     4 USOLI(LSOQ,10),USOLI(LH2SO3Q,10),
     5 USOLI(LHSQ,10),USOLI(LCOQ,10)
C     6 ,USOLI(LHQO,10),USOLI(LHQONO2,10)
      PRINT 192
      PRINT 183, USOL(LHO2NO2,1)*DENCOR,USOL(LHO2NO2,1)*DENCOR,
     2 USOL(LNO3,1)*DENCOR,USOL(LN2O5,1)*DENCOR,USOL(LN2O5,1)*DENCOR,
     3 USOL(LO2,1)*DENCOR,USOL(LSO,1)*DENCOR,USOL(LSO2,1)*DENCOR,
     4 USOL(LH2SO4,1)*DENCOR,USOL(LHSO,1)*DENCOR,USOL(LCO2,1)*DENCOR
C     5 ,USOL(LHO2,1),USOL(LHO2NO2,1)
      PRINT 183, USOL(LHO2NO2,10)*DENCOR,USOL(LHO2NO2,10)*DENCOR,
     2 USOL(LNO3,10)*DENCOR,USOL(LN2O5,10)*DENCOR,
     3 USOL(LN2O5,10)*DENCOR,USOL(LO2,10)*DENCOR,USOL(LSO,10)*DENCOR,
     4 USOL(LSO2,10)*DENCOR,USOL(LH2SO4,10)*DENCOR,
     5 USOL(LHSO,10)*DENCOR,USOL(LCO2,10)*DENCOR
C     5 ,USOL(LHO2,10),USOL(LHO2NO2,10)

C 
      PRINT 107 ! PRINT OUT NUMBER DENSITIES

      PRINT 185
      PRINT 182, DIZ(LH2CQ,1),DIZ(LQ,1),DIZ(LH2Q,1),
     * DIZ(LQH,1),DIZ(LHOQ,1),DIZ(LH2OQ,1),
     * DIZ(LO2Q,1),DIZ(LOQO,1),DIZ(LCQ,1)
      PRINT 182, DIZ(LH2CQ,10),DIZ(LQ,10),DIZ(LH2Q,10),
     * DIZ(LQH,10),DIZ(LHOQ,10),DIZ(LH2OQ,10),
     * DIZ(LO2Q,10),DIZ(LOQO,10),DIZ(LCQ,10)
      PRINT 186
      PRINT 182, SL(LH2CO,1),SL(LO,1),SL(LH2O,1),
     * SL(LOH,1),SL(LHO2,1),SL(LH2O2,1),
     * SL(LO3,1),SL(LO3,1),SL(LCO,1)
      PRINT 182, SL(LH2CO,10),SL(LO,10),SL(LH2O,10),
     * SL(LOH,10),SL(LHO2,10),SL(LH2O2,10),
     * SL(LO3,10),SL(LO3,10),SL(LCO,10)

      PRINT 187
      PRINT 182, DIZ(LCH3OQH,1),DIZ(LCH3QOH,1),DIZ(LCH3OQ,1),
     * DIZ(LCH3QO,1),DIZ(LN2Q,1),DIZ(LNQ,1),
     * DIZ(LNOQ,1),DIZ(LHNOQ,1),DIZ(LHNO2Q,1)
      PRINT 182, DIZ(LCH3OQH,10),DIZ(LCH3QOH,10),DIZ(LCH3OQ,10),
     * DIZ(LCH3QO,10),DIZ(LN2Q,10),DIZ(LNQ,10),
     * DIZ(LNOQ,10),DIZ(LHNOQ,10),DIZ(LHNO2Q,10)
      PRINT 188
      PRINT 182, SL(LCH3OOH,1),SL(LCH3OOH,1),SL(LCH3O2,1),
     * SL(LCH3O2,1),SL(LN2O,1),SL(LNO,1),
     * SL(LNO2,1),SL(LHNO2,1),SL(LHNO3,1)
      PRINT 182, SL(LCH3OOH,10),SL(LCH3OOH,10),SL(LCH3O2,10),
     * SL(LCH3O2,10),SL(LN2O,10),SL(LNO,10),
     * SL(LNO2,10),SL(LHNO2,10),SL(LHNO3,10)

      PRINT 189
      PRINT 183, DIZ(LHOQNO2,1)/2.,DIZ(LHO2NOQ,1),DIZ(LNO2Q,1),
     * DIZ(LN2O4Q,1),DIZ(LN2QO4,1),DIZ(LOQ,1),
     * DIZ(LSQ,1),DIZ(LSOQ,1),DIZ(LH2SO3Q,1)
     * ,DIZ(LHSQ,1),DIZ(LCOQ,1)
      PRINT 183, DIZ(LHOQNO2,10),DIZ(LHO2NOQ,10),DIZ(LNO2Q,10)
     * ,DIZ(LN2O4Q,10),DIZ(LN2QO4,10),DIZ(LOQ,10),
     * DIZ(LSQ,10),DIZ(LSOQ,10),DIZ(LH2SO3Q,10)
     * ,DIZ(LHSQ,10),DIZ(LCOQ,10)
      PRINT 192
      PRINT 183, SL(LHO2NO2,1),SL(LHO2NO2,1),SL(LNO3,1),
     * SL(LN2O5,1),SL(LN2O5,1),SL(LO2,1),
     * SL(LSO,1),SL(LSO2,1),SL(LH2SO4,1)
     * ,SL(LHSO,1),SL(LCO2,1)
      PRINT 183, SL(LHO2NO2,10),SL(LHO2NO2,10),SL(LNO3,10),
     * SL(LN2O5,10),SL(LN2O5,10),SL(LO2,10),
     * SL(LSO,10),SL(LSO2,10),SL(LH2SO4,10)
     * ,SL(LHSO,10),SL(LCO2,10)

C 
      CLOSE(90)


C-PL PRINT OUT FLUX FOR EVERY ISOTOPIC SPECIES IN EACH LAYER 07/2019
      DO 41 I=1,NZ
  41  ZF(I) = Z(I) + 0.5*DZ
C
      DO 43 K=1,NQI
      DO 42 I=1,NZ
  42  SLI(K,I) = USOLI(K,I)*DEN(I)
!C-PL Added molecular diffusion terms below because they were 
!     missing (7/25/19)
      DO 44 I=1,NZ-1
  44  FLUXCHI(K,I) = - DKI(K,I)*(USOLI(K,I+1) - USOLI(K,I))/DZ- 0.5*
     2  (HIZ(K,I)*DEN(I)*USOLI(K,I) + HIZ(K,I+1)*DEN(I+1)*USOLI(K,I+1))
  43  CONTINUE

      DO 45 K=1,NQI
      FLOWI(K) = FLUXCHI(K,1) - (YPI(K,1) - YLI(K,1)*SLI(K,1))*DZ     
      FUPI(K) = FLUXCHI(K,NZ1) + (YPI(K,NZ) - YLI(K,NZ)*SLI(K,NZ))*DZ
      CONI(K) = TPI(K) - TLI(K) + FLOWI(K) - FUPI(K)
  45  CONTINUE
!      PRINT 285, DKI(LOQ,1),USOLI(LOQ,1+1),USOLI(LOQ,1),DZ
C FOR WATER
      FLOWI(LH2Q) = FLUXCHI(LH2Q,11)
      CONI(LH2Q) = TPI(LH2Q) - TLI(LH2Q) + FLOWI(LH2Q) - FUPI(LH2Q)
      DO 46 I=1,10
  46  FLUXCHI(LH2Q,I) = 0.


C PRINT FLUX OUT
      PRINT 155
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED ISOTOPIC SPECIES')
      ZFL = 0.
      ZFT = ZF(NZ)
      IROW = 12
      LR = NQI/IROW + 1
      DO 47 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQI
      write(*,110) (ISPECI(K),K=K1,K2)
      write(*,122) ZFL,(FLOWI(K),K=K1,K2)
      DO 23 I=1,NZ,4
  23  write(*,122) ZF(I),(FLUXCHI(K,I),K=K1,K2)
      write(*,122) ZFT,(FUPI(K),K=K1,K2)
  47  CONTINUE
 110  FORMAT(/5X,'Z',8X,13(A8,1X))
 122  FORMAT(1X,1P13E9.2)
C **************************************************
C
C Calculated column-integrated reaction rates
      DO 120 L=1,NRI
      RATI(L) = 0.
      M = JCHEMI(1,L)
      K = JCHEMI(2,L)
      DO 120 J=1,NZ
  120 RATI(L) = RATI(L) + AI(L,J)*DIZ(M,J)*DIZ(K,J)*1.e5
     
      write(*,179)
  179 FORMAT(/1X,'INTEGRATED REACTION RATES'/)

      WRITE(*,181)
      IROW = 10
      LR = NRI/IROW + 1
      RL = FLOAT(NRI)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 170 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) THEN
        K2 = NRI
        WRITE(*,386) K1,(RATI(K),K=K1,K2),K2
  386   FORMAT(I3,2X,1P10E10.3,2X,I3)
        GO TO 170
      ENDIF
      WRITE(*,180) K1,(RATI(K),K=K1,K2),K2
  180 FORMAT(I3,2X,1P10E10.3,2X,I3)
  170 CONTINUE
      WRITE(*,181)
  181 FORMAT(9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,
     2    '8',9X,'9',8X,'10')

      GOTO 67
C
C ***********************************************
  25  PRINT 301,IERRI
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)
  20  PRINT 300,I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
C-AP END of the iterative loop   

                                       
  67  RETURN 

      END
