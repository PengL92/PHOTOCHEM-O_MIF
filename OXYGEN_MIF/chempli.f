      SUBROUTINE CHEMPLI(XPI,XLI,K)
!      INCLUDE "parameters_add_CO2.txt"
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/parNEQ_LDA.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
      INCLUDE '../INCLUDECHEM/parNF.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc' 
      INCLUDE 'INCLUDECHEM/isotope_name.inc'    

!      PARAMETER(NMAXI=70)
!      PARAMETER(NQI=10, NRI=99, NSPI=44, NSPI1=NSPI+1, NSPI2=NSPI+2)
c      PARAMETER(NMAXI=70)
c      PARAMETER(NQI=31, NRI=353, NSPI=78, NSPI1=NSPI+1, NSPI2=NSPI+2)

      DIMENSION XPI(NZ),XLI(NZ)
C
C      COMMON/ISOTOP/AI(NRI,NZ),ILOSSI(2,NSPI,NMAXI),
C     2  JCHEMI(5,NRI),NUMLI(NSPI),NUMPI(NSPI),TPI(NQI),TLI(NQI),
C     3  YP(NQI,NZ),YL(NQI,NZ),SRI(NQI),TLOSSI(NQI),PHIDEPI(NQI),
C     4 ISPECI(NSPI2),DIZ(NSPI2,NZ),IPRODI(NSPI,NMAXI),PSO4AER
C     5 ,LBOUNDI(NQI),FLOWI(NQI),FUPI(NQI),CONI(NQI)
C   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
C   FOR THE ISOTOPIC SPECIES USING THE INFORMATION IN THE MATRICES 
C   JCHEMI, ILOSSI, AND IPRODI
C   CALLED BY SUBROUTINE SULFUR.
C
 
      DO 1 I=1,NZ
      XPI(I) = 0.
   1  XLI(I) = 0. 
C
c      print *,'hello world'
c      print *,'K =',K 
c      print *,'AI(1,1) =',AI(1,1)
c      print *,'  NUMLI(K) =', NUMLI(k)
C   LOSS FREQUENCY XL
      NL = NUMLI(K)
c      print *, 'NL= ', NL
      DO 2 L=1,NL
      J = ILOSSI(1,K,L)
      M = ILOSSI(2,K,L)
      DO 2 I=1,NZ
   2  XLI(I) = XLI(I) + AI(J,I)*DIZ(M,I)
C

C   PRODUCTION RATE XP
      NP = NUMPI(K)
      DO 3 L=1,NP
      J = IPRODI(K,L)
      M = JCHEMI(1,J)
      N = JCHEMI(2,J)
      DO 3 I=1,NZ
   3  XPI(I) = XPI(I) + AI(J,I)*DIZ(M,I)*DIZ(N,I)
C
c      print *,'hello world' 

      RETURN
      END
C-AP
