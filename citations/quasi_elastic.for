C|------------------------------------------------------------------|
C|  A program to simulate kinematics of a quasi elastic scattering  |
C|  of high-energy protons on a cluster inside a nucleus.	    |
C|      AA -> 1 ==> BB + (A -> 1) =>  1prim + Aprim +  BB           | 
C|               (A is a cluster)                                   | 
C|  by L.Chulkov, GSI, 2001                                         |
C|------------------------------------------------------------------|
	IMPLICIT REAL*8(A-H,P-Z)
	COMMON/CONSTANTS/PI,R2D,UNIT,HBAR,ALPHA
        COMMON /RANDOM/ISEED
        open(20,file='quasi.out',status='unknown')
C--------     Constants are from  -------------------
C---------    Particle Physics  Booklet, July 2000 ---
                PI    = DACOS(-1.0D0)
                R2D   = 180.0D0/PI
                HBAR  = 197.326960D0
                ALPHA = 1.0D0/137.03599976
                UNIT  = 931.494013D0
C-----------------------------------------------
        MAX_STORY=200000   
	ISTORY  = 0
 1	ISTORY  = ISTORY+1
C:    DEP (1,A,B) - experimental uncertanties in momenta
	DEP1    = 1.0D-011
	DEPA    = 1.0D-011	
	DEPB    = 1.0D-011	
C	DEP1    = 5.0D-000
C	DEPA    = 8.3D-000	
C	DEPB    = 8.3D-000
C:	lower cut off 3-momentum	PTR_MIN (MeV/c)
	PTR_MIN =  0.0D0
	ISEED   = 1254
C-----------------------------------------------
C: The case is for  14Be + 1H -> 1H + alpha + 10He at 250 MeV/nucleon
c        TAA   = 250.0D0*14.d0
c        AMAA  =  14.042815522*UNIT
c        AMBBF =  10.052399713 *UNIT
c	Asterix =  gauss(1.0d0)	
c	amBB    =  AmbbF + Asterix
c        AMA  = 4.002603250*UNIT
c        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C: The case is for  14Be + 1H -> 1H + n + 13Be at 250 MeV/nucleon
c        TAA   = 250.0D0*14.d0
c        AMAA  =  14.042815522*UNIT
c        AMBBF =  13.036133834*UNIT
c	Asterix =  gauss(1.0d0)	
c	amBB    =  AmbbF + Asterix
c        AMA  = 1.008664923*UNIT
c        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C: The case is for  8He + 1H -> 1H + 1H + 7H at 250 MeV/nucleon
c        TAA  = 250.0D0*8.d0
c        AMAA = 8.033921838 *UNIT
c        AMBBF = (3.016049268+4*1.008664923)*UNIT
C	Asterix = gauss(4.0d0)	
C	ambb=AmbbF+dabs(Asterix)
C        AMA  = 1.007825032 *UNIT
C        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C: The case is for  8He + 1H -> 1H + 4He + 4n at 670 MeV/nucleon
C         TAA  = 240.0D0*8.d0
C         AMAA = 8.033921838 *UNIT
C         AMBBF = 4.*1.008664923*UNIT
c 	 Asterix = dabs(gauss(20.d0))	
C	 ambb=AmbbF+Asterix
C         AMA  = 4.002603250 *UNIT
C        AMA  = 6.018888072*UNIT 
C         Am1  = 1.007825032*UNIT
C-----------------------------------------------
C-----------------------------------------------
C: The case is for  8He + 1H -> 1H + 1n + 7He at 680 MeV/nucleon
C        TAA  = 82.3D0*8.d0
C        AMAA = 8.033921838 *UNIT
C        AMBBF = 7.028030527*UNIT
C	Asterix = gauss(0.0001d0)	
C	ambb=AmbbF+dabs(Asterix)
C        AMA  = 1.008664923 *UNIT
C        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C-----------------------------------------------
C: The case is for  8He + 1H -> 1H + 2n + 6He at 680 MeV/nucleon
c        TAA  = 82.3D0*8.d0
c        AMAA =   8.033921838 *UNIT
c        AMBBF =  2.*1.008664923*UNIT
C	Asterix = -DLOG(RNDM(ISEED))*10.
c        Asterix=0.0 
C	print *,asterix	
c	ambb=AmbbF+dabs(Asterix)
c        AMA  = 6.018888072 *UNIT
c        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C-----------------------------------------------
C: The case is for  4He + 1H -> 1H + 3He + n at 700 MeV/nucleon
C        TAA  = 700.0D0*4.d0
C        AMAA = 4.002603250*UNIT
C        AMBBF = 1.008664923*UNIT
C	AMBB=AMBBF
C        AMA  = 3.016029310*UNIT 
C	AM1  = 1.007825032*UNIT
C-----------------------------------------------
C: The case is for  6He + 1H -> 1H + 1n + 5He  at 700 MeV/nucleon
C        TAA  = 700.0D0*6.d0
C        AMAA = 6.018888072*UNIT 
C        AMBBF =  5.012223628*UNIT
C	 Asterix =  gauss(0.40d0)	
C	 amBB    =  AmbbF + Asterix
C        AMA  = 1.008664923*UNIT
C        AM1  = 1.007825032*UNIT
C-----------------------------------------------
C: The case is for  6He + 1H -> 4He+1H + 2n   at 700 MeV/nucleon
C: Kinetic energy
c          TAA  = 700.0D0*6.d0
C: Projectile
c          AMAA = 6.018888072*UNIT
C: rest of the projectile 
c          AMBBF =  1.008664923*UNIT*2.
c	  Asterix =  gauss(0.00040d0)	
c	  AMBB    =  AmbbF + Asterix
C: target 
c          AM1  = 1.007825032*UNIT
C: cluster
c         AMA  = 4.002603250*UNIT
C: iternal momentum distribution of a cluster
C----------------The case  12C+p-> p+p+11B -----------------------------
C         TAA  =  400.0D0*12.d0
C         AMAA =  12.00000000*UNIT
C         AMBBF = 11.009305466*UNIT
C	 Asterix = 0.	
c         AMBB    = AMBBF+Asterix
C         AMA  = 1.007825032*UNIT 
C         AM1  = 1.007825032*UNIT
C----------------The case  57Ni+p-> p+p+56Co -----------------------------
c         TAA  =  500.0D0*57.d0
c         AMAA =  56.939800489*UNIT
c         AMBBF = 55.939843937*UNIT
c         QQRAND =  RNDM(ISEED)
c         asterix=0.0
c         IF(QQRAND.le.0.5) asterix = 20.
c         if(ISTORY.le.20)  print  *,asterix,QQRAND          	
c         AMBB    = AMBBF+Asterix
c         AMA  = 1.007825032*UNIT 
c         AM1  = 1.007825032*UNITC
C----------------The case  17Ne+p-> p+p+16F -----------------------------
          TAA  =  500.0D0*17.d0
          AMAA =  17.017697565*UNIT
          AMBBF = 16.011465730*UNIT
          asterix=0.0    	
          AMBB    = AMBBF+Asterix
          AMA  = 1.007825032*UNIT 
          AM1  = 1.007825032*UNIT
C----------- GOLDHABER MODEL       ----------------------------------
             QQQQ = -(AMAA-AMA-AMBB)
             SIGMA=DSQRT(2.0D0*QQQQ*AMA*AMBB/(AMA+AMBB))   
        if(istory.eq.1) print *,' Goldhaber Sigma=', SIGMA 
C----------------------------------------------------------------------
C:           beam parameters
	PAA   = DSQRT(TAA*(TAA+2.D0*AMAA))
	EAA   = DSQRT(AMAA*AMAA+PAA*PAA)
	BETA  = -PAA/EAA
	GAMMA = 1.0D0/DSQRT(1.0D0-BETA*BETA)
	SFIRST= (EAA+AM1)**2-PAA**2
C---------------------------------------------------------------------
C-----  Internal momentum of a cluster (generated from a model) ------ 
  2	PAINXL = GAUSS(SIGMA)
	PAINYL = GAUSS(SIGMA)
	PAINZ  = GAUSS(SIGMA)
	PAIN   = DSQRT(PAINXL*PAINXL+PAINYL*PAINYL+PAINZ*PAINZ)
c	print *,painxl,painyl,painz
C-------------------------------------------------------------
C: off-shell mass of a cluster
        rrtt = AMAA*AMAA+AMBB**2-2.0D0*AMAA*
     1	          DSQRT(AMBB*AMBB+PAIN*PAIN)
	if(rrtt.le.0.0) go to 1
	AMAOFF = DSQRT(AMAA*AMAA+AMBB**2-2.0D0*AMAA*
     1	          DSQRT(AMBB*AMBB+PAIN*PAIN))
	EAIN=DSQRT(AMAOFF*AMAOFF+PAIN*PAIN)
C-------------------------------------------------------------
  	PBXL  = -PAINXL
	PBYL  = -PAINYL
	PBINZ = -PAINZ
	PBIN  = DSQRT(PBXL*PBXL+PBYL*PBYL+PBINZ*PBINZ)
	EBIN  = DSQRT(AMBB*AMBB+PBIN*PBIN)
C---------------------------------------------
C:  off-shell cluster momentum in lab.system (AMAOFF - off shell mass) 	
        CALL LORENTZ(GAMMA,BETA,EAIN,PAINZ,EAINL,PAINZL)
        CALL LORENTZ(GAMMA,BETA,EBIN,PBINZ,EBL,PBZL)
  	PAINL = DSQRT(PAINXL*PAINXL+PAINYL*PAINYL+PAINZL*PAINZL)
	PBL   = DSQRT(PBXL*PBXL+PBYL*PBYL+PBZL*PBZL)
C-----------------------------------------------	
C:     Now generare   << t >> from experimental dsigma/dt
C:     and calculate  << S >>  - Mandelstam variables
 10	 RUNNUM=RNDM(ISEED)
C----  When X-section is known
C         T = DSIGMADT(ISEED) 
C--------- isotropic distribution-----------------------
          T = -1500000.*RNDM(ISEED)
C	  PRINT *,' T: ',T 

C--------------------------------------------------------
	 IF(DSQRT(-T).LE.PTR_MIN) GO TO 10  
	 S = AMAOFF*AMAOFF+AM1*AM1+2.0D0*AM1*EAINL
C-------------------------------------------------
        CALL CENMASS(T,S,AMAOFF,AM1,AMA,AM1
     1		      ,EAPC,PAPC,E1PC,P1PC,THAC,TH1C)
	if(th1c.eq.1975.0d0) go to 1
C	PRINT *,' AMAOFF ',AMAOFF-AMA,PBIN 
C	PRINT *,' CM A',THAC,TH1C
C	PRINT *,' CM P',PAPC,P1PC
C	PRINT *,' CM E',EAPC,E1PC
C-------------- CLUSTER IN CM ------------------
	PHIA = 2.000*PI*RNDM(ISEED)
	PAXP = PAPC*DSIN(THAC)*DCOS(PHIA)
	PAYP = PAPC*DSIN(THAC)*DSIN(PHIA)
	PAZC = PAPC*DCOS(THAC)
C-------------- PROTON IN CM ------------------
	P1XP = -PAXP
	P1YP = -PAYP
	P1ZC = -PAZC
C------- Check CM calculations (OK)  --------------
C	PQQ=DSQRT(P1XP*P1XP+P1YP*P1YP+P1ZC*P1ZC)
C	PRINT *,P1PC,PQQ	
C-----------------------------------------------
	BETACM  = -PAINL/(EAINL+AM1)
	GAMMACM =  1.0D0/DSQRT(1.0D0-BETACM*BETACM)
C------     calculate relative to the direction ------  
C------     of the quasi-particle (cluster)     ------ 
         CALL LORENTZ(GAMMACM,BETACM,EAPC,PAZC,EAP,PAZP)
         CALL LORENTZ(GAMMACM,BETACM,E1PC,P1ZC,E1P,P1ZP)
C	 PRINT *,'PROTON',E1PC-AM1,P1ZC,E1P-AM1,P1ZP
C         PRINT *,'CLUSTER',EAPC-AMA,PAZC,EAP-AMA,PAZP 
C-----------------------------------------------
C	 PRINT *,PAZC,PAZP,P1ZC,P1ZP	
C-----------------------------------------------
C---------   Now rotate back to the beam direction ----------------
	  PAP  = DSQRT(PAXP*PAXP+PAYP*PAYP+PAZP*PAZP)
	  P1P  = DSQRT(P1XP*P1XP+P1YP*P1YP+P1ZP*P1ZP)
	  PBRL = DSQRT(PBXL*PBXL+PBYL*PBYL)
           CT = PAINZL/PAINL  
	   ST = DSQRT(1.0D0-CT*CT)
           CF = PAINXL/PAINL/ST
           SF = PAINYL/PAINL/ST
         CALL DREHUNG(PAXP,PAYP,PAZP,CT,ST,CF,SF,PAXL,PAYL,PAZL)
         CALL DREHUNG(P1XP,P1YP,P1ZP,CT,ST,CF,SF,P1XL,P1YL,P1ZL)
C********************** TEST FOR ROTATION *************************
C         CALL DREHUNG(0.0D0,0.0D0,PAINL,CT,ST,CF,SF,PQQX,PQQY,PQQZ)
C 	  PRINT *,' PQQQ',PQQX,PQQY,PQQZ
C	   PRINT *,PAINXL,PAINYL,PAINZL
C----------- SPOIL THE DATA AND ANALYSE THIS MIST -----------
	 PAXF = PAXL+GAUSS(DEPA)
	 PAYF = PAYL+GAUSS(DEPA)
	 PAZF = PAZL+GAUSS(DEPA)
	    P1XF = P1XL+GAUSS(DEP1)
	    P1YF = P1YL+GAUSS(DEP1)
	    P1ZF = P1ZL+GAUSS(DEP1)
C-----------  Randomize for 2p case  -----------
               CERBER = RNDM(ISEED)
          SAVEX=PAXF   
          SAVEY=PAYF  
          SAVEZ=PAZF   
          IF(CERBER.GE.0.5) GO TO 789
          PAXF=P1XF
          P1XF=SAVEX
          PAYF=P1YF
          P1YF=SAVEY
          PAZF=P1ZF
          P1ZF=SAVEZ
 789      CONTINUE 
C------------------------------------------------
	 PBXF = PBXL+GAUSS(DEPA)
	 PBYF = PBYL+GAUSS(DEPA)
	 PBZF = PBZL+GAUSS(DEPA)
 	 P1F = DSQRT(P1XF*P1XF+P1YF*P1YF+P1ZF*P1ZF)
	 PAF = DSQRT(PAXF*PAXF+PAYF*PAYF+PAZF*PAZF)
	 PBF = DSQRT(PBXF*PBXF+PBYF*PBYF+PBZF*PBZF)
	 E1F = DSQRT(AM1*AM1 + P1F*P1F ) 
	 EAF = DSQRT(AMA*AMA + PAF*PAF ) 
	 EBF = DSQRT(AMBB*AMBB + PBF*PBF ) 
C----------- NOW  ANALYSE THIS MIST ----------
         AMQ = DSQRT((EAA+AM1-E1F-EAF)**2-(PAA-PAZF-P1ZF)**2
     1     -(PAYF+P1YF)**2-(PAXF+P1XF)**2)
        CALL GETANGLE(P1XF,P1YF,P1ZF,TH1,PHI1)
	IF(PHI1.LE.0.0D0) PHI1=PHI1+2*PI
        CALL GETANGLE(PAXF,PAYF,PAZF,THA,PHIA)
        CALL GETANGLE(PBXF,PBYF,PBZF,THB,PHIB)
         AMQB = DSQRT((EBF+EAF)**2-(PAZF+PBZF)**2
     1     -(PAYF+PBYF)**2-(PAXF+PBXF)**2)
C-------------------- Missing Mass 
	    P3F=DSQRT(PAA*PAA+P1F*P1F-2.0*PAA*P1ZF)
	    A=AMAA*AMAA-2.0*EAA*E1F+2.0D0*PAA*P1ZF
	    B=(A+2.0d0*AM1*AM1)
	    QUATER=B+DSQRT(B*B+4.0d0*AM1*AM1*P3F*P3F-A*A)
	    ARAMISS=DSQRT(QUATER)-AMAA

c---------------------If cluster has numerical problems in mass (No problems!)
c	    AMAMISS=DSQRT(EAPC**2-PAPC**2)-AMA
c	    PRINT *,' AMAMISS IN CM',AMAMISS 
c	    AMAMISS=DSQRT(EAP**2-PAP**2)-AMA
c	    PRINT *,' AMAMISS IN LAB',AMAMISS 

C-------------------- Missing Mass 
C	ERES=EAF+EBF
C	PRES=PAZF+PAZF
C        CALL LORENTZ(GAMMA,-BETA,ERES,PRES,ERESF,PRESF)
C        CALL GETANGLE(PAXF+PBXF,PAYF+PBYF,PRESF,THBF,PHIBF)
C:  (PAXF,PAYF,PAZF) cluster (P1XF,P1YF,P1ZF) recoil 
C:  (PBXL,PBYL,PBZL) rest of the nucleus (spectator)
	WRITE(20,1010) PAXF,PAYF,PAZF,P1XF,P1YF,P1ZF,PBXL,PBYL,PBZL,
     1  AMAA-AMA-AMQ,E1F-AM1,TH1,PHI1,EAF-AMA,THA,PHIA,
     2 	EBF-AMBB,THB,PHIB,ARAMISS,T
	IJUMP=(ISTORY/10000)*10000
C	Write(40,*) T
	IF(IJUMP.EQ.ISTORY) PRINT *,'  Istory = ',ISTORY
C: Number of events ISTORY
	IF(ISTORY.LE.MAX_STORY) go to 1
 1010	 FORMAT(21E15.6)
	END
c+++++++++++++++++++++++++++++++++++++++++++++++
	DOUBLE PRECISION FUNCTION RNDM(I)
        DATA II/1/
	RNDM=RAND(II)
	IF(II.EQ.1) II=0
        RETURN
        END
C*********************************************************

C*********************************************************
C:         calculates all the variables in the center-of-mass
C:               when <<T>> and <<S>> are given
        SUBROUTINE CENMASS(T,S,AMA,AMB,AM1,AM2,E1,P1,E2,P2,A1,A2)
        IMPLICIT NONE
        REAL*8 T,S,AMA,AMB,AM1,AM2
        REAL*8 EA,EB,E1,E2,PA,PB,P1,P2,COA1,A1,A2
        REAL*8 SQRS,X,Y,Z,CINEMA,ERROR_CI
        REAL*8  PI,R2D,UNIT,HBAR,ALPHA
	COMMON/CONSTANTS/PI,R2D,UNIT,HBAR,ALPHA
C******* KINEMATICAL FUNCTION
        CINEMA(X,Y,Z)=X*X+Y*Y+Z*Z-2.0D0*(X*Y+Y*Z+Z*X)
        SQRS = DSQRT(S)
        X    = S
        Y    = AMA*AMA
        Z    = AMB*AMB
        PA   = DSQRT(CINEMA(X,Y,Z))/2.0D0/SQRS
        PB   = PA
        EA   = (S+Y-Z)/2.0D0/SQRS
        EB   = (S+Z-Y)/2.0D0/SQRS
        Y    = AM1*AM1
        Z    = AM2*AM2
	ERROR_CI = CINEMA(X,Y,Z)
	if(ERROR_CI.le.0.0d0) A2=1975.0d0
	if(ERROR_CI.le.0.0d0) return
        P1   = DSQRT(CINEMA(X,Y,Z))/2.0D0/SQRS
        P2   = P1
        E1   = (S+Y-Z)/2.0D0/SQRS
	PRINT *,' CM CLUSTER (P, P-codex): (',P1, DSQRT(E1**2-AM1**2),')'
        E2   = (S+Z-Y)/2.0D0/SQRS
        COA1 = (T-AMA*AMA-AM1*AM1+2.0D0*EA*E1)/(2.0D0*PA*P1)
	if(DABS(COA1).GE.1.0d0) A2=1975.0d0
	if(DABS(COA1).GE.1.0d0) RETURN
        A1   =  DACOS(COA1)
        A2   =  PI-A1
        RETURN
        END
C--------------------- from vector to angles ----------
        SUBROUTINE GETANGLE(PX,PY,PZ,TH,PHI)
        IMPLICIT REAL*8(A-H,O-Z)
	COMMON/CONSTANTS/PI,R2D,UNIT,HBAR,ALPHA
	P=DSQRT(PX*PX+PY*PY+PZ*PZ)
        EAP=DSQRT(AMA*AMA+PAP*PAP)
        CT=PZ/P
        TH=DACOS(CT)
        ST = DSQRT(PX*PX+PY*PY)/P
        CPHI = PX/P/ST
        SPHI = PY/P/ST
        PHI=DATAN2(SPHI,CPHI)
   	  RETURN
	  END
C*********************************************************
C:------------  Lorentz transformation ----------------
        SUBROUTINE LORENTZ(G,B,E,P,EP,PP)
        IMPLICIT NONE
        REAL*8 G,B,E,P,EP,PP
                EP  = G*E-G*B*P
                PP  =-G*B*E+G*P
        RETURN
        END
C*********************************************************
           REAL*8 FUNCTION   GAUSS(S)
C               GENERATE GAUSSIAN RANDOM NUMBERS
C               DIST(X)=EXP(-X**2/(2*S*S))/SQRT(2*PI*S*S)
C               WITH ZERO MEAN AND VARIANCE = S**2
C               ACCORDING TO BREND'S ALGOROTHM
C               ->INT.REV.NUCL.PHYS. VOL.4 1986 APPENDIX A
        REAL*8 S,RNDM         
        COMMON /RANDOM/ISEED
        DIMENSION D(32)
        DATA D/0.674489750,0.475859630,0.383771164,
     1  0.328611323,0.291142827,0.263684322,
     1  0.242508452,0.225567444,0.211634166,
     1  0.199924267,0.189910758,0.181225181,
     1  0.173601400,0.166841909,0.160796729,
     1  0.155349717,0.150409384,0.145902577,
     1  0.141770033,0.137963174,0.134441762,
     1  0.131172150,0.128125965,0.125279090,
     1  0.122610883,0.120103560,0.117741707,
     1  0.115511892,0.113402349,0.111402720,
     1  0.109503852,0.107697617/
        DATA U /0.0/
        A=0.0
        I=0
   10   U=U+U
        IF(U.LT.1.0) GO TO 20
        U=U-1.0
        I=I+1
        A=A-D(I)
        GO TO 10
   20   W=D(I+1)*U
        V=W*(0.5*W-A)
   30   U=RNDM(ISEED)
        IF(V.LE.U) GO TO 40
        V=RNDM(ISEED)
        IF(U.GT.V) GO TO 30
        U=(V-U)/(1.0-U)
        GO TO 20
   40   U=(U-V)/(1.0-V)
        U=U+U
        IF(U.LT.1.0) GO TO 50
        U=U-1.0
        GAUSS=S*(W-A)
        RETURN
   50   GAUSS=S*(A-W)
        RETURN
        end
c-------------------------------------------------
C-----  Simulate elastic scattering X-section  ---
C-----        on a free cluster      ------------- 
        REAL FUNCTION DSIGMADT(ISEED)
        IMPLICIT REAL*8 (A-H,O-Z)
	COMMON /SIGMA/ TT(800),SIG(800),PROB(800)
	DATA KEY/1234/
	IF(KEY.NE.1234) GO TO 1
        open(17,file='dsigma_dt.in',status='old')
C:------ prepare integrated distribution  ---------------------
	DO I=1,798
	READ(17,*) X2,X3
C*************!!!!!!!!!!!!! note here 5. is arbitrary chosen for a test case
	TT(I)  = X2
        SIG(I) = X3
	ENDDO
	    S   = 0.000
	PROB(1) = 0.0
	DO I=2,798
	PROB(I) = PROB(I-1)+SIG(I)*(TT(I)-TT(I-1))
	ENDDO    
	DO I=1,798
	PROB(I) = PROB(I)/PROB(798)
C	WRITE(50,*) TT(I),SIG(I),PROB(I)
	ENDDO
	            KEY = 0
C:------ end of preparation -------------------------
 1    	V=RNDM(ISEED)
	DO I=2,798
	INDE=I
	IF(V.GE.PROB(I-1).AND.V.LT.PROB(I)) GO TO 2
	ENDDO
        INDE=798
 2      CONTINUE
	DELTA=PROB(INDE)-PROB(INDE-1)
	EV=TT(INDE-1)+(PROB(INDE)-V)*(TT(INDE)-TT(INDE-1))/DELTA
     	DSIGMADT=-EV*1.0D+06
        if(dsigmadt.Ge.0.0) print *,INDE,v,DSIGMAT,TT(INDE),PROB(INDE) 
         RETURN
         END
C------------------------------------------------------------------
C----         Two consecutive rotations              --------------
C------ first around Z on <phi>, than around new X' on <theta> -----
      SUBROUTINE DREHUNG(PX,PY,PZ,CT,ST,CF,SF,PXP,PYP,PZP)
      IMPLICIT REAL*8(A-H,O-Z)	
	 PXP = PX*CT*CF-PY*SF+PZ*ST*CF	 
	 PYP = PX*CT*SF+PY*CF+PZ*ST*SF	 
	 PZP = -PX*ST+PZ*CT	 
      RETURN
      END
C----------------- Here we are ---------------------------------------
