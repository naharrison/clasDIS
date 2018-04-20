*72*********************************************************************
      SUBROUTINE DSTFU(KF,X,Q2,XDPQ)

       IMPLICIT NONE

C...Gives proton and neutron structure functions according to a few
C...different polarized parametrizations.
C...Note that what is coded is x times the probability distribution,
C...i.e. xdq(x,Q2) etc.

*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      COMMON /ARSTRF/ KFSAVE(2),XSAVE(2),XQ2SAV(2),
     +		XPQSAV(2,-6:6),XDPQSAV(2,-6:6)
      INTEGER KFSAVE
      REAL XSAVE,XQ2SAV,XPQSAV,XDPQSAV
      SAVE /ARSTRF/


      INTEGER KF,KFL,KFA
      REAL X,Q2,XPQ,XDPQ,XDPS
      DIMENSION XPQ(-6:6),XDPQ(-6:6)

C...Reset structure functions.
      DO 100 KFL=-6,6
      XDPQ(KFL)=0.
 100  XDPQSAV(1,KFL)=0.
      XSAVE(1)=X
      XQ2SAV(1)=Q2
      KFSAVE(1)=KF

C...Check x and particle species.
      IF(X.LE.0..OR.X.GE.1.) THEN
	WRITE(6,5000) X
	RETURN
      ENDIF
      KFA=IABS(KF)
      IF(KFA.NE.2112.AND.KFA.NE.2212) THEN
	WRITE(6,5100) KF
	RETURN
      ENDIF

*...  Put into output array  
      CALL PARTON(X,Q2,XPQ,XDPQ)


C...Check positivity and reset above maximum allowed flavour.
      DO 180 KFL=-6,6
      IF(IABS(KFL).GT.LST(12)) THEN 
        XDPQ(KFL)=0.
      ENDIF
 180  CONTINUE

C...  Isospin conjugation for neutron.
*HI>>      
      IF(KFA.EQ.2112.AND.LST(39).EQ.0) THEN
*     IF(KFA.EQ.2112) THEN
*HI<<         
	XDPS = XDPQ(1)
	XDPQ(1) = XDPQ(2)
	XDPQ(2) = XDPS
	XDPS = XDPQ(-1)
	XDPQ(-1) = XDPQ(-2)
	XDPQ(-2) = XDPS
      ENDIF

C...Charge conjugation for antiparticle.
      IF(KF.LT.0) THEN
	DO 170 KFL=1,6
	XDPS = XDPQ(KFL)
	XDPQ(KFL) = XDPQ(-KFL)
	XDPQ(-KFL) = XDPS
  170   CONTINUE
        ENDIF

      DO 120 KFL=-6,6
  120 XDPQSAV(1,KFL)=XDPQ(KFL)

C...Formats for error printouts.
 5000 FORMAT(' Error: x value outside physical range; x =',1P,E12.3)
 5100 FORMAT(' Error: illegal particle code for structure function;',
     &' KF =',I5)
 5200 FORMAT(' Error: unknown structure function; KF, library, set =',
     &3I5)

      RETURN
      END


*72*********************************************************************

      SUBROUTINE POLPAR(IFLAG,X,Q2,UVAL,DVAL,GLUE,QBAR,STR)

      IMPLICIT NONE

      INTEGER IFLAG
      DOUBLE PRECISION X,Q2,UVAL,DVAL,GLUE,QBAR,STR,AUX(5)

c ---- LO polarized parton distributions as described in
c ----    T. Gehrmann and W.J. Stirling:
c ----    "Polarized Parton Distributions of the Nucleon"
c ----    Durham preprint DTP/95/82
c ---- uval = u-ubar
c ---- dval = d-dbar 
c ---- glue 
c ---- qbar = ubar = dbar = 1/2 usea = 1/2 dsea
c ---- str = sbar = 1/2 strsea
 
      call rdarry(x,q2,aux,iflag)
      uval=aux(1)*(x**0.6d0*(1.d0-x)**3)
      dval=aux(2)*(x**0.75d0*(1.d0-x)**4)
      glue=aux(3)*(x**0.5d0*(1.d0-x)**5)
      qbar=aux(4)*(x**0.5d0*(1.d0-x)**6)
      str=aux(5)*(x**0.5d0*(1.d0-x)**6)
 
      return
      end

*72*********************************************************************
      subroutine rdarry(x,q2,aux,iflag)
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      data icount /0/
*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      
      dimension aux(5)
      common/pdist/arraya(151,20,6)

      nx=151
      ndata=nx-1
      nq2pts=20
      nq2inv=nq2pts-1
      q2sta=1.d0
      q2fin=1.d6
      ymin=5.d0
      xmin=10.d0**(-ymin)
      xmax=1.d0

      if (q2.lt.q2sta) then
         q2=q2sta
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
            WRITE(*,*) 'WARNING : rdarry : Q^2 set to minimal value !',Q2
         ENDIF   
      endif
      if (q2.gt.q2fin) then
         q2=q2fin
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
            WRITE(*,*) 'WARNING : rdarry : Q^2 set to maximal value !',Q2
         ENDIF
      endif
      if (x.lt.xmin) then
         x=xmin
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
            WRITE(*,*) 'WARNING : rdarry : X set to minimal value !',X
         ENDIF   
      endif
      if (x.gt.xmax) then
         x=xmax
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
            WRITE(*,*) 'WARNING : rdarry : X set to maximal value !',X
         ENDIF         
      endif
      
      y=dlog10(x)
      ram=(y+ymin)*ndata/ymin+1.d0
      iram=int(ram)
      fraci=ram-dble(iram)
      ram=dlog(q2/q2sta)*nq2inv/dlog(q2fin/q2sta)+1.d0
      jram=int(ram)
      fracj=ram-dble(jram)
      
      do i=1,5
         aux(i)=(arraya(iram,jram,i)*(1.d0-fraci)
     .        +arraya(iram+1,jram,i)*fraci)*(1.d0-fracj)+
     .        (arraya(iram,jram+1,i)*(1.d0-fraci)
     .        +arraya(iram+1,jram+1,i)*fraci)*fracj
      enddo
      
      return
      end

*72*********************************************************************
      subroutine polini
      implicit real*8(a-h,o-z)
      common/pdist/arraya(151,20,6)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      
      open(IPLST(3),file=CPOL,status='UNKNOWN')
      do i=1,20
         do j=1,151
            read(IPLST(3),901) arraya(j,i,1),arraya(j,i,2),
     &           arraya(j,i,3),arraya(j,i,4),arraya(j,i,5)
         enddo
      enddo
      close(IPLST(3))
      
 901  format(5f14.9)
      
      return
      end
*72*********************************************************************
      SUBROUTINE MRSEB
     &                  (X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C***************************************************************C
C								C
C     This is a package for the new MRS(A prime,G) parton       C
C     distributions. The minimum Q^2  value is 5 GeV^2 and the  C
C     x range is, as before 10^-5 < x < 1. MSbar factorization  C
C     is used. The package reads 2 grids, which are in separate C
C     files (A prime=for020.dat/ftn20, G=for021.dat/ftn21).     C  
C     Note that x times the parton distribution is returned,    C
C     Q is the scale in GeV,                                    C
C     and Lambda(MSbar,nf=4) = 231/255 MeV for A prime/G.       C
C								C
C	MODE=20 for MRS(A prime)                                C
C	MODE=21 for MRS(G)                                      C
C								C
C         The reference is :                                    C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         Phys. Lett. B354 (1995) 155-162                       C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE 20 GRID IS 0.00341            C
C         THE FIRST NUMBER IN THE 21 GRID IS 0.00269            C
C                                                               C
C         HI: changed SCALE to Q^2                              C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      DATA ICOUNT /0/
*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      Q2=SCALE
*     Q2STA= .5D0
      Q2STA = .625D0
      Q2FIN = 1310720.D0
      Xmin = 1D-5
      Xmax = 1D0

      if (q2.lt.q2sta) then
         q2=q2sta
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : Q^2 set to minimal value !',Q2
         ENDIF   
      endif
      if (q2.gt.q2fin) then
         q2=q2fin
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : Q^2 set to maximal value !',Q2
         ENDIF
      endif
      if (x.lt.xmin) then
         x=xmin
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : X set to minimal value !',X
         ENDIF   
      endif
      if (x.gt.xmax) then
         x=xmax
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : X set to maximal value !',X
         ENDIF         
      endif

      IF(MODE.EQ.20) 
     .   CALL STRC20(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.21) 
     .   CALL STRC21(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.10) THEN 
         IF(Q2.LT.0.625D0.OR.Q2.GT.1310720.D0) PRINT 99
      IF(Q2.GT.5D0) CALL STRC10(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     X        GLU)
      IF(Q2.LE.5D0) CALL STRC11(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     X        GLU)
 99      FORMAT('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
      ENDIF
*  check for positive definiteness

      IF(UPV.LT.0.0) UPV=0.0
      IF(DNV.LT.0.0) DNV=0.0
      IF(USEA.LT.0.0) USEA=0.0
      IF(DSEA.LT.0.0) DSEA=0.0
      IF(GLU.LT.0.0) GLU=0.0
      IF(STR.LT.0.0) STR=0.0
      IF(CHM.LT.0.0) CHM=0.0
      IF(BOT.LT.0.0) BOT=0.0

      RETURN
      END
*72*********************************************************************
      
      SUBROUTINE STRC20(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C     THIS IS THE NEW  "Aprime" FIT -- Feb 1995 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      save F
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
 
 
      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      OPEN(IPLST(2),FILE=CUNPOL,STATUS='OLD')
      DO 20 N=1,nx-1
         DO 20 M=1,19
            READ(IPLST(2),50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),
     &           F(5,N,M),F(7,N,M),F(6,N,M),F(8,N,M)
C     1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
            DO 25 I=1,8
 25         F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
 20   CONTINUE
      CLOSE(IPLST(2))
      DO 31 J=1,NTENTH-1
         XX(J)=DLOG10(XX(J))+1.1D0
         DO 31 I=1,8
            IF(I.EQ.7) GO TO 31
            DO 30 K=1,19
 30            F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)
     &              /DLOG(F(I,ntenth,K))
 31   CONTINUE
 50   FORMAT(8F10.5)
      DO 40 I=1,8
         DO 40 M=1,19
 40   F(I,nx,M)=0.D0
 10   CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
         G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .        + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
         IF(N.GE.ntenth) GOTO 65
         IF(I.EQ.7) GOTO 65
         FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
         G(I)=FAC**(G(I)/FAC)
 65      CONTINUE
         G(I)=G(I)*(1.D0-X)**N0(I)
 60   CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STRC21(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)

C     THIS IS THE NEW  "G" FIT -- Feb 1995 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      save F
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
 
 
      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      OPEN(IPLST(2),FILE='./pdf/MRSa1.dat',STATUS='OLD')
      DO 20 N=1,nx-1
         DO 20 M=1,19
            READ(IPLST(2),50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),
     &           F(5,N,M),F(7,N,M),F(6,N,M),F(8,N,M)
C     1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
            DO 25 I=1,8
 25         F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
 20   CONTINUE
      CLOSE(IPLST(2))
      DO 31 J=1,NTENTH-1
         XX(J)=DLOG10(XX(J))+1.1D0
         DO 31 I=1,8
            IF(I.EQ.7) GO TO 31
            DO 30 K=1,19
 30         F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
 50   FORMAT(8F10.5)
      DO 40 I=1,8
         DO 40 M=1,19
 40   F(I,nx,M)=0.D0
 10   CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
         G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     &        + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
         IF(N.GE.ntenth) GOTO 65
         IF(I.EQ.7) GOTO 65
         FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
         G(I)=FAC**(G(I)/FAC)
 65      CONTINUE
         G(I)=G(I)*(1.D0-X)**N0(I)
 60   CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STRC11(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)

C     THIS IS THE NEW  "A" FIT -- May 1994 -- low Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,8),G(8),XX(NX),N0(8)
      save F
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,0.625D0,5.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/

      xsave=x  ! don't let x be altered if it's out of range!!

      IF(INIT.NE.0) GOTO 10
      INIT=1
      OPEN(IPLST(3),FILE=CPOL,STATUS='OLD')
      DO 20 N=1,nx-1
         DO 20 M=1,7
            READ(IPLST(3),50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),
     &           F(5,N,M),F(7,N,M),F(6,N,M),F(8,N,M)
C     1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
            DO 25 I=1,8
  25        F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      CLOSE(IPLST(3))
      DO 31 J=1,NTENTH-1
         XX(J)=DLOG10(XX(J))+1.1D0
         DO 31 I=1,8
            IF(I.EQ.7.or.i.eq.5) GO TO 31
            DO 30 K=1,7
 30         F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)
     &              /DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,7
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)*2D0
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7.or.i.eq.5) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)

      x=xsave  !restore x

      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE STRC10(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)

C     THIS IS THE NEW  "A" FIT -- May 1994 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      save F
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
 
 
      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      OPEN(IPLST(2),FILE=CUNPOL,STATUS='OLD')
      DO 20 N=1,nx-1
         DO 20 M=1,19
            READ(IPLST(2),50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),
     &           F(5,N,M),F(7,N,M),F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
            DO 25 I=1,8
 25         F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
 20   CONTINUE
      CLOSE(IPLST(2))
      DO 31 J=1,NTENTH-1
         XX(J)=DLOG10(XX(J))+1.1D0
         DO 31 I=1,8
            IF(I.EQ.7) GO TO 31
            DO 30 K=1,19
 30         F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C============================================================================
C                CTEQ Parton Distribution Functions: Version 4.6
C                             June 21, 1996
C                   Modified: 10/17/96, 1/7/97, 1/15/97
C                             2/17/97, 2/21/97
C                   Last Modified on April 2, 1997
C
C   Ref[1]: "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF RECENT DEEP
C         INELASTIC SCATTERING AND INCLUSIVE JET DATA"
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C       Phys. Rev. D55, 1280 (1997)
C
C   Ref[2]: "CHARM PRODUCTION AND PARTON DISTRIBUTIONS"
C   By: H.L. Lai and W.K. Tung
C       MSU-HEP-61222, CTEQ-622, e-Print Archive: hep-ph/9701256
C       to appear in Z. Phys.
C
C   This package contains 13 sets of CTEQ4 PDFs. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C Ref[1]
C   1    CTEQ4M   Standard MSbar scheme   0.116     298   202    cteq4m.tbl
C   2    CTEQ4D   Standard DIS scheme     0.116     298   202    cteq4d.tbl
C   3    CTEQ4L   Leading Order           0.132     236   181    cteq4l.tbl
C   4    CTEQ4A1  Alpha_s series          0.110     215   140    cteq4a1.tbl
C   5    CTEQ4A2  Alpha_s series          0.113     254   169    cteq4a2.tbl
C   6    CTEQ4A3            ( same as CTEQ4M )
C   7    CTEQ4A4  Alpha_s series          0.119     346   239    cteq4a4.tbl
C   8    CTEQ4A5  Alpha_s series          0.122     401   282    cteq4a5.tbl
C   9    CTEQ4HJ  High Jet                0.116     303   206    cteq4hj.tbl
C   10   CTEQ4LQ  Low Q0                  0.114     261   174    cteq4lq.tbl
C ---------------------------------------------------------------------------
C Ref[2]
C   11   CTEQ4HQ  Heavy Quark             0.116     298   202    cteq4hq.tbl
C   12   CTEQ4HQ1 Heavy Quark:Q0=1,Mc=1.3 0.116     298   202    cteq4hq1.tbl
C        (Improved version of CTEQ4HQ, recommended)
C   13   CTEQ4F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=385)   cteq4f3.tbl
C   14   CTEQ4F4  Nf=4 FixedFlavorNumber  0.111     292   XXX    cteq4f4.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ(4HQ1) for which Q starts at a lower value of 0.7(1.0) GeV.  
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=5.0 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq4(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq4Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ4F3 has, by definition, only 3 flavors and gluon;
C              CTEQ4F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq4Pdf (Iparton, X, Q2)
      Implicit Double Precision (A-H,O-Z)

*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/


      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder
     > / XQrange / Qini, Qmax, Xmin

      Data Warn /.true./
      save Warn
      Q2STA = Qini
      Q2FIN = Qmax
      Xmax = 1D0
 

      if (q2.lt.q2sta) then
         q2=q2sta
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : CTEQ : Q^2 set to minimal value !',Q2
         ENDIF
      endif
      if (q2.gt.q2fin) then
         q2=q2fin
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : CTEQ : Q^2 set to maximal value !',Q2
         ENDIF
      endif
      if (x.lt.xmin) then
         x=xmin
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : CTEQ : X set to minimal value !',X
         ENDIF
      endif
      if (x.gt.xmax) then
         x=xmax
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : CTEQ : X set to maximal value !',X
         ENDIF
      endif
      
      Q=SQRT(Q2)

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq4Pdf: ', X
        Stop
      Endif
      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq4Pdf: ', Q
        Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq4Pdf: '
     &              , Iparton
         Endif
         Ctq4Pdf = 0D0
         Return
      Endif

      Ctq4Pdf = PartonX (Iparton, X, Q)
      if(Ctq4Pdf.lt.0.D0)  Ctq4Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X < Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q < Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                        Now interpolate in x for M1 Qs
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq4 (Iset)
      Implicit Double Precision (A-H,O-Z)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      Parameter (Isetmax=14)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq4m.tbl', 'cteq4d.tbl', 'cteq4l.tbl'
     > , 'cteq4a1.tbl', 'cteq4a2.tbl', 'cteq4m.tbl', 'cteq4a4.tbl'
     > , 'cteq4a5.tbl', 'cteq4hj.tbl', 'cteq4lq.tbl'
     > , 'cteq4hq.tbl', 'cteq4hq1.tbl', 'cteq4f3.tbl', 'cteq4f4.tbl' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
C    Rewrite file name to fit Pepsi distribution format
C
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq4 :', Iset
	    Stop
         Else
*            Tablefile=Flnm(Iset)
*            Tablefile='./pdf/' 
            Open(IPLST(2), File=CUNPOL, Status='OLD', Err=100)
            IU=IPLST(2)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     > 
      Stop
 101  Print*, Tablefile, ' cannot be opened '
      Print*, 'Please input the .tbl file:'
      Read (*,'(A)') Tablefile
      Goto 21
C                             ********************
      End

      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C*********************************************************************
      
      SUBROUTINE POLNLO(IFLAG,X,Q2,
     &                  DXDUV,DXDDV,DXDG,DXDUBAR,DXDDBAR,DXDS)
      IMPLICIT NONE      
      DOUBLE PRECISION AUX(6),XD,Q2D
      DOUBLE PRECISION DXDUV,DXDDV,DXDG,DXDUBAR,DXDDBAR,DXDS      
      INTEGER IFLAG
      REAL X,Q2


c ---- NLO polarized parton distributions as described in
c ---- T. Gehrmann and W.J. Stirling:
c ---- Polarized Parton Distributions of the Nucleon
c ---- Durham preprint DTP/95/82
c ----
c ---- DXDUV = u-ubar
c ---- DXDDV = d-dbar 
c ---- DXDG  = gluon
c ---- DXDUBAR = 1/2 usea
c ---- DXDDBAR = 1/2 dsea
c ---- DXDS = sbar = 1/2 strsea

      XD=DBLE(X)
      Q2D=DBLE(Q2)

      call rdarrynlo(XD,Q2D,AUX,IFLAG)
      DXDUV=aux(1)*(XD**0.6d0*(1.d0-XD)**3)
      DXDDV=aux(2)*(XD**0.75d0*(1.d0-XD)**4)
      DXDG=aux(3)*(XD**0.5d0*(1.d0-XD)**5)
      DXDUBAR=aux(4)*(XD**0.5d0*(1.d0-XD)**6)
      DXDDBAR=aux(5)*(XD**0.5d0*(1.d0-XD)**6)
      DXDS=aux(6)*(XD**0.5d0*(1.d0-XD)**6)
      
      return
      end


      subroutine rdarrynlo(x,q2,aux,iflag)
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension aux(6)
      common/pdist/arraya(151,20,6)
*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      
      nx=151
      ndata=nx-1
      nq2pts=20
      nq2inv=nq2pts-1
      q2sta=1.d0
      q2fin=1.d6
      ymin=5.d0
      xmin=10.d0**(-ymin)
      xmax=1.d0

      if (q2.lt.q2sta) then
         q2=q2sta
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : Q^2 set to minimal value !',Q2
         ENDIF   
      endif
      if (q2.gt.q2fin) then
         q2=q2fin
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : Q^2 set to maximal value !',Q2
         ENDIF
      endif
      if (x.lt.xmin) then
         x=xmin
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : X set to minimal value !',X
         ENDIF   
      endif
      if (x.gt.xmax) then
         x=xmax
         IF(CUT(1).lt.X.and.CUT(2).gt.X.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
         WRITE(*,*) 'WARNING : MRSEB : X set to maximal value !',X
         ENDIF         
      endif

      y=dlog10(x)
      ram=(y+ymin)*ndata/ymin+1.d0
      iram=int(ram)
      fraci=ram-dble(iram)
      ram=dlog(q2/q2sta)*nq2inv/dlog(q2fin/q2sta)+1.d0
      jram=int(ram)
      fracj=ram-dble(jram)
      
      do i=1,6
         aux(i)=(arraya(iram,jram,i)*(1.d0-fraci)
     .        +arraya(iram+1,jram,i)*fraci)*(1.d0-fracj)+
     .        (arraya(iram,jram+1,i)*(1.d0-fraci)
     .        +arraya(iram+1,jram+1,i)*fraci)*fracj
      enddo

      return
      end

C*********************************************************************
      subroutine nloini
      implicit real*8(a-h,o-z)
      common/pdist/arraya(151,20,6)
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      
      open(IPLST(3),file=CPOL,status='UNKNOWN')
      do i=1,20
         do j=1,151
            read(IPLST(3),901) arraya(j,i,1),arraya(j,i,2),
     &        arraya(j,i,3),arraya(j,i,4),arraya(j,i,5),arraya(j,i,6)
         enddo
      enddo
      close(IPLST(3))

 901  format(6f14.9)

      return
      end
*72*********************************************************************
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*    G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S     *
*                                                                 *
*                         1994 UPDATE                             *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.4   AND  1.E6                 *
*             X         BETWEEN  1.E-5  AND   1.                  *
*   LARGE-X REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION   *
*   IS NEGLIGIBLY SMALL, WERE EXCLUDED FROM THE FIT.              *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*      NLO :  LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131.                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*   IF NEEDED, HEAVY QUARK DENSITIES CAN BE TAKEN FROM THE 1991   *
*   GRV PARAMETRIZATION.                                          *
*                                                                 * 
*   NLO DISTRIBUTIONS ARE GIVEN IN MS-BAR FACTORIZATION SCHEME    *
*     (SUBROUTINE GRV94HO) AS WELL AS IN THE DIS SCHEME (GRV94DI),  *
*   THE LEADING ORDER PARAMETRIZATION IS PROVIDED BY "GRV94LO".   *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS : 
*
*    X   = MOMENTUM FRACTION 
*    Q2  = SCALE Q**2 IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION) :
*
*    UV  = U(VAL) = U - U(BAR)
*    DV  = D(VAL) = D - D(BAR)
*    DEL = D(BAR) - U(BAR)
*    UDB = U(BAR) + D(BAR) 
*    SB  = S = S(BAR)
*    GL  = GLUON
*
*...LO PARAMETRIZATION :
*
*72*********************************************************************      
      SUBROUTINE GRV94LO (X, Q2, UV, DV, DEL, UDB, SB, GL)
      IMPLICIT DOUBLE PRECISION (A - Z)
*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      INTEGER ICOUNT
      DATA ICOUNT/0/
      MU2  = 0.23
      LAM2 = 0.2322 * 0.2322
*HI>>
      IF(Q2.le.MU2) THEN
         Q2=MU2
         IF(CUT(5).lt.Q2.and.CUT(6).gt.Q2.AND.ICOUNT.lt.10) THEN
            ICOUNT = ICOUNT + 1
            WRITE(*,*)
     &'WARNING : GRV94LO predictions are not correct for your setup !'
         ENDIF
         S= 0.D0
         DS=S
         S2=S
         S3=S
      ELSE
         S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
         DS = DSQRT (S)
         S2 = S * S
         S3 = S2 * S
      ENDIF
*HI<<       
*...  UV :
      NU  =  2.284 + 0.802 * S + 0.055 * S2
      AKU =  0.590 - 0.024 * S
      BKU =  0.131 + 0.063 * S
      AU  = -0.449 - 0.138 * S - 0.076 * S2
      BU  =  0.213 + 2.669 * S - 0.728 * S2 
      CU  =  8.854 - 9.135 * S + 1.979 * S2
      DU  =  2.997 + 0.753 * S - 0.076 * S2 
      UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...  DV :
      ND  =  0.371 + 0.083 * S + 0.039 * S2 
      AKD =  0.376 
      BKD =  0.486 + 0.062 * S
      AD  = -0.509 + 3.310 * S - 1.248 * S2
      BD  =  12.41 - 10.52 * S + 2.267 * S2 
      CD  =  6.373 - 6.208 * S + 1.418 * S2
      DD  =  3.691 + 0.799 * S - 0.071 * S2 
      DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...  DEL :
      NE  =  0.082 + 0.014 * S + 0.008 * S2 
      AKE =  0.409 - 0.005 * S
      BKE =  0.799 + 0.071 * S
      AE  = -38.07 + 36.13 * S - 0.656 * S2
      BE  =  90.31 - 74.15 * S + 7.645 * S2 
      CE  =  0.0
      DE  =  7.486 + 1.217 * S - 0.159 * S2 
      DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...  UDB :
      ALX =  1.451
      BEX =  0.271 
      AKX =  0.410 - 0.232 * S 
      BKX =  0.534 - 0.457 * S
      AGX =  0.890 - 0.140 * S
      BGX = -0.981
      CX  =  0.320 + 0.683 * S
      DX  =  4.752 + 1.164 * S + 0.286 * S2
      EX  =  4.119 + 1.713 * S
      ESX =  0.682 + 2.978 * S
      UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...  SB :
      ALS =  0.914
      BES =  0.577
      AKS =  1.798 - 0.596 * S
      AS  = -5.548 + 3.669 * DS - 0.616 * S
      BS  =  18.92 - 16.73 * DS + 5.168 * S
      DST =  6.379 - 0.350 * S  + 0.142 * S2
      EST =  3.981 + 1.638 * S
      ESS =  6.402 
      SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...  GL :
      ALG =  0.524
      BEG =  1.088
      AKG =  1.742 - 0.930 * S
      BKG =                     - 0.399 * S2
      AG  =  7.486 - 2.185 * S
      BG  =  16.69 - 22.74 * S  + 5.779 * S2
      CG  = -25.59 + 29.71 * S  - 7.296 * S2
      DG  =  2.792 + 2.215 * S  + 0.422 * S2 - 0.104 * S3
      EG  =  0.807 + 2.005 * S  
      ESG =  3.841 + 0.316 * S 
      GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
      RETURN
      END
*
*...NLO PARAMETRIZATION (MS(BAR)) :
*
*72*********************************************************************      
      SUBROUTINE GRV94HO (X, Q2, UV, DV, DEL, UDB, SB, GL)

       IMPLICIT DOUBLE PRECISION (A - Z)
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
       DS = DSQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  1.304 + 0.863 * S 
       AKU =  0.558 - 0.020 * S
       BKU =          0.183 * S
       AU  = -0.113 + 0.283 * S - 0.321 * S2
       BU  =  6.843 - 5.089 * S + 2.647 * S2 - 0.527 * S3
       CU  =  7.771 - 10.09 * S + 2.630 * S2
       DU  =  3.315 + 1.145 * S - 0.583 * S2 + 0.154 * S3
       UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...DV :
       ND  =  0.102 - 0.017 * S + 0.005 * S2 
       AKD =  0.270 - 0.019 * S
       BKD =  0.260
       AD  =  2.393 + 6.228 * S - 0.881 * S2
       BD  =  46.06 + 4.673 * S - 14.98 * S2 + 1.331 * S3
       CD  =  17.83 - 53.47 * S + 21.24 * S2
       DD  =  4.081 + 0.976 * S - 0.485 * S2 + 0.152 * S3
       DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...DEL :
       NE  =  0.070 + 0.042 * S - 0.011 * S2 + 0.004 * S3 
       AKE =  0.409 - 0.007 * S
       BKE =  0.782 + 0.082 * S
       AE  = -29.65 + 26.49 * S + 5.429 * S2
       BE  =  90.20 - 74.97 * S + 4.526 * S2 
       CE  =  0.0
       DE  =  8.122 + 2.120 * S - 1.088 * S2 + 0.231 * S3
       DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...UDB :
       ALX =  0.877
       BEX =  0.561
       AKX =  0.275 
       BKX =  0.0
       AGX =  0.997
       BGX =  3.210 - 1.866 * S
       CX  =  7.300
       DX  =  9.010 + 0.896 * DS + 0.222 * S2
       EX  =  3.077 + 1.446 * S
       ESX =  3.173 - 2.445 * DS + 2.207 * S
       UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...SB :
       ALS =  0.756
       BES =  0.216
       AKS =  1.690 + 0.650 * DS - 0.922 * S
       AS  = -4.329 + 1.131 * S
       BS  =  9.568 - 1.744 * S
       DST =  9.377 + 1.088 * DS - 1.320 * S + 0.130 * S2
       EST =  3.031 + 1.639 * S
       ESS =  5.837 + 0.815 * S
       SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...GL :
       ALG =  1.014
       BEG =  1.738
       AKG =  1.724 + 0.157 * S
       BKG =  0.800 + 1.016 * S
       AG  =  7.517 - 2.547 * S
       BG  =  34.09 - 52.21 * DS + 17.47 * S
       CG  =  4.039 + 1.491 * S
       DG  =  3.404 + 0.830 * S
       EG  = -1.112 + 3.438 * S  - 0.302 * S2
       ESG =  3.256 - 0.436 * S 
       GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
*
*...NLO PARAMETRIZATION (DIS) :
*
*72*********************************************************************
      SUBROUTINE GRV94DI (X, Q2, UV, DV, DEL, UDB, SB, GL)
       IMPLICIT DOUBLE PRECISION (A - Z)
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
       DS = DSQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  2.484 + 0.116 * S + 0.093 * S2 
       AKU =  0.563 - 0.025 * S
       BKU =  0.054 + 0.154 * S
       AU  = -0.326 - 0.058 * S - 0.135 * S2
       BU  = -3.322 + 8.259 * S - 3.119 * S2 + 0.291 * S3
       CU  =  11.52 - 12.99 * S + 3.161 * S2
       DU  =  2.808 + 1.400 * S - 0.557 * S2 + 0.119 * S3
       UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...DV :
       ND  =  0.156 - 0.017 * S 
       AKD =  0.299 - 0.022 * S
       BKD =  0.259 - 0.015 * S
       AD  =  3.445 + 1.278 * S + 0.326 * S2
       BD  = -6.934 + 37.45 * S - 18.95 * S2 + 1.463 * S3
       CD  =  55.45 - 69.92 * S + 20.78 * S2
       DD  =  3.577 + 1.441 * S - 0.683 * S2 + 0.179 * S3
       DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...DEL :
       NE  =  0.099 + 0.019 * S + 0.002 * S2  
       AKE =  0.419 - 0.013 * S
       BKE =  1.064 - 0.038 * S
       AE  = -44.00 + 98.70 * S - 14.79 * S2 
       BE  =  28.59 - 40.94 * S - 13.66 * S2 + 2.523 * S3
       CE  =  84.57 - 108.8 * S + 31.52 * S2 
       DE  =  7.469 + 2.480 * S - 0.866 * S2 
       DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...UDB :
       ALX =  1.215
       BEX =  0.466
       AKX =  0.326 + 0.150 * S 
       BKX =  0.956 + 0.405 * S
       AGX =  0.272
       BGX =  3.794 - 2.359 * DS
       CX  =  2.014
       DX  =  7.941 + 0.534 * DS - 0.940 * S + 0.410 * S2
       EX  =  3.049 + 1.597 * S
       ESX =  4.396 - 4.594 * DS + 3.268 * S
       UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...SB :
       ALS =  0.175
       BES =  0.344
       AKS =  1.415 - 0.641 * DS 
       AS  =  0.580 - 9.763 * DS + 6.795 * S  - 0.558 * S2
       BS  =  5.617 + 5.709 * DS - 3.972 * S
       DST =  13.78 - 9.581 * S  + 5.370 * S2 - 0.996 * S3
       EST =  4.546 + 0.372 * S2
       ESS =  5.053 - 1.070 * S  + 0.805 * S2
       SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...GL :
       ALG =  1.258
       BEG =  1.846
       AKG =  2.423 
       BKG =  2.427 + 1.311 * S  - 0.153 * S2
       AG  =  25.09 - 7.935 * S
       BG  = -14.84 - 124.3 * DS + 72.18 * S
       CG  =  590.3 - 173.8 * S
       DG  =  5.196 + 1.857 * S
       EG  = -1.648 + 3.988 * S  - 0.432 * S2
       ESG =  3.232 - 0.542 * S 
       GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
*
*...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
*
       FUNCTION FV (X, N, AK, BK, A, B, C, D)
       IMPLICIT DOUBLE PRECISION (A - Z)
       DX = DSQRT (X)
       FV = N * X**AK * (1.+ A*X**BK + X * (B + C*DX)) * (1.- X)**D
       RETURN
       END
*
       FUNCTION FW (X, S, AL, BE, AK, BK, A, B, C, D, E, ES)
       IMPLICIT DOUBLE PRECISION (A - Z)
       LX = DLOG (1./X)
       FW = (X**AK * (A + X * (B + X*C)) * LX**BK + S**AL
     1      * DEXP (-E + DSQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
*
       FUNCTION FWS (X, S, AL, BE, AK, AG, B, D, E, ES)
       IMPLICIT DOUBLE PRECISION (A - Z)
       DX = DSQRT (X)
       LX = DLOG (1./X)
       FWS = S**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1       * DEXP (-E + DSQRT (ES * S**BE * LX))
       RETURN
       END


c********************************************************************
*********************************************************************
*                                                                   *
*    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    *
*                                                                   *
*         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         *
*                   DO-TH 95/13, RAL-TR-95-042                      *
*                  (TO APPEAR IN PHYS. REV. D)                      *       
*          PROBLEMS/QUESTIONS TO VOGELSANG@V2.RL.AC.UK              *
*                 OR  TO STRAT@HAL1.PHYSIK.UNI-DORTMUND.DE          *
*                                                                   *
*   INPUT:   ISET = number of the parton set :                      *
*              ISET = 1  'STANDARD' SCENARIO, NEXT-TO-LEADING ORDER *
*                        (DATA FILE 'STDNLO.GRID' UNIT=11, TO BE    *
*                         DEFINED BY THE USER )                     *
*              ISET = 2  'VALENCE' SCENARIO,  NEXT-TO-LEADING ORDER *
*                        (DATA FILE 'VALNLO.GRID' UNIT=22, TO BE    *
*                         DEFINED BY THE USER )                     *
*              ISET = 3  'STANDARD' SCENARIO, LEADING ORDER         *
*                        (DATA FILE 'STDLO.GRID' UNIT=33, TO BE     *
*                         DEFINED BY THE USER )                     *
*              ISET = 4  'VALENCE' SCENARIO,  LEADING ORDER         *
*                        (DATA FILE 'VALLO.GRID' UNIT=44, TO BE     *
*                         DEFINED BY THE USER )                     *
*                                                                   *
*            X  = Bjorken-x       (between  1.E-4  and  1)          *
*            Q2 = scale in GeV**2 (between  0.4  and   1.E4)        *
*             (for values outside the allowed range the program     *
*              writes a warning and extrapolates to the x and       *
*              Q2 values requested)                                 *
*                                                                   *
*   OUTPUT:  UV = x * ( DELTA u - DELTA u(bar) ),                   *
*            DV = x * ( DELTA d - DELTA d(bar) ),                   *        
*            QB = x * POLARIZED LIGHT SEA,                          *   
*                 QB= x * ( (DELTA U(BAR) + DELTA D(BAR))/2 )       * 
*            ST = x * DELTA STRANGE = x * DELTA STRANGE(BAR)        *     
*            GL = x * DELTA GLUON                                   *
*                                                                   *
*            always x times the distribution is returned            *
*                                                                   *
*            The sets are the result of a combined fit to           *
*            data for the spin asymmetries A_1 (p,n,d)              *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be zero when PARPOL is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
*********************************************************************
*
      SUBROUTINE PARPOL (ISET, X, Q2, UV, DV, QB, ST, GL)

      IMPLICIT NONE
      
      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      
      real*8 xb1,dfint,xb0,q2,uv,x,st,gl,dv,qb,xb,qs,
     &       parton,arrf,xt,xqbf,xuvf,xdvf,xgf,xsf
      integer n,iq,m,ix,iini,iset,npart,nx,nq,narg
     &        ,na
cywu<<
      PARAMETER (NPART=5, NX=42, NQ=21, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XQBF(NX,NQ), 
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      COMMON / INTINI / IINI
      SAVE XUVF, XDVF, XQBF, XSF, XGF, NA, ARRF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.4D0, 0.75D0, 1.0D0, 1.5D0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4 /
       DATA XB / 1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     2           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     3           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     4           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5           0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6           0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9, 1.0 /
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.1.0D-4) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C          GOTO 60
       ENDIF
       IF ( (Q2.LT.0.4D0) .OR. (Q2.GT.1.D4) ) THEN
cerwin           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
cerwin           write(6,*)'q2=',q2
C          GOTO 60
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    FILE - NO. = 11 FOR NLO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.040E-03 )     
*    FILE - NO. = 22 FOR NLO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 9.740E-04 )     
*    FILE - NO. = 33 FOR  LO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.731E-03 )     
*    FILE - NO. = 44 FOR  LO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.846E-03 )     
      IF (IINI.NE.0) GOTO 16
      OPEN(IPLST(3),FILE=CPOL,STATUS='OLD')
C
      DO 15 M = 1, NX-1 
         DO 15 N = 1, NQ
       READ(IPLST(3),90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M)
 15   CONTINUE
      CLOSE(IPLST(3))
 90   FORMAT (5(1PE10.3))
C
      IINI = 1
*.... ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0**0.7)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.6)
        XQBF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0**0.3) 
        XSF(IX,IQ)  = PARTON(4,IQ,IX) / (XB1**7 * XB0**0.3) 
        XGF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**10 * XB0**0.3)
  20  CONTINUE
        XUVF(NX,IQ) = 0.D0
        XDVF(NX,IQ) = 0.D0
        XQBF(NX,IQ) = 0.D0
        XSF(NX,IQ)  = 0.D0
        XGF(NX,IQ)  = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      UV = DFINT(NARG,XT,NA,ARRF,XUVF) * (1.D0-X)**3 * X**0.7
      DV = DFINT(NARG,XT,NA,ARRF,XDVF) * (1.D0-X)**4 * X**0.6 
      QB = DFINT(NARG,XT,NA,ARRF,XQBF) * (1.D0-X)**7 * X**0.3
      ST = DFINT(NARG,XT,NA,ARRF,XSF)  * (1.D0-X)**7 * X**0.3
      GL = DFINT(NARG,XT,NA,ARRF,XGF)  * (1.D0-X)**10 * X**0.3
 60   RETURN
      END
*72*********************************************************************
*     
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION DFINT(NARG,ARG,NENT,ENT,TABLE)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cerw
c     IMPLICIT INTEGER (I-N)
cerw
cywu>>
      DOUBLE PRECISION fac,d,ent,arg,dfint,table
      integer iadr,ifadr,j,jr,ncomb,ient,il,k,jb,nent,narg,
     &        ja,i,m,kd
cywu<<
  
      DIMENSION ARG(5),NENT(5),ENT(63),TABLE(882)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      DFINT=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      DFINT=DFINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									     C
C       POLARIZED PARTON DISTRIBUTIONS (NLO and LO)			     C
C       FROM CERN-TH/71-97 AND HEP-PH 9711440                                C
C	D. de Florian, R. Sassot and O. Sampayo                              C
C 									     C
C	MODE: 1 SET1 NLO (MSbar)					     C
C             2 SET2 NLO (MSbar)    					     C
C             3 SET3 NLO (MSbar) 					     C
C             4 SET1 LO              					     C
C             5 SET2 LO 						     C
C             6 SET3 LO							     C
C									     C
C                    Q2=Q^2						     C
C                    DUV :    X * U VALENCE DISTRIBUTION		     C
C                    DDV :    X * D VALENCE DISTRIBUTION 		     C
C                    DQBAR :  X * UBAR= X * DBAR DISTRIBUTION   	     C
C                    DSTR :   X * STR= X * STRBAR DISTRIBUTION		     C
C                    DGLU :   X * GLUON DISTRIBUTION 			     C
C                    G1P :    X * POLARIZED STRUCTURE FUNCTION (PROTON)	     C
C                    G1N :    X * POLARIZED STRUCTURE FUNCTION (NEUTRON)     C
C                    G1D :    X * POLARIZED STRUCTURE FUNCTION(DEUTERON)     C
C									     C
C       REMEMBER: ALWAYS X*DISTRIBUTION !!!				     C
C 	BEFORE CALLING THE SUBRUTINE `POLFIT` FOR THE FIRST TIME, THE        C
C       SUBROUTINE `INI` MUST BE CALLED (ONLY ONCE) TO READ THE GRIDS.       C
C              (CALL INI)						     C
C									     C
C	RANGE OF VALIDITY OF THE INTERPOLATION:  			     C
C					10**(-4)< X < 0.9		     C
C					1 < Q**2 < 5*10**4                   C
C									     C
C       IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE REPORT TO		     C
C        Daniel.de.Florian@cern.ch					     C
C									     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE POLFIT(MODE,X,Q2,DUV,DDV,DQBAR,DSTR,DGLU)
************************      
*HI    Speed optimized *
************************      
       IMPLICIT NONE

      COMMON/ XARRAY / XARRAY
      COMMON/ GRID / GU,GD,GL,GS,GG,GP,GN
      DOUBLE PRECISION XARRAY(102), GU(76,26), GD(76,26), GL(76,26),
     &     GS(76,26), GG(76,26), GP(76,26), GN(76,26)

      DOUBLE PRECISION XQ(2),PERINOLA,X3,X4,X5,X6,X7,X,Q2
      DOUBLE PRECISION DUV,DDV,DQBAR,DSTR,DGLU
      INTEGER MODE
 		
      XQ(1) = X
      XQ(2) = Q2
      X3=(1.D0-XQ(1))**3.D0
      X4=(1.D0-XQ(1))**4.D0
      X5=XQ(1)**0.5D0
      X6=XQ(1)**0.6D0
      X7=XQ(1)**0.7D0
      DUV = PERINOLA(XQ,GU) * X3* X6
      DDV = PERINOLA(XQ,GD) * X4 * X7 
      DQBAR = PERINOLA(XQ,GL) * X3 * X5
      DSTR = PERINOLA(XQ,GS)  * X3 * X5
      DGLU = PERINOLA(XQ,GG)  * X3 * X5
*HI      G1P = PERINOLA(XQ,GP)  * X3 * X5
*HI      G1N = PERINOLA(XQ,GN)  * X4 * X5
*HI      G1D=(G1P+G1N)*0.5D0*(1.D0-1.5D0*0.058D0)
      RETURN
      END
*72*********************************************************************
      DOUBLE PRECISION FUNCTION PERINOLA(ARG,TABLE)
       IMPLICIT NONE

      
      COMMON/ XARRAY / ENT
      DOUBLE PRECISION ENT(102)
      
      INTEGER NENT(2),NCOMB(7),IENT(7),NARG,KD,M,JA,I,JB,J,JR,IADR,IFADR
      INTEGER IL,K
      DOUBLE PRECISION ARG(2),TABLE(*),D(7),FAC

      NARG=2
      NENT(1)=76
      NENT(2)=26
      KD=1
      M=1
      JA=1
      DO 5 I=1,NARG
         NCOMB(I)=1
         JB=JA-1+NENT(I)
         DO 2 J=JA,JB
            IF (ARG(I).LE.ENT(J)) GO TO 3
 2       CONTINUE
         J=JB
 3       IF (J.NE.JA) GO TO 4
         J=J+1
 4       JR=J-1
         D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
         IENT(I)=J-JA
         KD=KD+IENT(I)*M
         M=M*NENT(I)
 5    JA=JB+1
      PERINOLA=0.D0
 10   FAC=1.D0
      IADR=KD
      IFADR=1
      DO 15 I=1,NARG
         IF (NCOMB(I).EQ.0) GO TO 12
         FAC=FAC*(1.D0-D(I))
         GO TO 15
 12      FAC=FAC*D(I)
         IADR=IADR-IFADR
 15   IFADR=IFADR*NENT(I)
      PERINOLA=PERINOLA+FAC*TABLE(IADR)
      IL=NARG
 40   IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
      DO 50  K=IL,NARG
 50   NCOMB(K)=1
      GO TO 10
 80   IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
*72*********************************************************************
      SUBROUTINE INIDEFLO
************************      
*HI    Speed optimized *
************************      
       IMPLICIT NONE

*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      


      INTEGER K,J,LX1,LQ
      DOUBLE PRECISION XX,QQ
      COMMON/ XARRAY / XARRAY
      COMMON/ GRID / GU,GD,GL,GS,GG,GP,GN
      DOUBLE PRECISION XARRAY(102), GU(76,26), GD(76,26), GL(76,26),
     &              GS(76,26), GG(76,26), GP(76,26), GN(76,26) 

      OPEN(IPLST(3),FILE=CPOL,STATUS='OLD')
      DO  K = 1, 76 
         DO  J = 1, 26
            READ(IPLST(3),40)GU(K,J),GD(K,J),GL(K,J),GS(K,J), 
     &                 GG(K,J),GP(K,J),GN(K,J)
 	END DO
      END DO
      CLOSE(IPLST(3))
 40   FORMAT (7(1PE15.7))

      DO LX1=0,25
         XX=10.D0**(-4.D0+LX1/25.D0*3.D0)
         XARRAY(LX1+1) = DLOG(XX)
      ENDDO
      DO LX1=1,50
         XX=0.1D0+LX1/50.D0*0.9D0
         XARRAY(LX1+26) = DLOG(XX)
      ENDDO
      DO LQ=0,25
         QQ=0.6D0*10**(5.D0*LQ/25.D0)
         XARRAY(LQ+77) = DLOG(QQ)
      ENDDO
 
      RETURN
      END 

*72*********************************************************************
      SUBROUTINE FITPARTON(X,Q2,XPQ,XDPQ)
      IMPLICIT NONE
      REAL X,Q2,XPQ(-6:6),XDPQ(-6:6)

      WRITE(*,*) 'FITPARTON : This is a user function'
      WRITE(*,*) 'Program will stop'
      STOP
      
      RETURN
      END

*72********************************************************************* 
      SUBROUTINE PARTON(X,Q2,XPQ,XDPQ)
      IMPLICIT NONE
C *********************************************************************  
C This subroutine gives the parton distributions (polarized 
C and unpolarized) for the polarized leptoproduction. The 
C parton distributions are optimized for polarized scattering.
C The following sets are included
C
C       LST(15) = 101  Schaefer, Phys. Lett. B 208,2 (1988) 175
C                      for comparison with older PEPSI versions
C       LST(15) = 102  free
C       LST(15) = 103  free
C       LST(15) = 104  Schaefer  et al hep-ph/9505306
C                      using Glueck et al. Z. Phys. C 53 (1992) 127
C       LST(15) = 105  Bartelski et al hep-ph/9502271 Set II(p,n)
C       LST(15) = 106  Bartelski et al hep-ph/9502271 Set II(P,n)
C       LST(15) = 107  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP  
C       LST(15) = 108  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP
C       LST(15) = 109  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP
C       LST(15) = 110  Gehrmann  et al hep-ph/9512406 Gluon A (LO)
C       LST(15) = 111  Gehrmann  et al hep-ph/9512406 Gluon B (LO)
C       LST(15) = 112  Gehrmann  et al hep-ph/9512406 Gluon C (LO)
C       LST(15) = 113  Gehrmann  et al hep-ph/9512406 Gluon A (LO) + (DGLAP)
C       LST(15) = 114  Gehrmann  et al hep-ph/9512406 Gluon B (LO) + (DGLAP)
C       LST(15) = 115  Gehrmann  et al hep-ph/9512406 Gluon C (LO) + (DGLAP)
C...............................................................................      
C       M. Glueck, E. Reya, M. Stratmann and W. Vogelsang,
C       DO-TH 95/13, RAL-TR-95-042                
C       LST(15) = 116  'standard' scenario, next-to-leading order
C       LST(15) = 117  'valence'  scenario, next-to-leading order 
C       LST(15) = 118  'standard' scenario, leading order
C       LST(15) = 119  'valence'  scenario, leading order
C       LST(15) = 120  Stanley J.Brodsky Nucl.Phys. B441(1995)
C       LST(15) = 121  S.Keler & J.F.Owens  Phys.Lett. B266(1991)
C                                         & Phys.Rev. D19(1994)1199
C       LST(15) = 124  D. de Florian et al., hep-ph/9711440  LO set 1      
C       LST(15) = 125  D. de Florian et al., hep-ph/9711440  LO set 2      
C       LST(15) = 126  D. de Florian et al., hep-ph/9711440  LO set 3      
C       LST(15) = 127  D. de Florian et al., hep-ph/9711440 NLO set 1      
C       LST(15) = 128  D. de Florian et al., hep-ph/9711440 NLO set 2      
C       LST(15) = 129  D. de Florian et al., hep-ph/9711440 NLO set 3      
C       LST(15) = 130  Fake sample :
C                      unpolarized  Gehrmann  et al hep-ph/9512406
C                      with Delta u(x) = 0.5 * u(x) and Delta d(x) = 0.
C       LST(15) = 131  Fake sample :
C                      unpolarized  Gehrmann  et al hep-ph/9512406
C                      with Delta d(x) = 0.5 * d(x) and Delta u(x) = 0.
C       LST(15) = 132  fit routine. (Is outside the official code.)
C
C   CTEQ4 collaboration: http: 
C       UNPOL: Low Q2 parametrization is the only one used here
C       POL:  BOGUS, du=0.5* u(x) dd=-0.3*d(x)  0.0 all else
C       LST(15) = 133  CTEQ4LQ  
C
C   MRS low Q2
C       UNPOL MRS 95 low Q2
C       POL:  BOGUS, du=0.5* u(x) dd=-0.3*d(x)  0.0 all else
C       LST(15)= 137 MRS low Q2
C       LST(15) = 144   grsv2000 hep-ph/0011215  LO standard scenario     
C       LST(15) = 145   grsv2000 hep-ph/0011215  LO valence scenario     
C       LST(15) = 146   grsv2000 hep-ph/0011215 NLO standard scenario     
C       LST(15) = 147   grsv2000 hep-ph/0011215 NLO valence scenario     
C..............................................................................
C       unpolarized distributions
C
C        LST(15) = 150   cteq5l                   LO 
C        LST(15) = 151   cteq5m                   NLO MSBAR
C        LST(15) = 152   cteq5m1                  NLO MSBAR (update)
C
C        LST(15) = 161  mrs99   cor01  central gluon, a_s    
C        LST(15) = 162  mrs99   cor02  higher gluon        
C        LST(15) = 163  mrs99   cor03  lower gluon       
C        LST(15) = 164  mrs99   cor04  lower a_s      
C        LST(15) = 165  mrs99   cor05  higher a_s     
C        LST(15) = 166  mrs99   cor06  quarks up        
C        LST(15) = 167  mrs99   cor07  quarks down            
C        LST(15) = 168  mrs99   cor08  strange up           
C        LST(15) = 169  mrs99   cor09  strange down          
C        LST(15) = 170  mrs99   cor10  charm up         
C        LST(15) = 171  mrs99   cor11  charm down       
C        LST(15) = 172  mrs99   cor12  larger d/u          
C
C        LST(15) = 173  cteq6l LO
C        LST(15) = 174  cteq6d DIS NLO
C        LST(15) = 175  cteq6m NLO MSBAR              
***********************************************************************

*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LEPTO2)
*
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),
     &                XLP,YLP,W2LP,Q2LP,ULP
      REAL CUT,PARL,XLP,YLP,W2LP,Q2LP,ULP
      INTEGER LST
      SAVE /LEPTOU/

      INTEGER IMXPDF
      PARAMETER (IMXPDF=40)
      COMMON /PEPADM/CPDFNAM(2,IMXPDF),IPDFNAM(2,IMXPDF),
     &       IPLST(10),CUNPOL,CPOL
      CHARACTER*256 CPDFNAM,CUNPOL,CPOL
      INTEGER IPLST,IPDFNAM
      SAVE /PEPADM/

**************************************************************
*      
*   IPLST(1) =  0 (default) : number of PDF warnings
*   IPLST(2) = 11 (default) : unit -1- for pdf files
*   IPLST(3) = 12 (default) : unit -2- for pdf files      
*
**************************************************************      

      
      REAL X,Q2,XPQ(-6:6),XDPQ(-6:6),ULMASS
      DOUBLE PRECISION DX,DQ2      
      REAL XUV, XDV, XS, XG  
      REAL XUB, XDB, XC, XBOT
      REAL XDUV,XDDV,XDS,XDG  
      REAL XDUB,XDDB,XDC

C...Variables for set 101
      REAL VALPHAU,VBETAU,GAMMA
      REAL VALPHAD,VBETAD
      REAL VGAMMA,VC3
      REAL VAU0,VAU1,VAD0,VAD1
      REAL FFU0,FFU1,FFD0,FFD1,FG1P,FG1N

C.. Variables for set 104
      REAL VS,VN,VKA,VSW
      REAL VALPHA,VBETA  
      REAL VKLEINA,VKLEINB
      REAL VA,VB,VC,VD,VE,VF
      REAL VA0,VA10
      REAL MNUCL,WWW2

      REAL XD,XU,XDD,XDU,XDDS,XDUS
      
C...Variables for sets 107-109
      DOUBLE PRECISION MRSUPV,MRSDNV,MRSDSEA,MRSUSEA
      DOUBLE PRECISION MRSSTR,MRSCHM,MRSBOT,MRSGLU
      INTEGER MODE

C...Variables for sets 122-129
      DOUBLE PRECISION DDXUV,DDXDV,DDXDEL,DDXUDB,DDXSB,DDXGL,
     &       DXDUV,DXDDV,DXDQBAR,DXDSTR,DXDGLU,
     &       DDXUS,DDXDS
C...Variables for sets 105 and 106
      REAL FDELTA,FM
      REAL VA1,VA2,VA3
      REAL VC1,VC2   

C...Variables for set 110-112
      REAL VETAU ,VETAD, VETAQB, VETAG
      REAL VNORMU,VNORMD,VNORMQB,VNORMG
      REAL VAU,   VAD,   VAQB,   VAG
      REAL VBU,   VBD,   VBQB,   VBG
      REAL VGU,   VGD,   VGQB,   VGG
      REAL VRU,   VRD,   VRQB,   VRG

C...Variables for set 113-115
      INTEGER IFLAG
      DOUBLE PRECISION DXDG,DXDUBAR,DXDDBAR,DXDS
      REAL   XUPDB,XDMUB   

      REAL F2N,F2P
      
C...Variables for sets 133-36
      DOUBLE PRECISION DDXPQ(-6:6)
      INTEGER INITCTEQ,CTEQSET
      DOUBLE PRECISION CTQ4PDF
      
C...Variables for set 116-119
      INTEGER ISET

   
C...Variables for set 120-121
      REAL Y,X1,SS,Z
      REAL AJ(4),BJ(4),CJ(4),RJ(4),W(3),XW(3),WW(3),AAA(6,4)
      REAL CONST0(6,4),CONST1(6,4),CONST2(6,4),CONST3(6,4),
     &     SUMR(2),QG(4),QJ(4)
      REAL AQ(4),BQ(4),CQ(4),DQ(4),ST1(4),ST2(4),
     &     QP(4),QM(4),ST3(4),ST4(4)
      DATA XW/-1.775E-02, -1.E-02, -2.2540E-03/,W/5.556E-03,
     &     8.889E-03, 5.556E-03/,WW/0.3872983, 3.1329002E-02,
     &     4.849505/,AJ/ -0.2, -0.2, -0.3, -0.4/,BJ/ 4.7, 3.7,
     &     8., 8./,CJ/ 0., 2.54, 0., 0./,RJ/-0.9582284, 2.068491,
     &     -0.3953224, 14.54266/,AQ/0.757, 3.784, 0.2897, 2./,
     &     BQ/-0.645, -3.672, -0.2637, -1.25/,CQ/3.23, 2.004,
     &     1., 2./,DQ/-3.118, -1.892, -0.9748, -1.25/,ST1/3., 3.,
     &     5., 4./,ST2/4., 4., 6., 5./,ST3/5., 5., 7., 6./,
     &     ST4/6., 6., 8., 7./      
      DATA CONST0/-0.335  , 3.614 , 0.8673 , 0.0 , 0.0 , 0.0  ,
     &     -0.1612 , 4.667 , 0.0    , 0.0 , 0.0 , 0.0  ,
     &     -1.0   , 7.278 , 0.0    , 0.0 , 0.0 , 0.909,
     &     -1.0   , 5.304 , 0.0    , 0.0 , 0.0 , 3.017/
      DATA CONST1/-0.1097 , .8395 ,-1.6637 , 1.1049, 0.0 , 0.0,
     &     -0.2092 , 0.7951,-1.0232 , 0.8616, 0.0 , 0.0,
     &     -0.3823,-0.7904,-1.6629,-0.0133,0.1211,-0.4023,
     &     -0.9342, 1.4654,-3.9141,9.0176,-5.9602,-4.7347/
      DATA CONST2/-0.002442,-0.02186,0.342,-0.2369,0.0,0.0,
     &     0.02657, 0.1081, 0.05799, 0.153,0.0,0.0,
     &     0.02766, 0.8108,0.5719,0.5299,-0.1739,0.0063,
     &     0.5454,  -1.4292,2.8445,-10.426,7.515,3.3594/
      DATA CONST3/ 0.0    , 0.0 , 0.0 , 0.0  , 0.0 , 0.0  ,
     &     0.0    , 0.0 , 0.0 , 0.0  , 0.0 , 0.0  ,
     &     0.0    , 0.0 , 0.0 , 0.0  , 0.0 , 0.0  ,
     &     -0.1668,0.7569,-0.8411,4.0983,-2.7329,-0.9443/
      DATA INITCTEQ/ 1 /
      DATA CTEQSET/10/
C.. Variables for helium amendment
      REAL PERCS,PERCSP,PERCD
      REAL ADMXS,ADMXSP,ADMXD
      REAL FACTS,FACTSP,FACTD
      REAL PARI11
C
C.. Variables for set 144-147
      DOUBLE PRECISION DXDD,DXDU,DUMMYG1P, DUMMYG1N
      DOUBLE PRECISION DXDUB,DXDDB

C.. Variables for set 150-152
      INTEGER IPARTON
      DOUBLE PRECISION CTQ5PDF

C.. Variable for set 173-175
C.. Use IPARTON from definition above 
      DOUBLE PRECISION CTQ6PDF

      INTEGER I,J,K,ICOUNT,ICOUN1,ICERR
      DATA ICOUNT /0/,ICOUN1/0/,ICERR/0/

      ICOUNT = ICOUNT + 1
      


C  Reset all partondensities
      DO I=-6,6
         XPQ(I)  = 0.
         XDPQ(I) = 0.
      ENDDO

C ***************************************************************************
C
C      set 101 Schaefer, Phys. Lett. B 208,2 (1988) 175
C
C ***************************************************************************


      IF(LST(15).EQ.101) THEN

C    unpolarized part  
            
         VALPHAU = 0.588
         VBETAU  = 2.69 
         XUV = GAMMA(VALPHAU)*GAMMA(VBETAU+1.)/GAMMA(VALPHAU+VBETAU+1.)
         XUV = 2./XUV*X**VALPHAU*(1.-X)**VBETAU
         
         VALPHAD = 1.03
         VBETAD  = 6.87
         XDV = GAMMA(VALPHAD)*GAMMA(VBETAD+1.)/GAMMA(VALPHAD+VBETAD+1.)
         XDV = 1./XDV*X**VALPHAD*(1.-X)**VBETAD
         
         VGAMMA = 14.6  
         VC3    = 0.0147
         XUB    = VC3*(1.+VGAMMA)*(1.-X)**VGAMMA
         
         XDB = XUB
         XS  = XUB
         XC  = 0.02*(1.-X)**19

         XPQ( 1)  = XDV + XDB
         XPQ(-1)  = XDB
         XPQ( 2)  = XUV + XUB
         XPQ(-2)  = XUB
         XPQ( 3)  = XS 
         XPQ(-3)  = XS 
         XPQ( 4)  = XC 
         XPQ(-4)  = XC 

C     polarized part

         VAU0 = 0.12
         VAU1 = 0.3*VAU0
         VAD0 = 7.0*VAU0
         VAD1 = 7.0*VAU1
         VALPHAU = 0.588
         VALPHAD = 1.03 
         
         FFU0 = 1./(1.+ VAU0*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFU1 = 1./(1.+ VAU1*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFD0 = 1./(1.+ VAD0*X**(-VALPHAD)*(1.-X)*(1.-X))
         FFD1 = 1./(1.+ VAD1*X**(-VALPHAD)*(1.-X)*(1.-X))
         FG1P =  0.5*(4*XUV*FFU0/9.
     &        - XDV*(4*FFU0+4*FFU1-FFD0+2*FFD1)/27.)
         FG1N =  0.5*( XUV*FFD0/9.
     &        - XDV*(-4*FFU0+8*FFU1+FFD0+FFD1)/27.)
         
         XDUV = 6*(4*FG1P-FG1N)/5.

         VAU0 = 0.12
         VAU1 = 0.3*VAU0
         VAD0 = 7*VAU0  
         VAD1 = 7*VAU1  
         
         VALPHAU = 0.588
         VALPHAD = 1.03 
         
         FFU0 = 1./(1.+ VAU0*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFU1 = 1./(1.+ VAU1*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFD0 = 1./(1.+ VAD0*X**(-VALPHAD)*(1.-X)*(1.-X))
         FFD1 = 1./(1.+ VAD1*X**(-VALPHAD)*(1.-X)*(1.-X))
         FG1P =  0.5*(4*XUV*FFU0/9.
     &        - XDV*( 4*FFU0+4*FFU1-FFD0+2*FFD1)/27.)
         FG1N =  0.5*(XUV*FFD0/9.
     &        - XDV*(-4*FFU0+8*FFU1+FFD0+FFD1)/27.)

         XDDV = 6*(4*FG1N-FG1P)/5.
         XDUB  = 0.
         XDDB  = 0.
         
         XDPQ( 1) = XDDV + XDDB
         XDPQ(-1) = XDDB
         XDPQ( 2) = XDUV + XDUB
         XDPQ(-2) = XDUB

      ELSEIF(LST(15).EQ.104) THEN      
C ***************************************************************************
C 
C      set  104  Schaefer  et al hep-ph/9505306
C 
C ***************************************************************************

C     unpolarized part
         IF (((Q2.LT.0.25 ).OR.(Q2.GT.1.E+8)).OR.
     &        (( X.LT.1.E-5).OR.( X.GT. 1.  )))    THEN
            ICOUN1 = ICOUN1 + 1
            IF(ICOUN1.lt.10) THEN
               WRITE(*,*) 'WARNING: (Q^2,x) exceeds allowed range'
            ENDIF
            IF(LST(22).eq.0.or.LST(25).eq.0) THEN 
               XPQ(2) = 1.
            ENDIF
            RETURN
         ENDIF
        
         MNUCL = ULMASS(2212)
         WWW2  = MNUCL*MNUCL + Q2*(1./X-1.) 
         
         VS  = LOG(LOG(Q2/(0.232*0.232))/LOG(0.25/(0.232*0.232)))
         
         VN  =  0.579 + 0.283*VS + 0.047*VS*VS
         VKA =  0.523 - 0.015*VS
         VA  =  2.22  - 0.59*VS - 0.27*VS*VS
         VB  =  5.95  - 6.19*VS + 1.55*VS*VS
         VD  =  3.57  + 0.94*VS - 0.16*VS*VS
         XDV = VN*X**VKA*(1.+VA*SQRT(X)+ VB*X)*(1.-X)**VD
         
         VN  =  0.663 + 0.191*VS - 0.041*VS*VS + 0.031*VS*VS*VS   
         VKA =  0.326
         VA  = -1.97  + 6.74*VS - 1.96*VS*VS
         VB  = 24.4   -20.7*VS  + 4.08*VS*VS
         VD  = 2.86   + 0.7*VS  - 0.02*VS*VS
         XUV = VN*X**VKA*(1.+VA*SQRT(X)+ VB*X)*(1.-X)**VD - XDV
         
         VALPHA  =  1.396
         VBETA   =  1.331
         VKLEINA =  0.412 - 0.171*VS
         VKLEINB =  0.566 - 0.496*VS
         VA      =  0.363
         VB      = -1.196
         VC      =  1.029 + 1.785*VS - 0.459*VS*VS
         VD      =  4.696 + 2.109*VS
         VE      =  3.838 + 1.944*VS
         VF      =  2.845
         XUB = ( X**VKLEINA*(VA+VB*X+VC*X*X)*(LOG(1./X))**VKLEINB
     &        +VS**VALPHA*EXP(-VE+SQRT(VF*VS**VBETA*LOG(1./X))))  
     &        *(1.-X)**VD
         XDB = XUB
         
         VSW = 0.
         VALPHA  =  0.803
         VBETA   =  0.563
         VKLEINA =  2.082 -  0.577*VS
         VA      = -3.055 +  1.024*VS**0.67
         VB      = 27.4   - 20.0*VS**0.154 
         VD      = 6.22
         VE      = 4.33 + 1.408*VS
         VF      = 8.27 - 0.437*VS
         
         IF (VS.GT.VSW) THEN
            XS = (VS-VSW)**VALPHA/(LOG(1./X)**VKLEINA)
     &           *(1. + VA*SQRT(X) + VB*X)*(1.-X)**VD
     &           *EXP(-VE + SQRT(VF*VS**VBETA*LOG(1./X)))
         ELSE
            XS = 0
         ENDIF
         
         VSW = 0.888
         VALPHA  =  1.01
         VBETA   =  0.37
         VKLEINA =  0.  
         VA      =  0.  
         VB      = 4.24  - 0.804*VS
         VD      = 3.46  + 1.076*VS
         VE      = 4.61  + 1.49*VS 
         VF      = 2.555 + 1.961*VS
         
*MM   >>    if w > 2* m_c = 3GeV
         IF (WWW2.GT.9.AND.VS.GT.VSW) THEN
            XC = (VS-VSW)**VALPHA/(LOG(1./X)**VKLEINA)
     &           *(1. + VA*SQRT(X) + VB*X)*(1.-X)**VD
     &           *EXP(-VE + SQRT(VF*VS**VBETA*LOG(1./X)))
         ELSE
            XC =0.
         ENDIF
         
         VSW = 1.351
         VALPHA  =  1.00
         VBETA   =  0.51
         VKLEINA =  0.  
         VA      =  0.  
         VB      = 1.848
         VD      = 2.929  + 1.396*VS
         VE      = 4.71   + 1.514*VS
         VF      = 4.02   + 1.239*VS
         
*MM   >>   if w > 2*m_b = 10 GeV
         IF (WWW2.GT.100.AND.VS.GT.VSW) THEN
            XBOT = (VS-VSW)**VALPHA/(LOG(1./X)**VKLEINA)
     &           *(1. + VA*SQRT(X) + VB*X)*(1.-X)**VD
     &           *EXP(-VE + SQRT(VF*VS**VBETA*LOG(1./X)))
         ELSE
            XBOT =0.
         ENDIF
         
         VALPHA  =  0.558
         VBETA   =  1.218
         VKLEINA =  1.00  - 0.17*VS
         VKLEINB =  0.
         VA      =          4.879*VS - 1.383*VS*VS
         VB      =  25.92 - 28.97*VS + 5.596*VS*VS
         VC      = -25.69 + 23.68*VS - 1.975*VS*VS
         VD      =  2.537 + 1.718*VS + 0.353*VS*VS
         VE      =  0.595 + 2.138*VS
         VF      =  4.066
         XG = ( X**VKLEINA*(VA+VB*X+VC*X*X)
     &        +VS**VALPHA*EXP(-VE+SQRT(VF*VS**VBETA*LOG(1./X))))
     &        *(1-X)**VD
         
         XPQ( 0)  = XG
         XPQ( 1)  = XDV + XDB
         XPQ(-1)  = XDB
         XPQ( 2)  = XUV + XUB
         XPQ(-2)  = XUB
         XPQ( 3)  = XS 
         XPQ(-3)  = XS 
         XPQ( 4)  = XC 
         XPQ(-4)  = XC 
         XPQ( 5)  = XBOT
         XPQ(-5)  = XBOT
         
C     polarized part
         
         VA0   =  0.225
         VA10  =  0.15 
         VALPHAU = 0.326
         VALPHAD = 0.505
         
         FFU0  = 1./(1.+VA0*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFU1  = 1./(1.+VA0*VA10*X**(-VALPHAU)*(1.-X)*(1.-X))
         FFD0  = 1./(1.+VA0*X**(-VALPHAD)*(1.-X)*(1.-X))
         FFD1  = 1./(1.+VA0*VA10*X**(-VALPHAD)*(1.-X)*(1.-X))
         
         XDUV  =  XUV*FFU0
     &        -XDV*(FFU0+FFU1)/3.
         XDDV  =  XDV*(FFD0-2*FFD1)/3.
         
         XDUB  = 0.
         XDDB  = 0.
         XDG  =  XUV*(1.-FFU0)/2.
     &        +XDV*(2*FFD1-FFD0+FFU1+FFU0-3.)/6.
         
         XDPQ( 0) = XDG
         XDPQ( 1) = XDDV + XDDB
         XDPQ(-1) = XDDB
         XDPQ( 2) = XDUV + XDUB
         XDPQ(-2) = XDUB
      
      ELSEIF((LST(15).EQ.105).OR.(LST(15).EQ.106)) THEN
C***************************************************************************
C     
C       set 105  Bartelski et al hep-ph/9502271 Set II(p,n)
C       set 106  Bartelski et al hep-ph/9502271 Set II(P,n)
C     
C***************************************************************************
 
         XUV =  1.996*X**(-0.462)*(1.-X)**3.96
     &        *(1. - 0.39*SQRT(X) + 5.13*X)*X
         XDV =  0.296*X**(-0.67)*(1.-X)**4.71
     &        *(1. + 5.03*SQRT(X) + 5.56*X)*X
         FDELTA =   0.099*X**(-0.462)*(1.-X)**9.27*(1.+25*X)
         FM     =   0.411*X**(-1.3  )*(1.-X)**9.27
     &        *(1. - 1.15*SQRT(X) + 15.6*X)  
         XUB  = 0.5*X*(0.392*FM - FDELTA)
         XDB  = 0.5*X*(0.392*FM + FDELTA)
         XS   = 0.5*X*0.196*FM
         XC   = 0.5*X*0.020*FM
         XBOT = 0.
         XG   = 0.775*X**(-1.3)*(1.-X)**5.3*(1.+5.2*X)*X
         
         XPQ( 0)  = XG
         XPQ( 1)  = XDV + XDB
         XPQ(-1)  = XDB
         XPQ( 2)  = XUV + XUB
         XPQ(-2)  = XUB
         XPQ( 3)  = XS 
         XPQ(-3)  = XS 
         XPQ( 4)  = XC 
         XPQ(-4)  = XC 
      
C     polarized part
         
         IF (LST(15).EQ.105) THEN
            VA1 =  0.924
            VA2 = -3.237
            VA3 = 11.3  
            XDUV = X**(-0.462)*(1.-X)**3.96*(VA1+VA2*SQRT(X)+VA3*X)*X
            
            VA1 = -0.029
            VA2 = -0.163
            VA3 = -1.644
            XDDV = X**(-0.670)*(1.-X)**4.71*(VA1+VA2*SQRT(X)+VA3*X)*X
            
            VC1 = 0.
            VC2 = -0.993
            VD  = -0.015
            FDELTA =   X**(-0.462)*(1.-X)**9.27*VD*(1.+25.*X)
            FM     =   X**(-0.8  )*(1.-X)**9.27
     &           *(VC1 + VC2*SQRT(X))

            XDUB = 0.5*X*(0.392*FM - FDELTA)            
            XDDB = 0.5*X*(0.392*FM + FDELTA)
            XDS = 0.5*X*0.196*FM
            XDC = 0.5*X*0.020*FM
            
         ELSEIF(LST(15).EQ.106) THEN
            
            VA1 =  0.929
            VA2 = -3.005
            VA3 = 10.24
            XDUV = X**(-0.462)*(1.-X)**3.96*(VA1+VA2*SQRT(X)+VA3*X)*X
            
            VA1 = -0.066
            VA2 = -0.064
            VA3 = -1.644
            XDDV = X**(-0.670)*(1.-X)**4.71*(VA1+VA2*SQRT(X)+VA3*X)*X
            
            VC1 = 0.
            VC2 = -0.861
            VD  = -0.013
            FDELTA =   X**(-0.462)*(1.-X)**9.27*VD*(1.+25.*X)
            FM     =   X**(-0.8  )*(1.-X)**9.27
     &           *(VC1+ VC2*SQRT(X))
            XDUB = 0.5*X*(0.392*FM - FDELTA)
            XDDB = 0.5*X*(0.392*FM + FDELTA)
            XDS = 0.5*X*0.196*FM
            XDC = 0.5*X*0.020*FM
         ENDIF

         XDPQ( 1) = XDDV + XDDB
         XDPQ(-1) = XDDB
         XDPQ( 2) = XDUV + XDUB
         XDPQ(-2) = XDUB
         XDPQ( 3) = XDS 
         XDPQ(-3) = XDS 
         XDPQ( 4) = XDC 
         XDPQ(-4) = XDC  
      ELSEIF((LST(15).GE.107).AND.(LST(15).LE.109)) THEN
C ***************************************************************************
C 
C       LST(15) = 107  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP  
C       LST(15) = 108  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP
C       LST(15) = 109  Gehrmann hep-ph/9512406 Gluon A (NLO) + DGLAP
C         
C ***************************************************************************
         
         IF (LST(15).EQ.107) IFLAG = 0
         IF (LST(15).EQ.108) IFLAG = 1
         IF (LST(15).EQ.109) IFLAG = 2
        
         MODE = 20       ! MSR(A prime) Parametrization


C  unpolarized part 
*     MRS(A prime)
        
         DX=DBLE(X)
         DQ2=DBLE(Q2)
         
         CALL MRSEB(DX,DQ2,MODE,MRSUPV,MRSDNV,MRSUSEA,
     &        MRSDSEA,MRSSTR,MRSCHM,MRSBOT,MRSGLU)

       
         XPQ( 0) = SNGL(MRSGLU)
         XPQ( 1) = SNGL(MRSDNV + MRSDSEA) 
         XPQ(-1) = SNGL(         MRSDSEA)
         XPQ( 2) = SNGL(MRSUPV + MRSUSEA)
         XPQ(-2) = SNGL(         MRSUSEA)
         XPQ( 3) = SNGL(MRSSTR)
         XPQ(-3) = SNGL(MRSSTR)
         XPQ( 4) = SNGL(MRSCHM)  
         XPQ(-4) = SNGL(MRSCHM)
         XPQ( 5) = SNGL(MRSBOT)   
         XPQ(-5) = SNGL(MRSBOT)
         
C     polarized part (evolved set of Gehrmann hep-ph/9512406 Gluon A,B,C (LO))
         CALL POLNLO(IFLAG,X,Q2,DXDUV,DXDDV,DXDG,DXDUBAR,DXDDBAR,DXDS)
         
         XDPQ( 0) = SNGL(DXDG)
         XDPQ( 1) = SNGL(DXDDV+DXDDBAR)
         XDPQ(-1) = SNGL(DXDDBAR)
         XDPQ( 2) = SNGL(DXDUV+DXDUBAR)
         XDPQ(-2) = SNGL(DXDUBAR)
         XDPQ( 3) = SNGL(DXDS) 
         XDPQ(-3) = SNGL(DXDS) 
         
C ***************************************************************************
C 
C       LST(15) = 110  Gehrmann  et al hep-ph/9512406 Gluon A (LO)
C       LST(15) = 111  Gehrmann  et al hep-ph/9512406 Gluon B (LO)
C       LST(15) = 112  Gehrmann  et al hep-ph/9512406 Gluon C (LO)
C
C       LST(15) = 130  Fake sample :
C                      unpolarized  Gehrmann  et al hep-ph/9512406
C                      with Delta u(x) = 0.5 * u(x) and Delta d(x) = 0.
C       LST(15) = 131  Fake sample :
C                      unpolarized  Gehrmann  et al hep-ph/9512406
C                      with Delta d(x) = 0.5 * d(x) and Delta u(x) = 0.
C
C ***************************************************************************



      ELSEIF((LST(15).GE.110).AND.(LST(15).LE.112) .or.
     &       (LST(15).GE.130).AND.(LST(15).LE.131) ) THEN

C      unpolarized part

         XUV = 3.221*X**0.564*(1.-X)**3.726
     &        *(1.-0.6889*X**0.200 + 2.254*X + 1.261*SQRT(X)*X)
         XDV = 0.507*X**0.376*(1.-X)**4.476
     &        *(1.+1.615 *X**0.553 + 3.651*X + 1.299*SQRT(X)*X)
         XUPDB = (X**0.158*(0.738 - 0.981*X + 1.063*X*X)
     &          * (-LOG(X))**0.037 
     &          + 0.00285*EXP(SQRT(-4.010*LOG(X))))*(1.-X)**6.356
         XDMUB = 0.1067 * X**0.4036 * (1.+ 0.4680*X**0.8763 + 
     &            X * 18.6760) * (1.- X)**8.6220
          XS = 0.0034*(-LOG(X))**(-1.15)
     &         *(1.-2.392*SQRT(X) + 7.094*X)*(1.-X)**6.166
     &         *EXP(SQRT(-6.719*LOG(X)))
          XG = (X**0.731*(5.110 - 1.204*X - 1.911*X*X)
     &        * (-LOG(X))**(-0.4718)
     &         + 0.0527*EXP(SQRT(-4.584*LOG(X))))*(1.-X)**5.566

         XPQ( 0)  = XG
         XPQ( 1)  = XDV + (XUPDB +  XDMUB)*0.5 
         XPQ(-1)  =       (XUPDB +  XDMUB)*0.5
         XPQ( 2)  = XUV + (XUPDB -  XDMUB)*0.5
         XPQ(-2)  =       (XUPDB -  XDMUB)*0.5
         XPQ( 3)  = XS 
         XPQ(-3)  = XS 

C  polarized part

C  Gluon A        
         IF(LST(15).EQ.110) THEN
            VETAU   =  0.823
            VETAD   = -0.303 
            VETAQB  = -0.0495 
            VETAG   =  1.9
            
            VAU     =  0.578
            VAD     =  0.666
            VAQB    =  0.520
            VAG     =  0.520
            
            VBU     =  3.73
            VBD     =  4.73
            VBQB    = 15.06 
            VBG     =  9.45
            
            VGU     =  9.38
            VGD     = 10.46
            VGQB    =  2.30
            VGG     =  0.0
            
            VRU     = -4.26 
            VRD     = -5.10  
            VRQB    = -2.00
            VRG     =  0.00
            
            VNORMU  =  (1.+ VGU*VAU/(VAU+VBU+1.))
     &           * GAMMA(VAU)*GAMMA(VBU+1.)/GAMMA(VAU+VBU+1.)
     &           + VRU*GAMMA(VAU+0.5)*GAMMA(VBU+1.)
     &           /GAMMA(VAU+VBU+1.5)
            VNORMD  =  (1.+ VGD*VAD/(VAD+VBD+1.))
     &           * GAMMA(VAD)*GAMMA(VBD+1.)/GAMMA(VAD+VBD+1.)
     &           + VRD*GAMMA(VAD+0.5)*GAMMA(VBD+1.)
     &           /GAMMA(VAD+VBD+1.5)
            VNORMQB  =  (1.+ VGQB*VAQB/(VAQB+VBQB+1.))
     &           * GAMMA(VAQB)*GAMMA(VBQB+1.)/GAMMA(VAQB+VBQB+1.)
     &           + VRQB*GAMMA(VAQB+0.5)*GAMMA(VBQB+1.)
     &           /GAMMA(VAQB+VBQB+1.5)
            VNORMG  =  (1.+ VGG*VAG/(VAG+VBG+1.))
     &           * GAMMA(VAG)*GAMMA(VBG+1.)/GAMMA(VAG+VBG+1.)
     &           + VRG*GAMMA(VAG+0.5)*GAMMA(VBG+1.)
     &           /GAMMA(VAG+VBG+1.5)
            VNORMU  = 1./VNORMU
            VNORMD  = 1./VNORMD
            VNORMQB = 1./VNORMQB
            VNORMG  = 1./VNORMG
            
            XDUV  = VETAU*VNORMU*X**VAU*(1.-X)**VBU*(1.+VGU*X
     &           + VRU*SQRT(X))
            XDDV  = VETAD*VNORMD*X**VAD*(1.-X)**VBD*(1.+VGD*X
     &           + VRD*SQRT(X))
            XDUB  = VETAQB*VNORMQB*X**VAQB*(1.-X)**VBQB*
     &           (1.+VGQB*X + VRQB*SQRT(X))       
            XDDB  = XDUB
            XDS   = XDUB
            XDG   = VETAG*VNORMG*X**VAG*(1.-X)**VBG*(1.+VGG*X
     &           + VRG*SQRT(X))

            XDPQ( 0) = XDG
            XDPQ( 1) = XDDV + XDDB
            XDPQ(-1) = XDDB
            XDPQ( 2) = XDUV + XDUB
            XDPQ(-2) = XDUB
            XDPQ( 3) = XDS 
            XDPQ(-3) = XDS 

C Gluon B
         ELSE IF(LST(15).EQ.111) THEN

            VETAU   =  0.823
            VETAD   = -0.303 
            VETAQB  = -0.0495 
            VETAG   =  1.9
            
            VAU     =  0.585
            VAD     =  0.662
            VAQB    =  0.524
            VAG     =  0.524
            
            VBU     =  3.73
            VBD     =  4.73
            VBQB    = 15.96 
            VBG     =  6.87
            
            VGU     =  9.31
            VGD     = 10.91
            VGQB    =  2.42
            VGG     =  1.
            
            VRU     = -4.28 
            VRD     = -5.09  
            VRQB    = -2.00
            VRG     = -2.00
           
            VNORMU  =  (1.+ VGU*VAU/(VAU+VBU+1.))
     &           * GAMMA(VAU)*GAMMA(VBU+1.)/GAMMA(VAU+VBU+1.)
     &             + VRU*GAMMA(VAU+0.5)*GAMMA(VBU+1.)
     &           /GAMMA(VAU+VBU+1.5)
            VNORMD  =  (1.+ VGD*VAD/(VAD+VBD+1.))
     &           * GAMMA(VAD)*GAMMA(VBD+1.)/GAMMA(VAD+VBD+1.)
     &           + VRD*GAMMA(VAD+0.5)*GAMMA(VBD+1.)
     &           /GAMMA(VAD+VBD+1.5)
            VNORMQB  =  (1.+ VGQB*VAQB/(VAQB+VBQB+1.))
     &           * GAMMA(VAQB)*GAMMA(VBQB+1.)/GAMMA(VAQB+VBQB+1.)
     &           + VRQB*GAMMA(VAQB+0.5)*GAMMA(VBQB+1.)
     &           /GAMMA(VAQB+VBQB+1.5)
            VNORMG  =  (1.+ VGG*VAG/(VAG+VBG+1.))
     &           * GAMMA(VAG)*GAMMA(VBG+1.)/GAMMA(VAG+VBG+1.)
     &           + VRG*GAMMA(VAG+0.5)*GAMMA(VBG+1.)
     &           /GAMMA(VAG+VBG+1.5)
            
            VNORMU  = 1./VNORMU
            VNORMD  = 1./VNORMD
            VNORMQB = 1./VNORMQB
            VNORMG  = 1./VNORMG
            
            XDUV  = VETAU*VNORMU*X**VAU*(1.-X)**VBU
     &           *(1.+VGU*X + VRU*SQRT(X))
            XDDV  = VETAD*VNORMD*X**VAD*(1.-X)**VBD
     &           *(1.+VGD*X + VRD*SQRT(X))
            XDUB  = VETAQB*VNORMQB*X**VAQB*(1.-X)**VBQB*
     &           (1.+VGQB*X + VRQB*SQRT(X))       
            XDDB  = XDUB
            XDS   = XDUB
            XDG   = VETAG*VNORMG*X**VAG*(1.-X)**VBG
     &           *(1.+VGG*X + VRG*SQRT(X))

            XDPQ( 0) = XDG
            XDPQ( 1) = XDDV + XDDB
            XDPQ(-1) = XDDB
            XDPQ( 2) = XDUV + XDUB
            XDPQ(-2) = XDUB
            XDPQ( 3) = XDS 
            XDPQ(-3) = XDS 

C Gluon C
         ELSEIF(LST(15).EQ.112) THEN

            VETAU   =  0. 823
            VETAD   = -0.303 
            VETAQB  = -0.0495 
            VETAG   =  1.9
            
            VAU     =  0.582
            VAD     =  0.660
            VAQB    =  0.456
            VAG     =  0.456
            
            VBU     =  3.73
            VBD     =  4.73
            VBQB    = 11.82 
            VBG     =  8.72
            
            VGU     =  9.50
            VGD     = 11.04
            VGQB    =  2.11
            VGG     =  0.0
            
            VRU     = -4.28 
            VRD     = -5.06  
            VRQB    = -1.95
            VRG     = -3.00

            VNORMU  =  (1.+ VGU*VAU/(VAU+VBU+1.))
     &           * GAMMA(VAU)*GAMMA(VBU+1.)/GAMMA(VAU+VBU+1.)
     &           + VRU*GAMMA(VAU+0.5)*GAMMA(VBU+1.)
     &           /GAMMA(VAU+VBU+1.5)
            VNORMD  =  (1.+ VGD*VAD/(VAD+VBD+1.))
     &           * GAMMA(VAD)*GAMMA(VBD+1.)/GAMMA(VAD+VBD+1.)
     &           + VRD*GAMMA(VAD+0.5)*GAMMA(VBD+1.)
     &           /GAMMA(VAD+VBD+1.5)
            VNORMQB  =  (1.+ VGQB*VAQB/(VAQB+VBQB+1.))
     &           * GAMMA(VAQB)*GAMMA(VBQB+1.)/GAMMA(VAQB+VBQB+1.)
     &           + VRQB*GAMMA(VAQB+0.5)*GAMMA(VBQB+1.)
     &           /GAMMA(VAQB+VBQB+1.5)            
            VNORMG  =  (1.+ VGG*VAG/(VAG+VBG+1.))
     &           * GAMMA(VAG)*GAMMA(VBG+1.)/GAMMA(VAG+VBG+1.)
     &           + VRG*GAMMA(VAG+0.5)*GAMMA(VBG+1.)
     &           /GAMMA(VAG+VBG+1.5)
            
            VNORMU  = 1./VNORMU
            VNORMD  = 1./VNORMD
            VNORMQB = 1./VNORMQB
            VNORMG  = 1./VNORMG
            
            XDUV  = VETAU*VNORMU*X**VAU*(1.-X)**VBU
     &           *(1.+VGU*X + VRU*SQRT(X))
            XDDV  = VETAD*VNORMD*X**VAD*(1.-X)**VBD
     &           *(1.+VGD*X + VRD*SQRT(X))
            XDUB  = VETAQB*VNORMQB*X**VAQB*(1.-X)**VBQB*
     &           (1.+VGQB*X + VRQB*SQRT(X))       
            XDDB  = XDUB
            XDS   = XDUB
            XDG   = VETAG*VNORMG*X**VAG*(1.-X)**VBG
     &           *(1.+VGG*X + VRG*SQRT(X))

            XDPQ( 0) = XDG
            XDPQ( 1) = XDDV + XDDB
            XDPQ(-1) = XDDB
            XDPQ( 2) = XDUV + XDUB
            XDPQ(-2) = XDUB
            XDPQ( 3) = XDS 
            XDPQ(-3) = XDS 

	ELSEIF(LST(15).EQ.130) THEN
C unpolarized  Gehrmann  et al hep-ph/9512406
C with Delta u(x) = 0.5 * u(x) and Delta d(x) = 0.
          XDPQ(2) = 0.5 * XPQ(2)

	ELSEIF(LST(15).EQ.131) THEN
C unpolarized  Gehrmann  et al hep-ph/9512406
C with Delta d(x) = 0.5 * d(x) and Delta u(x) = 0. 
          XDPQ(1) = 0.5 * XPQ(1)
            
         ENDIF

C***************************************************************************
C
C       LST(15) = 113  Gehrmann  et al hep-ph/9512406 Gluon A (LO) + (DGLAP)
C       LST(15) = 114  Gehrmann  et al hep-ph/9512406 Gluon B (LO) + (DGLAP)
C       LST(15) = 115  Gehrmann  et al hep-ph/9512406 Gluon C (LO) + (DGLAP)
C
C**************************************************************************

         ELSEIF ((LST(15).GE.113).AND.(LST(15).LE.115)) THEN 
 
         IF (LST(15).EQ.113) IFLAG = 0
         IF (LST(15).EQ.114) IFLAG = 1
         IF (LST(15).EQ.115) IFLAG = 2

C  unpolarized part 
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *

         DX=DBLE(X)
         DQ2=DBLE(Q2)
         CALL GRV94LO (DX,DQ2
     &          ,DDXUV,DDXDV,DDXDEL,DDXUDB,DDXSB,DDXGL)
                
         XPQ( 0) = SNGL(DDXGL)
         XPQ( 1) = SNGL(DDXDV  + (DDXUDB + DDXDEL)*0.5D0)
         XPQ(-1) = SNGL(         (DDXUDB + DDXDEL)*0.5D0)
         XPQ( 2) = SNGL(DDXUV  + (DDXUDB - DDXDEL)*0.5D0)
         XPQ(-2) = SNGL(         (DDXUDB - DDXDEL)*0.5D0)
         XPQ( 3) = SNGL(DDXSB) 
         XPQ(-3) = SNGL(DDXSB) 

C polarized part (evolved set of Gehrmann hep-ph/9512406 Gluon A,B,C (LO))
         CALL POLPAR(IFLAG,DX,DQ2
     &               ,DXDUV,DXDDV,DXDG,DXDQBAR,DXDS)

         XDPQ( 0) = SNGL(DXDG)
         XDPQ( 1) = SNGL(DXDDV + DXDQBAR)
         XDPQ(-1) = SNGL(DXDQBAR)
         XDPQ( 2) = SNGL(DXDUV + DXDQBAR)
         XDPQ(-2) = SNGL(DXDQBAR)
         XDPQ( 3) = SNGL(DXDS) 
         XDPQ(-3) = SNGL(DXDS) 


C ***************************************************************************
C
C       M. Glueck, E. Reya, M. Stratmann and W. Vogelsang,
C       DO-TH 95/13, RAL-TR-95-042
C
C           LST(15) = 116  'standard' scenario, next-to-leading order
C           LST(15) = 117  'valence'  scenario, next-to-leading order
C           LST(15) = 118  'standard' scenario, leading order
C           LST(15) = 119  'valence'  scenario, leading order
C
C ***************************************************************************

      ELSEIF ((LST(15).GE.116).AND.(LST(15).LE.119)) THEN
       IF (LST(15).EQ.116) ISET = 1
       IF (LST(15).EQ.117) ISET = 2
       IF (LST(15).EQ.118) ISET = 3
       IF (LST(15).EQ.119) ISET = 4

C  unpolarized part
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *
         

         DX=DBLE(X)
         DQ2=DBLE(Q2)
         CALL GRV94LO (DX,DQ2
     &          ,DDXUV,DDXDV,DDXDEL,DDXUDB,DDXSB,DDXGL)

         XPQ( 0) = SNGL(DDXGL)
         XPQ( 1) = SNGL(DDXDV  + (DDXUDB + DDXDEL)*0.5D0)
         XPQ(-1) = SNGL(         (DDXUDB + DDXDEL)*0.5D0)
         XPQ( 2) = SNGL(DDXUV  + (DDXUDB - DDXDEL)*0.5D0)
         XPQ(-2) = SNGL(         (DDXUDB - DDXDEL)*0.5D0)
         XPQ( 3) = SNGL(DDXSB)
         XPQ(-3) = SNGL(DDXSB)

C    Polarisierter Teil

      CALL  PARPOL (ISET,DX,DQ2,
     &      DXDUV, DXDDV, DXDQBAR, DXDS, DXDG)

      XDPQ( 0) = SNGL(DXDG)
      XDPQ( 1) = SNGL(DXDDV) + SNGL(DXDQBAR)
      XDPQ(-1) = SNGL(DXDQBAR)
      XDPQ( 2) = SNGL(DXDUV) + SNGL(DXDQBAR)
      XDPQ(-2) = SNGL(DXDQBAR)
      XDPQ( 3) = SNGL(DXDS)
      XDPQ(-3) = SNGL(DXDS)

C*************************************************************************

      ELSEIF(LST(15).EQ.120)THEN

*     This is polarization distribution of Stanley J.Brodsky et al.
*     Nucl.Phys. B441(1995)
         Y=1.-X
         X1=X**(-1.12)
         DO I = 1,4
            QP(I) = X1*(AQ(I)*Y**ST1(I)+BQ(I)*Y**ST2(I))
            QM(I) = X1*(CQ(I)*Y**ST3(I)+DQ(I)*Y**ST4(I))
            QG(I) = QP(I)+QM(I)
            QJ(I) = QP(I)-QM(I)
         ENDDO
         DO I = 1,2
            QG(I) = QG(I)-2.*QG(3)
            QJ(I) = QJ(I)-2.*QJ(3)
         ENDDO
         DO I = 1,4
            XPQ(I-4+INT(I/4.)*3) = QG(3)
            XDPQ(I-4+INT(I/4.)*3) = QJ(3)
         ENDDO
         
         XPQ(0) = QG(4)
         XPQ(1) = QG(1)
         XPQ(2) = QG(2)
         XDPQ(0) = QJ(4)
         XDPQ(1) = QJ(1)
         XDPQ(2) = QJ(2)
         
      ELSEIF(LST(15).EQ.121)THEN

*     This is polarization distribution of S.Keler & J.F.Owens
*     Phys.Lett. B266(1991)126 and Phys.Rev. D19(1994)1199
         
         IF(Q2.LE.4.) Q2 = 4.
         Z=ALOG(ALOG(Q2/WW(2))/WW(3))
         
         DO I = 1,6
            DO J = 1,4
               AAA(I,J) = CONST0(I,J) + CONST1(I,J)*Z +
     &              CONST2(I,J)*Z*Z + CONST3(I,J)*Z*Z*Z
            ENDDO
         ENDDO
         DO I = 1,2
            SUMR(I) = 0.
            DO 1 J = 1,3
               DO 1 K = 1,50
                  X1=XW(J)+K/50.
                  SS=1.+AAA(3,I)*X1+AAA(4,I)*X1*X1
                  SUMR(I)=SUMR(I)+X1**AAA(1,I)*(1.-X1)*
     &                 *AAA(2,I)*SS*W(J)
 1             CONTINUE
               AAA(6,I)=(5.-2.*I)/SUMR(I)
            ENDDO
            DO I = 1,4
               SS = 1.+AAA(3,I)*X+AAA(4,I)*X*X+AAA(5,I)*X*X*X    
               QG(I) = X**AAA(1,I)*(1.-X)**AAA(2,I)*SS*AAA(6,I)
            ENDDO
         SS=QG(2)
         QG(3)=QG(3)/6.
         QG(2)=QG(1)-SS
         QG(1)=SS
         DO I = 1,4
            SS = 1. + CJ(i)*X
            QJ(I) = X**AJ(I)*(1.-X)**BJ(I)*SS*RJ(I)
         ENDDO
         
         DO I = 1,4
            XPQ(I-4+INT(I/4.)*3) = QG(3)
            XDPQ(I-4+INT(I/4.)*3) = QJ(3)
         ENDDO         
         XPQ(0) = QG(4)
         XPQ(1) = QG(1)
         XPQ(2) = QG(2)
         XDPQ(0) = QJ(4)
         XDPQ(1) = QJ(1)
         XDPQ(2) = QJ(2)
C***********************************************************************
C       LST(15) = 124  D. de Florian et al., hep-ph/9711440  LO set 1      
C       LST(15) = 125  D. de Florian et al., hep-ph/9711440  LO set 2      
C       LST(15) = 126  D. de Florian et al., hep-ph/9711440  LO set 3      
C       LST(15) = 127  D. de Florian et al., hep-ph/9711440 NLO set 1      
C       LST(15) = 128  D. de Florian et al., hep-ph/9711440 NLO set 2      
C       LST(15) = 129  D. de Florian et al., hep-ph/9711440 NLO set 3      
C
C       using GRV 1995 as unpolarized reference set
C
C***********************************************************************
      ELSEIF ((LST(15).GE.124).AND.(LST(15).LE.129)) THEN    

c      IF(X.LE.0.0001)  THEN
c      WRITE(*,*) 'Warning: x exceeds minum value, x set to 0.0001001'
c      X=0.0001001
c      ENDIF
c      IF(X.GE.0.9)     THEN
c      WRITE(*,*) 'Warning: x exceeds maxinum value, x set to 0.89999'
c      X=0.89999
c      ENDIF
c      IF(Q2.LE.1.)     THEN
c      WRITE(*,*)'Warning: Q2 exceeds minimum value, Q2 set to 1.0001'
c      Q2 = 1.0001
c      ENDIF
c      IF(Q2.GE.50000.0) THEN
c      WRITE(*,*) 'Warning: Q2 exceeds maximum value, Q2 set to 49999.9'
c      Q2 = 49999.9
c      ENDIF 

      IF (LST(15).EQ.124) MODE=1
      IF (LST(15).EQ.125) MODE=2
      IF (LST(15).EQ.126) MODE=3
      IF (LST(15).EQ.127) MODE=4
      IF (LST(15).EQ.128) MODE=5
      IF (LST(15).EQ.129) MODE=6
      
C  unpolarized part 
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *

         DX=DBLE(X)
         DQ2=DBLE(Q2)
         IF ((LST(15).GE.124).AND.(LST(15).LE.126)) THEN
         CALL GRV94LO (DX,DQ2
     &          ,DDXUV,DDXDV,DDXDEL,DDXUDB,DDXSB,DDXGL)
         ELSEIF ((LST(15).GE.127).AND.(LST(15).LE.129)) THEN
            CALL GRV94HO (DX,DQ2
     &           ,DDXUV,DDXDV,DDXDEL,DDXUDB,DDXSB,DDXGL)
         ENDIF
 
         XPQ( 0) = SNGL(DDXGL)
         XPQ( 1) = SNGL(DDXDV  + (DDXUDB + DDXDEL)*0.5D0)
         XPQ(-1) = SNGL(         (DDXUDB + DDXDEL)*0.5D0)
         XPQ( 2) = SNGL(DDXUV  + (DDXUDB - DDXDEL)*0.5D0)
         XPQ(-2) = SNGL(         (DDXUDB - DDXDEL)*0.5D0)
         XPQ( 3) = SNGL(DDXSB) 
         XPQ(-3) = SNGL(DDXSB) 

*     polarized part

         CALL  POLFIT(MODE,DX,DQ2,DXDUV,DXDDV,DXDQBAR,DXDSTR,DXDGLU)
         XDPQ( 0) = SNGL(DXDGLU)
         XDPQ( 1) = SNGL(DXDDV) + SNGL(DXDQBAR)
         XDPQ(-1) =               SNGL(DXDQBAR) 
         XDPQ( 2) = SNGL(DXDUV) + SNGL(DXDQBAR)
         XDPQ(-2) =               SNGL(DXDQBAR)
         XDPQ( 3) =               SNGL(DXDSTR)
         XDPQ(-3) =               SNGL(DXDSTR)
      ELSEIF(LST(15).EQ.132)THEN         
*72*********************************************************************
         CALL FITPARTON(X,Q2,XPQ,XDPQ)
*72*********************************************************************
      ELSEIF ((LST(15).GE.133).AND.(LST(15).LE.136)) THEN

C  unpolarized part
*     "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF
*   RECENT DEEP INELASTIC SCATTERING AND INCLUSIVE JET DATA"
*   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, 
*    D. Soper W.K. Tung, H. Weerts
*    Phys. Rev. D55, 1280 (1997) 
*   ( low Q2  set )
*         

*     from data statement above
*       CTEQSET=10  !  Low Q2
         DX=DBLE(X)
         DQ2=DBLE(Q2)
         IF(INITCTEQ.eq.1) THEN
            INITCTEQ=0
            CALL SETCTQ4(CTEQSET)
         ENDIF

*   Note in CTEQ u and d indices are reversed WRT pepsi order:
*   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
*                            for (b, c, s, d, u, g, u_bar, ..., b_bar),


         DDXPQ(-3) = CTQ4PDF(-3,DX,DQ2)
         DDXPQ(-1) = CTQ4PDF(-2,DX,DQ2)
         DDXPQ(-2) = CTQ4PDF(-1,DX,DQ2)
         DDXPQ(0) =  CTQ4PDF(0,DX,DQ2)
         DDXPQ(2) =  CTQ4PDF(1,DX,DQ2)
         DDXPQ(1) =  CTQ4PDF(2,DX,DQ2)
         DDXPQ(3) =  CTQ4PDF(3,DX,DQ2)

* CTEQ routine double precision q(x).  Make single  xq(x)
         DO I=-3,3
            XPQ( I) = X*SNGL(DDXPQ(I))
         ENDDO
         
         XDPQ( 0) = 0.0
         XDPQ( 1) = -0.3*XPQ(1)
         XDPQ(-1) = 0.0
         XDPQ( 2) = 0.5*XPQ(2)
         XDPQ(-2) = 0.0
         XDPQ( 3) = 0.0
         XDPQ(-3) = 0.0

* require F2 non-zero
         IF(XPQ(2).EQ.0.0) XPQ(2)=0.000001

      ELSEIF((LST(15).GE.137).AND.(LST(15).LE.138)) THEN

         MODE = 10       ! MSR(A prime) Low Q2 Parametrization


C  unpolarized part 
*     MRS(A prime)
        
         DX=DBLE(X)
         DQ2=DBLE(Q2)

         CALL MRSEB(DX,DQ2,MODE,MRSUPV,MRSDNV,MRSUSEA,
     &        MRSDSEA,MRSSTR,MRSCHM,MRSBOT,MRSGLU)

       
         XPQ( 0) = SNGL(MRSGLU)
         XPQ( 1) = SNGL(MRSDNV + MRSDSEA) 
         XPQ(-1) = SNGL(         MRSDSEA)
         XPQ( 2) = SNGL(MRSUPV + MRSUSEA)
         XPQ(-2) = SNGL(         MRSUSEA)
         XPQ( 3) = SNGL(MRSSTR)
         XPQ(-3) = SNGL(MRSSTR)
         XPQ( 4) = SNGL(MRSCHM)  
         XPQ(-4) = SNGL(MRSCHM)
         XPQ( 5) = SNGL(MRSBOT)   
         XPQ(-5) = SNGL(MRSBOT)


         
         XDPQ( 0) = 0.0
         XDPQ( 1) = -0.3*XPQ(1)
         XDPQ(-1) = 0.0
         XDPQ( 2) = 0.5*XPQ(2)
         XDPQ(-2) = 0.0
         XDPQ( 3) = 0.0
         XDPQ(-3) = 0.0



C***********************************************************************
C       LST(15) = 144  grsv2000  hep-ph/0011215   LO standard scenario
C       LST(15) = 145  grsv2000  hep-ph/0011215   LO  valence scenario  
C       LST(15) = 146  grsv2000  hep-ph/0011215  NLO standard scenario
C       LST(15) = 147  grsv2000  hep-ph/0011215  NLO  valence scenario

C
C       using GRV 1998 as unpolarized reference set
C
C***********************************************************************
      ELSEIF ((LST(15).GE.144).AND.(LST(15).LE.147)) THEN    

      IF(X.LE.0.0001)  THEN
      WRITE(*,*) 'Warning: x exceeds minum value, x set to 0.0001001'
      X=0.0001001
      ENDIF
      IF(X.GE.1.0)     THEN
      WRITE(*,*) 'Warning: x exceeds maxinum value, x set to 0.99999'
      X=0.999999
      ENDIF
      IF(Q2.LE.0.8)     THEN
      WRITE(*,*)'Warning: Q2 exceeds minimum value, Q2 set to 0.80001'
      Q2 = 0.80001
      ENDIF
      IF(Q2.GE.1.0E6) THEN
      WRITE(*,*) 'Warning: Q2 exceeds maximum value, Q2 set to 999999.9'
      Q2 = 999999.9
      ENDIF 

      
C  unpolarized part 
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *

         DX=DBLE(X)
         DQ2=DBLE(Q2)
         IF ((LST(15).GE.144).AND.(LST(15).LE.145)) THEN
            ISET = 1 ! LO
         ELSEIF ((LST(15).GE.146).AND.(LST(15).LE.147)) THEN
            ISET = 2 ! NLO
         ENDIF
         CALL GRV98PA (ISET, DX, DQ2 
     &           ,DDXUV,DDXDV,DDXUS, DDXDS, DDXSB,DDXGL)

  
         XPQ( 0) = SNGL(DDXGL)
         XPQ( 1) = SNGL(DDXDV  + DDXDS)
         XPQ(-1) = SNGL(DDXDS)         
         XPQ( 2) = SNGL(DDXUV  + DDXUS)   
         XPQ(-2) = SNGL(DDXUS)         
         XPQ( 3) = SNGL(DDXSB) 
         XPQ(-3) = SNGL(DDXSB) 

*     polarized part


      IF (LST(15).EQ.144) ISET = 3
      IF (LST(15).EQ.145) ISET = 4
      IF (LST(15).EQ.146) ISET = 1
      IF (LST(15).EQ.147) ISET = 2



         CALL PARPOLNEW (ISET, DX, DQ2, 
     &          DXDU, DXDD, DXDUB, DXDDB, DXDSTR, DXDGLU, 
     &          DUMMYG1P, DUMMYG1N)
         XDPQ( 0) = SNGL(DXDGLU)
         XDPQ( 1) = SNGL(DXDD)
         XDPQ(-1) = SNGL(DXDDB)           
         XDPQ( 2) = SNGL(DXDU)
         XDPQ(-2) = SNGL(DXDUB)           
         XDPQ( 3) = SNGL(DXDSTR)           
         XDPQ(-3) = SNGL(DXDSTR)

C*************************************************************************
C        LST(15) = 150   cteq5l                   LO
C        LST(15) = 151   cteq5m                   NLO MSBAR
C        LST(15) = 152   cteq5m1                  NLO MSBAR (update)
C*************************************************************************
  
         ELSEIF((LST(15).GE.150).AND.(LST(15).LE.152)) THEN

         IF(X.LE.0.00001)  THEN
c$$$       WRITE(*,*) 'Warning: x exceeds minum value, x set to 0.0000101'
            X=0.0000101
         ENDIF
         IF(X.GE.1.0)     THEN
c$$$       WRITE(*,*) 'Warning: x exceeds maxinum value, x set to 0.99999'
            X=0.999999
         ENDIF
         IF(Q2.LE.1.0)     THEN
c$$$       WRITE(*,*)'Warning: Q2 exceeds minimum value, Q2 set to 1.00001'
            Q2 = 1.00001
         ENDIF
         IF(Q2.GE.1.0E4) THEN
c$$$       WRITE(*,*) 'Warning: Q2 exceeds maximum value, Q2 set to 9999.9'
            Q2 = 9999.9
         ENDIF 

        XPQ(-5) = X*SNGL( CTQ5PDF (-5, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-4) = X*SNGL( CTQ5PDF (-4, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-3) = X*SNGL( CTQ5PDF (-3, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-2) = X*SNGL( CTQ5PDF (-1, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-1) = X*SNGL( CTQ5PDF (-2, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 0) = X*SNGL( CTQ5PDF ( 0, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 1) = X*SNGL( CTQ5PDF ( 2, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 2) = X*SNGL( CTQ5PDF ( 1, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 3) = X*SNGL( CTQ5PDF ( 3, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 4) = X*SNGL( CTQ5PDF ( 4, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 5) = X*SNGL( CTQ5PDF ( 5, DBLE(X), DBLE(SQRT(Q2))))

C*****************************************************************************
C        LST(15) = 161  mrs99   cor01  central gluon, a_s    
C        LST(15) = 162  mrs99   cor02  higher gluon        
C        LST(15) = 163  mrs99   cor03  lower gluon       
C        LST(15) = 164  mrs99   cor04  lower a_s      
C        LST(15) = 165  mrs99   cor05  higher a_s     
C        LST(15) = 166  mrs99   cor06  quarks up        
C        LST(15) = 167  mrs99   cor07  quarks down            
C        LST(15) = 168  mrs99   cor08  strange up           
C        LST(15) = 169  mrs99   cor09  strange down          
C        LST(15) = 170  mrs99   cor10  charm up         
C        LST(15) = 171  mrs99   cor11  charm down       
C        LST(15) = 172  mrs99   cor12  larger d/u          
C****************************************************************************

           ELSEIF((LST(15).GE.161).AND.(LST(15).LE.172)) THEN

         IF(X.LE.0.00001)  THEN
c$$$       WRITE(*,*) 'Warning: x exceeds minum value, x set to 0.0000101'
            X=0.0000101
         ENDIF
         IF(X.GE.1.0)     THEN
c$$$       WRITE(*,*) 'Warning: x exceeds maxinum value, x set to 0.99999'
            X=0.999999
         ENDIF
         IF(Q2.LE.1.25)     THEN
c$$$        WRITE(*,*)'Warning: Q2 exceeds minimum value, Q2 set to 1.25001'
            Q2 = 1.25001
         ENDIF
        IF(Q2.GE.1.0E4) THEN
c$$$       WRITE(*,*) 'Warning: Q2 exceeds maximum value, ',
c$$$     &                          'Q2 set to 9999999.9'
            Q2 = 9999999.9
         ENDIF 

        MODE = LST(15)-160
        CALL MRS99(DBLE(X),DBLE(SQRT(Q2)),MODE,
     &       MRSUPV,MRSDNV,MRSUSEA,MRSDSEA,MRSSTR,MRSCHM,MRSBOT,
     &       MRSGLU) 

        XPQ(-5) = SNGL(MRSBOT)
        XPQ(-4) = SNGL(MRSCHM)
        XPQ(-3) = SNGL(MRSSTR)
        XPQ(-2) = SNGL(MRSUSEA)
        XPQ(-1) = SNGL(MRSDSEA)
        XPQ( 0) = SNGL(MRSGLU)
        XPQ( 1) = SNGL(MRSDSEA)+SNGL(MRSDNV)
        XPQ( 2) = SNGL(MRSUSEA)+SNGL(MRSUPV)
        XPQ( 3) = SNGL(MRSSTR)
        XPQ( 4) = SNGL(MRSCHM)
        XPQ( 5) = SNGL(MRSBOT)

C Added By Joe Seele for the inclusion of the CTEQ 6 PDF's

      ELSEIF((LST(15).GE.173).AND.(LST(15).LE.175)) THEN

         IF(X.LE.0.00001)  THEN
c$$$       WRITE(*,*) 'Warning: x exceeds minum value, x set to 0.0000101'
            X=0.0000101
         ENDIF
         IF(X.GE.1.0)     THEN
c$$$       WRITE(*,*) 'Warning: x exceeds maxinum value, x set to 0.99999'
            X=0.999999
         ENDIF
         IF(Q2.LE.1.0)     THEN
c$$$       WRITE(*,*)'Warning: Q2 exceeds minimum value, Q2 set to 1.00001'
            Q2 = 1.00001
         ENDIF
         IF(Q2.GE.1.0E4) THEN
c$$$       WRITE(*,*) 'Warning: Q2 exceeds maximum value, Q2 set to 9999.9'
            Q2 = 9999.9
         ENDIF 

        XPQ(-5) = X*SNGL( CTQ6PDF (-5, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-4) = X*SNGL( CTQ6PDF (-4, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-3) = X*SNGL( CTQ6PDF (-3, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-2) = X*SNGL( CTQ6PDF (-1, DBLE(X), DBLE(SQRT(Q2))))
        XPQ(-1) = X*SNGL( CTQ6PDF (-2, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 0) = X*SNGL( CTQ6PDF ( 0, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 1) = X*SNGL( CTQ6PDF ( 2, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 2) = X*SNGL( CTQ6PDF ( 1, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 3) = X*SNGL( CTQ6PDF ( 3, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 4) = X*SNGL( CTQ6PDF ( 4, DBLE(X), DBLE(SQRT(Q2))))
        XPQ( 5) = X*SNGL( CTQ6PDF ( 5, DBLE(X), DBLE(SQRT(Q2))))

C*****************************************************************************




C*************************************************************************
      ELSE
         WRITE(*,*) 'NOT KNOWN PARTON DENSITY FUNCTION: STOP',LST(15)
         STOP
      ENDIF
*72*********************************************************************
      DO I=-6,6
        IF(ABS(XDPQ(I)).gt.XPQ(I)) THEN
          IF(ICERR.lt.IPLST(1).and.Q2.gt.1.) THEN 
             WRITE(*,*)'Non proper parton set : ',XDPQ(I),XPQ(I),I,X,Q2
             ICERR=ICERR+1
          ENDIF
          IF(XDPQ(I).gt.0) THEN 
             XDPQ(I)=XPQ(I)
          ELSE
             XDPQ(I)=-XPQ(I)
          ENDIF
       ENDIF
      ENDDO      
*72*********************************************************************      
C   This Amendment is to simulate HeliumIII
      IF(LST(39).NE.0) THEN
         F2P = 4.*(XPQ(2)+XPQ(-2)+XPQ(4)+XPQ(-4))
     &       +XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3)
         IF(LST(36).eq.0) THEN 
            XDV=XPQ(1)-XPQ(-1)
            XUV=XPQ(2)-XPQ(-2)
            F2N = 4.*(XDV+2*XPQ(-2)+XPQ(4)+XPQ(-4))
     &          + XUV+2*XPQ(-1)+XPQ(3)+XPQ(-3)
         ELSE
            F2N = 4.*(XPQ(1)+XPQ(-1)+XPQ(4)+XPQ(-4))+
     &           XPQ(2)+XPQ(-2)+XPQ(3)+XPQ(-3)
         ENDIF
         
C     admitxurecoefficient for the unpolarized structure functions 
         PARI11 = 1./3.
 
         IF(LST(39).EQ.3) THEN
C..........................
C  This is a real helium 3.
C..........................
C Percentage of neutron and proton to the helium3 asymmetry
            PERCS   =  0.865
            PERCSP  =  0.
            PERCD   = -0.027
C admixture coefficients for the spin states
            ADMXS   = 1.0
            ADMXSP  = 0.
            ADMXD   = 0.
C factors for the spin state (number of polarized nucleons/total number
C of nucleons
            FACTS   = f2n/(f2n+2*f2p)
            FACTSP  = 0.
            FACTD   = 2*f2p/(f2n+2*f2p)
         ELSEIF(LST(39).EQ.2) THEN
C..........................
C  This is a real helium 3.
C..........................
C Percentage of the spin states
            PERCS   =  0.903
            PERCSP  =  0.014
            PERCD   = -0.083
C admixture coefficients for the spin states
            ADMXS   = 1.0
            ADMXSP  = 1./3.
            ADMXD   = 1./3.
C factors for the spin state (number of polarized nucleons/total number
C of nucleons
            FACTS   = f2n/(f2n+2*f2p)
            FACTSP  = 1./3.
            FACTD   = 1. 
C     admitxurecoefficient for the unpolarized structure functions 
            PARI11 = 1./3.
         ELSEIF(LST(39).EQ.1) THEN
C...........................
C  This is a ideal helium 3.
C...........................
C     Percentage of the spin states
            PERCS   = 1.0
            PERCSP  = 0.0
            PERCD   = 0.0
C     admixturecoefficients for the spin states
            ADMXS   = 1.
            ADMXSP  = 0.
            ADMXD   = 0.
C     factors for the spin state (number of polarized nucleons/total number
C     of nucleons
            FACTS   = f2n/(f2n+2*f2p)
            FACTSP  = 0.
            FACTD   = 0. 
C     admitxurecoefficient for the unpolarized structure functions
            PARI11 = 1./3.
         ENDIF
         
C...  For nuclear target, mix u- and d-valence distributions.
         IF(LST(36).EQ.0) THEN
C...  For nuclear target, mix u- and d-valence distributions.
            XPQ(1)=(1.-PARI11)*XDV+PARI11*XUV+XPQ(-1)
            XPQ(2)=(1.-PARI11)*XUV+PARI11*XDV+XPQ(-2)
C     Admixture in the polarized case         
            XDDV=XDPQ(1)-XDPQ(-1)
            XDUV=XDPQ(2)-XDPQ(-2)
            XDPQ(1)= ((1.-ADMXS )*XDDV+ADMXS *XDUV )
     &           * FACTS * PERCS
     &           +((1.-ADMXSP)*XDDV+ADMXSP*XDUV )
     &           * FACTSP* PERCSP 
     &           +((1.-ADMXD )*XDDV+ADMXD *XDUV )
     &           * FACTD * PERCD
            XDPQ(2)= ((1.-ADMXS )*XDUV+ADMXS *XDDV )
     &           * FACTS * PERCS
     &           +((1.-ADMXSP)*XDUV+ADMXSP*XDDV )
     &           * FACTSP* PERCSP
     &           +((1.-ADMXD )*XDUV+ADMXD *XDDV )
     &           * FACTD * PERCD
            DO I=-6,6
C...  For nuclear target, mix u- and d-valence distributions.
               IF(I.eq.1) THEN
* XDPQ(-1) has to be diluted here by (FACTS*PERCS+FACTSP...
                  XDPQ(1)= XDPQ(1) + XDPQ(-1)
               ELSEIF(I.eq.2) THEN
* XDPQ(-2) has to be diluted here by (FACTS*PERCS+FACTSP...                  
                  XDPQ(2)= XDPQ(2) + XDPQ(-2)
               ELSE
*     now for the sea quarks and gluon distributions                  
                  XDPQ(I)=XDPQ(I)*
     &                 (FACTS*PERCS+FACTSP*PERCSP+FACTD*PERCD)
               ENDIF
            ENDDO            
         ELSEIF(LST(36).EQ.1) THEN
            XD=XPQ(1)
            XU=XPQ(2)
            XPQ(1) =(1.-PARI11)*XD+PARI11*XU
            XPQ(2) =(1.-PARI11)*XU+PARI11*XD
            XD=XPQ(-1)
            XU=XPQ(-2)
            XPQ(-1)=(1.-PARI11)*XD+PARI11*XU
            XPQ(-2)=(1.-PARI11)*XU+PARI11*XD
C     Admixture in the polarized case
*.............................................
*           FOR u and d QUARKS
*.............................................         
            XDD=XDPQ(1)
            XDU=XDPQ(2)
            XDDS=XDPQ(-1)
            XDUS=XDPQ(-2)            
            DO I=-6,6
C...  For nuclear target, mix u- and d-valence distributions.
               IF(I.eq.1) THEN 
                  XDPQ(1)= ((1.-ADMXS )*XDD+ADMXS *XDU )
     &                 * FACTS * PERCS  
     &                 +((1.-ADMXSP)*XDD+ADMXSP*XDU )
     &                 * FACTSP* PERCSP 
     &                 +((1.-ADMXD )*XDD+ADMXD *XDU )
     &                 * FACTD * PERCD
               ELSEIF(I.eq.2) THEN 
                  XDPQ(2)= ((1.-ADMXS )*XDU+ADMXS *XDD )
     &                 * FACTS * PERCS
     &                 +((1.-ADMXSP)*XDU+ADMXSP*XDD )
     &                 * FACTSP* PERCSP
     &                 +((1.-ADMXD )*XDU+ADMXD *XDD )
     &                 * FACTD * PERCD
*     now for the sea quarks and gluon distributions
               ELSEIF(I.eq.-1) THEN 
                  XDPQ(-1)=((1.-ADMXS)*XDDS+ADMXS*XDUS)
     &                 * FACTS * PERCS  
     &                 +((1.-ADMXSP)*XDDS+ADMXSP*XDUS )
     &                 * FACTSP* PERCSP 
     &                 +((1.-ADMXD )*XDDS+ADMXD *XDUS )
     &                 * FACTD * PERCD
               ELSEIF(I.eq.-2) THEN
                  XDPQ(-2)=((1.-ADMXS)*XDUS+ADMXS*XDDS )
     &                 * FACTS * PERCS
     &                 +((1.-ADMXSP)*XDUS+ADMXSP*XDDS )
     &                 * FACTSP* PERCSP
     &                 +((1.-ADMXD )*XDUS+ADMXD *XDDS )
     &                 * FACTD * PERCD
               ELSE
*           now for the sea quarks and gluon distributions                  
                  XDPQ(I)=XDPQ(I)*
     &                 (FACTS*PERCS+FACTSP*PERCSP+FACTD*PERCD)
               ENDIF
            ENDDO
C...  For nuclear target, mix u- and d-valence distributions.
         ENDIF
      ENDIF
         
* polarisation enhancement factor 
      DO I=-6,6
         XDPQ(I) = FLOAT(LST(38))* XDPQ(I)
      ENDDO
      
      RETURN
      END
 



C************************************************************************
      FUNCTION GAMMA(XX)
      REAL GAMMAGAMMLN,XX
      INTEGER J
      REAL SER,STP,TMP,X,Y,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146E0,-86.50532032941677E0,
     & 24.01409824083091E0,-1.231739572450155E0,.1208650973866179E-2,
     & -.5395239384953E-5,2.5066282746310005E0/
      X=XX
      Y=X
      TMP=X+5.5E0
      TMP=(X+0.5E0)*LOG(TMP)-TMP
      SER=1.000000000190015E0
      DO J=1,6
         Y=Y+1.E0
         SER=SER+COF(J)/Y
      ENDDO
      GAMMLN=TMP+LOG(STP*SER/X)
      GAMMA = EXP(GAMMLN)
      RETURN
      END


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                   *
*     G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S      *
*                                                                   *
*                          1998 UPDATE                              *
*                                                                   *
*                  For a detailed explanation see                   *
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *
*                                                                   *
*   This package contains subroutines returning the light-parton    *
*   distributions in NLO (for the MSbar and DIS schemes) and LO;    * 
*   the respective light-parton, charm, and bottom contributions    *
*   to F2(electromagnetic); and the scale dependence of alpha_s.    *
*                                                                   *
*   The parton densities and F2 values are calculated from inter-   *
*   polation grids covering the regions                             *
*         Q^2/GeV^2  between   0.8   and  1.E6 ( 1.E4 for F2 )      *
*            x       between  1.E-9  and   1.                       *
*   Any call outside these regions stops the program execution.     *
*                                                                   *
*   At Q^2 = MZ^2, alpha_s reads  0.114 (0.125) in NLO (LO); the    *
*   heavy quark thresholds, QH^2 = mh^2, in the beta function are   *
*            mc = 1.4 GeV,  mb = 4.5 GeV,  mt = 175 GeV.            *
*   Note that the NLO alpha_s running is different from GRV(94).    * 
*                                                                   *
*    Questions, comments etc to:  avogt@physik.uni-wuerzburg.de     *
*                                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
      SUBROUTINE GRV98PA (ISET, X, Q2, UV, DV, US, DS, SS, GL)
*********************************************************************
*                                                                   *
*   THE PARTON ROUTINE.                                             *
*                                     __                            *
*   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         *
*                                                                   *
*   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       *
*            DS = d(bar),  SS = s = s(bar),  GL = gluon.            *
*            Always x times the distribution is returned.           *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIP / IINIP , and the     *
*            integer variable  IINIP  has always to be zero when    *
*            GRV98PA is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   *
*            (1+1809 lines with 6 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      COMMON / INTINIP / IINIP
      SAVE XUVF, XDVF, XDEF, XUDF, XSF, XGF, NA, ARRF
*
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8E0, 
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4, 2.2E4, 4.6E4,
     6           1.0E5, 2.2E5, 4.6E5, 
     7           1.E6 /
       DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9, 
     A           1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8, 
     B           1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7, 
     C           1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1           1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2           1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3           1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4           1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*
*...CHECK OF X AND Q2 VALUES : 
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
         STOP
      ENDIF
      IF (IINIP .NE. 0) GOTO 16
*
*...INITIALIZATION, IF REQUIRED :
*
*    SELECTION AND READING OF THE GRID : 
*    (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
        OPEN (11,FILE='./pdf/grv98lo.grid',STATUS='old')   ! 7.332E-05
      ELSE IF (ISET .EQ. 2) THEN
        OPEN (11,FILE='./pdf/grv98nlm.grid',STATUS='old')  ! 1.015E-04
      ELSE IF (ISET .EQ. 3) THEN
        OPEN (11,FILE='./pdf/grv98nld.grid',STATUS='old')  ! 1.238E-04
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'NO OR INVALID PARTON SET CHOICE')
        STOP
      END IF
      IINIP = 1
      READ(11,89) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(11,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1            PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M) 
  90  FORMAT (6(1PE10.3))
  15  CONTINUE
      CLOSE(11)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0V = XB(IX)**0.5 
        XB0S = XB(IX)**(-0.2) 
        XB1 = 1.-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
        XDEF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V) 
        XUDF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
  20  CONTINUE
        XUVF(NX,IQ) = 0.E0
        XDVF(NX,IQ) = 0.E0
        XDEF(NX,IQ) = 0.E0
        XUDF(NX,IQ) = 0.E0
        XSF(NX,IQ)  = 0.E0
        XGF(NX,IQ)  = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
*
*...CONTINUATION, IF INITIALIZATION WAS DONE PREVIOUSLY.
*
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      X1 = 1.- X
      XV = X**0.5
      XS = X**(-0.2)
      UV = FINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = FINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = FINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = FINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      US = 0.5 * (UD - DE)
      DS = 0.5 * (UD + DE)
      SS = FINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = FINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS 
*
 60   RETURN
      END
*
*
*

*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
*
*
*
*
      FUNCTION ALPHAS (Q2, NAORD)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE.                                            *
*                                                                   *
*   INPUT :  Q2    =  scale in GeV**2  (not too low, of course);    *
*            NAORD =  1 (LO),  2 (NLO).                             *
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) for use with the GRV(98) partons.      *  
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A - Z)
      INTEGER NF, K, I, NAORD
      DIMENSION LAMBDAL (3:6),  LAMBDAN (3:6), Q2THR (3)
*
*...HEAVY QUARK THRESHOLDS AND LAMBDA VALUES :
      DATA Q2THR   /  1.960,  20.25,  30625. /
      DATA LAMBDAL / 0.2041, 0.1750, 0.1320, 0.0665 /
      DATA LAMBDAN / 0.2994, 0.2460, 0.1677, 0.0678 /
*
*...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
      NF = 3
      DO 10 K = 1, 3
      IF (Q2 .GT. Q2THR (K)) THEN
         NF = NF + 1
      ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
         ALP  = 1./(B0 * DLOG (Q2/LAM2))
         GO TO 1
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
         WRITE (6,91)
  91     FORMAT ('INVALID CHOICE FOR ORDER IN ALPHA_S')
         STOP
       END IF
*
*...START VALUE FOR NLO ITERATION :
       LQ2 = DLOG (Q2 / LAM2)
       ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)
*
*...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE :
       DO 2 I = 1, 6
       XL  = DLOG (1./(B0*ALP) + B10)
       XLP = DLOG (1./(B0*ALP*1.01) + B10)
       XLM = DLOG (1./(B0*ALP*0.99) + B10)
       Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
       Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP
     1       + 1./ (B0*ALP*0.99) - B10 * XLP) / (0.02D0*ALP)
       ALP = ALP - Y/Y1
  2    CONTINUE
*
*...OUTPUT :
  1    ALPHAS = ALP
       RETURN
       END
*
*********************************************************************
*                                                                   *
*    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    *
*                                                                   *
*         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         *
*                        hep-ph/0011215                             *
*                                                                   *
*          PROBLEMS/QUESTIONS TO wvogelsang@bnl.gov                 *
*            OR TO marco.stratmann@physik.uni-regensburg.de         *
*                                                                   *
*   INPUT:   ISET = number of the parton set :                      *
*            ISET = 1  'STANDARD' SCENARIO, NEXT-TO-LEADING ORDER   *
*                      (MS-bar)                                     * 
*                      (DATA FILE 'std2000_nlo.grid' UNIT=11, TO BE *
*                       DEFINED BY THE USER )                       *
*            ISET = 2  'VALENCE' SCENARIO,  NEXT-TO-LEADING ORDER   *
*                      (MS-bar)                                     *   
*                      (DATA FILE 'val2000_nlo.grid' UNIT=22, TO BE *
*                       DEFINED BY THE USER )                       *
*            ISET = 3  'STANDARD' SCENARIO, LEADING ORDER           *
*                      (DATA FILE 'std2000_lo.grid' UNIT=33, TO BE  *
*                       DEFINED BY THE USER )                       *
*            ISET = 4  'VALENCE' SCENARIO,  LEADING ORDER           *
*                      (DATA FILE 'val2000_lo.grid' UNIT=44, TO BE  *
*                       DEFINED BY THE USER )                       *
*                                                                   *
*            X  = Bjorken-x       (between  1.E-4  and  1)          *
*            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        *
*                                                                   *
*   OUTPUT:  U = x * DELTA u                                        *
*            D = x * DELTA d                                        *        
*            UB = x * DELTA ubar                                    *   
*            DB = x * DELTA dbar                                    * 
*            ST = x * DELTA STRANGE                                 *     
*            GL = x * DELTA GLUON                                   *
*            G1P = g_1^proton                                       *
*            G1N = g_1^neutron                                      * 
*                                                                   *
*          (  For the parton distributions always x times           *
*                   the distribution is returned .                  *
*                 This is NOT the case for g1(p,n)  )               *
*                                                                   *
*            The sets are the result of a combined fit to           *
*            data for the spin asymmetries A_1 (p,n,d)              *
*                                                                   *
*            Note: No charm is included                             *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be zero when PARPOL is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
*********************************************************************
*
      SUBROUTINE PARPOLNEW (ISET, X, Q2, U, D, UB, DB, ST,GL,G1P,G1N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPART=8, NX=42, NQ=30, NARG=2)
      DIMENSION XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), XDBF(NX,NQ), 
     1          XSF(NX,NQ), XGF(NX,NQ), XG1P(NX,NQ), XG1N(NX,NQ),
     2          PARTON (NPART,NQ,NX-1), QS(NQ), XB(NX), XT(NARG), 
     3          NA(NARG), ARRF(NX+NQ) 
      COMMON / INTINI / IINI
      SAVE XUF, XDF, XUBF, XDBF, XSF, XGF, XG1P, XG1N, NA, ARRF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8D0, 1.0D0, 1.25d0, 1.5D0, 2.d0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 
     4           1.0D5, 1.8D5, 3.2D5, 5.8D5, 1.0D6  /
       DATA XB / 
     1           1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     2           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     3           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     4           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5           0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6           0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9, 1.0 /
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.1.0D-4) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
           STOP
c          GOTO 60
       ENDIF
       IF ( (Q2.LT.0.8D0) .OR. (Q2.GT.1.D6) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
           STOP
C          GOTO 60
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    FILE - NO. = 11 FOR NLO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.3478E-03 )
*    FILE - NO. = 22 FOR NLO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.5146E-05 )
*    FILE - NO. = 33 FOR  LO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 3.4686E-03 )     
*    FILE - NO. = 44 FOR  LO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 2.4395E-04 )
      IF (IINI.NE.0) GOTO 16
      IF (ISET.EQ.1) THEN
       IIREAD=11       
       OPEN(UNIT=11,FILE='./pdf/std2000_nlo_g1.grid',STATUS='OLD')
      ELSE IF (ISET.EQ.2) THEN
       IIREAD=22
       OPEN(UNIT=22,FILE='./pdf/val2000_nlo_g1.grid',STATUS='OLD')
      ELSE IF (ISET.EQ.3) THEN
       IIREAD=33       
       OPEN(UNIT=33,FILE='./pdf/std2000_lo_g1.grid',STATUS='OLD')
      ELSE IF (ISET.EQ.4) THEN
       IIREAD=44
       OPEN(UNIT=44,FILE='./pdf/val2000_lo_g1.grid',STATUS='OLD')
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE')
        GOTO 60
      END IF
C
       DO 15 M = 1, NX-1
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M)
  90   FORMAT (8(1PE12.4))
  15   CONTINUE
C
      IINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0)
        XDF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0)
        XUBF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**8 * XB0**0.5) 
        XDBF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**8 * XB0**0.5) 
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**8 * XB0**0.5) 
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0**2.)
        XG1P(IX,IQ)  = PARTON(7,IQ,IX) / XB1**3
        XG1N(IX,IQ)  = PARTON(8,IQ,IX) / XB1**3
  20  CONTINUE
        XUF(NX,IQ) = 0.D0
        XDF(NX,IQ) = 0.D0
        XUBF(NX,IQ) = 0.D0
        XDBF(NX,IQ) = 0.D0
        XSF(NX,IQ)  = 0.D0
        XGF(NX,IQ)  = 0.D0
        XG1P(NX,IQ)  = 0.D0
        XG1N(NX,IQ)  = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      U = DFINTNEW(NARG,XT,NA,ARRF,XUF) * (1.D0-X)**3 * X
      D = DFINTNEW(NARG,XT,NA,ARRF,XDF) * (1.D0-X)**4 * X
      UB = DFINTNEW(NARG,XT,NA,ARRF,XUBF) * (1.D0-X)**8 * X**0.5
      DB = DFINTNEW(NARG,XT,NA,ARRF,XDBF) * (1.D0-X)**8 * X**0.5
      ST = DFINTNEW(NARG,XT,NA,ARRF,XSF)  * (1.D0-X)**8 * X**0.5
      GL = DFINTNEW(NARG,XT,NA,ARRF,XGF)  * (1.D0-X)**5 * X**2.
      G1P = DFINTNEW(NARG,XT,NA,ARRF,XG1P)  * (1.D0-X)**3
      G1N = DFINTNEW(NARG,XT,NA,ARRF,XG1N)  * (1.D0-X)**3
 60   RETURN
      END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION DFINTNEW(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(5),NENT(5),ENT(73),TABLE(1200)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      DFINTNEW=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      DFINTNEW=DFINTNEW+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
