      SUBROUTINE UnFoldit(rawmat, calib, R,
     &                    unfmat, iVersion)
C A program for unfolding detector response effects from a continuos
C gamma ray spectrum. Modified version (Aug. 1995) which smoothes the 
C Compton contribution, which is to be subtracted.
C
C UCID - 30079
C Lawrence Livermore Laboratory
C William J. Oconnel September 1973
C
C UCRL - 9748
C A computer analysis for complex sodium iodide gamma spectra
c James F. Mollenauer   august 1961
C
C Oslo Cyclotron Laboratory version 1989,1993,1994,1995 M.G.
C See NIM paper of M. Guttormsen et al. 1996, in press
C        Raw    Raw spectrum
C        U      Unfolded spectrum
C        F      Folded spectrum, F=R(I,J)*U
C        a0,a1  Calibration coeff in units of keV and keV/ch
C        LEN    Channels per spectrum
C        mode   = -1 : Difference iteration (d)
C               = +1 : Ratio iteration      (r)
C        Iter   =Max number of iteration steps

CJEM  Allocation and f2py declarations of new variables:
      double precision rawmat(0:4095, 0:2047), unfmat(0:4095, 0:2047)
      double precision calib(2)
      double precision R(0:2047,0:2047)
Cf2py double precision intent(in) :: rawmat, calib, R
Cf2py double precision intent(out) :: unfmat

      INTEGER XDIM,YDIM,RDIM
      CHARACTER APP*4,waitname*8
CJEM      COMMON/Sp2Dim/rMAT(2,0:4095,0:2047),APP(2048),XDIM,YDIM
      COMMON/Sp2Dim/APP(2048),XDIM,YDIM
      COMMON/State/Istatus,ITYPE,IDEST,cal(2,2,2,3),Idim(2,2,2),
     & fname(2,2),comm(2,2)
      CHARACTER fname*8,comm*60
      COMMON/Sp1Dim/rSPEC(2,0:8191),MAXCH

CJEM      COMMON/response1/R(0:2047,0:2047),RDIM,a0,a1,FWHM,facFWHM
      COMMON/response1/RDIM,a0,a1,FWHM,facFWHM
      double precision a0, a1
      COMMON/response3/EffTot(0:2047),Fwhm1(0:2047),EffExp(0:2047)
      COMMON/response4/pf(0:2047),pc(0:2047),ps(0:2047),pd(0:2047),
     & pa(0:2047)

      INTEGER low,high,ChiLow,ChiHigh,lower,upper
      INTEGER sum,Fsum,sumCH,FsumCH
      DIMENSION Raw(0:2047),F(0:2047),U(0:2047),lower(0:2047),
     & upper(0:2047)
      DIMENSION Fsave(0:2047),Usave(0:2047)
      DIMENSION SoluF(200,0:2047),SoluP(200,7)
      DIMENSION rWAIT(0:8191)
      CHARACTER SCORE(200)*4,ANS*1
      CHARACTER modus*2,moil*2,moim*2,moir*2
      INTEGER COLORMAP(20),Colorc(0:19)
      COMMON /COLORMAP/ COLORMAP,Limit,Colorc
      REAL Limit(0:19)
      LOGICAL DISP
      INTEGER            IYAXIS,LDX,HDX,LDY,HDY,LOCH,HICH
      COMMON/DISPLA/DISP,IYAXIS,LDX,HDX,LDY,HDY,LDZ,HDZ,LOCH,HICH,LOCNT,
     & HICNT
      REAL LDZ,HDZ,LOCNT,HICNT
      DIMENSION u0(0:2047),su0(0:2047),v(0:2047),us(0:2047),ud(0:2047)
      DIMENSION ua(0:2047),sua(0:2047),w(0:2047),sw(0:2047),c(0:2047),
     & sc(0:2047)

CJEM       ISP=1
CJEM      IF(IDEST.EQ.1)ISP=2
      Emin  =30.          ! 30 keV lowest limit
      EminCh=100.         ! 100 keV lower limit for Chi2-test
      Iter  =33           ! Number of max. iterations

CJEM  Addition: Input of a0 and a1 from argument to function:
      a0 = calib(1)
      a1 = calib(2)  
      ITYPE = 2
CJEM      write(6,*)'ITYPE right after set =',ITYPE
      XDIM = 1600
      YDIM = 500

C Using rSPEC(IDEST,i) to display spectra. Putting back later
      DO i=0,MAXCH
        rWAIT(i)=rSPEC(IDEST,i)
        rSPEC(IDEST,i)=0
      ENDDO
      waitname=fname(2,IDEST)
      IF(ITYPE.GT.1)fname(2,IDEST)=fname(1,ISP)
      IF(ITYPE.EQ.1)fname(2,IDEST)=fname(2,ISP)

C Zeroing destination spectrum
CJEM      IF(ITYPE.GT.1)THEN                !matrix
        LEN=XDIM
CJEM        IF(RDIM.LT.XDIM)LEN=RDIM
        Iydim=YDIM
        DO I=0,4095
          DO J=0,2047
CJEM            rMAT(IDEST,I,J)=0
            unfmat(I,J)=0
          ENDDO
        ENDDO
CJEM      ELSE                              !singles
CJEM        LEN=MAXCH+1
CJEM        IF(RDIM.LT.LEN)LEN=RDIM
CJEM        Iydim=1
CJEM        DO I=0,8191
CJEM          rSPEC(IDEST,I)=0
CJEM        ENDDO
CJEM      ENDIF

C Lowest channel treated
      low=((Emin-a0)/a1)+0.5
      write(6,*)'low=',low
      IF(low.LT.0.OR.low.GE.LEN)low=0
      lowCh=((EminCh-a0)/a1)+0.5
      IF(lowCh.LT.0.OR.lowCh.GE.LEN)lowCh=0

C Defining upper border for the unfolding and chisq-test
      Ix1=LEN-1
      Ix2=LEN-1
      Iy1=0
      Iy2=Iydim-1
  
CJEM      IF(ITYPE.GT.1)THEN                !matrix
CJEM        WRITE(6,*)'Give upper limits for the unfolding. The boarder is'
CJEM        WRITE(6,*)'given by interpolation between (x1,y1) and (x2,y2)'
CJEM        WRITE(6,*)' '
CJEM        WRITE(6,*)'    (x2,y2) second point'
CJEM        WRITE(6,*)'xxxxxxx'
CJEM        WRITE(6,*)'xxxxxxxxx'
CJEM        WRITE(6,*)'xxxxxxxxxxx'
CJEM        WRITE(6,*)'xx matrix xxx'
CJEM        WRITE(6,*)'xxxxxxxxxxxxxxx'
CJEM        WRITE(6,*)'             (x1,y1) first point'
CJEM        WRITE(6,*)' '
CJEM
CJEM        WRITE(6,123)Ix1
CJEM 123    FORMAT(/'First point x1  <',I5,'>:',$)
CJEMCJEM        CALL READI(5,Ix1)
CJEM        WRITE(6,124)Iy1
CJEM 124    FORMAT( 'First point y1  <',I5,'>:',$)
CJEMCJEM        CALL READI(5,Iy1)
CJEM        WRITE(6,125)Ix2
CJEM 125    FORMAT( 'Second point x2 <',I5,'>:',$)
CJEMCJEM        CALL READI(5,Ix2)
CJEM        WRITE(6,126)Iy2
CJEM 126    FORMAT( 'Second point y2 <',I5,'>:',$)
CJEMCJEM        CALL READI(5,Iy2)
CJEM      ELSE                              !singles
CJEM        WRITE(6,12)Ix1
CJEM 12     FORMAT(/'Give upper channel for unfolding  <',I5,'>:',$)
CJEMCJEM        CALL READI(5,Ix1)
CJEM        Ix2=Ix1
CJEM      ENDIF
      IF(Istatus.NE.0)RETURN

CJEM  TODO make Ix1, Iy1, Ix2, Iy2 input arguments (tuple?)
      Ix1 = 100
      Iy1 = 10
      Ix2 = 800
      Iy2 = 499

      CF=0.
      Dx12=Ix2-Ix1
      Dy12=Iy2-Iy1
      IF(Iy2.NE.Iy1)CF=Dx12/Dy12
      DO J=0,Iydim-1
        upper(J)=Ix1-CF*(FLOAT(Iy1-J))+0.5
CJEM        write(6,*)'upper(',J,')=',upper(J)
        IF(upper(J).LT.low  )upper(J)=low
        IF(upper(J).GT.LEN-1)upper(J)=LEN-1
      ENDDO

C Setting lower limits for the chisquare test
CJEM      WRITE(6,130)
CJEM 130  FORMAT('  Give limits for the chisquare-test:',/,
CJEM     +       '  Opt. 1: Recommended for LaBr- and NaI-spll-',
CJEM     +/,     '          energy gammas above 2 MeV, we se at 500 keV.',
CJEM     +/,     '           Below, the limit is 1/4 of the Remember,'
CJEM     +/,     '          full-energy is taken from the up limit',
CJEM     +/,     '  Opt. 2: A fixed lower limit for the chi-ed'
CJEM     +/,     '  Opt. 3: Return and set proper upper limiing',/)
CJEM

CJEM  TODO Give Iopt as input argument (with default?)
      Iopt=1
CJEM      IF(Ix1.EQ.Ix2)Iopt=2
CJEM      WRITE(6,131)Iopt
CJEM 131  FORMAT(/'Option (1/2/3)                  <',I1,'>:',$)
CJEMCJEM      CALL READI(5,Iopt)
CJEM      IF(Istatus.NE.0)RETURN

      IF(Iopt.EQ.1)THEN
        DO J=0,Iydim-1
          Emax=a0+upper(J)*a1
          lower(J)=(500.-a0)/a1                     !500 keV IF E>2MeV
          IF(Emax.LT.2000)lower(J)=(Emax*.25-a0)/a1
          IF(lower(J).LT.low)lower(J)=low
        ENDDO
      ENDIF

      IF(Iopt.EQ.2)THEN
CJEM        WRITE(6,132)lowCh
CJEM 132    FORMAT(/'Lower channel for chi-test  <',I5,'>:',$)
CJEM        CALL READI(5,lowCh)
CJEM  TODO make lowCh input argument
        DO J=0,Iydim-1
          lower(J)=lowCh
        ENDDO
      ENDIF
      IF(Istatus.NE.0)RETURN

      IF(Iopt.NE.1.AND.Iopt.NE.2)THEN
        WRITE(6,*)'No unfolding performed'
        Istatus=1
        RETURN
      ENDIF

CJEM  TODO: Input argument ANS (boolean?) with default
      ANS='y'

CJEM      WRITE(6,133)ANS    
CJEM 133  FORMAT(/,'Include total detector efficiency (y/n) <',A1,'>:',$)
CJEM      CALL READA1(5,ANS)
      IF(Istatus.NE.0)RETURN
      IF(ANS.EQ.'y'.OR.ANS.EQ.'Y')THEN
        CALL ExpThres
      ENDIF

CJEM  TODO: iter as input argument
      iter=33
CJEM      WRITE(6,134)iter   
CJEM 134  FORMAT(/,'Number of iterations ( <200 )  <',I2,'>:',$)
CJEMCJEM      CALL READI(5,iter)
      IF(iter.GT.200)iter=200
      IF(Istatus.NE.0)RETURN

CJEM  TODO wfluc as input with default
      wfluc=0.2
CJEM      WRITE(6,*)' '
CJEM      WRITE(6,*)'The iteration is terminated when the folding of'
CJEM      WRITE(6,*)'the unfolded spectrum equals the raw spectrum. It'
CJEM      WRITE(6,*)'is however recommended to stop before the Chi2 is'
CJEM      WRITE(6,*)'at minimum. Thus, you can put a certain weight on'
CJEM      WRITE(6,*)'the importance not to have too strong fluctuations'
CJEM      WRITE(6,*)'in the final spectrum. We recommend a weight-factor'
CJEM      WRITE(6,*)'of 0.2 (valid range is 0.0 - 0.5)'
CJEM      WRITE(6,135)wfluc   
CJEM 135  FORMAT(/,'Weight on fluctuations <',F3.1,'>:',$)
CJEM      CALL READF(5,wfluc)
      IF(wfluc.GT.0.5)wfluc=0.2
      IF(Istatus.NE.0)RETURN



C Loop for spectrum J to be unfolded***************************
      write(6,*)'Iydim = ',Iydim
      DO J=0,Iydim-1
        high=   upper(J)
CJEM        write(6,*)'high=',high
        ChiHigh=upper(J)
        ChiLow= lower(J)

C Getting the raw spectrum Raw(I)
CJEM        write(6,*)'ITYPE =',ITYPE
        DO I=0,2047
CJEM          IF(ITYPE.GT.1)Raw(I)=rMAT(ISP,I,J)
CJEM          IF(ITYPE.GT.1)THEN
            Raw(I)=rawmat(I,J)
CJEM            IF(rawmat(I,J).gt.0)write(6,*)'Rawmat(',I,',',J,') =',
CJEM     &       rawmat(I,J)
CJEM            write(6,*)'Raw(',I,',',J,')=',Raw(I)
CJEM          ENDIF
CJEM          IF(ITYPE.EQ.1)Raw(I)=rSPEC(ISP,I)
        ENDDO
        sum=   0
        sumCH= 0
        DO I=low,high
          sum  =sum + Raw(I)
CJEM          write(6,*)'Hello from inside sum loop. Sum =',sum
        ENDDO
        write(6,*)'sum=',sum
        DO I=ChiLow,ChiHigh
          sumCH=sumCH+Raw(I)
        ENDDO

C Initialize parameters
        DO I=1,Iter
          DO JJ=0,LEN-1
            SoluF(I,JJ)=0
          ENDDO
          DO JJ=1,7
            SoluP(I,JJ)=0
          ENDDO
        ENDDO

        DO I=0,LEN-1
          U(I)    =0.0
          F(I)    =0.0
          Fsave(I)=0.0
          Usave(I)=0.0
        ENDDO

        ModeChange=0
        CHIsave=1000000000
        CHIold =1000000000
                              
        FlucRaw=Fluc(Raw,ChiLow,ChiHigh,a1)
       
C The experimental pulse height spectrum is the first
C approximation to the unfolded gamma spectrum U(I)
        DO I=low,high
          U(I)=Raw(I)
        ENDDO

C Iteration loop for spectrum J. Folding the unfolded spectrum
C We have to start the integration somewhat below ch I, due to the
C detector resolution. We take conservative 20% resolution
        mode=0                      !starts with no method
        iTry=0
        DO L=1,Iter
          DO I=low,high
            F(I)=0.0
            Klow=I*0.8
            DO K=Klow,high
              F(I)=F(I)+R(I,K)*U(K)
            ENDDO
          ENDDO

C Calculate chisquare between folded and raw spectrum
          CHISQ=0.
          Ichi=0
          DO I=ChiLow,ChiHigh
            sig2=Raw(I)
            IF(sig2.LT.4.)sig2=4.
            CHISQ=CHISQ+((F(I)-Raw(I))**2)/sig2
            Ichi=Ichi+1
          ENDDO
          IF(Ichi.LT.1)Ichi=1
          IF(Ichi.GT.1)Ichi=Ichi-1
          CHISQ=CHISQ/Ichi

C Compute sum of counts
          Fsum  =0.
          FsumCH=0.
          DO I=low,high
            Fsum=Fsum+F(I)
          ENDDO
          DO I=ChiLow,ChiHigh
            FsumCH=FsumCH+F(I)
          ENDDO

C           write(6,*)'Hello world'
 
          RelFluc=Fluc(U,ChiLow,ChiHigh,a1)/FlucRaw
          SoluP(L,1)=mode
          SoluP(L,2)=sum
          SoluP(L,3)=Fsum
          SoluP(L,4)=sumCH
          SoluP(L,5)=FsumCH
          SoluP(L,6)=CHISQ
          SoluP(L,7)=RelFluc
          DO JJ=low,high
            SoluF(L,JJ)=U(JJ)
          ENDDO
          Lmax=L

C Saves the best solution after at least 3 iterations
          IF(L.GT.3.AND.CHISQ.LT.CHIsave)THEN
            Lsave=L
            CHIsave=CHISQ
            DO I=low,high
              Usave(I)=U(I)
              Fsave(I)=F(I)
            ENDDO
          ENDIF

          IF(L.GT.3.AND.ABS( CHISQ-CHIold).LT.0.0006)       IterStop=1
          IF(L.GT.3.AND.ABS((CHISQ-CHIold)/CHISQ).LT.0.002) IterStop=1
          IF(iTry.GT.10.AND.CHISQ.GT.CHIold)                IterStop=1
          IF(iTry.GT.Iter/2.0 .AND.ModeChange.EQ.1)         IterStop=1
          iTry=iTry+1
          IF(IterStop.EQ.1)THEN
            IF(ModeChange.LT.10)THEN                   !Changing mode
              IterStop=0    !flags if going to end for one or another mode
              iTry    =0
              mode=mode*(-1)
              ModeChange=ModeChange+1
              DO I=low,high                            !Using the best solution
                F(I)=Fsave(I)                          !as initial function
                U(I)=Usave(I)
              ENDDO
            ELSE
              GO TO 100                                !End iteration          
            ENDIF
          ENDIF

C mode=-1 selects difference, mode=+1 ratio iteration
          IF(L.EQ.1)mode=-1                            !used for loop number 2
          IF(mode.EQ.-1) THEN
            DO I=low,high
              U(I)=U(I)+(Raw(I)-F(I))                  !difference mode
            ENDDO
          ELSE
            DO I=low,high
              IF(ABS(F(I)).GT.4)U(I)=U(I)*(Raw(I)/F(I))!ratio mode
            ENDDO
          ENDIF
          CHIold=CHISQ
        ENDDO

C Iteration loop for spectrum J ended

  100   CONTINUE

C Finding the best solution: It will be loop number Isc
        CALL SCORING(SoluP,Lsave,Lmax,SCORE,Isc,wfluc)
CJEM        write(6,*)'Scoring says best solution is',Isc

C Making compressed output in case of singles spectrum
        IF(ITYPE.EQ.1)THEN
          moil=' n'
          moim=' n'
          moir=' n'
          WRITE(6,26)
  26      FORMAT(30X,'  S U M M A R Y')
          WRITE(6,27)
          LLH=(L/3.)
          DO IL=1,LLH
            moil=' n'
            moim=' n'
            moir=' n'
            IM=IL+LLH
            IR=MIN0(IM+LLH,200)
            IF(SoluP(IL,1).EQ.-1)moil=' d'
            IF(SoluP(IL,1).EQ.+1)moil=' r'
            IF(SoluP(IM,1).EQ.-1)moim=' d'
            IF(SoluP(IM,1).EQ.+1)moim=' r'
            IF(SoluP(IR,1).EQ.-1)moir=' d'
            IF(SoluP(IR,1).EQ.+1)moir=' r'
            WRITE(6,28)IL,MOIL,SoluP(IL,6),SoluP(IL,7),SCORE(IL),
     +                 IM,MOIM,SoluP(IM,6),SoluP(IM,7),SCORE(IM),
     +                 IR,MOIR,SoluP(IR,6),SoluP(IR,7),SCORE(IR)
  27      FORMAT(' LP MD  CHISQ  Fluct SCR     LP MD  CHISQ  
     & Fluct SCR     LP MD  CHISQ  Fluct SCR')
  28        FORMAT(I3,A2,F8.2,F7.2,A4,2(I7,A2,F8.2,F7.2,A4))
          ENDDO
        ENDIF

        modus=' n'
        IF(SoluP(Isc,1).EQ.-1)modus=' d'
        IF(SoluP(Isc,1).EQ.+1)modus=' r'
        K2=SoluP(Isc,5)+0.5
        K3=SoluP(Isc,4)+0.5

        WRITE(6,30)J,modus,K2,K3,SoluP(Isc,6),SoluP(Isc,7)
  30    FORMAT(' Row:',I4,'  Mode:',A2,'  Area:',I9,'(',I9,')  
     & Chi:',F7.2,'  Fluct:',F6.2)

        IF(iVersion.EQ.1)THEN   !Dropping the Compton subtraction method
          DO i=0,LEN-1
            u(i)=0.0
          ENDDO
          DO i=low,high
            SoluF(Isc,i)=SoluF(Isc,i)   !April 2013 *pf(i)
            u(i)=SoluF(Isc,i)
          ENDDO
          GO TO 999  
        ENDIF



C***************************************************************(new begin)
C       New method: Compton Subtraction Method (Desember 1995/mg)
C       Reference: M. Guttormsen et al. NIM (1996), in press
C       The resolution in the unfolded spectrum u0 is about 0.87FWHM. 
C       Thus, it remains to smooth it with (1/10)*FWHM.
C       Then we deduce S - U (rawspectrum - unfolded) to get the 
C       Compton-background etc., that is a smooth function of
C       energy (except for single+double+511 keV peaks). 
C       Therefore, we smooth
C       this differanse with FWHM, and calculate Ufinal=S-smooth(S-U). 
C       Remember, the FWHM parameter contains only 50% of
C       the value given by user: FWHM = FWHMresp = (1/10)*FWHMexp
        DO i=0,LEN-1
          u0(i)   =0.0
          su0(i)  =0.0
          u(i)    =0.0
          v(i)    =0.0
          us(i)   =0.0
          ud(i)   =0.0           ! spectrum names as in NIM paper
          ua(i)   =0.0
          sua(i)  =0.0
          w(i)    =0.0
          sw(i)   =0.0
          c(i)    =0.0
          sc(i)   =0.0
        ENDDO

CJEM        CALL ERASE                 !erase the spectrum window
CJEM        icmap1=COLORMAP(1)         !saving this, and put back value at end
CJEM        ired  =COLORMAP(10)        !colors for spectra
CJEM        iblue =COLORMAP(1)
CJEM        igreen=COLORMAP(4)
CJEM        ITYPEwait=ITYPE            !have to simulate that we deal with
CJEM                                   !singles for the reason of displaying spectra
CJEM        ITYPE=1
CJEM        LOCH=0                     !display markers
CJEM        HICH=high

C Taking the best solution from earlier
        DO i=low,high
          u0(i)=SoluF(Isc,i)
        ENDDO

C Show spectrum u0
        DO i=low,high
          rSPEC(IDEST,i)=u0(i)
        ENDDO
        i1=1
        i2=2
        i3=0
CJEM        CALL SetMarker(0,2,0)
CJEM        COLORMAP(1)=iblue
CJEM        CALL DSPSP(i1,i2,i3,*333)
CJEM  333   CONTINUE

C Making us, ud and ua (single, double, annih) from unfolded spectrum
C Single and double escape peaks have the same shape as u0 (same FWHM).
        id511 =min0(high,INT(( 511./a1)+0.5))
        id1022=min0(high,INT((1022./a1)+0.5))
        DO i=high,id1022,-1
          us(i- id511)=u0(i)*ps(i)
        ENDDO

        DO i=high,id1022,-1
          ud(i-id1022)=u0(i)*pd(i)
        ENDDO

        i511= min0(high,INT(((511.-a0)/a1)+0.5))
        ua511=0.
        DO i=high,i511,-1
          ua511=ua511+u0(i)*pa(i)
        ENDDO
        Egami=(511.-a0)/a1
        il=Egami                          ! distributing counts on ch il and ih
        ih=il+1
        yl=(ih-Egami)*ua511
        yh=(Egami-il)*ua511
        IF(il.GE.0.AND.il.LE.high)ua(il)=ua(il)+yl
        IF(ih.GE.0.AND.ih.LE.high)ua(ih)=ua(ih)+yh
        ymax=AMAX1(yl,yh) !finding fwhm already present for 511 peak
        ymin=AMIN1(yl,yh)
        IF(ymax.GT.0)THEN
          w0=(ymin+ymax)/ymax
        ELSE
          w0=1.
        ENDIF
        factor=facFWHM*1.03
        m1=max0(INT((0.7*i511+0.5)-1),low)
        m2=min0(INT((1.3*i511+0.5)+1),high)
        CALL GaussSmoothing(ua,sua,m1,m2,factor,w0)! have to be smoothed since ua is
                                             ! a spike in two channels around i511
C Making the us + ud spectrum
        DO i=low,high
          w(i)=us(i)+ud(i)
        ENDDO
C Smoothing with additional 0.5*FWHM
        factor=1.
        w0=1.           
        CALL GaussSmoothing(w,sw,low,high,factor,w0)
C Adding the sua spectrum to get final sw spectrum
        DO i=low,high
          sw(i)=sw(i)+sua(i)
        ENDDO
C Showing sw
        DO i=low,high
          rSPEC(IDEST,i)=sw(i)
        ENDDO
        i1=1
        i2=2
        i3=0
CJEM        CALL SetMarker(0,2,0)
CJEM        COLORMAP(1)=igreen
CJEM        CALL DSPSP(i1,i2,i3,*555)
CJEM555     continue
C Smoothing the u0 spectrum
        factor=1.
        w0=1.
        CALL GaussSmoothing(u0,su0,low,high,factor,w0)
C Multiplying down with pf
        DO i=low,high
          su0(i)=su0(i)*pf(i)        !April 2013 *pf(i)
        ENDDO
C Making the v spectrum
        DO i=low,high
          v(i)=su0(i)+sw(i)
        ENDDO
C Making the Compton c spectrum
        DO i=low,high
          c(i)=Raw(i)-v(i)              !c is Compton-background
        ENDDO

C Showing Compton spectrum c
        i1=1
        i2=2
        i3=0
        DO i=low,high
          rSPEC(IDEST,i)=c(i)
        ENDDO
        COLORMAP(1)=ired
CJEM        CALL SetMarker(0,2,0)
CJEM        CALL DSPSP(i1,i2,i3,*666)
CJEM  666   CONTINUE

C Smoothing Compton with 50% of FWHM
        factor=facFWHM/2.
        w0=1.
        CALL GaussSmoothing(c,sc,low,high,factor,w0)   !sc is smoothed Compton-background
C Showing original raw spectrum 
        i1=2
        i2=2
        i3=0
        DO i=low,high
          rSPEC(IDEST,i)=Raw(i)
        ENDDO
        COLORMAP(1)=ired
CJEM        CALL SetMarker(0,2,0)
CJEM        CALL DSPSP(i1,i2,i3,*777)
CJEM  777   CONTINUE
C Showing smoothed Compton sc+single+double+ann
        i1=2
        i2=2
        i3=0
        DO i=low,high
          rSPEC(IDEST,i)=sc(i)+sw(i)
        ENDDO
        COLORMAP(1)=igreen
CJEM        CALL DSPSP(i1,i2,i3,*888)
CJEM  888   CONTINUE

        DO i=low,high
          u(i)=Raw(i)-sc(i)-sw(i)         !u is the unfolded spectrum
        ENDDO

C Using the photo/total probability pf(i) to make u(i) contain all counts
        DO i=low,high
          IF(pf(i).GT.0.)THEN
            u(i)=u(i)/pf(i)                 !April 2013 /pf(i)
          ELSE
            u(i)=0.
          ENDIF
CJEM          if(u(i).gt.0)write(6,*)'J=',J,' u(',i,')=',u(i)
        ENDDO

C Showing final u/pf
        DO i=low,high
          rSPEC(IDEST,i)=u(i)
        ENDDO
        COLORMAP(1)=iblue
CJEM        CALL SetMarker(0,2,0)
CJEM        CALL DSPSP(i1,i2,i3,*999)

C Putting back what was in rSPEC(IDEST,i)
        DO i=0,MAXCH
          rSPEC(IDEST,i)=rWAIT(i)
        ENDDO
        ITYPE=ITYPEwait
C Putting back color
        COLORMAP(1)=icmap1

 999    CONTINUE


C********************************************************(new ended)

C Correcting for detector efficiency as function of gamma-energy
        negstat=-1
        DO i=0,high
CJEM          rMAT(IDEST,I,J)=0
          unfmat(I,J)=0
CJEM          rSPEC(IDEST,I)=0
          IF(u(i).GT.0)negstat=+1
          IF(i.lt.low)u(i)=0
          IF(negstat.EQ.-1)u(i)=0 !remove neg. counts in the first channels
CJEM  TODO: What is the ANS variable here? Is it needed? Something with old/new version?
CJEM          IF(ANS.EQ.'y'.OR.ANS.EQ.'Y')THEN
            IF(u(i).GT.0)write(6,*)'Hello from inside last loop, i =',i
            effi=EffTot(I)*EffExp(I)
CJEM            IF(ITYPE.GT.1.AND.effi.GT.0.)rMAT(IDEST,I,J)=u(i)/effi
CJEM removed ITYPE>1 test            IF(ITYPE.GT.1.AND.effi.GT.0)unfmat(I,J)=u(i)/effi
            IF(effi.GT.0)unfmat(I,J)=u(i)/effi
CJEM            IF(ITYPE.EQ.1.AND.effi.GT.0.)rSPEC(IDEST,I)=u(i)/effi
CJEM          ELSE
CJEM            IF(ITYPE.GT.1)rMAT(IDEST,I,J)=u(i)
CJEM removed ITYPE>1 test        IF(ITYPE.GT.1)unfmat(I,J)=u(i)
CJEM            unfmat(I,J)=u(i)
CJEM            IF(ITYPE.EQ.1)rSPEC(IDEST,I)=u(i)
CJEM          ENDIF
        ENDDO

      ENDDO                                !J-loop for all spectra ended
      fname(2,IDEST)=waitname              !Putting back old name
      END














      SUBROUTINE SCORING(SoluP,Lsave,Lmax,SCORE,Isc,wfluc)
c Routine to calculate the score of each solution. The score-function,
c which should be as small as possible, is simply a weighted sum of
c chisquare and fluctuations

      DIMENSION SoluP(200,7),SC(200)
      CHARACTER SCORE(200)*4

      Isc=0
      DO I=1,200
        SCORE(I)='    '
        SC(200)=0
      ENDDO
      IF(wfluc.LT.0.OR.wfluc.GT.0.5)wfluc=0.2
      wchi=1.-wfluc

C Calculating the score-function
      THRES=0.3*SoluP(LSAVE,6)
      DO I=1,Lmax
        CHI=SoluP(I,6)
        FLU=SoluP(I,7)
        IF(FLU.LT.1.45)FLU=1.45
        SC(I)=wchi*CHI + wfluc*FLU
        EPS=SoluP(I,6)-SoluP(LSAVE,6)
        IF(EPS.GT.THRES)SC(I)=SC(I)+1000
      ENDDO

C Finding the favorit solution
      SCR=1000000
      DO I=1,Lmax
        IF(SC(I).LT.SCR)THEN
          SCR=SC(I)
          Isc=I
        ENDIF
      ENDDO
C Marking the solutions
      THRES1=0.1*THRES
      THRES2=0.3*THRES
      THRES3=0.7*THRES
      DO I=1,Lmax
        EPS=SC(I)-SC(Isc)
        IF(EPS.LT.THRES3)SCORE(I)=' *  '
        IF(EPS.LT.THRES2)SCORE(I)=' ** '
        IF(EPS.LT.THRES1)SCORE(I)=' ***'
      ENDDO

      SCORE(LSAVE)='<   '
      IF(LSAVE.EQ.Isc)SCORE(LSAVE)='<***'
      RETURN
      END















      SUBROUTINE GaussSmoothing(x,y,l1,l2,factor,w0)
C Folding the gamma-detector with a Gauss distribution with sig=FWHM/2.35
C The response matrix is made with FWHMresp=FWHMexp/facFWHM. Thus, the resolution
C of the detector is FWHMexp=factor*FWHMresp, using factor = facFWHM
C The w0 parameter is the resolution already present due to bin-width of channels
C It takes values between 1 and 2
      INTEGER RDIM
      COMMON/response1/R(0:2047,0:2047),RDIM,a0,a1,FWHM,facFWHM
      COMMON/response3/EffTot(0:2047),Fwhm1(0:2047),EffExp(0:2047)
      DIMENSION x(0:2047),y(0:2047)
      W=0
      DO i=0,2047
        y(i)=0.
      ENDDO
      IF(w0.LT.0.OR.w0.GT.2)THEN
        write(6,*)'w0 = ',w0,' is out of range, changed to w0 = 1.0'
        w0=1.
      ENDIF
      w0=w0*ABS(a1)        !from channels to energy
      DO I=l1,l2
        E=a0+I*a1
        Wtot=Fwhm1(I)*(factor*FWHM/100.)*E
        xx=(Wtot*Wtot)-(w0*w0)
        IF(xx.GT.0         )W=SQRT(xx)
        IF(W .LT.ABS(a1)/5.)W=ABS(a1)/5.  

C Finding integration limits. Going 3*sigma to each side
        Kmin=((E-a0-(6.*W/2.35))/a1)+0.5
        Kmax=((E-a0+(6.*W/2.35))/a1)+0.5
        IF(Kmin.LT.l1)Kmin=l1
        IF(Kmax.LT.l1)Kmax=l1
        IF(Kmax.GT.l2)Kmax=l2
        yK=0                     !used as normalization (should be =1)
        DO K=Kmin,Kmax
          EE=a0+K*a1
          yK=yK+GAUSS(E,EE,W)*a1
        ENDDO
        IF(yK.LE.0)yK=10000000
        DO K=Kmin,Kmax
          EE=a0+K*a1
          y(K)=y(K)+((x(I)*GAUSS(E,EE,W)*a1)/yK)
        ENDDO
      ENDDO
      END














      SUBROUTINE ExpThres
      INTEGER RDIM
      COMMON/State/Istatus,ITYPE,IDEST,cal(2,2,2,3),Idim(2,2,2),
     & fname(2,2),comm(2,2)
      CHARACTER fname*8,comm*60, ANS*1
      COMMON/response1/R(0:2047,0:2047),RDIM,a0,a1,FWHM,facFWHM
      COMMON/response3/EffTot(0:2047),Fwhm1(0:2047),EffExp(0:2047)
      DIMENSION EgamD(10),EffD(10),x(10),y(10)
      DATA EgamD/30.,80.,122.,183.,244.,294.,344.,562.,779.,1000./
      DATA EffD/ 0.0,0.0,0.0,0.06,0.44,0.60,0.87,0.99,1.00,1.000/
c      DATA EffD/ 0.0,0.1,0.51,0.66,0.78,0.85,0.89,0.99,1.00,1.000/  !Replaced April 2013

      ID=10
CJEM      write(6,*)' '
CJEM      write(6,*)'The efficiency at low energy (< 1000 ko be given.'
CJEM      write(6,*)'It depends on various experimental con thresholds'
CJEM      write(6,*)'on ADCs, gamma-absorber (2 mm Cu), timhe program'
CJEM      write(6,*)'always assumes Eff = 0.0 at Egam < 30 f = 1.00 at'
CJEM      write(6,*)'Egam > 1000 keV. However, in between ties you can'
CJEM      write(6,*)'define a new shape of the discrimination.'
CJEM      write(6,*)' '

C Want to save EffD values as default. Therefore, defines x and y to work with.
      DO i=1,ID
        x(i)=EgamD(i)
        y(i)=EffD(i)
      ENDDO

C Finding x as a function of channels with energy calibration a0 and a1
 9999   DO I=0,2047
        E=a0+FLOAT(I)*a1
        I1=1                      !finding interpolation points (I1 and I2)
        DO ii=1,ID
          IF(E.GT.x(ii))I1=ii
        ENDDO
        I2=I1+1
        IF(I1.EQ.ID) THEN
          I1=ID-1
          I2=ID
        ENDIF
        EffExp(i)=y(I1)+(y(I2)-y(I1))*(E-x(I1))/(x(I2)-x(I1))
        IF(EffExp(i).LE.0.)EffExp(i)=0.
        IF(E.LT. 30.)EffExp(i)=0.
        IF(E.GE.1408.)EffExp(i)=1.                           
      ENDDO

C Displaying discriminator function
CJEM      CALL DrawStars

C Should function be changed?
      ANS='n'
CJEM      WRITE(6,1)ANS
CJEM   1  FORMAT(/'Do you want to change the discriminator threshold 
CJEM     & <',A1,'>:',$)
CJEM      CALL READA1(5,ANS)
      IF(Istatus.NE.0)RETURN
      IF(ANS.EQ.'y'.OR.ANS.EQ.'Y')THEN
        DO i=2,ID-1
          WRITE(6,2)x(i),y(i)
   2      FORMAT('Give efficiency at ',F6.1,' keV     <',F4.2,'>:',$)
CJEM          CALL READF(5,y(i))
          IF(Istatus.NE.0)RETURN
        ENDDO
        GO TO 9999
      ELSE

C Writting data to file:resp.dat
        Idev=27
        OPEN (Idev,FILE='resp.dat',ACCESS='APPEND',IOSTAT=IOS)
        IF(IOS.EQ.0)THEN
          WRITE(Idev,*)'Details on efficiency up to 1000 keV:'
          WRITE(Idev,22)
          WRITE(Idev,24)
22        FORMAT('    Egam  EffTot  EffExp  EffTot*EffExp')
24        FORMAT('=======================================')
25        FORMAT(   F8.1,  F8.3,    F8.3,     F8.3      )
          icha=((1000.-a0)/a1)+0.5
          iStep=1
          IF(icha.GT.100)iStep=icha/100
          DO i=0,icha,istep
            E=a0+FLOAT(i)*a1
            IF(E.GE.0)WRITE(Idev,25)E,EffTot(i),EffExp(i),EffTot(i)*
     &       EffExp(i)
          ENDDO
          CLOSE(Idev)
        ENDIF
      ENDIF
      RETURN
      END








      FUNCTION Fluc(F,ChiLow,ChiHigh,a1)
C Calculates fluctuations in a spectrum F between ChiLow and ChiHigh
C by (F(i)-Faverage)/Faverage. The average is taken over fwhm=0.12
C at 662 keV, and so on (Oslo 4/3-1988 /M. Guttormsen)
      DIMENSION F(0:2047)
      INTEGER ChiLow,ChiHigh
      Fluc=0.0
      nChannels=0
      i1=ChiLow
      Egam=(ChiLow+1)*a1
      iFW=0.12*(SQRT(662.*Egam))/a1
      IF(iFW.LT.2)iFW=2
      Egam=(ChiLow+iFW/2)*a1
      iFW=0.12*(SQRT(662.*Egam))/a1
      IF(iFW.LT.2)iFW=2
      i2=i1+iFW

      DO WHILE(i2.LE.ChiHigh)
        Sum=0
        DO I=i1,i2
          Sum=Sum+F(I)
        ENDDO
        Average=Sum/(i2-i1+1)
        IF(Average.LT.2.)Average=2.
        DO I=i1,i2
          nChannels=nChannels+1
          Fluc=Fluc+ABS(F(I)-Average)/Average
        ENDDO
        i1=i2+1
        Egam=(i1+iFW/2.)*a1
        iFW=0.12*(SQRT(662.*Egam))/a1
        IF(iFW.LT.2)iFW=2
        i2=i1+iFW
      ENDDO

      Fluc=Fluc/nChannels
      RETURN
      END








      FUNCTION GAUSS(E,EE,W)
C Calculates a Gaussian distribution with centroid at E and half-width W
      SIG=W/2.35
      A=-((EE-E)**2)/(2.*SIG**2)
      GAUSS= (1./(SIG*SQRT(2.*3.141592)))*EXP(A)
      END