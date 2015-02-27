; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE CO-TUNE SEQUENCES
; OBJECTIVE: TO DERIVE THE NON-TUNABLE TRANSMISSION PROFILE
; WE FIT FOR THE FSRs OF NB AND WB AT THE SAME TIME
;----------------------------------------------------------------------
;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser2_FSR,draw

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0        ; Fe I 6173 line central wavelength in air

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 256              ; number of rows
ny          = 256              ; number of columns
factor      = 100.
nseq        = 10               ; number of positions in the co-tune sequence
nl0         = 26               ; number of wavelength positions for the dye laser
nparam      = 12               ; number of filter parameters we fit for with the cotune: 4 phases, 4 contrasts, laser intensity
anglim      = 930
xcenter     = nx/2
ycenter     = nx/2
distance    = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;------------------------------------------------------

FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.17076972d0;0.172457d0         ; for the narrow-band Michelson
FSR[1]      = 0.34158512d0;0.344242d0         ; for the broad-band  Michelson
FSR[2]      = 0.693d0            ; for E1
FSR[3]      = 1.405d0            ; for E2
FSR[4]      = 2.779d0            ; for E3
FSR[5]      = 5.682d0            ; for E4
FSR[6]      = 11.354d0           ; for E5


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;------------------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin'
blocker0     = transmission/100.d0;INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
lam          = wavelength*10.d0-lamref
;blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.d0-lamref,lam) ; I center the profile
blocker1     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.61328-1.8-3.75-lamref,lam)
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.61328-lamref,lam) ; I center the profile


; IF INTENSITY METER IS BEHIND A BLOCKER FILTER
;-------------------------------------------------------

; blocker0[*]  = 1.d0 ;!!!!!!!!!!!!!! WARNING !!!!!
; blocker1[*]  = 1.d0

; MEASURED INTENSITIES
;-------------------------------------------------------

list    = STRARR(nl0)
expo    = FLTARR(nl0)
lam0    = FLTARR(nl0)
relI    = FLTARR(nl0)


; CALMODE
;GOTO,OBSMODE
;list[0 ]= 'listLaser070222_78702'
;list[1 -1]= 'listLaser070222_78728'
list[2 -2]= 'listLaser070222_78767'
list[3 -2]= 'listLaser070222_78832'
list[4 -2]= 'listLaser070222_78899'
list[5 -2]= 'listLaser070222_78925'
list[6 -2]= 'listLaser070222_78951'
list[7 -2]= 'listLaser070222_78977'
list[8 -2]= 'listLaser070222_79003'
list[9 -2]= 'listLaser070222_79029'
list[10-2]= 'listLaser070222_79107'
list[11-2]= 'listLaser070222_79133'
list[12-2]= 'listLaser070222_79159'
list[13-2]= 'listLaser070222_79185' ; weird shapes
;list[1-24]= 'listLaser070222_79263'
list[14-2]= 'listLaser070222_79379'
list[15-2]= 'listLaser070222_79418'
list[16-2]= 'listLaser070222_79457'
list[17-2]= 'listLaser070222_79496'
list[18-2]= 'listLaser070222_79535'
list[19-2]= 'listLaser070222_79561'
list[20-2]= 'listLaser070222_79587'
list[21-2]= 'listLaser070222_79613'
list[22-2]= 'listLaser070222_79665'
list[23-2]= 'listLaser070222_79717'
list[24-2]= 'listLaser070222_79743'
list[25-2]= 'listLaser070222_79795'
list[26-2]= 'listLaser070222_79821'
list[27-2]= 'listLaser070222_79847'

; WAVELENGTHS FROM THE WAVEMETER FILE
;lam0[0 ] = 6173.403d0;6173.4027
;lam0[1 -1] = 6173.150d0;6173.1539
lam0[2 -2] = 6172.9005d0;6172.9064
lam0[3 -2] = 6172.647d0;6172.6565 ; wavemeter has a couple of jumps
lam0[4 -2] = 6172.399d0;6172.4097
lam0[5 -2] = 6172.150d0;6172.1544 ; laser drifts a little bit (<1 mA)
lam0[6 -2] = 6171.900d0;6171.9081
lam0[7 -2] = 6171.400d0;6171.4038
lam0[8 -2] = 6170.898d0;6170.9121 ; laser drifts a little bit (<1 mA)
lam0[9 -2] = 6170.399d0;6170.4111
lam0[10-2] = 6169.899d0;6169.9048
lam0[11-2] = 6169.399d0;6169.4091
lam0[12-2] = 6168.899d0;6168.9073
lam0[13-2] = 6168.399d0;6168.4114
;lam0[1-24] = 6173.341d0;6173.3396
lam0[14-2] = 6173.339d0
lam0[15-2] = 6173.650d0;6173.6569 ; laser drifts slowly to 6173.649
lam0[16-2] = 6173.899d0;6173.9069 ; wavemeter has several jumps
lam0[17-2] = 6174.149d0;6174.1503
lam0[18-2] = 6174.399d0;6174.4071
lam0[19-2] = 6174.649d0;6174.6563 ; laser drifts a little bit (<1 mA)
lam0[20-2] = 6174.900d0;6174.9044
lam0[21-2] = 6175.399d0;6175.4104
lam0[22-2] = 6175.906d0;6175.9062
lam0[23-2] = 6176.400d0;6176.3933
lam0[24-2] = 6176.900d0;6176.8926
lam0[25-2] = 6177.4005d0;6177.3982 ; laser drifts a little bit (<1 mA)
lam0[26-2] = 6177.9005d0;6177.8900 ; laser drifts a little bit (<1 mA)
lam0[27-2] = 6178.400d0;6178.3897
lam0     = lam0-lamref

;relI[0 -1] = 764.233d0
;relI[1 -1] = 673.793d0
relI[2 -2] = 684.916d0
relI[3 -2] = 564.503d0
relI[4 -2] = 441.862d0
relI[5 -2] = 628.823d0
relI[6 -2] = 480.376d0
relI[7 -2] = 542.084d0
relI[8 -2] = 640.077d0
relI[9 -2] = 505.523d0
relI[10-2] = 1907.35d0
relI[11-2] = 1618.34d0
relI[12-2] = 1513.63d0
relI[13-2] = 1826.68d0
;relI[1-24] = 539.917d0
relI[14-2] = 484.695d0
relI[15-2] = 524.810d0
relI[16-2] = 407.981d0
relI[17-2] = 929.936d0
relI[18-2] = 1202.92d0
relI[19-2] = 1347.16d0
relI[20-2] = 1583.57d0
relI[21-2] = 1023.45d0
relI[22-2] = 1675.66d0
relI[23-2] = 1136.39d0
relI[24-2] = 1109.03d0
relI[25-2] = 1346.94d0
relI[26-2] = 1355.20d0
relI[27-2] = 1329.41d0


; OBSMODE
OBSMODE:
GOTO,OBSMODE2
list[0] = 'listLaser070223_78715'
list[1] = 'listLaser070223_78754'
;list[2] = 'listLaser070223_78793' ; wavelength might have jumped at the end of the detune
list[2] = 'listLaser070223_78819' ; wavelength jumped a couple of times
list[3] = 'listLaser070223_78845'
list[4] = 'listLaser070223_78912'
list[5] = 'listLaser070223_78938'
list[6] = 'listLaser070223_78964' ; wavelength had a very short jump
list[7] = 'listLaser070223_78990'
list[8] = 'listLaser070223_79016'
;list[9] = 'listLaser070223_79042'
;list[9] = 'listLaser070223_79068'
list[9] = 'listLaser070223_79094'
list[10]= 'listLaser070223_79120'
list[11]= 'listLaser070223_79146'
list[12]= 'listLaser070223_79198'
list[13]= 'listLaser070223_79237' ; wavemeter shows big jumps
;list[13]= 'listLaser070223_79366' ; etalon A and logs disagree for wavelength
list[14]= 'listLaser070223_79405' ; maybe beginning of sequence was at bad wavelength
list[15]= 'listLaser070223_79431' ; a bit dark
;list[15]= 'listLaser070223_79444' ; a few short wavelength jumps
list[16]= 'listLaser070223_79483'
list[17]= 'listLaser070223_79522'
list[18]= 'listLaser070223_79548'
list[19]= 'listLaser070223_79574'
list[20]= 'listLaser070223_79600'
;list[21]= 'listLaser070223_79626'
list[21]= 'listLaser070223_79652'
list[22]= 'listLaser070223_79704'
list[23]= 'listLaser070223_79730'
;list[24]= 'listLaser070223_79756'
list[24]= 'listLaser070223_79782' ; wavelength may have jumped a bit
list[25]= 'listLaser070223_79808'
list[26]= 'listLaser070223_79834'
;list[27]= 'listLaser070223_79860'


lam0[0] =  6173.150d0
lam0[1] =  6172.900d0
lam0[2] =  6172.648d0
lam0[3] =  6172.400d0
lam0[4] =  6172.150d0
lam0[5] =  6171.900d0
lam0[6] =  6171.400d0
lam0[7] =  6170.899d0
lam0[8] =  6170.399d0
lam0[9] =  6169.899d0
lam0[10]=  6169.399d0 
lam0[11]=  6168.900d0
lam0[12]=  6168.399d0 
lam0[13]=  6173.342d0 
lam0[14]=  6173.650d0
lam0[15]=  6173.900d0
lam0[16]=  6174.150d0 
lam0[17]=  6174.400d0
lam0[18]=  6174.650d0
lam0[19]=  6174.900d0
lam0[20]=  6175.400d0
lam0[21]=  6175.907d0
lam0[22]=  6176.400d0
lam0[23]=  6176.900d0
lam0[24]=  6177.400d0
lam0[25]=  6177.901d0
lam0[26]=  6178.400d0
;lam0[27]= 6173.3006
lam0    = lam0-lamref

relI[0 ]= 679.495d0;673.637d0;
relI[1 ]= 669.347d0;669.347d0;
relI[2 ]= 668.128d0;668.128d0
relI[3 ]= 480.968d0;474.200d0;
relI[4 ]= 603.670d0;607.348d0;
relI[5 ]= 481.708d0;481.708d0
relI[6 ]= 542.603d0;543.051d0;
relI[7 ]= 633.874d0;633.334d0;
relI[8 ]= 489.055d0;502.073d0;
relI[9 ]= 1817.60d0;1817.60d0
relI[10]= 1592.77d0;1610.84d0;
relI[11]= 1511.54d0;1528.19d0;
relI[12]= 1886.16d0;1905.35d0;
relI[13]= 551.192d0;552.00d0;
;relI[13]= 486.d0;488.530d0
relI[14]= 524.292d0;524.292d0
relI[15]= 408.035;405.018d0;
;relI[15]= 407.607d0;407.607d0
relI[16]= 956.450d0;945.622d0;
relI[17]= 1192.77d0;1183.50d0;
relI[18]= 1333.62d0;1333.62d0
relI[19]= 1555.47d0;1555.47d0
relI[20]= 1041.86d0;1041.86d0
relI[21]= 1670.22d0;1670.22d0
relI[22]= 1138.39d0;1138.39d0
relI[23]= 1102.86d0;1114.67d0;
relI[24]= 1337.77d0;1337.77d0
relI[25]= 1299.48d0;1328.69d0;
relI[26]= 1396.10d0;1396.10d0
;relI[27]= 495.603;482.80d0;
OBSMODE2:

Inten    = DBLARR(nx,ny,nseq*nl0) ; measured output intensities (measured on a HMI CCD) with the cotune
;SET_PLOT,'x'
;WINDOW,0,RETAIN=2,xsize=800,ysize=1000
FOR i=0,nl0-1 DO BEGIN
    RESTORE,'SEQUENCE_'+STRTRIM(list[i],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
    Inten[*,*,nseq*i:nseq*(i+1)-1] = imx[*,*,2:11]/factor
    expo[i] = FLOAT(exposuretime)
ENDFOR

; WE NORMALIZE THE INTENSITIES
FOR i=0,nl0-1 DO BEGIN
    FOR j=0,nseq-1 DO BEGIN
        Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/relI[i]/expo[i]*MEAN(expo)*MEAN(relI)
    ENDFOR
ENDFOR

FOR i=0,nl0-1 DO Inten[*,*,nseq*i:nseq*(i+1)-1] = Inten[*,*,nseq*i:nseq*(i+1)-1]/INTERPOL(blocker1,lam,lam0[i])

nx          = 32;64
ny          = 32;64
Inten       = REBIN(Inten,nx,ny,nseq*nl0)
distance    = REBIN(distance,nx,nx)

; VARIABLE DEFINITION
;--------------------------------------------------------

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
FSRg        = DBLARR(nx,ny,3)
Residual    = DBLARR(nseq*nl0) ; residual for the least-squares fit
dIntendB3   = DBLARR(nseq*nl0) ; derivatives to compute the Jacobian matrix
dIntendB4   = DBLARR(nseq*nl0)
dIntendB5   = DBLARR(nseq*nl0)
dIntendB6   = DBLARR(nseq*nl0)
dIntendPhi3 = DBLARR(nseq*nl0)
dIntendPhi4 = DBLARR(nseq*nl0)
dIntendPhi5 = DBLARR(nseq*nl0)
dIntendPhi6 = DBLARR(nseq*nl0)
dIntendInt  = DBLARR(nseq*nl0)
dIntendFSR0 = DBLARR(nseq*nl0)
dIntendFSR1 = DBLARR(nseq*nl0)
dIntendFSR2 = DBLARR(nseq*nl0)
lateconv    = INTARR(nx,ny)    ; to discard spurious results
maxsteps    = 150
Jac         = DBLARR(nparam,nseq*nl0) ; Jacobian matrix of the least-squares fit


; DEFINITION OF THE CO-TUNE SEQUENCE
; co-tuning in the range [-2 FSR[E1],+2 FSR[E1]]
; TABLE PROVIDED BY JESPER
;----------------------------------------------------------

tuningJS      = DBLARR(nseq,3)  ; co-tune sequence of Jesper

tuningJS[0,*] = [ 7, 53, 78]
tuningJS[1,*] = [13, 41, 54]
tuningJS[2,*] = [19, 29, 90]
tuningJS[3,*] = [25, 77, 66]
tuningJS[4,*] = [31, 65,102]
tuningJS[5,*] = [37, 53, 78]
tuningJS[6,*] = [43, 41, 54]
tuningJS[7,*] = [49, 29, 90]
tuningJS[8,*] = [55, 77, 66]
tuningJS[9,*] = [61, 65,102]

tuningJS      = REVERSE(tuningJS,2)
FOR i=0,nseq-1 DO tuningJS[i,*]=REFORM(tuningJS[i,*])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]

IF draw EQ 1 THEN GOTO,draw


; GUESS VALUES
;----------------------------------------------------

;RESTORE,'RESULTS/RESULTS_listLaser070220_74548_256.BIN' ;OBSMODE (WARNING CANCEL nx)
;RESTORE,'RESULTS/RESULTS_listLaser070223_79336_256.BIN' ; OBSMODE
;RESTORE,'RESULTS/RESULTS_listLaser070222_79916_256.BIN' ; CALMODE
RESTORE,'RESULTS/RESULTS_listLaser070222_78530_256.BIN' ; CALMODE
;RESTORE,'RESULTS/RESULTS_listLaser070222_79306_256.BIN' ; CALMODE

temp=REFORM(Phig0[*,*,0]) ; FOR SERIES 78530
a=WHERE(temp GT 150.d0*!pi/180.d0)
temp[a]=temp[a]-2.d0*!dpi
Phig0[*,*,0]=temp
a=WHERE(FINITE(Bg0) EQ 0)
IF (a[0] NE -1) THEN Bg0[a]  =1.d0
a=WHERE(FINITE(Phig0) EQ 0)
IF (a[0] NE -1) THEN Phig0[a]=0.d0

nx  = 32
ny  = 32
Bg0 = REBIN(Bg0,nx,nx,3)
Phig0=REBIN(Phig0,nx,nx,3)

Bg[*,*,0:2]   = Bg0[*,*,*]         ; contrasts and phases of the tuning elements
Phig[*,*,0:2] = Phig0[*,*,*]

; FOR 78530 WITH nx=1
;Phig[0,0,0:2]=[-179.55462,1.1009741,-158.49477]*!pi/180.d0
;Bg[0,0,0:2] = [0.97110963,0.98628050,0.96947014]

FOR jjj=0,ny-1 DO BEGIN

    FOR iii=0,nx-1 DO BEGIN

   ; BEGINNING OF THE LEAST-SQUARES FIT
   ;-----------------------------------


        IF(distance[iii,jjj] LE anglim) THEN BEGIN

            Bg[iii,jjj,3:6]  = [0.95d0,0.959d0,0.975d0,0.92d0]
            Phig[iii,jjj,3:6]= [-12.0,-4.57,-5.43,-20.98]*!dpi/180.d0
            FSRg[iii,jjj,*]  = FSR[0:2]

            ; guess of the initial intensity
            sum = 0.d0
            FOR j=0,nseq-1 DO BEGIN
                FOR ii=0,nl0-1 DO BEGIN
                    profileg  = INTERPOL(blocker0,lam,lam0[ii])
                    FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSRg[iii,jjj,i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0
                    sum = sum+profileg
                ENDFOR
            ENDFOR
            thresh       = TOTAL(REBIN(Inten,1,1,nl0*nseq))/sum;TOTAL(Inten[iii,jjj,*])/sum
            I0g[iii,jjj] = thresh
            converg = 0
            jj      = 0

            WHILE converg EQ 0 DO BEGIN

                FOR j=0,nseq-1 DO BEGIN
                    FOR ii=0,nl0-1 DO BEGIN

                        profileg  = INTERPOL(blocker0,lam,lam0[ii])
                        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSRg[iii,jjj,i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                        profilegf = profileg
                        FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0

                        Residual[ii*nseq+j]  = Inten[iii,jjj,ii*nseq+j] - profileg*I0g[iii,jjj]
            
                        dIntendB3[ii*nseq+j] = profilegf*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB4[ii*nseq+j] = profilegf*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4])*(1.d0+Bg[iii,jjj,3]$
                       *COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB5[ii*nseq+j] = profilegf*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))$
                       *(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendB6[ii*nseq+j] = profilegf*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6])*(1.d0+Bg[iii,jjj,4]$
                       *COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]*COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))/16.d0*I0g[iii,jjj]

                        dIntendPhi3[ii*nseq+j]= profilegf*(-Bg[iii,jjj,3]*SIN(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendPhi4[ii*nseq+j]= profilegf*(-Bg[iii,jjj,4]*SIN(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))$
                       *(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendPhi5[ii*nseq+j]= profilegf*(-Bg[iii,jjj,5]*SIN(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,3]$
                       *COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))*(1.d0+Bg[iii,jjj,6]*COS(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))/16.d0*I0g[iii,jjj]

                        dIntendPhi6[ii*nseq+j]= profilegf*(-Bg[iii,jjj,6]*SIN(dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]))$
                       *(1.d0+Bg[iii,jjj,4]*COS(dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]))*(1.d0+Bg[iii,jjj,5]$
                       *COS(dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]))*(1.d0+Bg[iii,jjj,3]*COS(dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]))/16.d0*I0g[iii,jjj]

                        dIntendInt[ii*nseq+j] = profileg

                        dIntendFSR0[ii*nseq+j]= INTERPOL(blocker0,lam,lam0[ii])*(Bg[iii,jjj,0]/FSRg[iii,jjj,0]^2.d0*dpi*lam0[ii]*SIN(dpi/FSRg[iii,jjj,0]*lam0[ii]+Phig[iii,jjj,0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[iii,jjj,1]*COS(dpi/FSRg[iii,jjj,1]*lam0[ii]+Phig[iii,jjj,1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[iii,jjj,2]*COS(dpi/FSRg[iii,jjj,2]*lam0[ii]+Phig[iii,jjj,2]+tuningJS[j,2]))/2.d0*I0g[iii,jjj]
                        FOR i=3,6 DO dIntendFSR0[ii*nseq+j] = dIntendFSR0[ii*nseq+j] * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0

                        dIntendFSR1[ii*nseq+j]=INTERPOL(blocker0,lam,lam0[ii])*(Bg[iii,jjj,1]/FSRg[iii,jjj,1]^2.d0*dpi*lam0[ii]*SIN(dpi/FSRg[iii,jjj,1]*lam0[ii]+Phig[iii,jjj,1]+tuningJS[j,1]))/2.d0*(1.d0+Bg[iii,jjj,0]*COS(dpi/FSRg[iii,jjj,0]*lam0[ii]+Phig[iii,jjj,0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[iii,jjj,2]*COS(dpi/FSRg[iii,jjj,2]*lam0[ii]+Phig[iii,jjj,2]+tuningJS[j,2]))/2.d0*I0g[iii,jjj]
                        FOR i=3,6 DO dIntendFSR1[ii*nseq+j] = dIntendFSR1[ii*nseq+j] * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0

                        dIntendFSR2[ii*nseq+j]=INTERPOL(blocker0,lam,lam0[ii])*(Bg[iii,jjj,2]/FSRg[iii,jjj,2]^2.d0*dpi*lam0[ii]*SIN(dpi/FSRg[iii,jjj,2]*lam0[ii]+Phig[iii,jjj,2]+tuningJS[j,2]))/2.d0*(1.d0+Bg[iii,jjj,0]*COS(dpi/FSRg[iii,jjj,0]*lam0[ii]+Phig[iii,jjj,0]+tuningJS[j,0]))/2.d0*(1.d0+Bg[iii,jjj,1]*COS(dpi/FSRg[iii,jjj,1]*lam0[ii]+Phig[iii,jjj,1]+tuningJS[j,1]))/2.d0*I0g[iii,jjj]
                        FOR i=3,6 DO dIntendFSR2[ii*nseq+j] = dIntendFSR2[ii*nseq+j] * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0

                    ENDFOR
                ENDFOR



              ; Computation of the Jacobian matrix
                FOR i= 0,nseq*nl0-1 DO BEGIN
                    Jac[0,i]  = dIntendB3[i]
                    Jac[1,i]  = dIntendB4[i]
                    Jac[2,i]  = dIntendB5[i]
                    Jac[3,i]  = dIntendB6[i]
                    Jac[4,i]  = dIntendPhi3[i]
                    Jac[5,i]  = dIntendPhi4[i]
                    Jac[6,i]  = dIntendPhi5[i]
                    Jac[7,i]  = dIntendPhi6[i]
                    Jac[8,i]  = dIntendInt[i]
                    Jac[9,i]  = dIntendFSR0[i]
                    Jac[10,i] = dIntendFSR1[i]
                    Jac[11,i] = dIntendFSR2[i]
                ENDFOR

                LA_SVD,Jac,W,U,V,/DOUBLE ; Singular Value Decomposition using the LAPACK algorithm

                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

                err = TRANSPOSE(Dx)/[REFORM(Bg[iii,jjj,3]),REFORM(Bg[iii,jjj,4]),REFORM(Bg[iii,jjj,5]),REFORM(Bg[iii,jjj,6])$
,REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj],REFORM(FSRg[iii,jjj,0]),REFORM(FSRg[iii,jjj,1]),REFORM(FSRg[iii,jjj,2])]

                err = MAX(ABS(err))

                IF(jj EQ maxsteps-1 AND err GT 1.d-7 )      THEN          converg         = 2
                IF(err LE 1.d-7)                            THEN          converg         = 1

                Bg[iii,jjj,3]   = Bg[iii,jjj,3]  +Dx[0]
                IF(Bg[iii,jjj,3] LT 0.d0 OR Bg[iii,jjj,3] GT 10.4d0) THEN Bg[iii,jjj,3]   = 0.95
                Bg[iii,jjj,4]   = Bg[iii,jjj,4]  +Dx[1]
                IF(Bg[iii,jjj,4] LT 0.d0 OR Bg[iii,jjj,4] GT 10.4d0) THEN Bg[iii,jjj,4]   = 0.96
                Bg[iii,jjj,5]   = Bg[iii,jjj,5]  +Dx[2]
                IF(Bg[iii,jjj,5] LT 0.d0 OR Bg[iii,jjj,5] GT 10.4d0) THEN Bg[iii,jjj,5]   = 0.97
                Bg[iii,jjj,6]   = Bg[iii,jjj,6]  +Dx[3]
                IF(Bg[iii,jjj,6] LT 0.d0 OR Bg[iii,jjj,6] GT 10.4d0) THEN Bg[iii,jjj,6]   = 0.92
                Phig[iii,jjj,3] = Phig[iii,jjj,3]+Dx[4]
                IF(ABS(Phig[iii,jjj,3]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,3] = 0.d0
                Phig[iii,jjj,4] = Phig[iii,jjj,4]+Dx[5]
                IF(ABS(Phig[iii,jjj,4]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,4] = 0.d0
                Phig[iii,jjj,5] = Phig[iii,jjj,5]+Dx[6]
                IF(ABS(Phig[iii,jjj,5]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,5] = 0.d0
                Phig[iii,jjj,6] = Phig[iii,jjj,6]+Dx[7]
                IF(ABS(Phig[iii,jjj,6]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,6] = 0.d0
                I0g[iii,jjj]    = I0g[iii,jjj]   +Dx[8]
                IF(I0g[iii,jjj] LT 0.2*thresh OR I0g[iii,jjj] GT 4.0*thresh) THEN I0g[iii,jjj] = thresh
                FSRg[iii,jjj,0] = FSRg[iii,jjj,0]+Dx[9]
                IF(ABS(FSRg[iii,jjj,0]-FSR[0])/FSR[0] GT 0.5)        THEN FSRg[iii,jjj,0] = FSR[0]
                FSRg[iii,jjj,1] = FSRg[iii,jjj,1]+Dx[10]
                IF(ABS(FSRg[iii,jjj,1]-FSR[1])/FSR[1] GT 0.5)        THEN FSRg[iii,jjj,1] = FSR[1]
                FSRg[iii,jjj,2] = FSRg[iii,jjj,2]+Dx[11]
                IF(ABS(FSRg[iii,jjj,2]-FSR[2])/FSR[2] GT 0.5)        THEN FSRg[iii,jjj,2] = FSR[2]

                lateconv[iii,jjj] = converg
                jj = jj+1

            ENDWHILE

            PRINT,'PIXEL',iii,jjj
            PRINT,'PARAMETERS='
            PRINT,'CONTRASTS'
            FOR i=3,6 DO PRINT,'Bgi =',Bg[iii,jjj,i]
            PRINT,'PHASES'
            Phig[iii,jjj,*] = Phig[iii,jjj,*]*180.d0/!dpi  ; in degrees
            FOR i=3,6 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'LASER INTENSITY'
            PRINT,'I0g=',I0g[iii,jjj]
            PRINT,'FSRs=',FSRg[iii,jjj,0],FSRg[iii,jjj,1],FSRg[iii,jjj,2]
            PRINT,'lateconv=',converg

        ENDIF

    ENDFOR
ENDFOR

SAVE,Bg,Phig,I0g,lateconv,FSRg,FILE='RESULTS/RESULTS2_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'ELAPSED TIME=',SYSTIME(1)-time0

draw:

;-----------------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS2_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a=WHERE(distance GT anglim,COMPLEMENT=b) ; depends on the size of the target ;963 obsmode; 925 calmode

IF a[0] NE -1 THEN BEGIN
    FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & temp[a]=-10000.d0 & phig[*,*,i]=temp[*,*] & temp=REFORM(Bg[*,*,i]) & temp[a]=-10000.d0 & Bg[*,*,i]=temp[*,*] & phig[0:1,*,i]=-10000.d0 & Bg[0:1,*,i]=-10000.d0 & ENDFOR
    a = WHERE(FINITE(Phig) EQ 0)
    IF(a[0] NE -1) THEN Phig[a]= -10000.d0
    a = WHERE(FINITE(Bg) EQ 0)
    IF(a[0] NE -1) THEN Bg[a]  = -10000.d0
ENDIF

; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(4)
moyb=moy
moyf=DBLARR(3)
FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & a=WHERE(temp NE -10000.d0 AND temp NE 0.d0,COMPLEMENT=b) & moy[i-3]=MEAN(temp[a]) & IF b[0] NE -1 THEN temp[b]=moy[i-3] & Phig[*,*,i]=temp & temp=REFORM(Bg[*,*,i]) & a=WHERE(temp NE -10000.d0 AND temp ne 0.d0,COMPLEMENT=b) & moyb[i-3]=MEAN(temp[a]) & IF b[0] NE -1 THEN temp[b]=moyb[i-3] &  Bg[*,*,i]=temp & ENDFOR

FOR i=0,2 DO BEGIN &  temp=REFORM(FSRg[*,*,i]) & a=where(lateconv eq 1) & moyf[i]=MEAN(temp[a]) & ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
IF nx NE 1 THEN BEGIN
    !P.MULTI=[0,2,2]
    device,file='yo.ps',bits=24,xoffset=-0.5,yoffset=1,xsize=22.5,ysize=19,/color
    LOADCT,4
    tvim,phig[*,*,3],/scale,tit='!17E2',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)'
    tvim,phig[*,*,4],/scale,tit='!17E3',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)'
    tvim,phig[*,*,5],/scale,tit='!17E4',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)'
    tvim,phig[*,*,6],/scale,tit='!17E5',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)'
    tvim,Bg[*,*,3]  ,/scale,tit='!17E2',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,3]),MAX(Bg[*,*,3])]
    tvim,Bg[*,*,4]  ,/scale,tit='!17E3',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,4]),MAX(Bg[*,*,4])]
    tvim,Bg[*,*,5]  ,/scale,tit='!17E4',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,5]),MAX(Bg[*,*,5])]
    tvim,Bg[*,*,6]  ,/scale,tit='!17E5',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,6]),MAX(Bg[*,*,6])]
    DEVICE,/close
ENDIF
PRINT,'PHASE AVERAGES'
FOR i=0,3 DO PRINT,moy[i]
PRINT,'CONTRAST AVERAGES'
FOR i=0,3 DO PRINT,moyb[i]
FOR i=0,2 DO PRINT,'FSR AVERAGES',moyf[i]


; RECONSTRUCTION OF THE SEQUENCE
Intenrecons = FLTARR(nx,ny,nseq*nl0)
a=WHERE(distance GT anglim,COMPLEMENT=b)
IF a[0] NE -1 THEN distance[a] = 0.0
IF b[0] NE -1 THEN distance[b] = 1.0
FOR i=0,nseq*nl0-1 DO Inten[*,*,i] = Inten[*,*,i]*distance

FOR jjj=0,ny-1 DO BEGIN
    PRINT,jjj
    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] NE 0.0) THEN BEGIN
            
            FOR j=0,nseq-1 DO BEGIN
                FOR ii=0,nl0-1 DO BEGIN
                    
                    profileg  = INTERPOL(blocker0,lam,lam0[ii])
                    FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSRg[iii,jjj,i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.+tuningJS[j,i]))/2.d0
                    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.))/2.d0
                    
                    Intenrecons[iii,jjj,ii*nseq+j]  = profileg*I0g[iii,jjj]
                    
                ENDFOR
            ENDFOR
            
        ENDIF
        
    ENDFOR
ENDFOR


!P.MULTI=0
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=20,ysize=15,xoffset=0,yoffset=0,/color
LOADCT,3

a = WHERE(lateconv EQ 1)
b = WHERE(lateconv EQ 2)
IF b[0] NE -1 THEN BEGIN
    FOR i=0,nseq*nl0-1 DO BEGIN 
        temp = REFORM(Intenrecons[*,*,i])
        temp[b] = MEAN(temp[a])
        Intenrecons[*,*,i] = temp
    ENDFOR
ENDIF

PLOT,REBIN(Inten,1,1,nseq*nl0),xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5,yst=1;,/ylog
OPLOT,REBIN(Intenrecons,1,1,nseq*nl0),color=140;,thick=2,linestyle=3
device,/close
PRINT,'RESIDUAL=',SQRT(TOTAL( (REFORM(REBIN(Inten,1,1,nseq*nl0))-REFORM(REBIN(Intenrecons,1,1,nseq*nl0)))^2.0 ))/TOTAL(REFORM(REBIN(Intenrecons,1,1,nseq*nl0)))

a           = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0
DEVICE,file='yo2.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,bits=24,/color
LOADCT,3
!P.multi=[0,2,5]
FOR i=0,nseq*nl0-1 DO BEGIN & temp=REFORM((inten[*,*,i]-intenrecons[*,*,i])/inten[*,*,i])*distance & a=WHERE(FINITE(temp) EQ 0) & IF(a[0] NE -1) THEN temp[a]=0.0 & TVIM,temp,/scale,barwidth=0.5,range=[-0.1,0.1] & ENDFOR
DEVICE,/CLOSE

SET_PLOT,'x'
!P.MULTI=0
READ,pause

END

