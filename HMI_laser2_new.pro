; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE CO-TUNE SEQUENCES
; OBJECTIVE: TO DERIVE THE NON-TUNABLE TRANSMISSION PROFILE

; LIKE HMI_LASER2.PRO BUT THE TRANSMITTANCE OF LYOT ELEMENTS CHANGED
; SO THAT CONTRASTS ARE NOT CORRELATED TO LASER INTENSITY ANYMORE


;----------------------------------------------------------------------


;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser2_new,draw


fitE5       = 1 ;fit for contrasts (0=no, 1=yes) ?

;SET_PLOT,'X'
;WINDOW,0,RETAIN=2,xsize=900,ysize=900
;!P.MULTI    = [0,2,2]
time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0        ; Fe I 6173 line central wavelength in air

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 256               ; number of rows
ny          = 256               ; number of columns
factor      = 10000.
nseq        = 27               ; number of positions in the co-tune sequence
nparam      = 10                ; number of filter parameters we fit for with the cotune: 4 phases, 4 contrasts, laser intensity
IF (fitE5 EQ 0) THEN nparam = 6
anglim      = 1020.01; 971 IN OBSMODE AND 1020 IN CALMODE
xcenter     = (nx-1.)/2.
ycenter     = (nx-1.)/2.
ratio       = 4096./FLOAT(nx)
distance    = FLTARR(nx,nx)
for iii=0,nx-1 do for jjj=0,nx-1 do distance[iii,jjj]=sqrt( (iii-xcenter)^2.d0+(jjj-ycenter)^2.d0 )*ratio*0.5



; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;------------------------------------------------------

; VALUES OF JUNE 2009
FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.17094240d0;0.1710098d0            ; for the narrow-band Michelson
FSR[1]      = 0.34192317d0;0.3421506d0            ; for the broad-band  Michelson
FSR[2]      = 0.69348347d0;0.6943613d0            ; for E1
FSR[3]      = 1.407d0;                            ; for E2
FSR[4]      = 2.779d0;                            ; for E3
FSR[5]      = 5.682d0                             ; for E4
FSR[6]      = 11.354d0                            ; for E5


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;------------------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
;RESTORE,'frontwindow.bin' ; THIS FRONT WINDOW OF THE S/N 3?
; IN OCTOBER 2007 THE FRON T WINDOW WAS S/N 3:
transmission = DBLARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker0     = transmission/100.d0
lam          = wavelength*10.d0-lamref
q            = READFITS('blocker11.fits')
; MAKE SURE THE BLOCKER CENTER IS CORRECT (DEPENDS ON TEMPERATURE, AND
; VACUUM VS. AIR?
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lamref,lam) ; I center the profile (VALUE OF JUNE 2009)

; MEASURED INTENSITIES
;-------------------------------------------------------

; OBSMODE

;list    = [710420,710600,710720,710900,711140,711260,711380,711500,711620,711740,711860,711980,712100,712220,712408,712528,712708,712828,712948,713068,713188,713308,713428,713608,713728,713848,713968,714088,714208]
;list    = [710600,710720,710900,711140,711260,711380,711500,711620,711740,711860,711980,712100,712220,712408,712528,712708,712828,712948,713068,713188,713308,713428,713608,713728,713848,713968,714088,714208]; 28 detune sequences (USE THIS ONE)
;list    = [710600,710720,710900,711140,711260,711380,711500,711620,711740,711860,711980,712100,712220,712408,712708,712828,712948,713068,713188,713308,713428,713608,713728,713848,713968,714088,714208]
;list    = [710720,710900,711140,711260,711380,711500,711620,711740,711860,711980,712100,712220,712408,712708,712828,712948,713068,713188,713308,713428,713608,713728,713848,713968,714088,714208]



;CALMODE
list    = [710660,710780,710960,711200,711320,711440,711560,711680,711800,711920,712040,712160,712280,712468,712588,712768,712888,713008,713128,713248,713368,713488,713668,713788,713908,714028,714148,714268]; 28 detune sequences (USE THIS ONE)
;list    = [710660,710780,710960,711200,711320,711440,711560,711680,711800,711920,712040,712160,712280,712468,712768,712888,713008,713128,713248,713368,713488,713668,713788,713908,714028,714148,714268]
;list    = [710780,710960,711200,711320,711440,711560,711680,711800,711920,712040,712160,712280,712468,712588,712768,712888,713008,713128,713248,713368,713488,713668,713788,713908,714028,714148,714268]
;list    = [710660,710780,711200,711320,711440,711560,711680,711800,711920,712040,712160,712280,712468,712588,712888,713008,713128,713248,713368,713488,713668,713788,713908,714028,714148,714268]

nl0     = N_ELEMENTS(list)              ; number of wavelength positions for the dye laser
lam0    = DBLARR(nl0)

;OBSMODE
;allwavelength= [[710420,6173.3197d0],[710600,6173.4588d0],[710720,6173.6989d0],[710900,6173.8433d0],[711020,6174.0972d0],[711140,6174.2943d0],[711260,6174.6390d0],[711380,6174.8988d0],[711500,6175.8092d0],[711620,6176.0035d0],[711740,6176.3030d0],[711860,6176.5034d0],[711980,6178.5969d0],[712100,6178.9934d0],[712220,6175.4254d0],[712408,6173.2952d0],[712528,6173.1101d0],[712708,6172.9077d0],[712828,6172.7339d0],[712948,6172.5079d0],[713068,6172.3051d0],[713188,6172.0072d0],[713308,6171.7053d0],[713428,6170.8120d0],[713608,6170.6120d0],[713728,6170.3063d0],[713848,6168.0095d0],[713968,6167.6062d0],[714088,6171.3046d0],[714208,6173.6128d0],[710540,6173.4585d0],[710840,6173.8434d0],[712648,6172.9073d0] ]

;CALMODE
allwavelength= [[710480,6173.3317d0],[710660,6173.4582d0],[710780,6173.6985d0],[710960,6173.8483d0],[711080,6174.0966d0],[711200,6174.2943d0],[711320,6174.6383d0],[711440,6174.8986d0],[711560,6175.8087d0],[711680,6176.0031d0],[711800,6176.3027d0],[711920,6176.5035d0],[712040,6178.5968d0],[712160,6178.9928d0],[712280,6175.4247d0],[712468,6173.2949d0],[712588,6173.1096d0],[712768,6172.9072d0],[712888,6172.7329d0],[713008,6172.5080d0],[713128,6172.3052d0],[713248,6172.0071d0],[713368,6171.7048d0],[713488,6170.8117d0],[713668,6170.6113d0],[713788,6170.3056d0],[713908,6168.0089d0],[714028,6167.6054d0],[714148,6171.3044d0],[714268,6173.6119d0] ] ; USING AN AVERAGE OF ALL THE LASER DETUNES

;allwavelength=[[710480d0,6173.3321d0],[710660d0,6173.4585d0],[710780d0,6173.6989d0],[710960d0,6173.8486d0],[711080d0,6174.0969d0],[711200d0,6174.2947d0],[711320d0,6174.6386d0],[711440d0,6174.8989d0],[711560d0,6175.8091d0],[711680d0,6176.0035d0],[711800d0,6176.3030d0],[711920d0,6176.5039d0],[712040d0,6178.5972d0],[712160d0,6178.9931d0],[712280d0,6175.4251d0],[712468d0,6173.2953d0],[712588d0,6173.1099d0],[712768d0,6172.9076d0],[712888d0,6172.7332d0],[713008d0,6172.5083d0],[713128d0,6172.3055d0],[713248d0,6172.0074d0],[713368d0,6171.7052d0],[713488d0,6170.8120d0],[713668d0,6170.6117d0],[713788d0,6170.3060d0],[713908d0,6168.0092d0],[714028d0,6167.6057d0],[714148d0,6171.3047d0],[714268d0,6173.6123d0]] ; USING THE RESULT OF HMI_SUN1.PRO ON 707090 WITH A FIRST ESTIMATE OF THE NON-TUNABLE ELEMENT PROFILES (CAREFUL CONTAMINATION BY CONVECTIVE BLUESHIFT AND LIMB EFFECT)


allwavelength[1,*]=allwavelength[1,*]-lamref

FOR i=0,nl0-1 DO BEGIN
    a = WHERE(allwavelength[0,*] EQ list[i])
    lam0[i]=allwavelength[1,a[0]]
ENDFOR

FOR i=0,nl0-1 DO PRINT,list[i],lam0[i]+lamref
READ,pause

nx2     = 256
ny2     = 256
Inten   = DBLARR(nx2,ny2,nseq*nl0) ; measured output intensities (measured on a HMI CCD) with the cotune
FSN     = LONARR(nseq*nl0);

front   = 1 ; front=1 means use the front camera images (0= side camera)
obsmode = 0 ;=0 for calmode, =1 for obsmode

; OPEN DETUNE SEQUENCES TREATED WITH
; CPT/wavelength_dependence_test_cotunes.pro (FLATFIELDED AND DARK FRAME REMOVED)
FOR i=0,nl0-1 DO BEGIN
    RESTORE,'CPT/SEQUENCE_listLaser071015_'+STRTRIM(list[i],1)+'_'+STRTRIM(STRING(LONG(nx2)),1)+'.BIN'
    IF(front EQ 1) THEN BEGIN
        Inten[*,*,nseq*i:nseq*(i+1)-1] = imf[*,*,2:28]/factor
        FSN[nseq*i:nseq*(i+1)-1] = list[i]+FINDGEN(nseq)*2+4
    ENDIF ELSE BEGIN
        Inten[*,*,nseq*i:nseq*(i+1)-1] = ims[*,*,2:28]/factor
        FSN[nseq*i:nseq*(i+1)-1] = list[i]+FINDGEN(nseq)*2+5
    ENDELSE
ENDFOR

 RESTORE,'CORRECTION_INT_710420.BIN' ; includes all the sequences
;OBTAINED BY intensity_correction.pro RECOMPUTE EACH TIME !!!!
;RESTORE,'CORRECTION_INTENSITY_710420.BIN' ;OBSMODE
;RESTORE,'CORRECTION_INTENSITY_710600.BIN' ;OBSMODE WITH FIRST DETUNE REMOVED
;RESTORE,'CORRECTION_INTENSITY_710720.BIN' ;OBSMODE WITH FIRST 2 DETUNES REMOVED
;RESTORE,'CORRECTION_INTENSITY_710480.BIN' ;CALMODE
;RESTORE,'CORRECTION_INTENSITY_710660.BIN' ;CALMODE WITH FIRST DETUNE REMOVED

INTENSITYdetunesf=INTENSITYdetunesf-0.03156d-9 ; background intensity
INTENSITYdetuness=INTENSITYdetuness-0.03156d-9 ; background intensity



; WE NORMALIZE THE INTENSITIES
FOR i=0,nl0-1 DO BEGIN
;moncul=DBLARR(nseq)
;moncul2=moncul
    FOR j=0,nseq-1 DO BEGIN
        IF(front EQ 1) THEN BEGIN
            a=WHERE(FSNdetunesf EQ FSN[nseq*i+j])
            Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/INTENSITYdetunesf[a[0]]*MEAN(INTENSITYdetunesf)/EXPOSUREf[a[0]]*MEAN(EXPOSUREf)
            ;moncul[j]=INTENSITYdetunesf[a[0]]
            ;moncul2[j]=EXPOSUREf[a[0]]
        ENDIF ELSE BEGIN
            a=WHERE(FSNdetuness EQ FSN[nseq*i+j])
            Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/INTENSITYdetuness[a[0]]*MEAN(INTENSITYdetuness)/EXPOSUREs[a[0]]*MEAN(EXPOSUREs)
        ENDELSE
    ENDFOR
;Inten[*,*,i*nseq:i*nseq+nseq-1]=Inten[*,*,i*nseq:i*nseq+nseq-1]/MEAN(moncul)*MEAN(INTENSITYdetunesf)/MEAN(moncul2)*MEAN(EXPOSUREf)
;plot,moncul
;read,pause
ENDFOR

IF(nx NE nx2) THEN Inten   = REBIN(Inten,nx,ny,nseq*nl0)
; FOR PLOT ONLY:
inten=rebin(inten,1,1,756)
inten2=fltarr(28)
for i=0,27 do inten2[i]=total(inten[0,0,i*27:(i+1)*27-1])
read,pause


; VARIABLE DEFINITION
;--------------------------------------------------------

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
threshold   = DBLARR(nx,ny)
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
dIntendThr  = DBLARR(nseq*nl0)
errorbar    = DBLARR(nseq*nl0)
lateconv    = INTARR(nx,ny)    ; to discard spurious results
maxsteps    = 150;555
Jac         = DBLARR(nparam,nseq*nl0) ; Jacobian matrix of the least-squares fit
contrasts   = 0.99d0


; DEFINITION OF THE CO-TUNE SEQUENCE
; co-tuning in the range [-2 FSR[E1],+2 FSR[E1]]
; TABLE PROVIDED BY JESPER
;----------------------------------------------------------

tuningJS       = DBLARR(nseq,3)
tuningJS[0 ,*] = [         0.d0,          0.d0,          0.d0]
tuningJS[1 ,*] = [        80.d0,          0.d0,          0.d0]
tuningJS[2 ,*] = [       160.d0,          0.d0,          0.d0]
tuningJS[3 ,*] = [         0.d0,         80.d0,          0.d0]
tuningJS[4 ,*] = [        80.d0,         80.d0,          0.d0]
tuningJS[5 ,*] = [       160.d0,         80.d0,          0.d0]
tuningJS[6 ,*] = [         0.d0,        160.d0,          0.d0]
tuningJS[7 ,*] = [        80.d0,        160.d0,          0.d0]
tuningJS[8 ,*] = [       160.d0,        160.d0,          0.d0]
tuningJS[9 ,*] = [         0.d0,          0.d0,         80.d0]
tuningJS[10,*] = [        80.d0,          0.d0,         80.d0]
tuningJS[11,*] = [       160.d0,          0.d0,         80.d0]
tuningJS[12,*] = [         0.d0,         80.d0,         80.d0] 
tuningJS[13,*] = [        80.d0,         80.d0,         80.d0]
tuningJS[14,*] = [       160.d0,         80.d0,         80.d0]
tuningJS[15,*] = [         0.d0,        160.d0,         80.d0]
tuningJS[16,*] = [        80.d0,        160.d0,         80.d0]
tuningJS[17,*] = [       160.d0,        160.d0,         80.d0]
tuningJS[18,*] = [         0.d0,          0.d0,        160.d0]
tuningJS[19,*] = [        80.d0,          0.d0,        160.d0]
tuningJS[20,*] = [       160.d0,          0.d0,        160.d0]
tuningJS[21,*] = [         0.d0,         80.d0,        160.d0]
tuningJS[22,*] = [        80.d0,         80.d0,        160.d0]
tuningJS[23,*] = [       160.d0,         80.d0,        160.d0]
tuningJS[24,*] = [         0.d0,        160.d0,        160.d0]
tuningJS[25,*] = [        80.d0,        160.d0,        160.d0]
tuningJS[26,*] = [       160.d0,        160.d0,        160.d0]

FOR i=0,nseq-1 DO tuningJS[i,*]  = tuningJS[i,*] * [6.d0,6.d0,-6.d0]*!dpi/180.d0


IF draw EQ 1 THEN GOTO,draw

; PHASE AND CONTRAST OF TUNING ELEMENTS (OBTAINED WITH
; CPT/CPT_laser.pro)
;-----------------------------------------------------

IF( OBSMODE EQ 1) THEN BEGIN
; FROM MEDIAN VALUE OF ALL LASER DETUNES IN OBSMODE
meanNB=-113.287d0
meanWB= 5.26063d0
meanE1=-138.755d0
 RESTORE,'CPT/CPT_laser_front_710600.BIN' ; OBSMODE
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-112.86541d0-meanNB,6.0024801d0-meanWB,-138.18929d0-meanE1]*!dpi/180.d0
 p1=PHIG0
 B1=BG0
 RESTORE,'CPT/CPT_laser_side_710600.BIN' ; OBSMODE'
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-113.46681d0-meanNB,5.9007082d0-meanWB,-138.12755d0-meanE1]*!dpi/180.d0
 p2=PHIG0
 B2=BG0
 RESTORE,'CPT/CPT_laser_front_712408.BIN' ; OBSMODE
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-112.84284d0-meanNB,5.8180826d0-meanWB,-138.89860d0-meanE1]*!dpi/180.d0
 p3=PHIG0
 B3=BG0
 RESTORE,'CPT/CPT_laser_side_712408.BIN' ; OBSMODE
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-113.90070d0-meanNB,6.0767232d0-meanWB,-138.86148d0-meanE1]*!dpi/180.d0
 p4=PHIG0
 B4=BG0
ENDIF ELSE BEGIN
; FROM MEDIAN VALUE OF ALL LASER DETUNES IN CALMODE
meanNB=-112.853d0
meanWB= 6.20089d0
meanE1=-135.252d0
; FROM HMI_SUN1.PRO WITH FSN=707090 AND FIRST ESTIMATE OF NON-TUNABLE
; (CAREFUL THERE IS CONTAMINATION BY THE CONVECTIVE BLUESHIFT AND LIMB EFFECT)
; ELEMENT PROFILES (FRONT CAMERA ONLY)
;meanNB=-114.852d0
;meanWB= 10.186d0
;meanE1=-139.284d0

 RESTORE,'CPT/CPT_laser_front_710660.BIN' ; CALMODE (1025 arcsec, wavelength=6173.4585, 256x256)
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-112.08587d0-meanNB,7.3293818d0-meanWB,-136.15046d0-meanE1]*!dpi/180.d0
 p1=PHIG0
 B1=BG0
 RESTORE,'CPT/CPT_laser_side_710660.BIN' ; CALMODE (1025 arcsec, wavelength=6173.4585, 256x256)
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-112.85749d0-meanNB,7.4471143d0-meanWB,-135.61035d0-meanE1]*!dpi/180.d0
 p2=PHIG0
 B2=BG0
 RESTORE,'CPT/CPT_laser_front_712468.BIN' ; CALMODE (1025 arcsec, wavelength=6173.2953, 256x256)
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-113.10764d0-meanNB,7.6432744d0-meanWB,-134.67950d0-meanE1]*!dpi/180.d0
 p3=PHIG0
 B3=BG0
 RESTORE,'CPT/CPT_laser_side_712468.BIN' ; CALMODE (1025 arcsec, wavelength=6173.2953, 256x256)
 FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]-[-113.69024d0-meanNB,7.4421384d0-meanWB,-134.63083d0-meanE1]*!dpi/180.d0
 p4=PHIG0
 B4=BG0
ENDELSE
 phig0=(p1+p2+p3+p4)/4.d0
 bg0  =(b1+b2+b3+b4)/4.d0

;RESTORE,'CPT/CPT_laser_front_712408_NEWFSR.BIN' ; OBSMODE NEW FSRs
;RESTORE,'CPT/CPT_laser_front_712468_NEWFSR.BIN' ; CALMODE NEW FSRs
;RESTORE,'CPT/CPT_laser_front_AVERAGE_NEWFSR.BIN' ; created by average_phase_contrast.pro



a=WHERE(FINITE(Bg0) EQ 0)
IF (a[0] NE -1) THEN Bg0[a]  =1.d0
a=WHERE(FINITE(Phig0) EQ 0)
IF (a[0] NE -1) THEN Phig0[a]=0.d0

IF(nx NE nx2) THEN BEGIN
    Bg0 = REBIN(Bg0,nx,nx,3)
    Phig0=REBIN(Phig0,nx,nx,3)
ENDIF

Bg[*,*,0:2]   = Bg0[*,*,*]         ; contrasts and phases of the tuning elements
Phig[*,*,0:2] = Phig0[*,*,*]

;----------------------------------------------------

!p.multi=[0,2,2]

FOR jjj=0,ny-1 DO BEGIN

    TVIM,phig[*,*,3]
    TVIM,phig[*,*,4]
    TVIM,phig[*,*,5]
    TVIM,phig[*,*,6]

    FOR iii=0,nx-1 DO BEGIN

   ; BEGINNING OF THE LEAST-SQUARES FIT
   ;-----------------------------------


        IF(distance[iii,jjj] LE anglim) THEN BEGIN

            Phig[iii,jjj,3:6]= [-6.,0.,-5.,1.0]/180.d0*!dpi
            Bg[iii,jjj,3:6]  = [.96,.96,.96,.96]
            IF (fitE5 EQ 0) THEN BEGIN
                Phig[iii,jjj,3:6]= [0.0,0.0,0.0,0.0]/180.d0*!dpi 
                Bg[iii,jjj,3:6]  = [0.99d0,0.99d0,0.99d0,0.99d0] ; that means contrasts=0.98
            ENDIF

            thresh       = 7.;1.4;sum/nl0
            I0g[iii,jjj] = thresh

            converg = 0
            jj      = 0

            WHILE converg EQ 0 DO BEGIN

                FOR j=0,nseq-1 DO BEGIN
                    FOR ii=0,nl0-1 DO BEGIN

                        profileg  = INTERPOL(blocker0,lam,lam0[ii])
                        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                        profilegf = profileg
                        FOR i=3,6 DO profileg = profileg * (1.d0-Bg[iii,jjj,i]*SIN(!dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]/2.d0)^2.d0)
                        errorbar[ii*nseq+j]  = 1.d0;rmsINTENSITY[ii]*Inten[iii,jjj,ii*nseq+j]/MEAN(INTENSITYdetunes[ii*nseq:(ii+1)*nseq-1])+.001d0
                        errorbar[ii*nseq+j]  = SQRT(1.d0+256.d0*Inten[iii,jjj,ii*nseq+j])
                        Residual[ii*nseq+j]  = (Inten[iii,jjj,ii*nseq+j] - (profileg*I0g[iii,jjj]+threshold[iii,jjj]))/errorbar[ii*nseq+j]
            
                        IF(fitE5 EQ 1) THEN BEGIN
                            dIntendB3[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB4[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]$
                       *SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB5[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB6[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]
                        ENDIF

                        dIntendPhi3[ii*nseq+j]= profilegf*(-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)*COS(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi4[ii*nseq+j]= profilegf*(-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)*COS(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0))$
                       *(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi5[ii*nseq+j]= profilegf*(-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)*COS(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]$
                       *SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi6[ii*nseq+j]= profilegf*(-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)*COS(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendInt[ii*nseq+j] = profileg/errorbar[ii*nseq+j]

                        dIntendThr[ii*nseq+j] = 1.d0/errorbar[ii*nseq+j] ;/errorbar NOT NEEDED HERE !!!!

                    ENDFOR
                ENDFOR

              ; Computation of the Jacobian matrix
                FOR i= 0,nseq*nl0-1 DO BEGIN
                    IF(fitE5 NE 0) THEN BEGIN
                        Jac[0,i]  = dIntendB3[i]
                        Jac[1,i]  = dIntendB4[i]
                        Jac[2,i]  = dIntendB5[i]
                        Jac[3,i]  = dIntendB6[i]
                        Jac[4,i]  = dIntendPhi3[i]
                        Jac[5,i]  = dIntendPhi4[i]
                        Jac[6,i]  = dIntendPhi5[i]
                        Jac[7,i]  = dIntendPhi6[i]
                        Jac[8,i]  = dIntendInt[i]
                        Jac[9,i]  = dIntendThr[i]
                    ENDIF ELSE BEGIN
                        Jac[0,i]  = dIntendPhi3[i]
                        Jac[1,i]  = dIntendPhi4[i]
                        Jac[2,i]  = dIntendPhi5[i]
                        Jac[3,i]  = dIntendPhi6[i]
                        Jac[4,i]  = dIntendInt[i]
                        Jac[5,i]  = dIntendThr[i]
                    ENDELSE
                ENDFOR               

                TJac    = TRANSPOSE(Jac)
                Hessian = TJac##Jac;+lambda*In

               ;LA_SVD,Jac,W,U,V,/DOUBLE ; Singular Value Decomposition using the LAPACK algorithm
               ;LA_SVD,Hessian,W,U,V,/DOUBLE 
                Hessiani=LA_INVERT(Hessian,/DOUBLE)
                
               ;Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
               ;Dx     = (V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U))##TJac##Residual
                Dx     =  Hessiani##TJac##Residual

                IF (fitE5 NE 0) THEN err = TRANSPOSE(Dx)/[REFORM(Bg[iii,jjj,3]),REFORM(Bg[iii,jjj,4]),REFORM(Bg[iii,jjj,5]),REFORM(Bg[iii,jjj,6]),REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj],threshold[iii,jjj]] ELSE err = TRANSPOSE(Dx)/[REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj],threshold[iii,jjj]]

                err = MAX(ABS(err))

                IF(jj EQ maxsteps-1 AND err GT 1.d-3 )      THEN          converg         = 2 ;d-7
                IF(err LE 1.d-3)                            THEN          converg         = 1

                IF( fitE5 NE 0) THEN BEGIN
                    Bg[iii,jjj,3]   = Bg[iii,jjj,3]  +Dx[0]
                    IF(Bg[iii,jjj,3] LT 0.d0 OR Bg[iii,jjj,3] GT 10.4d0) THEN Bg[iii,jjj,3]   = contrasts ;10.4
                    Bg[iii,jjj,4]   = Bg[iii,jjj,4]  +Dx[1]
                    IF(Bg[iii,jjj,4] LT 0.d0 OR Bg[iii,jjj,4] GT 10.4d0) THEN Bg[iii,jjj,4]   = contrasts
                    Bg[iii,jjj,5]   = Bg[iii,jjj,5]  +Dx[2]
                    IF(Bg[iii,jjj,5] LT 0.d0 OR Bg[iii,jjj,5] GT 10.4d0) THEN Bg[iii,jjj,5]   = contrasts
                    Bg[iii,jjj,6]   = Bg[iii,jjj,6]  +Dx[3]
                    IF(Bg[iii,jjj,6] LT 0.d0 OR Bg[iii,jjj,6] GT 10.4d0) THEN Bg[iii,jjj,6]   = contrasts
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
                    threshold[iii,jjj]=threshold[iii,jjj]+Dx[9]
                    IF(threshold[iii,jjj] LT -0.5) THEN threshold[iii,jjj] = 0.0d0
                ENDIF ELSE BEGIN
                    Phig[iii,jjj,3] = Phig[iii,jjj,3]+Dx[0]
                    IF(ABS(Phig[iii,jjj,3]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,3] = 0.d0
                    Phig[iii,jjj,4] = Phig[iii,jjj,4]+Dx[1]
                    IF(ABS(Phig[iii,jjj,4]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,4] = 0.d0
                    Phig[iii,jjj,5] = Phig[iii,jjj,5]+Dx[2]
                    IF(ABS(Phig[iii,jjj,5]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,5] = 0.d0
                    Phig[iii,jjj,6] = Phig[iii,jjj,6]+Dx[3]
                    IF(ABS(Phig[iii,jjj,6]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,6] = 0.d0
                    I0g[iii,jjj]    = I0g[iii,jjj]   +Dx[4]
                    IF(I0g[iii,jjj] LT 0.2*thresh OR I0g[iii,jjj] GT 4.0*thresh) THEN I0g[iii,jjj] = thresh
                    threshold[iii,jjj]=threshold[iii,jjj]+Dx[5]
                    IF(threshold[iii,jjj] LT -0.5) THEN threshold[iii,jjj] = 0.0d0
                ENDELSE


                lateconv[iii,jjj] = converg
                jj = jj+1

            ENDWHILE

            PRINT,' '
            PRINT,'PIXEL',iii,jjj
            PRINT,'PARAMETERS='
            PRINT,'CONTRASTS'
            Bg[iii,jjj,3:6]=Bg[iii,jjj,3:6]/(2.d0-Bg[iii,jjj,3:6])
            FOR i=3,6 DO PRINT,'Bgi =',Bg[iii,jjj,i]
            PRINT,'PHASES'
            Phig[iii,jjj,*] = Phig[iii,jjj,*]*180.d0/!dpi  ; in degrees
            FOR i=3,6 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'LASER INTENSITY'
            I0g[iii,jjj]=I0g[iii,jjj]*2.d0/(1.d0+Bg[iii,jjj,3])*2.d0/(1.d0+Bg[iii,jjj,4])*2.d0/(1.d0+Bg[iii,jjj,5])*2.d0/(1.d0+Bg[iii,jjj,6])
            PRINT,'I0g=',I0g[iii,jjj]
            PRINT,'threshold=',threshold[iii,jjj]
            PRINT,'lateconv=',converg            
            
        ENDIF

    ENDFOR
ENDFOR

;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS2b_May08_710600_nocontrast_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_FRONT.BIN'
;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN'
 SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_FRONT.BIN'

;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE_BLOCKER+2.BIN'
;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_FRONT_BLOCKER+2.BIN'
;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_FRONT_BLOCKER+2.BIN'
;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN' ;FRONT2 AND SIDE2 are based on the phases of the tunable elements obtained with HMI_sun1.pro. while SIDE and FRONT are based on those obtained from averaging the laser detunes



;SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN'


PRINT,'ELAPSED TIME=',SYSTIME(1)-time0

draw:

;-----------------------------------------------------------------------------------------------

;RESTORE,'RESULTS/RESULTS2_May08_710600_nocontrast_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
RESTORE,'RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_FRONT.BIN'

a=WHERE(distance GT anglim OR lateconv NE 1,COMPLEMENT=b)

FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & temp[a]=-10000.d0 & phig[*,*,i]=temp[*,*] & temp=REFORM(Bg[*,*,i]) & temp[a]=-10000.d0 & Bg[*,*,i]=temp[*,*] & phig[0:1,*,i]=-10000.d0 & Bg[0:1,*,i]=-10000.d0 & ENDFOR
I0g[a]=-10000.d0
threshold[a]=-10000.d0

a = WHERE(FINITE(Phig) EQ 0)
IF(a[0] NE -1) THEN Phig[a]= -10000.d0
a = WHERE(FINITE(Bg) EQ 0)
IF(a[0] NE -1) THEN Bg[a]  = -10000.d0
a = WHERE(FINITE(I0g) EQ 0)
IF(a[0] NE -1) THEN I0g[a]  = -10000.d0
a = WHERE(FINITE(threshold) EQ 0)
IF(a[0] NE -1) THEN threshold[a]  = -10000.d0


; ADD 180 DEGREES TO PHASES < -180
FOR i=3,6 DO BEGIN
temp = REFORM(Phig[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        Phig[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        Phig[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(4)
moyb=moy
FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy[i-3]=MEAN(temp[a]) & temp[b]=moy[i-3] & Phig[*,*,i]=temp & temp=REFORM(Bg[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb[i-3]=MEAN(temp[a]) & temp[b]=moyb[i-3] & Bg[*,*,i]=temp  & ENDFOR
temp=REFORM(I0g[*,*])
a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
moyc=MEAN(temp[a])
temp[b] = 0.d0
I0g = temp
temp=REFORM(threshold[*,*])
a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
moyd=MEAN(temp[a])


; WE PLOT THE RESULT
SET_PLOT,'ps'
!P.MULTI=[0,2,2]
device,file='yo.ps',bits=24,xoffset=0.5,yoffset=1,xsize=21,ysize=19,/color
LOADCT,4
tvim,phig[*,*,3],/scale,tit='!17E2 '+STRTRIM(STRING(moy[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,3]),MAX(phig[*,*,3])]
tvim,phig[*,*,4],/scale,tit='!17E3 '+STRTRIM(STRING(moy[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,4]),MAX(phig[*,*,4])]
tvim,phig[*,*,5],/scale,tit='!17E4 '+STRTRIM(STRING(moy[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,5]),MAX(phig[*,*,5])]
tvim,phig[*,*,6],/scale,tit='!17E5 '+STRTRIM(STRING(moy[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,6]),MAX(phig[*,*,6])]
tvim,Bg[*,*,3]  ,/scale,tit='!17E2 '+STRTRIM(STRING(moyb[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,3]),MAX(Bg[*,*,3])]
tvim,Bg[*,*,4]  ,/scale,tit='!17E3 '+STRTRIM(STRING(moyb[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,4]),MAX(Bg[*,*,4])]
tvim,Bg[*,*,5]  ,/scale,tit='!17E4 '+STRTRIM(STRING(moyb[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,5]),MAX(Bg[*,*,5])]
tvim,Bg[*,*,6]  ,/scale,tit='!17E5 '+STRTRIM(STRING(moyb[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,6]),MAX(Bg[*,*,6])]
tvim,I0g,/scale,tit='!17Intensity',barwidth=0.5,range=[MIN(I0g),MAX(I0g)]
DEVICE,/close
PRINT,'PHASE AVERAGES'
FOR i=0,3 DO PRINT,moy[i]
PRINT,'CONTRAST AVERAGES'
FOR i=0,3 DO PRINT,moyb[i]
PRINT,'I0g AVERAGE'
print,moyc
PRINT,'THRESHOLD AVERAGE'
print,moyd

; RECONSTRUCTION OF THE SEQUENCE
Intenrecons = FLTARR(nx,ny,nseq*nl0)
a=WHERE(distance GT anglim,COMPLEMENT=b)
distance[a] = 0.0
distance[b] = 1.0
FOR i=0,nseq*nl0-1 DO Inten[*,*,i] = Inten[*,*,i]*distance

FOR jjj=0,ny-1 DO BEGIN
    PRINT,jjj
    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] EQ 1.0) THEN BEGIN
            
            FOR j=0,nseq-1 DO BEGIN
                FOR ii=0,nl0-1 DO BEGIN
                    
                    profileg  = INTERPOL(blocker0,lam,lam0[ii],/QUADRATIC)
                    FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!dpi/180.d0+tuningJS[j,i]))/2.d0
                    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!dpi/180.d0))/2.d0
                    
                    Intenrecons[iii,jjj,ii*nseq+j]  = profileg*I0g[iii,jjj]+threshold[iii,jjj]
                    
                ENDFOR
            ENDFOR
            
        ENDIF
        
    ENDFOR
ENDFOR

; WE GET RID OF SOME CRAP IN INTENRECONS. DUNNO WHY IT'S HERE IN THE
; FIRST PLACE (CONCERNS ONLY ABOUT 0.02% OF PIXELS)
; IF NOT CORRECTED, THEY SCREW THE CALCULATION OF RESIDUAL
a=WHERE(Intenrecons Eq -10000.0,na)
Intenrecons[a]=0.0

!P.MULTI=[0,1,2]
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=20,ysize=26,xoffset=0,yoffset=0,/color
LOADCT,3
temp=REBIN(Inten,1,1,nseq*nl0)
temp2=REBIN(Intenrecons,1,1,nseq*nl0)
PLOT,temp,xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5,yst=1,thick=4,xrange=[0,375]
OPLOT,temp2,color=150
PLOT,temp,xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5,yst=1,thick=4,xrange=[375,nseq*nl0]
OPLOT,temp2,color=150
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten,1,1,nseq*nl0))-REFORM(REBIN(Intenrecons,1,1,nseq*nl0)))^2.0 )


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

