; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW HOW TO CO-TUNE HMI
;
; ASSUMPTIONS: the transmission profile of the non-tunable part of the
; filter comes from measurements performed in-air by R. Shine in sunlight
; with the Lyot filter (excluding the blocker) outside the instrument
;
; THE SOLAR LINE IS MODELED BY A GAUSSIAN PROFILE:
; I(l) = Ic - Ic*d*EXP(-l^2/vl^2)
; WHERE Ic IS THE CONTINUUM, d IS THE LINEDEPTH
; AND vl IS THE LINEWIDTH
; WE FIT FOR 3 PARAMETERS FOR THE LINE PROFILE:
; Ic, depth=Ic*d, and vl
;
; THIS NONLINEAR LEAST-SQUARES PROBLEM IS SOLVED BY A GAUSS-NEWTON ALGORITHM
;
; version 1.4 December 5, 2006
; last modified in June 2009 for the detune FSN 707090-707149
;
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
; NB: be aware of the fact that there might be a systematic error on the
; determination of the phases due to the lack of
; absolute reference in the wavelength scale (we assume that we are
; centered on 6173.3433 A)
;
;----------------------------------------------------------------------

PRO HMI_sun1,draw



lamref      = 6173.3433d0 ; reference central wavelength for the FeI 6173 line IN AIR
time0       = SYSTIME(1); variable used to determine the elapsed time
dpi         = 2.d0*!dpi


; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

nseq        = 27                  ; number of filtergrams in the detune sequence
nx          = 256                 ; number of rows (we rebin the filtergrams)
ny          = nx                  ; number of columns
anglim      = 965.                ; maximum radius of the solar disk in arcseconds
factor      = 10000.              ; we divide the intensity by this value to make the fit easier
nx2         = 256                 ; to rebin the images

; MEASURED INTENSITIES
; these are the input data: the list must contains
; the nseq+2 names of the filtergrams to be analyzed
; (+2 because of the dark frames)
;------------------------------------------------

;list ='listSun060207_230329' ; in Calmode  VIGNETTING
;list ='listSun060208_180841' ; in Calmode  VIGNETTING
;list ='listSun060222_204505' ; in Obsmode
;list ='listSun060224_225806' ; in Calmode  NO APPARENT PROBLEM
;list ='listSun060224_000059' ; in Calmode  CLOUDS ?
;list ='listSun060304_232326' ; in Calmode  CLOUDS ?
;list ='listSun060616_205620' ; in Calmode
;list ='listSun060616_211450' ; in Calmode
;list ='listSun060616_224250' ; in Calmode
;list ='listSun060621_170559' ; in Calmode
;list ='listSun060621_171740' ; in Obsmode
;list ='listSun060623_230811' ; in Obsmode
;list ='listSun060703_205709' ; in Calmode
;list ='listSun060714_221512' ; in Calmode WOBBLE 36 positions
;list ='listSun060714_223907' ; in Calmode WOBBLE 36 positions (with zero position shifted by 20 steps)
;list='listSun060714_231213'  ; in Calmode WOBBLE 36 positions
;list='listSun060714_230106'  ; in Calmode
;list='listSun061031_205938'  ; in Calmode
;list='listSun061031_200735'  ; in Obsmode
;list='listSun070114_42753'   ; in Calmode
;list='listSun070114_42723'   ; in Calmode
;list='listSun070217_66787'   ; in Calmode CRAP !
;list='listSun070217_66926'   ; in Calmode
;list='listSun070217_66956'   ; in Calmode
;list='listSun070217_66675'   ; in Calmode
;list='listSun070219_70777'   ; in Calmode
;list='listSun070219_71079'   ; in Calmode
;list='listSun070219_71161'   ; in Calmode
;list='listSun070219_72027'   ; in Calmode, HMI_TS12_OVN_LYOT=29.8 degrees C
;list='listSun070219_72087'   ; in Calmode
;list='listSun070227_89915'   ; in Calmode ; vignetting at bottom, HMI_TS12_OVN_LYOT=29.8 degrees C
;list='listSun070227_89945'   ; in Calmode
;list='listSun070302_102330'  ; in Calmode ; vignetting at bottom
;list='listSun070302_102360'  ; in Calmode ; vignetting at bottom
;list='listSun070302_102391'  ; in Calmode ; vignetting at bottom
;list='listSun070302_103436'  ; in Calmode ; vignetting at top
;list='listSun070819_450863'  ; in Obsmode
;list='listSun070819_452443'  ; in Obsmode
;list='listSun070819_453647'  ; in Obsmode
;list='listSun070819_452115'  ; in Obsmode
;list='listSun070903_539794'  ; in Calmode  ; VACUUM CALIBRATION SEPTEMBER 2007 
;list='listSun070903_539881'  ; in Obsmode  ; vignetting at bottom right
;list='listSun070903_539949'  ; in Calmode  ; vignetting at bottom right
;list='listSun070903_541538'  ; in Obsmode
;list='listSun070903_541606'  ; in Calmode
;list='listSun070903_543082'  ; in Calmode
;list='listSun070903_543150'  ; in Obsmode
;list='listSun070905_548149'  ; in Calmode
;list='listSun070905_553018'  ; in Calmode
;list='listSun070905_553086'  ; in Obsmode
;list='listSun070905_553154'  ; in Obsmode
;list='listSun070905_553222'  ; in Calmode
;list='listSun070905_554777'  ; in Calmode
;list='listSun070905_554845'  ; in Obsmode
 list='listSun071014_707090'  ; in Calmode, vacuum, oven should be at 30 degrees
;list='listSun071023_719470'  ; in Calmode, air, oven not at stabilized temperature
;list='listSun071023_744572'  ; in Calmode, air
;list='listSun071023_744648'  ; in Calmode, air



 IF list EQ 'listSun070114_42753' THEN lam0 = 370.58260 ; in m/s for 01/14/2007 at 23:46:30 UT
 IF list EQ 'listSun070114_42723' THEN lam0 = 363.74918 ; in m/s for 01/14/2007 at 23:39:30 UT
 IF list EQ 'listSun070217_66787' THEN lam0 = 475.63660 ; in m/s for 02/17/2007 at 21:39:30 UT
 IF list EQ 'listSun070217_66926' THEN lam0 = 562.51958 ; in m/s for 02/17/2007 at 22:41:00 UT
 IF list EQ 'listSun070217_66956' THEN lam0 = 570.25230 ; in m/s for 02/17/2007 at 22:47:00 UT
 IF list EQ 'listSun070217_66675' THEN lam0 = 350.49657 ; in m/s for 02/17/2007 at 20:19:00 UT
 IF list EQ 'listSun070219_70777' THEN lam0 = 115.14830 ; in m/s for 02/19/2007 at 17:26:00 UT
 IF list EQ 'listSun070219_71079' THEN lam0 = 258.14143 ; in m/s for 02/19/2007 at 19:12:00 UT
 IF list EQ 'listSun070219_71161' THEN lam0 = 321.68256 ; in m/s for 02/19/2007 at 19:53:00 UT
 IF list EQ 'listSun070219_72027' THEN lam0 = 676.89487 ; in m/s for 02/20/2007 at 00:14:00 UT
 IF list EQ 'listSun070219_72087' THEN lam0 = 684.77375 ; in m/s for 02/20/2007 at 00:24:00 UT
 IF list EQ 'listSun070227_89915' THEN lam0 = 211.24707 ; in m/s for 02/27/2007 at 18:09:00 UT
 IF list EQ 'listSun070227_89945' THEN lam0 = 217.98412 ; in m/s for 02/27/2007 at 18:14:00 UT
 IF list EQ 'listSun070302_102330'THEN lam0 = 207.25771 ; in m/s for 03/02/2007 at 17:55:00 UT
 IF list EQ 'listSun070302_102360'THEN lam0 = 213.74955 ; in m/s for 03/02/2007 at 18:00:00 UT
 IF list EQ 'listSun070302_102391'THEN lam0 = 251.24194 ; in m/s for 03/02/2007 at 18:27:23 UT
 IF list EQ 'listSun070302_103436'THEN lam0 = 648.81597 ; in m/s for 03/02/2007 at 22:50:00 UT
 IF list EQ 'listSun070819_450863'THEN lam0 =-526.55052 ; in m/s for 08/19/2007 at 18:09:00 UT
 IF list EQ 'listSun070819_452443'THEN lam0 =-43.822087 ; in m/s for 08/19/2007 at 23:56:00 UT
 IF list EQ 'listSun070819_453647'THEN lam0 = 2.7723159 ; in m/s for 08/20/2007 at 01:05:00 UT
 IF list EQ 'listSun070819_452115'THEN lam0 =-421.67158 ; in m/s for 08/19/2007 at 19:21:00 UT
 IF list EQ 'listSun070903_539794'THEN lam0 =-612.55602 ; in m/s for 09/03/2007 at 18:02:04 UT
 IF list EQ 'listSun070903_539881'THEN lam0 =-596.94765 ; in m/s for 09/03/2007 at 18:13:44 UT
 IF list EQ 'listSun070903_539949'THEN lam0 =-585.83781 ; in m/s for 09/03/2007 at 18:21:34 UT
 IF list EQ 'listSun070903_541538'THEN lam0 =-255.40686 ; in m/s for 09/03/2007 at 21:55:54 UT
 IF list EQ 'listSun070903_541606'THEN lam0 =-245.26494 ; in m/s for 09/03/2007 at 22:03:04 UT
 IF list EQ 'listSun070903_543082'THEN lam0 =-150.82231 ; in m/s for 09/03/2007 at 23:18:24 UT
 IF list EQ 'listSun070903_543150'THEN lam0 =-143.77344 ; in m/s for 09/03/2007 at 23:25:34 UT
 IF list EQ 'listSun070905_548149'THEN lam0 =-637.13439 ; in m/s for 09/04/2007 at 17:46:54 UT
 IF list EQ 'listSun070905_553018'THEN lam0 =-729.19216 ; in m/s for 09/05/2007 at 16:28:54 UT
 IF list EQ 'listSun070905_553086'THEN lam0 =-714.99916 ; in m/s for 09/05/2007 at 16:43:4 UT
 IF list EQ 'listSun070905_553154'THEN lam0 =-703.55026 ; in m/s for 09/05/2007 at 16:54:34 UT
 IF list EQ 'listSun070905_553222'THEN lam0 =-695.39378 ; in m/s for 09/05/2007 at 17:02:04 UT
 IF list EQ 'listSun070905_554777'THEN lam0 =-447.59522 ; in m/s for 09/05/2007 at 19:56:24 UT
 IF list EQ 'listSun070905_554845'THEN lam0 =-434.84500 ; in m/s for 09/05/2007 at 20:04:24 UT
 IF list EQ 'listSun071014_707090'THEN lam0 =-181.42364 ; in m/s for 10/14/2007 at 23:51:17 UT
 IF list EQ 'listSun071023_719470'THEN lam0 =-637.91286 ; in m/s for 10/23/2007 at 18:09:00 UT
 IF list EQ 'listSun071023_744572'THEN lam0 =-787.22122 ; in m/s for 10/24/2007 at 15:53:39 UT
 IF list EQ 'listSun071023_744648'THEN lam0 =-773.91412 ; in m/s for 10/24/2007 at 16:09:28 UT


; THE FOLLOWING LINE LOCATES THE CENTER OF THE SOLAR DISK
; ON THE FILTERGRAMS. IT HAS TO BE ADJUSTED:
;center  = [70,58] (for February data in 128*128)
;center   = [126,129]
center   = [(nx-1.)/2.,(ny-1.)/2.]
distance = FLTARR(nx,nx)
FOR i=0,nx-1 do for j=0,nx-1 do distance[i,j]=SQRT((i-center[0])^2.d0+(j-center[1])^2.d0)*0.5d0*4096.d0/double(nx) ; distance in arcseconds from the image center

; READ DETUNE SEQUENCE
; the detune sequence, imx, must have been produced using the codes
; readimages.pro and clean.pro from J. Schou, and must have been saved
; in the idl binary format with the name:
; 'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
; where list contains the list of filtergrams to read and nx is the
; rebinned size of the filtergrams. The binary file must also contain
; the variable lam0, which is obtained from my code sunearthvel.pro
;---------------------------------------------------------------------------------


RESTORE,'CPT/SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN' ; FOR 707090, THE DATA HAVE ALREADY BEEN FLATFIELDED
dlamdv  = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam0    = dlamdv*lam0         ; positive velocities increase the wavelength

;Inten   = imx[*,*,1:27]/factor ; imx is after dark removal, imx[0] and imx[27] are the dark frames
;Inten   = imf[*,*,2:28]/factor ; front camera
Inten   = ims[*,*,2:28]/factor ; side camera


;Inten=imx/factor
;Inten=Inten/factor

IF (nx2 EQ 1) THEN BEGIN
    a = WHERE(distance LE anglim, COMPLEMENT=b)
    temp = FLTARR(nseq)
    FOR i=0,nseq-1 DO BEGIN
        temp2  = REFORM(Inten[*,*,i])
        temp[i]= MEAN(temp2[a])
    ENDFOR
    Inten=FLTARR(1,1,nseq)
    Inten[0,0,*] = temp
ENDIF


wobble:
;nseq=36
;inten=imx

;lam0=0.0018532850 ; for listSun060714_221512
;lam0=0.0024861747 ; for listSun060714_223907
;lam0=0.0032872231 ; for listSun060714_231213
;Inten[*,*,18:26] = Inten[*,*,18:26] * 0.875

; TO CORRECT FOR THE I-RIPPLE ON THE JUNE 2006 DATA ONLY
;RESTORE,'correction.ps'
;correction=correction/MEAN(correction)
;inten[0,0,*]=inten[0,0,*]/correction


; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------

;contrasts   = [0.93,0.96,0.97,0.93]
;phase       = [-6.7,+0.08,-1.45,-5.0]*!dpi/180.d0 ;-6.7
;contrasts   = [0.99,0.968,0.983,0.99]
;phase       = [-4.1,+1.17,-3.6,1.23]*!dpi/180.d0 ;-6.7

;RESTORE,'RESULTS/RESULTS_June09_710660_CAL_256_SIDE.BIN' ;NON-TUNABLE ELEMENTS PHASE+CONTRASTS ; phases of tunable elements and wavelengths obtained from average over laser detune sequences
RESTORE,'RESULTS/RESULTS_June09_710660_CAL_256_SIDE2.BIN' ;NON-TUNABLE ELEMENTS PHASE+CONTRASTS ; phases of tunable elements and wavelengths obtained from HMI_sun1.pro
phase       = PHIG[*,*,3:6]*!dpi/180.d0
contrasts   = BG[*,*,3:6]
IF(nx2 EQ 1) THEN BEGIN
    temp=FLTARR(4)
    FOR i=0,3 DO BEGIN
        temp2  = REFORM(phase[*,*,i])
        temp[i]= MEAN(temp2[a])
    ENDFOR
    phase=FLTARR(1,1,4)
    phase[0,0,*]=temp
    temp=FLTARR(4)
    FOR i=0,3 DO BEGIN
        temp2  = REFORM(contrasts[*,*,i])
        temp[i]= MEAN(temp2[a])
    ENDFOR
    contrasts=FLTARR(1,1,4)
    contrasts[0,0,*]=temp
    nx=nx2
    ny=nx2
    distance=FLTARR(1,1)
ENDIF

RESTORE,'CPT/CPT_laser_side_710660.BIN' ; for contrasts of tunable elements, obtained by CPT_laser1.pro, SIDE CAMERA (256x256)
B3=BG0
RESTORE,'CPT/CPT_laser_side_712468.BIN'
B4=BG0

Bg          = (B3+B4)/2.d0      ; contrasts of the tunable elements


Phig        = DBLARR(nx,ny,3)     ; relative phases
continuum   = DBLARR(nx,ny)       ; relative sun continuum intensity
linewidth   = DBLARR(nx,nx)       ; solar linewidth
linedepth   = DBLARR(nx,nx)       ; solar linedepth
mu0         = 0.d0                ; regularization parameter for the least-squares fit. OPTIONAL: DEPENDS ON THE NOISE LEVEL
nparam      = 6                   ; number of parameters we fit for
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelengths in Angstrom RELATIVE TO lamref
Residual    = DBLARR(nseq)        ; Residual of the least-squares fit
dIntendI0   = DBLARR(nseq)        ; derivatives to compute the Jacobian matrix
dIntendIc   = DBLARR(nseq)
dIntendw    = DBLARR(nseq)
dIntendPhi0 = DBLARR(nseq)
dIntendPhi1 = DBLARR(nseq)
dIntendPhi2 = DBLARR(nseq)
dlinedIc    = 1.d0
tuning      = DBLARR(3,nseq)      ; tuning sequence
maxsteps    = 105                 ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(nx,ny)       ; to discard spurious results

depth0      = 0.52                ; guess for the solar linedepth
width0      = 0.073               ; guess for the solar linewidth
thresh      = 75000./factor       ; guess of the continuum intensity

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172-0.0010576d0;0.1710098d0;0.17076972d0;0.172457; for the narrow-band Michelson
FSR[1]      = 0.344-0.00207683d0;0.3421506d0;0.34158512d0;0.344242; for the broad-band  Michelson
FSR[2]      = 0.693+0.000483467d0;0.6943613d0                ; for E1
FSR[3]      = 1.407d0;1.405d0       ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

FSR         = dpi/FSR

; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
; format:            NB MICHELSON   WB MICHELSON     E1 LYOT
;-----------------------------------------------------------

IF(nseq EQ 27) THEN BEGIN
; FOR DETUNE SEQUENCES
    tuning[*,0]  = [         0.d0,          0.d0,          0.d0]
    tuning[*,1]  = [        80.d0,          0.d0,          0.d0]
    tuning[*,2]  = [       160.d0,          0.d0,          0.d0]
    tuning[*,3]  = [         0.d0,         80.d0,          0.d0]
    tuning[*,4]  = [        80.d0,         80.d0,          0.d0]
    tuning[*,5]  = [       160.d0,         80.d0,          0.d0]
    tuning[*,6]  = [         0.d0,        160.d0,          0.d0]
    tuning[*,7]  = [        80.d0,        160.d0,          0.d0]
    tuning[*,8]  = [       160.d0,        160.d0,          0.d0]
    tuning[*,9]  = [         0.d0,          0.d0,         80.d0]
    tuning[*,10] = [        80.d0,          0.d0,         80.d0]
    tuning[*,11] = [       160.d0,          0.d0,         80.d0]
    tuning[*,12] = [         0.d0,         80.d0,         80.d0] 
    tuning[*,13] = [        80.d0,         80.d0,         80.d0]
    tuning[*,14] = [       160.d0,         80.d0,         80.d0]
    tuning[*,15] = [         0.d0,        160.d0,         80.d0]
    tuning[*,16] = [        80.d0,        160.d0,         80.d0]
    tuning[*,17] = [       160.d0,        160.d0,         80.d0]
    tuning[*,18] = [         0.d0,          0.d0,        160.d0]
    tuning[*,19] = [        80.d0,          0.d0,        160.d0]
    tuning[*,20] = [       160.d0,          0.d0,        160.d0]
    tuning[*,21] = [         0.d0,         80.d0,        160.d0]
    tuning[*,22] = [        80.d0,         80.d0,        160.d0]
    tuning[*,23] = [       160.d0,         80.d0,        160.d0]
    tuning[*,24] = [         0.d0,        160.d0,        160.d0]
    tuning[*,25] = [        80.d0,        160.d0,        160.d0]
    tuning[*,26] = [       160.d0,        160.d0,        160.d0]
    FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0
ENDIF ELSE BEGIN
    
; FOR WOBBLE SEQUENCE 2
;tuning[*,0]  = [   0,   20,   20]
;tuning[*,1]  = [  10,   20,   20]
;tuning[*,2]  = [  20,   20,   20]
;tuning[*,3]  = [  30,   20,   20]
;tuning[*,4]  = [  40,   20,   20]
;tuning[*,5]  = [  50,   20,   20]
;tuning[*,6]  = [  60,   20,   20]
;tuning[*,7]  = [  70,   20,   20]
;tuning[*,8]  = [  80,   20,   20]
;tuning[*,9]  = [  90,   20,   20]
;tuning[*,10] = [ 100,   20,   20]
;tuning[*,11] = [ 110,   20,   20]
;tuning[*,12] = [  20,   0,   20]
;tuning[*,13] = [  20,  10,   20]
;tuning[*,14] = [  20,  20,   20]
;tuning[*,15] = [  20,  30,   20]
;tuning[*,16] = [  20,  40,   20]
;tuning[*,17] = [  20,  50,   20]
;tuning[*,18] = [  20,  60,   20]
;tuning[*,19] = [  20,  70,   20]
;tuning[*,20] = [  20,  80,   20]
;tuning[*,21] = [  20,  90,   20]
;tuning[*,22] = [  20, 100,   20]
;tuning[*,23] = [  20, 110,   20]
;tuning[*,24] = [  20,   20,   0]
;tuning[*,25] = [  20,   20,  10]
;tuning[*,26] = [  20,   20,  20]
;tuning[*,27] = [  20,   20,  30]
;tuning[*,28] = [  20,   20,  40]
;tuning[*,29] = [  20,   20,  50]
;tuning[*,30] = [  20,   20,  60]
;tuning[*,31] = [  20,   20,  70]
;tuning[*,32] = [  20,   20,  80]
;tuning[*,33] = [  20,   20,  90]
;tuning[*,34] = [  20,   20, 100]
;tuning[*,35] = [  20,   20, 110]
; FOR WOBBLE SEQUENCE 1
    tuning[*,0]  = [   0,   0,   0]
    tuning[*,1]  = [  10,   0,   0]
    tuning[*,2]  = [  20,   0,   0]
    tuning[*,3]  = [  30,   0,   0]
    tuning[*,4]  = [  40,   0,   0]
    tuning[*,5]  = [  50,   0,   0]
    tuning[*,6]  = [  60,   0,   0]
    tuning[*,7]  = [  70,   0,   0]
    tuning[*,8]  = [  80,   0,   0]
    tuning[*,9]  = [  90,   0,   0]
    tuning[*,10] = [ 100,   0,   0]
    tuning[*,11] = [ 110,   0,   0]
    tuning[*,12] = [  0,   0,   0]
    tuning[*,13] = [  0,  10,   0]
    tuning[*,14] = [  0,  20,   0]
    tuning[*,15] = [  0,  30,   0]
    tuning[*,16] = [  0,  40,   0]
    tuning[*,17] = [  0,  50,   0]
    tuning[*,18] = [  0,  60,   0]
    tuning[*,19] = [  0,  70,   0]
    tuning[*,20] = [  0,  80,   0]
    tuning[*,21] = [  0,  90,   0]
    tuning[*,22] = [  0, 100,   0]
    tuning[*,23] = [  0, 110,   0]
    tuning[*,24] = [  0,   0,   0]
    tuning[*,25] = [  0,   0,  10]
    tuning[*,26] = [  0,   0,  20]
    tuning[*,27] = [  0,   0,  30]
    tuning[*,28] = [  0,   0,  40]
    tuning[*,29] = [  0,   0,  50]
    tuning[*,30] = [  0,   0,  60]
    tuning[*,31] = [  0,   0,  70]
    tuning[*,32] = [  0,   0,  80]
    tuning[*,33] = [  0,   0,  90]
    tuning[*,34] = [  0,   0, 100]
    tuning[*,35] = [  0,   0, 110]
    
    tuning       = REVERSE(tuning,1)
    FOR i=0,nseq-1 DO tuning[*,i]=REFORM(tuning[*,i])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]
;for wobble sequence 2:
;FOR i=0,nseq-1 DO tuning[*,i]=tuning[*,i]+[-15.,-15.,0.]*!dpi/180.d0 ; when I increase polarizer by 20 steps, I change NB and WB by -15 degrees

ENDELSE


; NON-TUNABLE FILTER PROFILE  PROVIDED BY LOCKHEED-MARTIN
;-----------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
;RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover

; AFTER JULY 2007, FRONT WINDOW S/N 1 WAS REPLACED BY S/N 3
transmission = FLTARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

q            = READFITS('blocker11.fits') ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0    = blocker0*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lamref,lam) ; MODIF OF JUNE 2009

;RESTORE,'Lyot.bin' ; transmission profile of E2-E5 obtained in sunlight on Oct 16, 2006 (see LyotShine_Sun.pro) from lyotsunsetSprofiles2006oct16.fits
;RESTORE,'Lyot2.bin' ; transmission profile of E2-E5 obtained in lamp light on Oct 16, 2006 (see LyotShine_Lamp.pro) from lyotLampsetBprofiles2006oct16.fits
;blocker0     = INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lamref,lam) ;,/LSQUADRATIC) ; uncentered blocker (correction of blocker11.fits because of f8 instead of collimated light)
;blocker0     = blocker0*INTERPOL(data,lresp-0.0075,lam) ; contains the blocker+Lyot+front window profile

;RESTORE,'nontunable.bin'



IF draw EQ 1 THEN GOTO,draw ; if draw=1 we do not want to fit for the phases, we just want to plot an existing result

; BEGINNING OF THE LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=900
!P.MULTI=[0,2,3]

FOR jjj=0,ny-1 DO BEGIN

    TVIM,Phig[*,*,0],/scale 
    TVIM,Phig[*,*,1],/scale
    TVIM,Phig[*,*,2],/scale
    TVIM,linewidth  ,/scale,range=[0.059,0.092]
    TVIM,linedepth  ,/scale
    TVIM,continuum  ,/scale,range=[0.7*thresh,1.3*thresh]

    FOR iii=0,nx-1 DO BEGIN
   
        IF(distance[iii,jjj] LE anglim) THEN BEGIN

           ; GUESS VALUES
           ;------------------------------------------------

            Icg             = thresh                       ; estimate of the solar continuum
            fdepthg         = depth0*thresh                ; depth of the solar line according to Stenflo & Lindegren (1977)
            fwidthg         = width0                       ; width of the solar line according to Stenflo & Lindegren (1977)
            Phig[iii,jjj,*] = [-114.0,7.,-139.]*!dpi/180./FSR[0:2]
            
            converg         = 0
            jj              = 0
           ;mu              = mu0
            
            profilef        = blocker0 ; we use the averaged blocker+front window profile
            FOR i=3,6 DO profilef = profilef * (1.d0+contrasts[iii,jjj,i-3]*COS(FSR[i]*lam+phase[iii,jjj,i-3]))/2.d0
            ;profilef=INTERPOL(nontunable,wavelength,lam)

            WHILE (converg EQ 0) DO BEGIN
                
                templine= EXP(-(lam-lam0)^2.d0/fwidthg^2.d0)
                line    = Icg -fdepthg*templine
                dlinedI0=     -        templine
                dlinedw =     -fdepthg*templine*(lam-lam0)^2.d0*2.d0/fwidthg^3.d0
                
                FOR j=0,nseq-1 DO BEGIN
                    
                    profileg      = 0.125d0 * profilef * (1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j])) * (1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j])) * (1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))
                  
                    Residual[j]   = Inten[iii,jjj,j] - TOTAL(line*profileg)*dlam
                    
                    dIntendI0[j]  = TOTAL(dlinedI0*profileg)*dlam
                    dIntendw [j]  = TOTAL(dlinedw *profileg)*dlam
                    dIntendIc[j]  = TOTAL(dlinedIc*profileg)*dlam
                    dIntendPhi0[j]= TOTAL(line*(-Bg[iii,jjj,0]*FSR[0]*SIN(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef)*dlam
                    dIntendPhi1[j]= TOTAL(line*(-Bg[iii,jjj,1]*FSR[1]*SIN(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef)*dlam
                    dIntendPhi2[j]= TOTAL(line*(-Bg[iii,jjj,2]*FSR[2]*SIN(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))/8.d0*profilef)*dlam
                ENDFOR
                

              ; Jacobian matrix
                Jac  = DBLARR(nparam,nseq)
                FOR i= 0,nseq-1 DO BEGIN
                    
                    Jac[0,i]  = dIntendI0[i]
                    Jac[1,i]  = dIntendIc[i]
                    Jac[2,i]  = dIntendw[i]
                    Jac[3,i]  = dIntendPhi0[i]
                    Jac[4,i]  = dIntendPhi1[i]
                    Jac[5,i]  = dIntendPhi2[i]
                    
                ENDFOR

                LA_SVD,Jac,W,U,V,/DOUBLE
                
              ; regularization:
              ; filter = W*W/(W*W+mu^2.d0) ; mu is the regularization parameter
              ; Dx     = V##DIAG_MATRIX(1.d0/W*filter)##TRANSPOSE(U)##Residual
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
                
                temp   = TRANSPOSE(Dx)/[fdepthg,Icg,fwidthg,REFORM(Phig[iii,jjj,0]),REFORM(Phig[iii,jjj,1]),REFORM(Phig[iii,jjj,2])]

                err[jj]= MAX(ABS(temp))
               ;PRINT,temp,err[jj]

              ; IF(FINITE(err[jj]) EQ 0 AND jj GT 0)          THEN mu = 2.5d0*mu ELSE mu = mu0
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7 )    THEN converg = 2 ; no convergence 
                IF(err[jj] LE 1.d-7)                          THEN converg = 1 ; convergence
                
                IF(converg EQ 2) THEN BEGIN
                    j = WHERE(err EQ MIN(err))
                    fdepthg = history[j[0],0]
                    Icg     = history[j[0],1]
                    fwidthg = history[j[0],2]
                    Phig[iii,jjj,0] = history[j[0],3]
                    Phig[iii,jjj,1] = history[j[0],4]
                    Phig[iii,jjj,2] = history[j[0],5]
                ENDIF
                
                IF(converg NE 2) THEN BEGIN

                    fdepthg = fdepthg+Dx[0]
                    IF(fdepthg LT 0.2d0*thresh*depth0 OR fdepthg GT 3.0d0*depth0*thresh) THEN fdepthg = depth0*thresh
                    Icg     = Icg    +Dx[1]
                    IF(Icg LT 0.2d0*thresh OR Icg GT 3.0d0*thresh) THEN Icg = thresh
                    fwidthg = fwidthg+Dx[2]
                    IF(fwidthg GT 5.d0*width0 OR fwidthg LT 0.3d0*width0) THEN fwidthg = width0
                    Phig[iii,jjj,0] = Phig[iii,jjj,0]+Dx[3]
                    IF(ABS(Phig[iii,jjj,0]) GT 1.2d0*!dpi/FSR[0]) THEN Phig[iii,jjj,0] = 0.
                    Phig[iii,jjj,1] = Phig[iii,jjj,1]+Dx[4]
                    IF(ABS(Phig[iii,jjj,1]) GT 1.2d0*!dpi/FSR[1]) THEN Phig[iii,jjj,1] = 0. 
                    Phig[iii,jjj,2] = Phig[iii,jjj,2]+Dx[5]
                    IF(ABS(Phig[iii,jjj,2]) GT 1.2d0*!dpi/FSR[2]) THEN Phig[iii,jjj,2] = 0.
                    
                ENDIF
                
                history[jj,*]     = [fdepthg,Icg,fwidthg,Phig[iii,jjj,0],Phig[iii,jjj,1],Phig[iii,jjj,2]]
                lateconv[iii,jjj] = converg
                
                jj = jj+1
                
            ENDWHILE
            Phig[iii,jjj,*]    = Phig[iii,jjj,*]*FSR[*]*180.d0/!dpi ; converted into degree units
            linewidth[iii,jjj] = fwidthg
            linedepth[iii,jjj] = fdepthg
            continuum[iii,jjj] = Icg
            PRINT,'PARAMETERS='
            FOR i=0,2 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'Ic=',Icg
            PRINT,'fdepthg=',fdepthg
            PRINT,'fwidthg=',fwidthg
            PRINT,'lateconv=',converg
            PRINT,iii,jjj

        ENDIF
        
    ENDFOR
ENDFOR

;SAVE,Phig,Bg,linewidth,linedepth,continuum,lateconv,lam0,blocker0,lam,FILE='RESULTS/RESULTS_'+STRTRIM(list,1)+'_side_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
SAVE,Phig,Bg,linewidth,linedepth,continuum,lateconv,lam0,blocker0,lam,FILE='RESULTS/RESULTS_'+STRTRIM(list,1)+'_side_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:

; COMPUTE THE AVERAGE PHASES
;--------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS_'+STRTRIM(list,1)+'_side_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a = WHERE(lateconv EQ 1,COMPLEMENT=aa)
IF(a[0] EQ -1) THEN BEGIN
    PRINT,'NO EARLY CONVERGENCES: PROBLEM'
    STOP
ENDIF

PRINT,'AVERAGED PHASES:'
tab = DBLARR(nx,ny)
moy = DBLARR(3)
dis = moy
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

PRINT,'ELAPSED TIME',SYSTIME(1)-time0

a=WHERE(distance LE anglim,COMPLEMENT=b)
distance[a]=1.d0
IF(b[0] NE -1) THEN distance[b]=0.d0

; TO SUPRESS THE POTENTIAL STRIPES AT THE LEFT AND UPPER EDGE
;distance[0:7,*]=0.0
;distance[*,120:nx-1]=0.0
;distance[*,0:2]=0.0
;distance[251:nx-1,*]=0.0

FOR i=0,2 DO Phig[*,*,i]=Phig[*,*,i]*distance
TVIM,distance

PRINT,'CORRECTED AVERAGED PHASES'
a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 1,COMPLEMENT=aa)
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

SET_PLOT,'ps'
IF(nx GT 1) THEN BEGIN
    !P.MULTI=[0,2,3]
    LOADCT,4
    DEVICE,file='yo.ps',/color,bits=24,xoffset=-0.5,yoffset=0,xsize=22,ysize=27
    titl=STRARR(3)
    titl[0]='Narrow-Band Michelson'
    titl[1]='Broad-Band Michelson'
    titl[2]='Lyot Element E1'
    FOR i=0,2 DO TVIM,Phig[*,*,i],/scale,range=[moy[i]-2.5*dis[i],moy[i]+2.5*dis[i]],tit='!17'+titl[i],stit='!17Phase (in degrees)',barwidth=0.5,pcharsize=1.5,xtit='Pixel number',ytit='Pixel number'
    LOADCT,3
    TVIM,continuum*distance,/scale,stit='!17Continuum intensity',barwidth=0.5,pcharsize=1.5,range=[MIN(continuum[a]),MAX(continuum[a])],xtit='Pixel number',ytit='Pixel number'
    TVIM,linewidth*distance,/scale,stit='!17Linewidth (A)',barwidth=0.5,pcharsize=1.5,range=[MIN(linewidth[a]),MAX(linewidth[a])],xtit='Pixel number',ytit='Pixel number'
    TVIM,(linedepth/continuum)*distance,/scale,stit='!17Linedepth',barwidth=0.5,pcharsize=1.5,range=[MIN(linedepth[a]/continuum[a]),MAX(linedepth[a]/continuum[a])],xtit='Pixel number',ytit='Pixel number'
    DEVICE,/close
    !P.MULTI=0
ENDIF


; RECONSTRUCTION OF THE DETUNE SEQUENCE
;--------------------------------------------------------------------------------------


FOR i=0,2 DO Phig[*,*,i] = Phig[*,*,i]/FSR[i]/180.d0*!dpi ; convert from degrees to radians/FSR

FOR i=0,nseq-1 DO Inten[*,*,i] = Inten[*,*,i]*distance
Inten0       = Inten
a            = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0

Inten2  = FLTARR(nx,nx,nseq) 
FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,nx-1 DO BEGIN
        IF(distance[i,j] NE 0.0) THEN BEGIN

            profilef     = blocker0 ; we use the averaged blocker+front window profile
            FOR k=3,6 DO profilef = profilef * (1.d0+contrasts[i,j,k-3]*COS(FSR[k]*lam+phase[i,j,k-3]))/2.d0
            
            templine     = EXP(-(lam-lam0)^2.d0/linewidth[i,j]^2.d0)
            line         = continuum[i,j] -linedepth[i,j]*templine
                        
            FOR k=0,nseq-1 DO BEGIN
                profileg      = 0.125d0 * profilef * (1.d0+Bg[i,j,0]*COS(FSR[0]*(lam+Phig[i,j,0])+tuning[0,k])) * (1.d0+Bg[i,j,1]*COS(FSR[1]*(lam+Phig[i,j,1])+tuning[1,k])) * (1.d0+Bg[i,j,2]*COS(FSR[2]*(lam+Phig[i,j,2])+tuning[2,k]))
                Inten2[i,j,k] = TOTAL(line*profileg)*dlam
            ENDFOR

        ENDIF
    ENDFOR
ENDFOR

read,pause

; WE PLOT THE MAPS OF THE RELATIVE RECONSTRUCTION ERROR
!P.MULTI=[0,2,3]
loadct,4
DEVICE,file='yo2.ps',xoffset=-0.8,yoffset=0,xsize=22,ysize=27,/color,bits=24

FOR i=0,nseq-1 DO BEGIN
    temp  = REFORM((Inten2[*,*,i]-Inten[*,*,i])/Inten[*,*,i])*distance
    a=WHERE(FINITE(temp) EQ 0)
    IF(a[0] NE -1) THEN temp[a]=0.0
    tempa = MEAN(temp[WHERE(temp NE 0.0)])
    TVIM,temp,range=[-0.03+tempa,0.03+tempa],/scale,tit=STRING(i+1),stit='Relative error',pcharsize=1.5
ENDFOR


FOR i=0,nseq-1 DO BEGIN
    temp = REFORM((Inten2[*,*,i]-Inten[*,*,i])/Inten[*,*,i])
    aa   = WHERE(temp LT 0.2 AND temp GT -0.2,na,COMPLEMENT=b)
    hist = histogram(temp[aa],binsize=0.001,min=-0.2,max=0.2)
    plot,FINDGEN(N_ELEMENTS(hist))*0.001-0.2,hist/FLOAT(na),psym=10,tit='!7r!17='+STRING(SIGMA(temp[aa]))+'/!7l!17='+STRING(MEAN(temp[aa])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5
    PRINT,MIN(temp[aa]),MAX(temp[aa]),SIGMA(temp[aa]),MEAN(temp[aa])
ENDFOR
DEVICE,/CLOSE


!P.MULTI=0
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=20,ysize=14,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,REBIN(Inten0,1,1,nseq),xst=1,tit='!17',xtit='position number',ytit='intensity (a.u.)',charsize=1.5
OPLOT,REBIN(Inten2,1,1,nseq),linestyle=2,color=180
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten0,1,1,nseq))-REFORM(REBIN(Inten2,1,1,nseq)))^2.0 )/SQRT(TOTAL(REBIN(Inten2,1,1,nseq)^2.0))


SET_PLOT,'X'
READ,pause

END
