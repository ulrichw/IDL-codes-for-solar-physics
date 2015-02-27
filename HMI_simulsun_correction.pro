; PROGRAM TO SIMULATE THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW THE ACCURACY WE CAN REACH AS A FUNCTION OF
; THE NUMBER OF POSITIONS IN THE DETUNE
;
; PROGRAM TO PRODUCE DATA FOR HMI_sun1_correction.pro
; THAT TRIES TO CORRECT FOR THE WRONG HCM POSITIONS
;
; ASSUMPTIONS: contrasts and phase of the non-tunable part are 
;              respectively equal to contrasts and 0 
;
; THE SOLAR LINE IS MODELED BY A GAUSSIAN PROFILE:
; I(l) = Ic - Ic*d*EXP(-l^2/vl^2)
; WHERE Ic IS THE CONTINUUM, d IS THE LINEDEPTH
; AND vl IS THE LINEWIDTH
; WE FIT FOR 3 PARAMETERS FOR THE LINE PROFILE:
; Ic, depth=Ic*d, and vl
;
; THE PROGRAM REMOVES THE BAD EXPOSURES AND OVERSCANS
;
; version 1.3 June 12, 2006
;
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
; NB: be aware of the fact that there might be a systematic error on the
; determination of the solar Doppler velocity due to the lack of
; absolute reference in the wavelength scale (we assume that we are
; centered on 6173.3433 A)
;
; NB: IT IS VERY IMPORTANT TO HAVE A GOOD ESTIMATE OF THE THRESH
; I.E. THE EXPOSURE TIME, BECAUSE THAT IMPACTS ON THE PARAMETERS
;
;----------------------------------------------------------------------

PRO HMI_simulsun_correction,seed



lamref      = 6173.3433;reference central wavelength for the FeI 6173 line
                      ; IN AIR
                      ; used to center the averaged blocker
                      ; and front window transmission profiles
                      ; value from Stenflo & Lindegren (1977)
                      ; NB: IF THIS VALUE IS NOT THE ACTUAL ONE, THEN
                      ; WE HAVE AN AVERAGE VELOCITY ON THE SOLAR DISK
                      ; DIFFERENT FROM 0. TO BE CORRECTED
                      ; I have at least an uncertainty of
                      ; +/- 7 mA on this value

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelengths in Angstrom RELATIVE TO lamref

; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

nseq        = 27

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172d0            ; for the narrow-band Michelson
FSR[1]      = 0.344d0            ; for the broad-band  Michelson
FSR[2]      = 0.693d0;0.702d0               ; for E1
FSR[3]      = 1.405d0               ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

FSR         = dpi/FSR

; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN
;-----------------------------------------------

RESTORE,'frontwindow.bin'
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits') 
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.61328-lamref,lam);,/LSQUADRATIC)
;res          = poly_fit(lam,blocker0,2,yfit=smoothblocker) ; "ideal" blocker+front window


; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------

Plpos        = DBLARR(3,nseq) ; tuning sequence
tuning       = Plpos
sign         = tuning
Pldel        = tuning
Plpos[*,0]  = [         0.d0,          0.d0,          0.d0]
Plpos[*,1]  = [        80.d0,          0.d0,          0.d0]
Plpos[*,2]  = [       160.d0,          0.d0,          0.d0]
Plpos[*,3]  = [         0.d0,         80.d0,          0.d0]
Plpos[*,4]  = [        80.d0,         80.d0,          0.d0]
Plpos[*,5]  = [       160.d0,         80.d0,          0.d0]
Plpos[*,6]  = [         0.d0,        160.d0,          0.d0]
Plpos[*,7]  = [        80.d0,        160.d0,          0.d0]
Plpos[*,8]  = [       160.d0,        160.d0,          0.d0]
Plpos[*,9]  = [         0.d0,          0.d0,         80.d0]
Plpos[*,10] = [        80.d0,          0.d0,         80.d0]
Plpos[*,11] = [       160.d0,          0.d0,         80.d0]
Plpos[*,12] = [         0.d0,         80.d0,         80.d0] 
Plpos[*,13] = [        80.d0,         80.d0,         80.d0]
Plpos[*,14] = [       160.d0,         80.d0,         80.d0]
Plpos[*,15] = [         0.d0,        160.d0,         80.d0]
Plpos[*,16] = [        80.d0,        160.d0,         80.d0]
Plpos[*,17] = [       160.d0,        160.d0,         80.d0]
Plpos[*,18] = [         0.d0,          0.d0,        160.d0]
Plpos[*,19] = [        80.d0,          0.d0,        160.d0]
Plpos[*,20] = [       160.d0,          0.d0,        160.d0]
Plpos[*,21] = [         0.d0,         80.d0,        160.d0]
Plpos[*,22] = [        80.d0,         80.d0,        160.d0]
Plpos[*,23] = [       160.d0,         80.d0,        160.d0]
Plpos[*,24] = [         0.d0,        160.d0,        160.d0]
Plpos[*,25] = [        80.d0,        160.d0,        160.d0]
Plpos[*,26] = [       160.d0,        160.d0,        160.d0]

sign[*,0] = [     -1.00000  ,   -1.00000   ,   -1.00000]
sign[*,1] = [     -1.00000  ,   -1.00000   ,   -1.00000]
sign[*,2] = [     -1.00000  ,   -1.00000   ,   -1.00000]
sign[*,3] = [     -1.00000  ,   -1.00000   ,   -1.00000]
sign[*,4] = [     1.00000  ,   -1.00000   ,   -1.00000]
sign[*,5] = [     1.00000  ,   -1.00000   ,   -1.00000]
sign[*,6] = [     1.00000  ,   -1.00000   ,   -1.00000]
sign[*,7] = [     1.00000  ,   1.00000   ,  -1.00000]
sign[*,8] = [     1.00000  ,   1.00000   ,  -1.00000]
sign[*,9] = [     1.00000  ,   1.00000   ,  -1.00000]
sign[*,10] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,11] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,12] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,13] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,14] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,15] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,16] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,17] = [    -1.00000  ,   1.00000   ,  -1.00000]
sign[*,18] = [    -1.00000  ,   1.00000   ,   1.00000]
sign[*,19] = [    -1.00000  ,   1.00000   ,   1.00000]
sign[*,20] = [    -1.00000  ,   1.00000   ,  -1.00000]
sign[*,21] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,22] = [     1.00000  ,  -1.00000   ,   1.00000]
sign[*,23] = [     1.00000  ,  -1.00000   ,  -1.00000]
sign[*,24] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,25] = [     1.00000  ,   1.00000   ,   1.00000]
sign[*,26] = [     1.00000  ,   1.00000   ,   1.00000]

Pldel[*,0] = [      3123.00  ,    3978.00   ,   4000.00]
Pldel[*,1] = [      3785.00  ,    3978.00   ,   4000.00]
Pldel[*,2] = [      3411.00  ,    3978.00   ,   4000.00]
Pldel[*,3] = [      3379.00  ,    3872.00   ,   4000.00]
Pldel[*,4] = [      3784.00  ,    3872.00   ,   4000.00]
Pldel[*,5] = [      3155.00  ,    3872.00   ,   4000.00]
Pldel[*,6] = [      3251.00  ,    3827.00   ,   4000.00]
Pldel[*,7] = [      3785.00  ,    3827.00   ,   4000.00]
Pldel[*,8] = [      3283.00  ,    3827.00   ,   4000.00]
Pldel[*,9] = [      2995.00  ,    3967.00   ,   3955.00]
Pldel[*,10] = [      3784.00  ,    3967.00   ,   3955.00]
Pldel[*,11] = [      3219.00  ,    3967.00   ,   3955.00]
Pldel[*,12] = [      3123.00  ,    4000.00   ,   3955.00]
Pldel[*,13] = [      3785.00  ,    4000.00   ,   3955.00]
Pldel[*,14] = [      3251.00  ,    4000.00   ,   3955.00]
Pldel[*,15] = [      3379.00  ,    4083.00   ,   3955.00]
Pldel[*,16] = [      3784.00  ,    4083.00   ,   3955.00]
Pldel[*,17] = [      3507.00  ,    4083.00   ,   3955.00]
Pldel[*,18] = [      3251.00  ,    3973.00   ,   3955.00]
Pldel[*,19] = [      3785.00  ,    3973.00   ,   3955.00]
Pldel[*,20] = [      3379.00  ,    3973.00   ,   3955.00]
Pldel[*,21] = [      2995.00  ,    3936.00   ,   3955.00]
Pldel[*,22] = [      3784.00  ,    3936.00   ,   3955.00]
Pldel[*,23] = [      3123.00  ,    3936.00   ,   3955.00]
Pldel[*,24] = [      3123.00  ,    3955.00   ,   3955.00]
Pldel[*,25] = [      3528.00  ,    3955.00   ,   3955.00]
Pldel[*,26] = [      3251.00  ,    3955.00   ,   3955.00]



; SIMULATED INTENSITIES
;------------------------------------------------

; INITIAL PARAMETERS
B            = [0.97,0.985,0.97,0.92,0.96,0.98,0.96]
Phi          = [-12.0,14.56,13.20,-13.00,-8.4,-24.7,-75.1]*!dpi/180./FSR[*]
Ic           = 120000.0
fdepth       = 0.4*Ic
fwidth       = 0.06851d0
lam0         = 0.d0
line         = Ic - fdepth*EXP(-(lam -lam0)^2.d0/fwidth^2.d0)
Inten0       = FLTARR(nseq)
thcm         = 5000.0 
chcm         = 0.0004d0
thcm2        = 5000.0 
chcm2        = 0.00025d0
thcm3        = 5000.0 
chcm3        = 0.000634d0

FOR j=0,nseq-1 DO BEGIN
    tuning[0,j]=Plpos[0,j]+(Pldel[0,j]-thcm) *sign[0,j]*chcm
    tuning[1,j]=Plpos[1,j]+(Pldel[1,j]-thcm2)*sign[1,j]*chcm2
    tuning[2,j]=Plpos[2,j]+(Pldel[2,j]-thcm3)*sign[2,j]*chcm3
ENDFOR

FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0

FOR j=0,nseq-1 DO BEGIN
    profile = blocker0
    FOR i=0,2 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam+Phi[i])+tuning[i,j]))/2.d0
    FOR i=3,6 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam+Phi[i])) )/2.d0
    Inten0  [j] = TOTAL(line*profile)*dlam
ENDFOR

SET_PLOT,'x'
!P.MULTI=0
WINDOW,0,RETAIN=2
PLOT,inten0,xst=1
; ADD NOISE
;Inten0 = Inten0 + 0.00282843*MEAN(Inten0)*RANDOMN(seed,nseq)
OPLOT,inten0,color=180
SAVE,Plpos,Pldel,sign,Inten0,lam0,file='test.bin'

END
