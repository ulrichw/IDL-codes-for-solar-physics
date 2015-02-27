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
; I ALSO TRY TO FIT FOR THE FSRs OF NB, WB, and E1
; USING 2 DETUNES TAKEN AT DIFFERENT TIME OF THE DAY
;
; version 1.4 December 5, 2006
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

PRO HMI_sun1_FSR,draw



lamref      = 6173.3433 ; reference central wavelength for the FeI 6173 line IN AIR
time0       = SYSTIME(1); variable used to determine the elapsed time
dpi         = 2.d0*!dpi


; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

nseq        = 27               ; number of filtergrams in the detune sequence
nx          = 256                 ; number of rows (we rebin the filtergrams)
ny          = 256                 ; number of columns
anglim      = 2048.;970.                ; maximum radius of the solar disk in arcseconds
factor      = 10000.              ; we divide the intensity by this value to make the fit easier

; MEASURED INTENSITIES
; these are the input data: the list must contains
; the nseq+2 names of the filtergrams to be analyzed
; (+2 because of the dark frames)
;------------------------------------------------
;list1='listSun070217_66956'   ; in Calmode
;list2='listSun070217_66675'   ; in Calmode
;IF list1 EQ 'listSun070217_66956' THEN lam1 = 570.25230 ; in m/s for 02/17/2007 at 22:47:00 UT
;IF list2 EQ 'listSun070217_66675' THEN lam2 = 350.49657 ; in m/s for 02/17/2007 at 20:19:00 UT

list1='listSun070302_102330'
list2='listSun070302_103436'
IF list1 EQ 'listSun070302_102330'THEN lam1 = 207.25771 ; in m/s for 03/02/2007 at 17:55:00 UT
IF list2 EQ 'listSun070302_103436'THEN lam2 = 648.81597 ; in m/s for 03/02/2007 at 22:50:00 UT


; THE FOLLOWING LINE LOCATES THE CENTER OF THE SOLAR DISK
; ON THE FILTERGRAMS. IT HAS TO BE ADJUSTED:
center   = [nx/2,ny/2]
distance = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; READ DETUNE SEQUENCE
; the detune sequence, imx, must have been produced using the codes
; readimages.pro and clean.pro from J. Schou, and must have been saved
; in the idl binary format with the name:
; 'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
; where list contains the list of filtergrams to read and nx is the
; rebinned size of the filtergrams. The binary file must also contain
; the variable lam0, which is obtained from my code sunearthvel.pro
;---------------------------------------------------------------------------------

Inten   = FLTARR(nx,nx,2*nseq)

RESTORE,'SEQUENCE_'+STRTRIM(list1,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
dlamdv  = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam1    = dlamdv*lam1         ; positive velocities increase the wavelength
Inten[*,*,0:nseq-1]        = imx[*,*,2:28]/factor ; for data after February 2007 obtained with readimages2.pro
RESTORE,'SEQUENCE_'+STRTRIM(list2,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
dlamdv  = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam2    = dlamdv*lam2         ; positive velocities increase the wavelength
Inten[*,*,nseq:2*nseq-1]   = imx[*,*,2:28]/factor ; for data after February 2007 obtained with readimages2.pro

wobble:
;nseq=36
;inten=imx

nx=1
ny=1
Inten1=REBIN(inten[*,*,0:nseq-1],1,1,nseq)
Inten2=REBIN(inten[*,*,nseq:2*nseq-1],1,1,nseq)
ricardo=MEAN(Inten1)
Inten = FLTARR(nx,ny,2*nseq)
Inten[0,0,*] = [REFORM(Inten1),REFORM(Inten2/MEAN(Inten2))*ricardo]

;nx=64
;ny=64
;Inten=REBIN(inten,nx,ny,nseq)
;distance=rebin(distance,nx,ny)

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


Bg          = DBLARR(nx,ny,3)     ; contrasts of the tunable elements
Phig        = DBLARR(nx,ny,3)     ; relative phases
FSRg        = DBLARR(nx,ny,3)     ; FSR of NB and WB
continuum   = DBLARR(nx,ny)       ; relative sun continuum intensity
linewidth   = DBLARR(nx,nx)       ; solar linewidth
linedepth   = DBLARR(nx,nx)       ; solar linedepth
mu0         = 0.d0                ; regularization parameter for the least-squares fit. OPTIONAL: DEPENDS ON THE NOISE LEVEL
nparam      = 8;9                   ; number of parameters we fit for
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelengths in Angstrom RELATIVE TO lamref
Residual    = DBLARR(2*nseq)        ; Residual of the least-squares fit
dIntendI0   = DBLARR(2*nseq)        ; derivatives to compute the Jacobian matrix
dIntendIc   = DBLARR(2*nseq)
dIntendw    = DBLARR(2*nseq)
dIntendPhi0 = DBLARR(2*nseq)
dIntendPhi1 = DBLARR(2*nseq)
dIntendPhi2 = DBLARR(2*nseq)
dIntendFSR0 = DBLARR(2*nseq)        ; NB
dIntendFSR1 = DBLARR(2*nseq)        ; WB
;dIntendFSR2 = DBLARR(2*nseq)        ; WB
dlinedIc    = 1.d0
tuning      = DBLARR(3,nseq)      ; tuning sequence
maxsteps    = 105                 ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(nx,ny)       ; to discard spurious results

contrasts   = 0.97d0              ; guess for Bg
depth0      = 0.52                ; guess for the solar linedepth
width0      = 0.073               ; guess for the solar linewidth
thresh      = 75000./factor;41000./factor       ; guess of the continuum intensity

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.693d0;0.690d0  ; for E1
FSR[3]      = 1.405d0               ; for E2
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
    
; NON-TUNABLE FILTER PROFILE  PROVIDED BY LOCKHEED-MARTIN
;-----------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

q            = READFITS('blocker11.fits') ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
                                          ; this profile, obtained with a Cary, is centered on 6169.8 A, but measures made by Andover
                                          ; with a spectrograph show the center on 6172 A instead of 6173.3433. Also
                                          ; we shift the profile by 1.3433 A
;q[*,1]       = SHIFT(q[*,1],2.75/(q[1,0]-q[0,0])) ; careful, the shift depends on the focal distance
;blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-lamref,lam) ; slightly uncentered blocker
;blocker0     = blocker0 *INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam)  ;centered blocker
 blocker0     = blocker0*INTERPOL(q[*,1]/100.d0,q[*,0]+2.d0-lamref,lam) ;,/LSQUADRATIC) ; uncentered blocker (correction of blocker11.fits because of f8 instead of collimated light)

;RESTORE,'Lyot.bin' ; transmission profile of E2-E5 obtained in sunlight on Oct 16, 2006 (see LyotShine_Sun.pro) from lyotsunsetSprofiles2006oct16.fits
;RESTORE,'Lyot2.bin' ; transmission profile of E2-E5 obtained in lamp light on Oct 16, 2006 (see LyotShine_Lamp.pro) from lyotLampsetBprofiles2006oct16.fits
;blocker0     = INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lamref,lam) ;,/LSQUADRATIC) ; uncentered blocker (correction of blocker11.fits because of f8 instead of collimated light)
;blocker0     = blocker0*INTERPOL(data,lresp-0.0075,lam) ; contains the blocker+Lyot+front window profile


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
            fdepthg         = depth0*thresh;0.62d0*Icg     ; depth of the solar line according to Stenflo & Lindegren (1977)
            fwidthg         = width0;2880.d0*dlamdv        ; width of the solar line according to Stenflo & Lindegren (1977)
            Bg[iii,jjj,*]   = contrasts
            Phig[iii,jjj,*] = [-113.86,29.97,-143.05]*!dpi/180./FSR[0:2]
            FSRg[iii,jjj,*] = FSR[0:2]
          
            converg         = 0
            jj              = 0
           ;mu              = mu0
            
            profilef        = blocker0 ; we use the averaged blocker+front window profile
            FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0
            
            WHILE (converg EQ 0) DO BEGIN
                
                FOR j=0,2*nseq-1 DO BEGIN 
                
                    IF j LE nseq-1 THEN BEGIN
                        templine= EXP(-(lam-lam1)^2.d0/fwidthg^2.d0)
                        line    = Icg -fdepthg*templine
                        dlinedI0=     -        templine
                        dlinedw =     -fdepthg*templine*(lam-lam1)^2.d0*2.d0/fwidthg^3.d0
                    ENDIF ELSE BEGIN
                        templine= EXP(-(lam-lam2)^2.d0/fwidthg^2.d0)
                        line    = Icg -fdepthg*templine
                        dlinedI0=     -        templine
                        dlinedw =     -fdepthg*templine*(lam-lam2)^2.d0*2.d0/fwidthg^3.d0
                    ENDELSE
    
                    profileg      = 0.125d0 * profilef * (1.d0+Bg[iii,jjj,0]*COS(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq])) * (1.d0+Bg[iii,jjj,1]*COS(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq])) * (1.d0+Bg[iii,jjj,2]*COS(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))
                  
                    Residual[j]   = Inten[iii,jjj,j] - TOTAL(line*profileg)*dlam
                    
                    dIntendI0[j]  = TOTAL(dlinedI0*profileg)*dlam
                    dIntendw [j]  = TOTAL(dlinedw *profileg)*dlam
                    dIntendIc[j]  = TOTAL(dlinedIc*profileg)*dlam
                    dIntendPhi0[j]= TOTAL(line*(-Bg[iii,jjj,0]*FSRg[iii,jjj,0]*SIN(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq]))*(1.d0+Bg[iii,jjj,1]*COS(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq]))*(1.d0+Bg[iii,jjj,2]*COS(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))/8.d0*profilef)*dlam
                    dIntendPhi1[j]= TOTAL(line*(-Bg[iii,jjj,1]*FSRg[iii,jjj,1]*SIN(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq]))*(1.d0+Bg[iii,jjj,0]*COS(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq]))*(1.d0+Bg[iii,jjj,2]*COS(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))/8.d0*profilef)*dlam
                    dIntendPhi2[j]= TOTAL(line*(-Bg[iii,jjj,2]*FSRg[iii,jjj,2]*SIN(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))*(1.d0+Bg[iii,jjj,0]*COS(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq]))*(1.d0+Bg[iii,jjj,1]*COS(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq]))/8.d0*profilef)*dlam
                    dIntendFSR0[j]= TOTAL(line*(-Bg[iii,jjj,0])*(lam+Phig[iii,jjj,0])*SIN(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq])*(1.d0+Bg[iii,jjj,1]*COS(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq])) * (1.d0+Bg[iii,jjj,2]*COS(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))/8.d0*profilef)*dlam
                    dIntendFSR1[j]= TOTAL(line*(-Bg[iii,jjj,1])*(lam+Phig[iii,jjj,1])*SIN(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq])*(1.d0+Bg[iii,jjj,0]*COS(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq])) * (1.d0+Bg[iii,jjj,2]*COS(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq]))/8.d0*profilef)*dlam
                   ;dIntendFSR2[j]= TOTAL(line*(-Bg[iii,jjj,2])*(lam+Phig[iii,jjj,2])*SIN(FSRg[iii,jjj,2]*(lam+Phig[iii,jjj,2])+tuning[2,j MOD nseq])*(1.d0+Bg[iii,jjj,0]*COS(FSRg[iii,jjj,0]*(lam+Phig[iii,jjj,0])+tuning[0,j MOD nseq])) * (1.d0+Bg[iii,jjj,1]*COS(FSRg[iii,jjj,1]*(lam+Phig[iii,jjj,1])+tuning[1,j MOD nseq]))/8.d0*profilef)*dlam
                ENDFOR
                



              ; Jacobian matrix
                Jac  = DBLARR(nparam,2*nseq)
                FOR i= 0,2*nseq-1 DO BEGIN
                    
                    Jac[0,i]  = dIntendI0[i]
                    Jac[1,i]  = dIntendIc[i]
                    Jac[2,i]  = dIntendw[i]
                    Jac[3,i]  = dIntendPhi0[i]
                    Jac[4,i]  = dIntendPhi1[i]
                    Jac[5,i]  = dIntendPhi2[i]
                    Jac[6,i]  = dIntendFSR0[i]
                    Jac[7,i]  = dIntendFSR1[i]
                   ;Jac[8,i]  = dIntendFSR2[i]

                ENDFOR

                LA_SVD,Jac,W,U,V,/DOUBLE
                
              ; regularization:
              ; filter = W*W/(W*W+mu^2.d0) ; mu is the regularization parameter
              ; Dx     = V##DIAG_MATRIX(1.d0/W*filter)##TRANSPOSE(U)##Residual
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
                
                temp   = TRANSPOSE(Dx)/[fdepthg,Icg,fwidthg,REFORM(Phig[iii,jjj,0]),REFORM(Phig[iii,jjj,1]),REFORM(Phig[iii,jjj,2]),REFORM(FSRg[iii,jjj,0]),REFORM(FSRg[iii,jjj,1])];,REFORM(FSRg[iii,jjj,2])]
                err[jj]= MAX(ABS(temp))
                
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
                    FSRg[iii,jjj,0] = history[j[0],6]
                    FSRg[iii,jjj,1] = history[j[0],7]
                    ;FSRg[iii,jjj,2] = history[j[0],8]
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
                    FSRg[iii,jjj,0] = FSRg[iii,jjj,0] +Dx[6]
                    IF(ABS(FSRg[iii,jjj,0]-FSR[0])/FSR[0] GT .1) THEN FSRg[iii,jjj,0]=FSR[0]
                    FSRg[iii,jjj,1] = FSRg[iii,jjj,1] +Dx[7]
                    IF(ABS(FSRg[iii,jjj,1]-FSR[1])/FSR[1] GT .1) THEN FSRg[iii,jjj,1]=FSR[1]
                    ;FSRg[iii,jjj,2] = FSRg[iii,jjj,2] +Dx[8]
                    ;IF(ABS(FSRg[iii,jjj,2]-FSR[2])/FSR[2] GT .1) THEN FSRg[iii,jjj,2]=FSR[2]
                ENDIF
                
                history[jj,*]     = [fdepthg,Icg,fwidthg,Phig[iii,jjj,0],Phig[iii,jjj,1],Phig[iii,jjj,2],FSRg[iii,jjj,0],FSRg[iii,jjj,1]];,FSRg[iii,jjj,2]]
                lateconv[iii,jjj] = converg
                
                jj = jj+1
                
            ENDWHILE
            FOR i=0,2 DO Phig[iii,jjj,i]    = Phig[iii,jjj,i]*FSRg[iii,jjj,i]*180.d0/!dpi ; converted into degree units 
            linewidth[iii,jjj] = fwidthg
            linedepth[iii,jjj] = fdepthg
            continuum[iii,jjj] = Icg
            PRINT,'PARAMETERS='
            FOR i=0,2 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'Ic=',Icg
            PRINT,'fdepthg=',fdepthg
            PRINT,'fwidthg=',fwidthg
            PRINT,'FSRs=',dpi/FSRg[iii,jjj,0],dpi/FSRg[iii,jjj,1],dpi/FSRg[iii,jjj,2]
            PRINT,'lateconv=',converg
            PRINT,iii,jjj

        ENDIF
        
    ENDFOR
ENDFOR

SAVE,Phig,Bg,FSRg,linewidth,linedepth,continuum,lateconv,lam1,lam2,blocker0,lam,FILE='RESULTS/RESULTS_FSR_'+STRTRIM(list1,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;SAVE,Phig,linewidth,linedepth,continuum,lateconv,lam0,FILE='RESULTS/RESULTS_simulation.bin'

draw:

; COMPUTE THE AVERAGE PHASES
;--------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS_FSR_'+STRTRIM(list1,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;RESTORE,'RESULTS/RESULTS_simulation.bin'

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
    FOR i=0,2 DO TVIM,Phig[*,*,i],/scale,range=[moy[i]-2.5*dis[i],moy[i]+2.5*dis[i]],tit='!17'+titl[i],stit='!17Phase (in degrees)',barwidth=0.5,pcharsize=0.75,xtit='Pixel number',ytit='Pixel number'
    LOADCT,3
    TVIM,continuum*distance,/scale,stit='!17Continuum intensity',barwidth=0.5,pcharsize=0.75,range=[MIN(continuum[a]),MAX(continuum[a])],xtit='Pixel number',ytit='Pixel number'
    TVIM,linewidth*distance,/scale,stit='!17Linewidth (A)',barwidth=0.5,pcharsize=0.75,range=[MIN(linewidth[a]),MAX(linewidth[a])],xtit='Pixel number',ytit='Pixel number'
    TVIM,(linedepth/continuum)*distance,/scale,stit='!17Linedepth',barwidth=0.5,pcharsize=0.75,range=[MIN(linedepth[a]/continuum[a]),MAX(linedepth[a]/continuum[a])],xtit='Pixel number',ytit='Pixel number'
    DEVICE,/close
    !P.MULTI=0
ENDIF


; RECONSTRUCTION OF THE DETUNE SEQUENCE
;--------------------------------------------------------------------------------------


FOR i=0,2 DO Phig[*,*,i] = Phig[*,*,i]/FSRg[*,*,i]/180.d0*!dpi ; convert from degrees to radians/FSR

FOR i=0,2*nseq-1 DO Inten[*,*,i] = Inten[*,*,i]*distance
Inten0       = Inten
a            = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0

profilef     = blocker0            ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0

Inten2  = FLTARR(nx,nx,2*nseq) 
FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,nx-1 DO BEGIN
        IF(distance[i,j] NE 0.0) THEN BEGIN

                        
            FOR k=0,2*nseq-1 DO BEGIN

                IF k LE nseq-1 THEN BEGIN
                    templine     = EXP(-(lam-lam1)^2.d0/linewidth[i,j]^2.d0)
                    line         = continuum[i,j] -linedepth[i,j]*templine
                ENDIF ELSE BEGIN
                    templine     = EXP(-(lam-lam2)^2.d0/linewidth[i,j]^2.d0)
                    line         = continuum[i,j] -linedepth[i,j]*templine
                ENDELSE


                profileg      = 0.125d0 * profilef * (1.d0+Bg[i,j,0]*COS(FSRg[i,j,0]*(lam+Phig[i,j,0])+tuning[0,k MOD nseq])) * (1.d0+Bg[i,j,1]*COS(FSRg[i,j,1]*(lam+Phig[i,j,1])+tuning[1,k MOD nseq])) * (1.d0+Bg[i,j,2]*COS(FSR[2]*(lam+Phig[i,j,2])+tuning[2,k MOD nseq]))
                Inten2[i,j,k] = TOTAL(line*profileg)*dlam
            ENDFOR

        ENDIF
    ENDFOR
ENDFOR


; WE PLOT THE MAPS OF THE RELATIVE RECONSTRUCTION ERROR
!P.MULTI=[0,2,3]
loadct,4
DEVICE,file='yo2.ps',xoffset=-0.8,yoffset=0,xsize=22,ysize=27,/color,bits=24

FOR i=0,2*nseq-1 DO BEGIN
    temp  = REFORM((Inten2[*,*,i]-Inten[*,*,i])/Inten[*,*,i])*distance
    a=WHERE(FINITE(temp) EQ 0)
    IF(a[0] NE -1) THEN temp[a]=0.0
    tempa = MEAN(temp[WHERE(temp NE 0.0)])
    TVIM,temp,range=[-0.03+tempa,0.03+tempa],/scale,tit=STRING(i+1),stit='Relative error',pcharsize=1.5
ENDFOR


FOR i=0,2*nseq-1 DO BEGIN
    temp = REFORM((Inten2[*,*,i]-Inten[*,*,i])/Inten[*,*,i])
    aa   = WHERE(temp LT 0.2 AND temp GT -0.2,na,COMPLEMENT=b)
    hist = histogram(temp[aa],binsize=0.001,min=-0.2,max=0.2)
    plot,FINDGEN(N_ELEMENTS(hist))*0.001-0.2,hist/FLOAT(na),psym=10,tit='!7r!17='+STRING(SIGMA(temp[aa]))+'/!7l!17='+STRING(MEAN(temp[aa])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5
    PRINT,MIN(temp[aa]),MAX(temp[aa]),SIGMA(temp[aa]),MEAN(temp[aa])
ENDFOR
DEVICE,/CLOSE


!P.MULTI=0
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=15,ysize=12,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,REBIN(Inten0,1,1,2*nseq),xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5
OPLOT,REBIN(Inten2,1,1,2*nseq),linestyle=2,color=180
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten0,1,1,2*nseq))-REFORM(REBIN(Inten2,1,1,2*nseq)))^2.0 )


SET_PLOT,'X'
READ,pause

END
