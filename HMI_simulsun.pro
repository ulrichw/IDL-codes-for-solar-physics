; PROGRAM TO SIMULATE THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW THE ACCURACY WE CAN REACH AS A FUNCTION OF
; THE NUMBER OF POSITIONS IN THE DETUNE
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

PRO HMI_simulsun,seed



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

nseq        = 27;21

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.7039;0.702d0               ; for E1
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
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam);,/LSQUADRATIC)
res          = poly_fit(lam,blocker0,2,yfit=smoothblocker) ; "ideal" blocker+front window


; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------

IF nseq EQ 27 THEN BEGIN
    tuning       = DBLARR(3,nseq) ; tuning sequence
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
ENDIF ELSE BEGIN
    tuning       = DBLARR(3,nseq) ; tuning sequence
    tuning[*,0]  = [         0.d0,          0.d0,          0.d0]
    tuning[*,1]  = [        80.d0,          0.d0,          0.d0]
    tuning[*,2]  = [       160.d0,          0.d0,          0.d0]
    tuning[*,3]  = [         0.d0,         80.d0,          0.d0]
    tuning[*,4]  = [        80.d0,         80.d0,          0.d0]
   ;tuning[*,5]  = [       160.d0,         80.d0,          0.d0]
    tuning[*,5]  = [         0.d0,        160.d0,          0.d0]
   ;tuning[*,7]  = [        80.d0,        160.d0,          0.d0]
    tuning[*,6]  = [       160.d0,        160.d0,          0.d0]
    tuning[*,7]  = [         0.d0,          0.d0,         80.d0]
    tuning[*,8] = [        80.d0,          0.d0,         80.d0]
   ;tuning[*,11] = [       160.d0,          0.d0,         80.d0]
    tuning[*,9] = [         0.d0,         80.d0,         80.d0] 
    tuning[*,10] = [        80.d0,         80.d0,         80.d0]
    tuning[*,11] = [       160.d0,         80.d0,         80.d0]
   ;tuning[*,15] = [         0.d0,        160.d0,         80.d0]
    tuning[*,12] = [        80.d0,        160.d0,         80.d0]
    tuning[*,13] = [       160.d0,        160.d0,         80.d0]
    tuning[*,14] = [         0.d0,          0.d0,        160.d0]
   ;tuning[*,19] = [        80.d0,          0.d0,        160.d0]
    tuning[*,15] = [       160.d0,          0.d0,        160.d0]
   ;tuning[*,21] = [         0.d0,         80.d0,        160.d0]
    tuning[*,16] = [        80.d0,         80.d0,        160.d0]
    tuning[*,17] = [       160.d0,         80.d0,        160.d0]
    tuning[*,18] = [         0.d0,        160.d0,        160.d0]
    tuning[*,19] = [        80.d0,        160.d0,        160.d0]
    tuning[*,20] = [       160.d0,        160.d0,        160.d0]
ENDELSE



    FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0


; SIMULATED INTENSITIES
;------------------------------------------------

; INITIAL PARAMETERS
B            = [0.85,0.91,0.875,0.923,0.865,0.9454,0.9023]
Phi          = [-15.,23.,12.,5.,-11.,21.34,-14.56]*!dpi/180./FSR[*]
Ic           = 51974.024
fdepth       = 24043.43
fwidth       = 0.06851d0
lam0         = 0.d0
line         = Ic - fdepth*EXP(-(lam -lam0)^2.d0/fwidth^2.d0)
Inten0       = FLTARR(nseq)

FOR j=0,nseq-1 DO BEGIN
    profile = blocker0
    FOR i=0,2 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam+Phi[i])+tuning[i,j]))/2.d0
    FOR i=3,6 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam+Phi[i])) )/2.d0
    Inten0  [j] = TOTAL(line*profile)*dlam
ENDFOR

; on the dark frame /tmp20/schou/hmi060224b/i_060224_225806.fit
; after the linear trend is removed, the dark current is fitted by
; A[2]*exp(-(x-A[0])^2.0/2.0/A[1]^2.0), where A[1] = 84.86 (with
; NO averaging, i.e. in 4096*4096)
; and A[1]=14.00 (with averaging in 256*256)
; the dark current follows a Gaussian distribution in very good approximation
; for this detune sequence in Calmode and Sunlight, Ic=51974.024,
; fwidth=0.068510673, and fdepthg=24043.43

; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------

Bg          = DBLARR(3)     ; contrasts of the tunable elements
Phig        = DBLARR(3)     ; relative phases
continuum   = DBLARR(1)       ; relative sun intensity
linewidth   = DBLARR(1)       ; solar linewidth
linedepth   = DBLARR(1)       ; solar linedepth
mu0         = 0.d0 ;0.02          ; regularization parameter for the least-squares fit. OPTIONAL: DEPENDS ON THE NOISE LEVEL
nparam      = 6                   ; number of parameters we fit for
nrepet      = 2000
error       = FLTARR(nrepet,nparam)
Residual    = DBLARR(nseq)        ; Residual of the least-squares fit
dIntendI0   = DBLARR(nseq)        ; derivatives to compute the Jacobian matrix
dIntendIc   = DBLARR(nseq)
dIntendw    = DBLARR(nseq)
dIntendPhi0 = DBLARR(nseq)
dIntendPhi1 = DBLARR(nseq)
dIntendPhi2 = DBLARR(nseq)
dlinedIc    = 1.d0
maxsteps    = 105                 ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(1)       ; to discard spurious results


photons     = 680000.d0           ; electrons per seconds detected by a pixel of the CCD
                                  ; estimate by Jesper Schou (Radiometric Analysis document)
exposure    = 0.155;0.2d0         ; integration time
gain        = 19.d0               ; CCD inverse gain
eqwidth     = 0.045d0             ; equivalent width of the blocker + window + Lyot + Michelson profile
                                  ; when contrasts=1 and phases=0
;eqwidth    = 0.0571d0            ; approximate value if no front window !!!!!!!!!!!!!!!!!!!!!!!
thresh      = photons/gain*exposure/eqwidth*0.55 ;*0.55 absorption by the heliostat ?
PRINT,'ESTIMATE OF Ic=',thresh

;!!! NB: CHANGE THE ABOVE VALUES (heliostat absorption and eqwidth) IF CODE DOES NOT CONVERGE !!!!
; NB integral(blocker+front_window+Lyot+Michelson * dlam) = 0.05d0 si
; B=1 et phi=0
contrasts   = 0.9d0

; LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=900
!P.MULTI=[0,2,3]


; GUESS
;------------------------

FOR repet=0,nrepet-1 DO BEGIN
    PRINT,repet

; ADDITION OF PHOTON NOISE AND DARK CURRENT
    Inten   = Inten0 + 14.00*RANDOMN(seed,nseq)

    
    depth0  = 0.49
    width0  = 0.071
    Icg     = thresh            ; estimate of the solar continuum
    fdepthg = depth0*thresh ;0.62d0*Icg     ; depth of the solar line according to Stenflo & Lindegren (1977)
    fwidthg = width0 ;2880.d0*dlamdv ; width of the solar line according to Stenflo & Lindegren (1977)
    Bg[*]   = contrasts
    Phig[*] = 0.0
    
    converg = 0
    jj      = 0
    
    profilef= smoothblocker ; we use the averaged blocker+front window profile
    FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0
    
    WHILE (converg EQ 0) DO BEGIN
        
        templine= EXP(-(lam-lam0)^2.d0/fwidthg^2.d0)
        line    = Icg -fdepthg*templine
        dlinedI0=     -        templine
        dlinedw =     -fdepthg*templine*(lam-lam0)^2.d0*2.d0/fwidthg^3.d0
        
        FOR j=0,nseq-1 DO BEGIN
            
            profileg      = 0.125d0 * profilef * (1.d0+Bg[0]*COS(FSR[0]*(lam+Phig[0])+tuning[0,j])) * (1.d0+Bg[1]*COS(FSR[1]*(lam+Phig[1])+tuning[1,j])) * (1.d0+Bg[2]*COS(FSR[2]*(lam+Phig[2])+tuning[2,j]))
            
            Residual[j]   = Inten[j] - TOTAL(line*profileg)*dlam
            
            dIntendI0[j]  = TOTAL(dlinedI0*profileg)*dlam
            dIntendw [j]  = TOTAL(dlinedw *profileg)*dlam
            dIntendIc[j]  = TOTAL(dlinedIc*profileg)*dlam
            dIntendPhi0[j]= TOTAL(line*(-Bg[0]*FSR[0]*SIN(FSR[0]*(lam+Phig[0])+tuning[0,j]))*(1.d0+Bg[1]*COS(FSR[1]*(lam+Phig[1])+tuning[1,j]))*(1.d0+Bg[2]*COS(FSR[2]*(lam+Phig[2])+tuning[2,j]))/8.d0*profilef)*dlam
            dIntendPhi1[j]= TOTAL(line*(-Bg[1]*FSR[1]*SIN(FSR[1]*(lam+Phig[1])+tuning[1,j]))*(1.d0+Bg[0]*COS(FSR[0]*(lam+Phig[0])+tuning[0,j]))*(1.d0+Bg[2]*COS(FSR[2]*(lam+Phig[2])+tuning[2,j]))/8.d0*profilef)*dlam
            dIntendPhi2[j]= TOTAL(line*(-Bg[2]*FSR[2]*SIN(FSR[2]*(lam+Phig[2])+tuning[2,j]))*(1.d0+Bg[0]*COS(FSR[0]*(lam+Phig[0])+tuning[0,j]))*(1.d0+Bg[1]*COS(FSR[1]*(lam+Phig[1])+tuning[1,j]))/8.d0*profilef)*dlam
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
        
        temp   = TRANSPOSE(Dx)/[fdepthg,Icg,fwidthg,Phig[0],Phig[1],Phig[2]]
                                ;err[jj]= TOTAL(temp^2.d0)/DOUBLE(nparam)
        err[jj]= MAX(ABS(temp))
        
                                ; IF(FINITE(err[jj]) EQ 0 AND jj GT 0)          THEN mu = 2.5d0*mu ELSE mu = mu0
        IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-8 )    THEN converg = 2 ; 1.d-14 
        IF(err[jj] LE 1.d-8)                          THEN converg = 1
        
        IF(converg EQ 2) THEN BEGIN
            j = WHERE(err EQ MIN(err))
            fdepthg = history[j[0],0]
            Icg     = history[j[0],1]
            fwidthg = history[j[0],2]
            Phig[0] = history[j[0],3]
            Phig[1] = history[j[0],4]
            Phig[2] = history[j[0],5]
        ENDIF
        
        IF(converg NE 2) THEN BEGIN
            
            Icg     = Icg    +Dx[1]
            IF(Icg LT 0.2d0*thresh OR Icg GT 3.0d0*thresh) THEN Icg = thresh
            fdepthg = fdepthg+Dx[0]
            IF(fdepthg LT 0.2d0*thresh*depth0 OR fdepthg GT 3.0d0*depth0*thresh) THEN fdepthg = depth0*thresh
            fwidthg = fwidthg+Dx[2]
            IF(fwidthg GT 5.d0*width0 OR fwidthg LT 0.3d0*width0) THEN fwidthg = width0
            Phig[0] = Phig[0]+Dx[3]
            IF(ABS(Phig[0]) GT 1.25d0*!dpi/FSR[0]) THEN Phig[0] = 0.
            Phig[1] = Phig[1]+Dx[4]
            IF(ABS(Phig[1]) GT 1.25d0*!dpi/FSR[1]) THEN Phig[1] = 0. 
            Phig[2] = Phig[2]+Dx[5]
            IF(ABS(Phig[2]) GT 1.25d0*!dpi/FSR[2]) THEN Phig[2] = 0.
                    
        ENDIF
        
        history[jj,*]     = [fdepthg,Icg,fwidthg,Phig[0],Phig[1],Phig[2]]
        lateconv          = converg
        
                                ;PRINT,jj,err[jj]
        jj = jj+1
        
    ENDWHILE
    
    error[repet,0]    = (Phig[0]-Phi[0])*FSR[0]*180.d0/!dpi ; converted into degree units
    error[repet,1]    = (Phig[1]-Phi[1])*FSR[1]*180.d0/!dpi ; converted into degree units
    error[repet,2]    = (Phig[2]-Phi[2])*FSR[2]*180.d0/!dpi ; converted into degree units
    error[repet,3]    = (fwidthg-fwidth)/fwidth
    error[repet,4]    = (fdepthg-fdepth)/fdepth
    error[repet,5]    = (Icg-Ic)/Ic

    plot,error[*,0],charsize=1.5,tit=STRING(SIGMA(error[*,0]))
    plot,error[*,1],charsize=1.5,tit=STRING(SIGMA(error[*,1]))
    plot,error[*,2],charsize=1.5,tit=STRING(SIGMA(error[*,2]))
    plot,error[*,3],charsize=1.5,tit=STRING(SIGMA(error[*,3]))
    plot,error[*,4],charsize=1.5,tit=STRING(SIGMA(error[*,4]))
    plot,error[*,5],charsize=1.5,tit=STRING(SIGMA(error[*,5]))

ENDFOR

SET_PLOT,'ps'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=21,ysize=26
    plot,error[*,0],charsize=1.5,tit=STRING(SIGMA(error[*,0]))+STRING(MEAN(error[*,0]))
    plot,error[*,1],charsize=1.5,tit=STRING(SIGMA(error[*,1]))+STRING(MEAN(error[*,1]))
    plot,error[*,2],charsize=1.5,tit=STRING(SIGMA(error[*,2]))+STRING(MEAN(error[*,2]))
    plot,error[*,3],charsize=1.5,tit=STRING(SIGMA(error[*,3]))+STRING(MEAN(error[*,3]))
    plot,error[*,4],charsize=1.5,tit=STRING(SIGMA(error[*,4]))+STRING(MEAN(error[*,4]))
    plot,error[*,5],charsize=1.5,tit=STRING(SIGMA(error[*,5]))+STRING(MEAN(error[*,5]))
DEVICE,/CLOSE
SET_PLOT,'x'

READ,pause

END
