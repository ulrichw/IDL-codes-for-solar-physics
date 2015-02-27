; PROGRAM TO SIMULATE THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW THE ACCURACY WE CAN REACH AS A FUNCTION OF
; THE NUMBER OF POSITIONS IN THE TUNING SEQUENCE
; PLUS: WE DO NOT FIT FOR DEFECTS IN THE TUNING HALF-WAVE PLATES
; UNLIKE HMI_simullasermore.pro
;
; version 1.3 November 26, 2006
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

PRO HMI_simullasermore2,seed



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
dlam        = 0.001d0          ; resolution in wavelength
lam0        = 0.d0

; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

ntun        = 2
nseq        = 23;27;ntun*9

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.7039d0              ; for E1
FSR[3]      = 1.405d0               ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

FSR         = dpi/FSR

; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN
;-----------------------------------------------

RESTORE,'frontwindow.bin'
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam0)
q            = READFITS('blocker11.fits') 
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam0)


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
   ;tuning[*,3]  = [         0.d0,         80.d0,          0.d0]
    tuning[*,3]  = [        80.d0,         80.d0,          0.d0]
    tuning[*,4]  = [       160.d0,         80.d0,          0.d0]
    tuning[*,5]  = [         0.d0,        160.d0,          0.d0]
    tuning[*,6]  = [        80.d0,        160.d0,          0.d0]
    tuning[*,7]  = [       160.d0,        160.d0,          0.d0]
    tuning[*,8]  = [         0.d0,          0.d0,         80.d0]
   ;tuning[*,9] =  [        80.d0,          0.d0,         80.d0]
    tuning[*,9]  = [       160.d0,          0.d0,         80.d0]
    tuning[*,10] = [         0.d0,         80.d0,         80.d0]
    tuning[*,11] = [        80.d0,         80.d0,         80.d0]
    tuning[*,12] = [       160.d0,         80.d0,         80.d0]
    tuning[*,13] = [         0.d0,        160.d0,         80.d0]
    tuning[*,14] = [        80.d0,        160.d0,         80.d0]
    tuning[*,15] = [       160.d0,        160.d0,         80.d0]
   ;tuning[*,18] = [         0.d0,          0.d0,        160.d0]
    tuning[*,16] = [        80.d0,          0.d0,        160.d0]
    tuning[*,17] = [       160.d0,          0.d0,        160.d0]
    tuning[*,18] = [         0.d0,         80.d0,        160.d0]
    tuning[*,19] = [        80.d0,         80.d0,        160.d0]
    tuning[*,20] = [       160.d0,         80.d0,        160.d0]
    tuning[*,21] = [         0.d0,        160.d0,        160.d0]
    tuning[*,22] = [        80.d0,        160.d0,        160.d0]
   ;tuning[*,22] = [       160.d0,        160.d0,        160.d0]   

    ;tuning = DBLARR(3,ntun*9) ; tuning sequence
    ;FOR i=0,ntun-1 DO tuning[*,i+0*ntun] = [i*24.,0.d0  ,  0.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+1*ntun] = [i*24.,0.d0  , 80.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+2*ntun] = [i*24.,0.d0  ,160.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+3*ntun] = [i*24.,80.d0 ,  0.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+4*ntun] = [i*24.,80.d0 , 80.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+5*ntun] = [i*24.,80.d0 ,160.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+6*ntun] = [i*24.,160.d0,  0.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+7*ntun] = [i*24.,160.d0, 80.d0]
    ;FOR i=0,ntun-1 DO tuning[*,i+8*ntun] = [i*24.,160.d0,160.d0]
ENDELSE



    FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0]*!dpi/180.d0


; SIMULATED INTENSITIES
;------------------------------------------------

; INITIAL PARAMETERS
B            = [0.95,0.95,0.95,0.923,0.865,0.9454,0.9023]
Phi          = [ 0.,-120.,-139.,5.,-11.,21.34,-14.56]*!dpi/180./FSR[*]
Ic           = 5.657d6
Inten0       = FLTARR(nseq)
gamma=[1.d0,1.d0,1.d0]

FOR j=0,nseq-1 DO BEGIN
    profile = blocker0
    FOR i=0,2 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam0+Phi[i])+tuning[i,j]))/2.d0
    FOR i=3,6 DO profile = profile * (1.d0+B[i]*COS(FSR[i]*(lam0+Phi[i])) )/2.d0
    Inten0  [j] = Ic*profile*dlam
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
nparam      = 6                   ; number of parameters we fit for
nrepet      = 500
error       = FLTARR(nrepet,nparam)
Residual    = DBLARR(nseq)        ; Residual of the least-squares fit
dIntendPhi0 = DBLARR(nseq)
dIntendPhi1 = DBLARR(nseq)
dIntendPhi2 = DBLARR(nseq)
dIntendB0   = DBLARR(nseq)
dIntendB1   = DBLARR(nseq)
dIntendB2   = DBLARR(nseq)
maxsteps    = 105                 ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(1)       ; to discard spurious results



contrasts   = 0.9d0

; LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=900

; GUESS
;------------------------

FOR repet=0,nrepet-1 DO BEGIN
   ;PRINT,repet

; ADDITION OF PHOTON NOISE AND DARK CURRENT
    noiselev= 10.0;78.62  
    Inten   = Inten0+ noiselev*RANDOMN(seed,nseq)

    Bg[*]   = contrasts
    Phig[*] = 0.	
    converg = 0
    jj      = 0
    
    profilef= blocker0 ; we use the averaged blocker+front window profile
    FOR i=3,6 DO profilef = profilef * (1.d0+B[i]*COS(FSR[i]*(lam0+Phi[i])))/2.d0

    WHILE (converg EQ 0) DO BEGIN
               
        FOR j=0,nseq-1 DO BEGIN
            
            NB = (1.d0+Bg[0]*COS(FSR[0]*(lam0+Phig[0])+tuning[0,j]))/2.d0
            WB = (1.d0+Bg[1]*COS(FSR[1]*(lam0+Phig[1])+tuning[1,j]))/2.d0
            E1 = (1.d0+Bg[2]*COS(FSR[2]*(lam0+Phig[2])+tuning[2,j]))/2.d0
            
            Residual[j]      =  Inten[j] - Ic*profilef*NB*WB*E1*dlam

            dIntendPhi0[j]   =  (-Bg[0]*FSR[0]*SIN(FSR[0]*(lam0+Phig[0])+tuning[0,j])/2.d0) * WB * E1 * profilef * Ic * dlam
            dIntendPhi1[j]   =  (-Bg[1]*FSR[1]*SIN(FSR[1]*(lam0+Phig[1])+tuning[1,j])/2.d0) * NB * E1 * profilef * Ic * dlam
            dIntendPhi2[j]   =  (-Bg[2]*FSR[2]*SIN(FSR[2]*(lam0+Phig[2])+tuning[2,j])/2.d0) * NB * WB * profilef * Ic * dlam

            dIntendB0[j]     =  COS(FSR[0]*(lam0+Phig[0])+tuning[0,j])/2.d0 * WB * E1 * profilef * Ic * dlam
            dIntendB1[j]     =  COS(FSR[1]*(lam0+Phig[1])+tuning[1,j])/2.d0 * NB * E1 * profilef * Ic * dlam
            dIntendB2[j]     =  COS(FSR[2]*(lam0+Phig[2])+tuning[2,j])/2.d0 * NB * WB * profilef * Ic * dlam

        ENDFOR
        
        
                                ; Jacobian matrix
        Jac  = DBLARR(nparam,nseq)
        FOR i= 0,nseq-1 DO BEGIN
            
            Jac[0,i]  = dIntendPhi0[i]
            Jac[1,i]  = dIntendPhi1[i]
            Jac[2,i]  = dIntendPhi2[i]
            Jac[3,i]  = dIntendB0[i]
            Jac[4,i]  = dIntendB1[i]
            Jac[5,i]  = dIntendB2[i]

        ENDFOR
        
        LA_SVD,Jac,W,U,V,/DOUBLE
        
        Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
        
        temp   = TRANSPOSE(Dx)/[Phig[0],Phig[1],Phig[2],Bg[0],Bg[1],Bg[2]];,gammag[0]];,alphag[1],alphag[2]]
        err[jj]= MAX(ABS(temp))
        
        IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-10 )    THEN converg = 2 ; 1.d-14 
        IF(err[jj] LE 1.d-10)                          THEN converg = 1
        
        IF(converg EQ 2) THEN BEGIN
            j = WHERE(err EQ MIN(err))
            Phig[0] = history[j[0],0]
            Phig[1] = history[j[0],1]
            Phig[2] = history[j[0],2]
            Bg[0]   = history[j[0],3]
            Bg[1]   = history[j[0],4]
            Bg[2]   = history[j[0],5]
        ENDIF
        
        IF(converg NE 2) THEN BEGIN
            
            Phig[0] = Phig[0]+Dx[0]
            IF(ABS(Phig[0]) GT 1.25d0*!dpi/FSR[0]) THEN Phig[0] = 0.
            Phig[1] = Phig[1]+Dx[1]
            IF(ABS(Phig[1]) GT 1.25d0*!dpi/FSR[1]) THEN Phig[1] = 0. 
            Phig[2] = Phig[2]+Dx[2]
            IF(ABS(Phig[2]) GT 1.25d0*!dpi/FSR[2]) THEN Phig[2] = 0.
            Bg[0] = Bg[0] +Dx[3]
            IF(Bg[0] GT 1.5d0 OR Bg[0] LT 0.d0) THEN Bg[0] = contrasts
            Bg[1] = Bg[1] +Dx[4]
            IF(Bg[1] GT 1.5d0 OR Bg[1] LT 0.d0) THEN Bg[1] = contrasts
            Bg[2] = Bg[2] +Dx[5]
            IF(Bg[2] GT 1.5d0 OR Bg[2] LT 0.d0) THEN Bg[2] = contrasts

        ENDIF
        
        history[jj,*]     = [Phig[0],Phig[1],Phig[2],Bg[0],Bg[1],Bg[2]]
        lateconv          = converg
        jj = jj+1
        
    ENDWHILE

   ;PRINT,'convergence',lateconv
    error[repet,0]    = (Phig[0]-Phi[0])*FSR[0]*180.d0/!dpi ; converted into degree units
    error[repet,1]    = (Phig[1]-Phi[1])*FSR[1]*180.d0/!dpi ; converted into degree units
    error[repet,2]    = (Phig[2]-Phi[2])*FSR[2]*180.d0/!dpi ; converted into degree units
    error[repet,3]    = (Bg[0]-B[0])/B[0]
    error[repet,4]    = (Bg[1]-B[1])/B[1]
    error[repet,5]    = (Bg[2]-B[2])/B[2]

    IF(repet MOD 20 EQ 0) THEN BEGIN
        !P.MULTI=[0,3,3]
        plot,error[*,0],charsize=1.5,tit=STRING(SIGMA(error[*,0]))
        plot,error[*,1],charsize=1.5,tit=STRING(SIGMA(error[*,1]))
        plot,error[*,2],charsize=1.5,tit=STRING(SIGMA(error[*,2]))
        plot,error[*,3],charsize=1.5,tit=STRING(SIGMA(error[*,3]))
        plot,error[*,4],charsize=1.5,tit=STRING(SIGMA(error[*,4]))
        plot,error[*,5],charsize=1.5,tit=STRING(SIGMA(error[*,5]))

    ENDIF

ENDFOR

SET_PLOT,'ps'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=21,ysize=26
!P.MULTI=[0,3,3]
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
