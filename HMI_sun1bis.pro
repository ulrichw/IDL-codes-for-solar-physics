; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW HOW TO CO-TUNE HMI
; UNLIKE HMI_sun1.pro WE ALSO FIT FOR THE FRINGES ON THE FRONT WINDOW+BLOCKER FILTER
; WITH nfreq FREQUENCIES
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
; NB: thresh, width0, depth0 are guess values that need
; to be adjusted otherwise the fit might be bad !!!
;
; ver 1.4 December 5, 2006
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
;
;----------------------------------------------------------------------

PRO HMI_sun1bis,draw,freq
; freq must be an array of nfreq thicknesses in meters


lamref      = 6173.3433 ;reference central wavelength for the FeI 6173 line IN AIR
time0       = SYSTIME(1); variable used to determine the elapsed time
dpi         = 2.d0*!dpi

; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

nseq        = 27                  ; number of filtergrams in the detune sequence
nx          = 256                 ; number of rows (we rebin the filtergrams)
ny          = 256                 ; number of columns
anglim      = 990.;980.d0         ; size of solar disk in arcseconds
factor      = 10000.              ; we divide the intensity by this value to make the fit easier

; MEASURED INTENSITIES
;------------------------------------------------


;list ='listSun060207_230329' ; in Calmode
;list ='listSun060208_180841' ; in Calmode
;list ='listSun060222_204505' ; in Obsmode
;list ='listSun060224_225806' ; in Calmode
;list ='listSun060224_000059' ; in Calmode
;list ='listSun060304_232326' ; in Calmode
;list ='listSun060616_205620' ; in Calmode
;list ='listSun060616_211450' ; in Calmode
;list ='listSun060616_224250' ; in Calmode
;list='listSun061031_205938'  ; in Calmode
;list='listSun070114_42753'   ; in Calmode
 list='listSun070114_42723'   ; in Calmode
 IF list EQ 'listSun070114_42753' THEN lam0 = 370.58260 ; in m/s for 01/14/2007 at 23:46:30 UT
 IF list EQ 'listSun070114_42723' THEN lam0 = 363.74918 ; in m/s for 01/14/2007 at 23:39:30 UT

; THE FOLLOWING LINE LOCATES THE CENTER OF THE SOLAR DISK
; ON THE FILTERGRAMS. IT HAS TO BE ADJUSTED:
;center=[70,58]
center   = [nx/2,ny/2]
;center   = [126,129]
distance = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; READ DETUNE SEQUENCE AND GET RID OF OVERSCANS AND BAD EXPOSURES
;---------------------------------------------------------------------------------

RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
dlamdv      = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam0        = dlamdv*lam0         ; positive velocities increase the wavelength

Inten = imx[*,*,1:27]/factor

;nx=64
;ny=64
;nx=1
;ny=1
;Inten = REBIN(Inten,nx,ny,27)
;distance=rebin(distance,nx,ny)

; FREQUENCIES TO FIT FOR ON THE FRONT WINDOW + BLOCKER
;-------------------------------------------------------

;freq        = ([.003d0,.005d0,.006d0,.009d0]*1.d10*2.d0*1.518569d0+lamref+lam0)/(lamref+lam0)^2.d0; dominant frequency in A^{-1} of the front window + blocker fringes
;freq        = ([.004d0,.005d0,.006d0,.009d0]*1.d10*2.d0*1.518569d0+lamref+lam0)/(lamref+lam0)^2.d0; dominant frequency in A^{-1} of the front  WARNING, ONLY FOR TESTING PURPOSE !!!!!
freq        = (freq*1.d10*2.d0*1.518569d0+lamref+lam0)/(lamref+lam0)^2.d0 ; angular frequency
nfreq       = N_ELEMENTS(freq)


; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

Bg          = DBLARR(nx,ny,3)     ; contrasts of the tunable elements
Phig        = DBLARR(nx,ny,3)     ; relative phases
continuum   = DBLARR(nx,ny)       ; relative sun intensity
linewidth   = DBLARR(nx,nx)       ; solar linewidth
linedepth   = DBLARR(nx,nx)       ; solar linedepth
proxy1      = DBLARR(nx,nx,nfreq) ; Cosine coefficient
proxy2      = DBLARR(nx,nx,nfreq) ; Sine coefficient
mu0         = 0.d0                ; regularization parameter for the least-squares fit. OPTIONAL: DEPENDS ON THE NOISE LEVEL
nparam      = 6+nfreq*2           ; number of parameters we fit for
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelengths in Angstrom
Residual    = DBLARR(nseq)        ; Residual of the least-squares fit
dIntendI0   = DBLARR(nseq)        ; derivatives to compute the Jacobian matrix
dIntendIc   = DBLARR(nseq)
dIntendw    = DBLARR(nseq)
dIntendPhi0 = DBLARR(nseq)
dIntendPhi1 = DBLARR(nseq)
dIntendPhi2 = DBLARR(nseq)
dproxy1     = DBLARR(nseq,nfreq)
dproxy2     = DBLARR(nseq,nfreq)
dlinedIc    = 1.d0
tuning      = DBLARR(3,nseq)      ; tuning sequence
maxsteps    = 105                  ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(nx,ny)       ; to discard spurious results

; GUESS VALUES
contrasts   = 0.98d0
depth0      = 0.53
width0      = 0.071
thresh      = 41000./factor
prox0       = FLTARR(nfreq)
prox1       = prox0
prox0[0]    = 0.0;-0.0087;-0.014127391
prox0[1]    = 0.0;-0.00147; 0.031075751
prox0[2]    = 0.0;-0.00277;-0.026653074
IF nfreq GE 4 THEN prox0[3]    = 0.0;-0.0044137; 0.0035977823
prox1[0]    = 0.0;-0.00786; 0.018048106
prox1[1]    = 0.0;-0.006036;-0.0046438219
prox1[2]    = 0.0;-0.004409; 0.014074757
IF nfreq GE 4 THEN prox1[3]    = 0.0;-0.015568; 0.0026015436


; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)      ; FSR in Angstrom
FSR[0]      = 0.172457d0     ; for the narrow-band Michelson
FSR[1]      = 0.344242d0     ; for the broad-band  Michelson
FSR[2]      = 0.690d0        ; for E1
FSR[3]      = 1.405d0        ; for E2
FSR[4]      = 2.779d0        ; for E3
FSR[5]      = 5.682d0        ; for E4
FSR[6]      = 11.354d0       ; for E5

FSR         = dpi/FSR

; DEFINITION OF THE DETUNE SEQUENCE
; the detune sequence is given in steps for the HCMs
; it is then converted in actual angles
; the NB tuning waveplate rotates in the opposite
; direction to the WB and E1 waveplates
;                NB MICHELSON   WB MICHELSON         E1 LYOT
;-----------------------------------------------------------

tuning[*,0]  = [         0.d0,          0.d0 ,         0.d0]
tuning[*,1]  = [        80.d0,          0.d0 ,         0.d0]
tuning[*,2]  = [       160.d0,          0.d0 ,         0.d0]
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


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN
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

;blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam);,/LSQUADRATIC) ; centered blocker
;blocker0    = INTERPOL(q[*,1]/100.d0,q[*,0]+2.08916d0-lamref,lam) ; uncentered blocker (correction of blocker11.fits because of f8 instead of f18)
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0-lamref,lam) ; uncentered blocker
;RESTORE,'Lyot.bin' ; transmission profile of E2-E5 obtained in sunlight on Oct 16, 2006 (see LyotShine_Sun.pro) from lyotsunsetSprofiles2006oct16.fits
;blocker0    = blocker0*INTERPOL(F/MAX(F),l,lam) ; contains the blocker+Lyot+front window profile

IF draw EQ 1 THEN GOTO,draw ; if draw=1 we do not want to fit for the phases, we just want to plot an existing result


; BEGINNING OF THE LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=950

Phig0         = [-30.7,-90.7,-142.4]/FSR[*]/180.d0*!dpi ; GUESS VALUES
proxy1[*,*,0] = prox0[0]
proxy1[*,*,1] = prox0[1]
proxy1[*,*,2] = prox0[2]
IF nfreq GE 4 THEN proxy1[*,*,3] = prox0[3]
proxy2[*,*,0] = prox1[0]
proxy2[*,*,1] = prox1[1]
proxy2[*,*,2] = prox1[2]
IF nfreq GE 4 THEN proxy2[*,*,3] = prox1[3]
Bg[*,*,0]     = contrasts
Bg[*,*,1]     = contrasts
Bg[*,*,2]     = contrasts
Phig[*,*,0]   = Phig0[0]
Phig[*,*,1]   = Phig0[1]
Phig[*,*,2]   = Phig0[2]


FOR jjj=0,ny-1 DO BEGIN

    !P.MULTI=[0,3,3]
    TVIM,Phig[*,*,0],/scale,pcharsize=1.5
    TVIM,Phig[*,*,1],/scale,pcharsize=1.5
    TVIM,Phig[*,*,2],/scale,pcharsize=1.5
    FOR i=0,nfreq-1 DO TVIM,proxy1[*,*,i],/scale,pcharsize=1.5
    TVIM,linewidth,/scale,pcharsize=1.5,range=[0.059,0.092]
    TVIM,continuum,/scale,pcharsize=1.5,range=[0.7*thresh,1.3*thresh]


    FOR iii=0,nx-1 DO BEGIN
   
        IF(distance[iii,jjj] LE anglim) THEN BEGIN


           ; GUESS VALUES
           ;------------------------------------------------

            Icg     = thresh        ; estimate of the solar continuum
            fdepthg = depth0*thresh ; depth of the solar line according to Stenflo & Lindegren (1977)
            fwidthg = width0        ; width of the solar line according to Stenflo & Lindegren (1977)            
            
            converg = 0
            jj      = 0
           ;mu      = mu0
            
            profilef= blocker0 ; we use the averaged blocker+front window profile
            FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0
            
            WHILE (converg EQ 0) DO BEGIN
                
                templine= EXP(-(lam-lam0)^2.d0/fwidthg^2.d0)
                line    = Icg -fdepthg*templine
                dlinedI0=     -        templine
                dlinedw =     -fdepthg*templine*(lam-lam0)^2.d0*2.d0/fwidthg^3.d0

                ripples = DBLARR(nlam)+1.d0
                FOR i=0,nfreq-1 DO ripples = ripples + proxy1[iii,jjj,i]*COS(dpi*lam*freq[i]) + proxy2[iii,jjj,i]*SIN(dpi*lam*freq[i])
                
                FOR j=0,nseq-1 DO BEGIN
                  
                    profileg      = 0.125d0 * profilef * (1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j])) * (1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j])) * (1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))
                  
                    Residual [j]  = Inten[iii,jjj,j] - TOTAL(line*profileg*ripples)*dlam
                    
                    dIntendI0[j]  = TOTAL(dlinedI0*profileg*ripples)*dlam
                    dIntendw [j]  = TOTAL(dlinedw *profileg*ripples)*dlam
                    dIntendIc[j]  = TOTAL(dlinedIc*profileg*ripples)*dlam
                    dIntendPhi0[j]= TOTAL(line*(-Bg[iii,jjj,0]*FSR[0]*SIN(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef*ripples)*dlam
                    dIntendPhi1[j]= TOTAL(line*(-Bg[iii,jjj,1]*FSR[1]*SIN(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef*ripples)*dlam
                    dIntendPhi2[j]= TOTAL(line*(-Bg[iii,jjj,2]*FSR[2]*SIN(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))/8.d0*profilef*ripples)*dlam
                    FOR i=0,nfreq-1 DO dproxy1[j,i]    = TOTAL(line*profileg *  COS(dpi*lam*freq[i]) )*dlam
                    FOR i=0,nfreq-1 DO dproxy2[j,i]    = TOTAL(line*profileg *  SIN(dpi*lam*freq[i]) )*dlam
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
                    FOR j=0,nfreq-1 DO BEGIN
                        Jac[6+j*2,i]  = dproxy1[i,j]
                        Jac[7+j*2,i]  = dproxy2[i,j]
                    ENDFOR
                ENDFOR
                
                LA_SVD,Jac,W,U,V,/DOUBLE               
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
                
                temp   = TRANSPOSE(Dx)/[fdepthg,Icg,fwidthg,REFORM(Phig[iii,jjj,0]),REFORM(Phig[iii,jjj,1]),REFORM(Phig[iii,jjj,2]),REFORM(proxy1[iii,jjj,0]),REFORM(proxy2[iii,jjj,0]),REFORM(proxy1[iii,jjj,1]),REFORM(proxy2[iii,jjj,1]),REFORM(proxy1[iii,jjj,2]),REFORM(proxy2[iii,jjj,2]),REFORM(proxy1[iii,jjj,3]),REFORM(proxy2[iii,jjj,3])];,REFORM(proxy1[iii,jjj,4]),REFORM(proxy2[iii,jjj,4])]
                err[jj]= MAX(ABS(temp))
                
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7 )    THEN converg = 2 ; no convergence
                IF(err[jj] LE 1.d-7)                          THEN converg = 1 ; convergence
                
                IF(converg EQ 2) THEN BEGIN
                    j = WHERE(err EQ MIN(err))
                    fdepthg           = history[j[0],0]
                    Icg               = history[j[0],1]
                    fwidthg           = history[j[0],2]
                    Phig[iii,jjj,0]   = history[j[0],3]
                    Phig[iii,jjj,1]   = history[j[0],4]
                    Phig[iii,jjj,2]   = history[j[0],5]
                    FOR i=0,nfreq-1 DO BEGIN
                        proxy1[iii,jjj,i] = history[j[0],6+2*i]
                        proxy2[iii,jjj,i] = history[j[0],7+2*i]
                    ENDFOR
                ENDIF
                
                IF(converg NE 2) THEN BEGIN
                    fdepthg = fdepthg+Dx[0]
                    IF(fdepthg LT 0.2d0*thresh*depth0 OR fdepthg GT 3.0d0*depth0*thresh) THEN fdepthg = depth0*thresh
                    Icg     = Icg    +Dx[1]
                    IF(Icg LT 0.3d0*thresh OR Icg GT 3.0d0*thresh) THEN Icg = thresh
                    fwidthg = fwidthg+Dx[2]
                    IF(fwidthg GT 3.5d0*width0 OR fwidthg LT 0.3d0*width0) THEN fwidthg = width0
                    Phig[iii,jjj,0] = Phig[iii,jjj,0]+Dx[3]
                    IF(ABS(Phig[iii,jjj,0]) GT 1.2d0*!dpi/FSR[0]) THEN Phig[iii,jjj,0] = Phig0[0] 
                    Phig[iii,jjj,1] = Phig[iii,jjj,1]+Dx[4]
                    IF(ABS(Phig[iii,jjj,1]) GT 1.2d0*!dpi/FSR[1]) THEN Phig[iii,jjj,1] = Phig0[1]
                    Phig[iii,jjj,2] = Phig[iii,jjj,2]+Dx[5]
                    IF(ABS(Phig[iii,jjj,2]) GT 1.2d0*!dpi/FSR[2]) THEN Phig[iii,jjj,2] = Phig0[2]
                    FOR i=0,nfreq-1 DO BEGIN
                        proxy1[iii,jjj,i] = proxy1[iii,jjj,i]+Dx[6+2*i]
                        IF(ABS(proxy1[iii,jjj,i]) GT 2.d0) THEN proxy1[iii,jjj,i] = prox0[i]
                        proxy2[iii,jjj,i] = proxy2[iii,jjj,i]+Dx[7+2*i]
                        IF(ABS(proxy2[iii,jjj,i]) GT 2.d0) THEN proxy2[iii,jjj,i] = prox1[i]
                    ENDFOR
                ENDIF

                history[jj,*]     = [fdepthg,Icg,fwidthg,Phig[iii,jjj,0],Phig[iii,jjj,1],Phig[iii,jjj,2],proxy1[iii,jjj,0],proxy2[iii,jjj,0],proxy1[iii,jjj,1],proxy2[iii,jjj,1],proxy1[iii,jjj,2],proxy2[iii,jjj,2],proxy1[iii,jjj,3],proxy2[iii,jjj,3]];,proxy1[iii,jjj,4],proxy2[iii,jjj,4]]
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
            PRINT,'cosine=',proxy1[iii,jjj,0],proxy1[iii,jjj,1],proxy1[iii,jjj,2]
            PRINT,'sine=',proxy2[iii,jjj,0],proxy2[iii,jjj,1],proxy2[iii,jjj,2]
            PRINT,iii,jjj

        ENDIF
        
    ENDFOR
ENDFOR

SAVE,Phig,linewidth,linedepth,continuum,lateconv,lam0,proxy1,proxy2,freq,maxsteps,contrasts,blocker0,lam,FILE='RESULTS/RESULTSbis_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:

; DRAWING OF THE RESULTS
;-----------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTSbis_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a = WHERE(lateconv EQ 1,COMPLEMENT=aa)
IF(a[0] EQ -1) THEN BEGIN
    a = WHERE(lateconv EQ 2,COMPLEMENT=aa)  ; !!!!WARNING!!!!
    PRINT,'NO EARLY CONVERGENCES: PROBLEM'
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

; TO SUPPRESS THE STRIPES AT THE LEFT AND UPPER EDGE
;distance[0:7,*]=0.0
;distance[*,120:nx-1]=0.0
IF nx EQ 256 THEN BEGIN
    distance[*,0:2]=0.0
    distance[251:nx-1,*]=0.0
ENDIF

FOR i=0,2 DO Phig[*,*,i]=Phig[*,*,i]*distance

PRINT,'CORRECTED AVERAGED PHASES'
a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 1,COMPLEMENT=aa)
IF(a[0] EQ -1) THEN a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 2,COMPLEMENT=aa) ; !!!!WARNING!!!!
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

LOADCT,4
SET_PLOT,'ps'
IF(nx GT 1) THEN BEGIN
;GOTO,recons
!p.multi=[0,2,3]
DEVICE,file='yo.ps',/color,bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27
titl=STRARR(3)
titl[0]='Narrow-Band Michelson'
titl[1]='Broad-Band Michelson'
titl[2]='Lyot Element E1'
FOR i=0,2 DO TVIM,Phig[*,*,i],/scale,range=[moy[i]-2.5*dis[i],moy[i]+2.5*dis[i]],tit='!17'+titl[i]+STRING(moy[i]),stit='!17Phase (in degrees)'
FOR i=0,nfreq-1 DO TVIM,proxy1[*,*,i]*distance,/scale,stit='!17Cosine coefficient',range=[-0.03+MEAN(proxy1[nx/2-10:nx/2+10,nx/2-10:nx/2+10,i]),0.03+MEAN(proxy1[nx/2-10:nx/2+10,nx/2-10:nx/2+10,i])],tit=STRING(freq[i])
FOR i=0,nfreq-1 DO TVIM,proxy2[*,*,i]*distance,/scale,stit='!17Sine coefficient',range=[-0.03+MEAN(proxy2[nx/2-10:nx/2+10,nx/2-10:nx/2+10,i]),0.03+MEAN(proxy2[nx/2-10:nx/2+10,nx/2-10:nx/2+10,i])]
TVIM,continuum*distance,/scale,stit='!17Continuum intensity'
TVIM,linewidth*distance,/scale,stit='!17Linewidth (A)',barwidth=0.5,pcharsize=0.75,xtit='Pixel number',ytit='Pixel number',range=[0.062,0.08];,range=[MIN(linewidth[a]),MAX(linewidth[a])]
TVIM,(linedepth/continuum)*distance,/scale,stit='!17Linedepth',barwidth=0.5,pcharsize=0.75,range=[MIN(linedepth[a]/continuum[a]),MAX(linedepth[a]/continuum[a])],xtit='Pixel number',ytit='Pixel number'
DEVICE,/close
;recons:
ENDIF

; RECONSTRUCTION OF THE DETUNE SEQUENCE
;--------------------------------------------------------------------------------------


FOR i=0,2 DO Phig[*,*,i] = Phig[*,*,i]/FSR[i]/180.d0*!dpi ; convert from degrees to radians/FSR

FOR i=0,nseq-1 DO Inten[*,*,i] = Inten[*,*,i]*distance
Inten0 = Inten

profilef     = blocker0  ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0 ;!!!!!CONTRASTS!!!!

Inten2  = FLTARR(nx,nx,nseq) 
FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,nx-1 DO BEGIN
        IF(distance[i,j] NE 0.0) THEN BEGIN
            templine     = EXP(-(lam-lam0)^2.d0/linewidth[i,j]^2.d0)
            line         = continuum[i,j] -linedepth[i,j]*templine
            ripples = DBLARR(nlam)+1.d0
            FOR k=0,nfreq-1 DO ripples = ripples + proxy1[i,j,k]*COS(dpi*lam*freq[k]) + proxy2[i,j,k]*SIN(dpi*lam*freq[k])
            
            FOR k=0,nseq-1 DO BEGIN
                profileg      = 0.125d0 * profilef * (1.d0+Bg[i,j,0]*COS(FSR[0]*(lam+Phig[i,j,0])+tuning[0,k])) * (1.d0+Bg[i,j,1]*COS(FSR[1]*(lam+Phig[i,j,1])+tuning[1,k])) * (1.d0+Bg[i,j,2]*COS(FSR[2]*(lam+Phig[i,j,2])+tuning[2,k]))
                Inten2[i,j,k] = TOTAL(line*profileg*ripples)*dlam
            ENDFOR

        ENDIF
    ENDFOR
ENDFOR

; WE PLOT THE MAPS OF THE RELATIVE RECONSTRUCTION ERROR
!P.MULTI=[0,2,3]
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
    aa   = WHERE(temp LT 0.05 AND temp GT -0.05,na,COMPLEMENT=b)
    hist = histogram(temp[aa],binsize=0.001,min=-0.05,max=0.05)
    plot,FINDGEN(N_ELEMENTS(hist))*0.001-0.05,hist/FLOAT(na),psym=10,tit='!7r!17='+STRING(SIGMA(temp[aa]))+'/!7l!17='+STRING(MEAN(temp[aa])),ytit='Histogram in percentage',xtit='Relative error',charsize=1.5
    PRINT,MIN(temp[aa]),MAX(temp[aa]),SIGMA(temp[aa]),MEAN(temp[aa])
ENDFOR
DEVICE,/CLOSE

; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
!P.MULTI=0
residual=TOTAL( (REFORM(REBIN(Inten0,1,1,nseq))-REFORM(REBIN(Inten2,1,1,nseq)))^2.0 )
device,file='yo3.ps',xsize=15,ysize=12,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,REBIN(Inten0,1,1,nseq),xst=1,tit='!17'+STRING(residual),xtit='position number',ytit='intensity',charsize=1.5
OPLOT,REBIN(Inten2,1,1,nseq),linestyle=2,color=180
device,/close
PRINT,'RESIDUAL=',residual


SET_PLOT,'X'
READ,pause

END
