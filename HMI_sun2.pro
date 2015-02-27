; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO OBTAIN A FIRST DOPPLERGRAM (IN OBSMODE WHILE THE PHASES OF THE
; TUNABLE ELEMENTS WERE OBTAINED IN CALMODE, HENCE SLIGHT INCOHERENCE)
;
; ASSUMPTIONS: we assume contrasts and phase of the non-tunable part 
; to be equal to respectively contrasts and 0 
;
; THE SOLAR LINE IS MODELED BY A GAUSSIAN PROFILE:
; I(l) = Ic - Ic*d EXP(-l^2/vl^2)
; WHERE Ic IS THE CONTINUUM, d IS THE LINEDEPTH
; AND vl IS THE LINEWIDTH
; WE FIT FOR 3 PARAMETERS FOR THE LINE PROFILE:
; Ic, depth=Ic*d, and vl
;
; THIS PROGRAM MENTIONS THE PRESENCE OF BAD EXPOSURES AND OVERSCANS
;
; ver 1.4 February 20, 2007
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
;----------------------------------------------------------------------

PRO HMI_sun2,draw

lamref      = 6173.3433d0;reference central wavelength for the FeI 6173 line
                      ; IN AIR
                      ; used to center the averaged blocker
                      ; and front window transmission profiles
                      ; value from Stenflo & Lindegren (1977)
                      ; NB: IF THIS VALUE IS NOT THE ACTUAL ONE, THEN
                      ; WE HAVE AN AVERAGE VELOCITY ON THE SOLAR DISK
                      ; DIFFERENT FROM 0. TO BE CORRECTED
time0       = SYSTIME(1)
dpi         = 2.d0*!dpi

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------
nseq        = 27                 ; number of positions in the sequence
nx          = 256                ; number of rows
ny          = 256                ; number of columns
anglim      = 920. ; less than 950. because image not centered, unlike the phase maps in calmode...
factor      = 10000.              ; we divide the intensity by this value to make the fit easier
;center     = [nx/2,ny/2]
center      = [130,124]
distance    = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; MEASURED INTENSITIES
;--------------------------------------------

;list    = 'listSun070217_66986'    ; in Obsmode
;list    = 'listSun070217_66705'    ; CRAP !
;list     = 'listSun070217_66817'    ; in Obsmode
;list='listSun070903_539881'  ; in Obsmode  ; vignetting at bottom right
list='listSun070903_543150'
IF list EQ 'listSun070217_66986' THEN lam0 = 576.58489 ; in m/s for 02/17/2007 at 22:52:00 UT
IF list EQ 'listSun070217_66705' THEN lam0 = 361.55206 ; in m/s for 02/17/2007 at 20:26:00 UT
IF list EQ 'listSun070217_66817' THEN lam0 = 484.54499 ; in m/s for 02/17/2007 at 21:45:30 UT
IF list EQ 'listSun070903_539881'THEN lam0 =-596.94765 ; in m/s for 09/03/2007 at 18:13:44 UT
IF list EQ 'listSun070903_543150'THEN lam0 =-143.77344 ; in m/s for 09/03/2007 at 23:25:34 UT
RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
dlamdv  = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam01   = dlamdv*lam0 
Inten   = imx[*,*,2:28]/factor ; for data after February 2007 obtained with readimages2.pro

; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------

Bg          = DBLARR(nx,ny,3)    ; contrasts of the tunable elements
Phig        = DBLARR(nx,ny,3)    ; relative phases
lam0g       = DBLARR(nx,ny)      ; Doppler shift
continuum   = DBLARR(nx,ny)      ; solar continuum, saved for the throughput test
linewidth   = DBLARR(nx,ny)
linedepth   = DBLARR(nx,ny)
mu0         = 0.d0   ;0.02       ; regularization parameter for the least-squares fit. usual value. OPTIONAL
nparam      = 4                  ; number of parameters we fit for
nlam        = 2250               ; number of wavelengths
dlam        = 3.6d0/1.d3         ; resolution in wavelength
lam         = dlam*(DINDGEN(nlam)-(nlam-1)/2.)
Residual    = DBLARR(nseq)       ; residuals of the least-squares fit
dIntendlam0 = DBLARR(nseq)       ; derivatives to compute the Jacobian matrix
dIntendI0   = DBLARR(nseq)
dIntendIc   = DBLARR(nseq)
dIntendw    = DBLARR(nseq)
dlinedIc    = 1.d0
tuning      = DBLARR(3,nseq)     ; tuning sequence
maxsteps    = 105                ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)   ; error estimate at each step

contrasts   = [0.92,0.96,0.98,0.96]
phase       = [-13.00,-8.4,-24.7,-45.1]*!dpi/180.d0
depth0      = 0.52               ; guess for the solar linedepth
width0      = 0.062              ; guess for the solar linewidth
thresh      = 100000./factor     ; guess of the continuum intensity
lateconv    = INTARR(nx,ny)       ; to discard spurious results


; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)    ; FSR in Angstrom
FSR[0]      = 0.172d0;0.172457d0     ; for the narrow-band Michelson
FSR[1]      = 0.344d0;0.344242d0     ; for the broad-band  Michelson
FSR[2]      = 0.693d0;0.690d0;0.7039;0.702        ; for E1
FSR[3]      = 1.407d0        ; for E2
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

IF nseq EQ 27 THEN BEGIN
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
ENDIF ELSE BEGIN   

; CO-TUNE SEQUENCE (22 POSITIONS, JUNE 2006)


    tuning[*, 0]  = [ 53,72,44]
    tuning[*, 1]  = [ 41,66,47]
    tuning[*, 2]  = [ 29,60,50]
    tuning[*, 3]  = [ 77,54,53]
    tuning[*, 4]  = [ 65,48,56]
    tuning[*, 5]  = [ 53,42,59]
    tuning[*, 6]  = [ 41,96,62]
    tuning[*, 7]  = [ 29,90,65]
    tuning[*, 8]  = [ 77,84,68]
    tuning[*, 9]  = [ 65,78,71]
    tuning[*,10]  = [ 53,72,74]
    tuning[*,11]  = [ 41,66,77]
    tuning[*,12]  = [ 29,60,80]
    tuning[*,13]  = [ 77,54,83]
    tuning[*,14]  = [ 65,48,86]
    tuning[*,15]  = [ 53,42,89]
    tuning[*,16]  = [ 41,96,92]
    tuning[*,17]  = [ 29,90,95]
    tuning[*,18]  = [ 77,84,98]
    tuning[*,19]  = [ 65,78,101]
    FOR i=0,nseq-1 DO tuning[*,i]=REFORM(tuning[*,i])*6.d0*!dpi/180.d0*[1.d0,1.d0,-1.d0]
ENDELSE


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVDED BY LOCKHEED-MARTIN
;-----------------------------------------------

; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile obtained from the file
; ANDV9601_27336_Final_1-13.csv provided by Rock Bush, e-mail 01/03/2006
; related to the front window with the serial number 27336 form Andover
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')           ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
;blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8d0,lam);,/LSQUADRATIC)
blocker0      = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.61328-lamref,lam);,/LSQUADRATIC) ; uncentered blocker (correction of blocker11.fits because of f8 instead of collimated light)

IF draw EQ 1 THEN GOTO,draw

; RELATIVE PHASES OBTAINED BY HMI_sun1.pro
;------------------------------------------------

;RESTORE,'RESULTS/RESULTS_listSun060224_225806_256.BIN'; phases in calmode
;RESTORE,'RESULTS/RESULTS_listSun060616_211450_128.BIN'
;RESTORE,'RESULTS/RESULTS_listLaser060622_221251_128.BIN'
;RESTORE,'RESULTS/RESULTS_listSun070217_66956_256.BIN' ; !WARNING THE LAM0 OF THIS FILE DELETES THE PREVIOUS LAM0, IDEM FOR BLOCKER0
;RESTORE,'RESULTS/RESULTS_listSun070217_66675_256.BIN'
RESTORE,'RESULTS/RESULTS_listSun070903_539794_256.BIN' ; PROVIDES BLOCKER0, LAM0, PHIG, AND BG
FOR j=0,ny-1 DO FOR i=0,nx-1 DO Phig[i,j,*] = Phig[i,j,*] / FSR[*] / 180.d0 * !dpi

continuum = DBLARR(nx,ny)
linewidth = DBLARR(nx,ny)
linedepth = DBLARR(nx,ny)
lateconv  = DBLARR(nx,ny)

; LEAST-SQUARES ALGORITHM
;-------------------------------------------------

; ESTIMATE OF THE FIXED TRANSMISSION PROFILE
profilef= blocker0  ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts[i-3]*COS(FSR[i]*lam+phase[i-3]))/2.d0

WINDOW,0,RETAIN=2

FOR jjj=0,ny-1 DO BEGIN

    TVIM,(lam0g-lam01)/dlamdv,/scale

    FOR iii=0,nx-1 DO BEGIN

        

        IF(distance[iii,jjj] LE anglim) THEN BEGIN

          ; GUESS PARAMETERS

            Icg      = thresh
            fdepthg  = depth0*thresh  ; depth of the solar line according to Stenflo & Lindegren (1977)
            fwidthg  = width0         ; width of the solar line according to Stenflo & Lindegren (1977)
            lam0g[iii,jjj] = lam01    ; average solar Doppler velocity due to Earth movement
                     
          ; BEGINNING OF DETUNE SEQUENCE
   
            converg  = 0
            jj       = 0
           ;mu       = mu0

            WHILE (converg EQ 0) DO BEGIN
                
                templine= EXP(-(lam-lam0g[iii,jjj])^2.d0/fwidthg^2.d0)
                line    = Icg -       fdepthg             *templine
                dlinedla=     -2.d0 * fdepthg/fwidthg^2.d0*templine*(lam-lam0g[iii,jjj])
                dlinedI0=     -                            templine
                dlinedw =     -2.d0 * fdepthg/fwidthg^3.d0*templine*(lam-lam0g[iii,jjj])^2.d0
                
                FOR j=0,nseq-1 DO BEGIN
                                        
                    profileg      = 0.125d0 * profilef * (1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j])) * (1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j])) * (1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))
                    Residual [j]  = Inten[iii,jjj,j] - TOTAL(profileg*line)*dlam
                    
                    dIntendI0[j]  = TOTAL(dlinedI0*profileg)*dlam
                    dIntendw [j]  = TOTAL(dlinedw *profileg)*dlam
                    dIntendIc[j]  = TOTAL(dlinedIc*profileg)*dlam
                    dIntendlam0[j]= TOTAL(dlinedla*profileg)*dlam
                ENDFOR
                
              ; Jacobian matrix
                
                Jac  = DBLARR(nparam,nseq)
                FOR i= 0,nseq-1 DO BEGIN
                    
                    Jac[0,i]  = dIntendI0[i]
                    Jac[1,i]  = dIntendIc[i]
                    Jac[2,i]  = dIntendw[i]
                    Jac[3,i]  = dIntendlam0[i]
                    
                ENDFOR
                
                LA_SVD,Jac,W,U,V,/DOUBLE
    
              ; regularization:
              ; filter = W*W/(W*W+mu^2.d0) ; mu is the regularization parameter
              ; Dx     = V##DIAG_MATRIX(1.d0/W*filter)##TRANSPOSE(U)##Residual
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

                temp   = TRANSPOSE(Dx)/[fdepthg,Icg,fwidthg,REFORM(lam0g[iii,jjj])]
               ;err[jj]= TOTAL(temp^2.d0)/DOUBLE(nparam)
                err[jj]= MAX(ABS(temp))

              ; IF(FINITE(err[jj]) EQ 0 AND jj GT 0)          THEN mu = 2.5d0*mu ELSE mu = mu0
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7)     THEN converg = 2 
                IF(err[jj] LE 1.d-7)                          THEN converg = 1
                
                IF(converg EQ 2) THEN BEGIN
                    j       = WHERE(err EQ MIN(err))
                    fdepthg = history[j[0],0]
                    Icg     = history[j[0],1]
                    fwidthg = history[j[0],2]
                    lam0g[iii,jjj] = history[j[0],3]
                ENDIF
                
                IF(converg NE 2) THEN BEGIN
                    fdepthg = fdepthg + Dx[0]
                    IF(fdepthg LT 0.2d0*thresh*depth0 OR fdepthg GT 3.0d0*depth0*thresh) THEN fdepthg = depth0*thresh
                    Icg     = Icg     + Dx[1]
                    IF(Icg LT 0.2d0*thresh OR Icg GT 3.d0*thresh)                        THEN Icg     = thresh
                    fwidthg = fwidthg + Dx[2]
                    IF(fwidthg GT 5.d0*width0 OR fwidthg LT 0.3d0*width0)                THEN fwidthg = width0
                    lam0g[iii,jjj]   = lam0g[iii,jjj] + Dx[3]
                    IF(ABS(lam0g[iii,jjj]-lam01) GT  .350d0)                             THEN lam0g[iii,jjj] =  lam01 ;206 mA equals +/- 10 km/s
                ENDIF
                
                history[jj,*] = [fdepthg,Icg,fwidthg,lam0g[iii,jjj]]
                lateconv[iii,jjj] = converg
               ;PRINT,jj,err[jj]
                jj = jj+1
        
            ENDWHILE
            PRINT,'PARAMETERS='
            PRINT,'Doppler Shift=',lam0g[iii,jjj]/dlamdv
            PRINT,'Icg=',Icg
            PRINT,'fdepthg=',fdepthg
            PRINT,'widthg=',fwidthg
            PRINT,'lateconv=',converg
            PRINT,iii,jjj
            continuum[iii,jjj] = Icg
            linewidth[iii,jjj] = fwidthg
            linedepth[iii,jjj] = fdepthg
            
        ENDIF

    ENDFOR
ENDFOR

lam0g = lam0g/dlamdv       ; we convert Doppler shifts into velocities
lam01 = lam01/dlamdv
SAVE,lam01,lam0g,continuum,linewidth,linedepth,lateconv,FILE='RESULTS/RESULTS_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:

RESTORE,'RESULTS/RESULTS_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a=WHERE(distance LE anglim,COMPLEMENT=b)
distance[a]=1.d0
distance[b]=0.d0
lam0g[*,*]=(lam0g[*,*]-lam01)*distance
linewidth=linewidth*distance
linedepth=linedepth*distance
continuum=continuum*distance
a=WHERE(lam0g NE 0.d0)
PRINT,'AVERAGE VELOCITY=',MEAN(lam0g[a])

; TO TRACE THE ROTATION AXIS
;a=where(ABS(lam0g-621.37)/621.37 LT 0.01,complement=b)
;lam0g[a] = -2500.

LOADCT,3
SET_PLOT,'ps'
DEVICE,file='yo.ps',/color,bits=24,xoffset=-0.5,yoffset=0,xsize=26,ysize=20
!P.MULTI=[0,2,2]
TVIM,lam0g[*,*],/scale,tit='!17',stit='!17Doppler Velocity (m/s)',barwidth=0.5;,lam0/dlamdv
TVIM,linewidth,/scale,tit='!17',stit='Linewidth (A)',barwidth=0.5,range=[MIN(linewidth[a]),0.09];MAX(linewidth[a])]
TVIM,linedepth/continuum,/scale,tit='!17',stit='Linedepth',barwidth=0.5,range=[MIN(linedepth[a]/continuum[a]),MAX(linedepth[a]/continuum[a])]
TVIM,continuum,/scale,tit='!17',stit='Intensity',barwidth=0.5,range=[MIN(continuum[a]),MAX(continuum[a])]
DEVICE,/close
SET_PLOT,'x'

READ,pause

END

