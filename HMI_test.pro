; TEST ANOTHER MODEL FOR THE SOLAR LINE

; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE SUN AS A SOURCE, AND FOR A DE-TUNE SEQUENCE
; OBJECTIVE: TO KNOW HOW TO CO-TUNE HMI
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

PRO HMI_test,draw



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


; VARIOUS PARAMETER AND VARIABLE DEFINITIONS
;--------------------------------------------

nseq0       = 27
nx          = 128                 ; number of rows (we rebin the filtergrams)
ny          = 128                 ; number of columns

; MEASURED INTENSITIES
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
 list ='listSun060621_170559' ; in Calmode
;list ='listSun060621_171740' ; in Obsmode

;center  = [70,58] ; CHANGE MANUALLY !!!!!!!!!!!!!!!!!!!!!!
 center   = [nx/2,ny/2]
distance = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; READ DETUNE SEQUENCE AND GET RID OF OVERSCANS AND BAD EXPOSURES
;---------------------------------------------------------------------------------


;lam0        = 0.0
;imx         = READFILE(nx,ny,lam0,list) ; SUBROUTINE WRITTEN BY S. COUVIDAT 
RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
medc        = MEDIAN(overscan)
wgood1      = WHERE(overscan LE medc+20 ,COMPLEMENT=wbad)
PRINT,'OVERSCAN',wbad
medc2       = MEDIAN(exposure)
wgood2      = WHERE(exposure LE medc2+100,COMPLEMENT=wbad)
PRINT,'BAD EXPOSURES',wbad
dlamdv      = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam0        = dlamdv*lam0         ; positive velocities increase the wavelength

wgood       = WHERE(exposure LE medc2+100 AND overscan LE medc+20)
IF(wgood[0] EQ 0) THEN wgood=wgood[1:N_ELEMENTS(wgood)-1] ELSE PRINT,'WARNING: PROBLEM WITH THE FIRST DARK' ; 1st image is a dark frame
IF(wgood[N_ELEMENTS(wgood)-1] EQ nseq0+1) THEN wgood=wgood[0:N_ELEMENTS(wgood)-2] ELSE PRINT,'WARNING: PROBLEM WITH THE SECOND DARK'; 2nde image is a dark frame
wgood       = wgood-1

;wgood=INDGEN(18)

nseq        = N_ELEMENTS(wgood)
Inten       = imx[*,*,wgood]

; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------


Bg          = DBLARR(nx,ny,3)     ; contrasts of the tunable elements
Phig        = DBLARR(nx,ny,3)     ; relative phases
continuum   = DBLARR(nx,ny)       ; relative sun intensity
linewidth   = DBLARR(nx,nx)       ; solar linewidth
linedepth   = DBLARR(nx,nx)       ; solar linedepth
mu0         = 0.d0 ;0.02          ; regularization parameter for the least-squares fit. OPTIONAL: DEPENDS ON THE NOISE LEVEL
nparam      = 4                   ; number of parameters we fit for
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
tuning      = DBLARR(3,nseq0)     ; tuning sequence
maxsteps    = 105                  ; maximum steps allowed for the iterative Least-Squares algorithm
history     = DBLARR(maxsteps,nparam) ; we save the temporary parameters here
    err     = DBLARR(maxsteps)    ; error estimate at each step
lateconv    = INTARR(nx,ny)       ; to discard spurious results


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
contrasts   = 0.98d0
anglim      = 980;1024.0;980.d0              ; in arcseconds

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

tuning       = tuning[*,wgood]


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

blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam);,/LSQUADRATIC)


; MODEL OF THE IRON LINE
;q            = READFITS('iron_line_model.fits',h)
;q0           = [lam[0],REFORM(q[0,*]),lam[nlam-1]]
;q1           = [1.0,REFORM(q[1,*]),1.0]
;line0        = INTERPOL(q1,q0-lamref,lam)
;dlinedIc     = line0

; ATLAS PROFILE: OBSERVED SOLAR LINE
RESTORE,'atlas_profile.sav'
lb1=lb1-6175.025  ; vacuum solar line
lb1=[lam[0],lb1,lam[nlam-1]]
pr1=[1.0,pr1,1.0]
line0=INTERPOL(pr1,lb1,lam)
dlinedIc     = line0

IF draw EQ 1 THEN GOTO,draw

; LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=900
!P.MULTI=[0,2,2]


; GUESS
;------------------------

FOR jjj=0,ny-1 DO BEGIN

    TVIM,Phig[*,*,0],/scale
    TVIM,Phig[*,*,1],/scale
    TVIM,Phig[*,*,2],/scale
    TVIM,continuum  ,/scale

    FOR iii=0,nx-1 DO BEGIN
   
        IF(distance[iii,jjj] LE anglim) THEN BEGIN


           ; GUESS VALUES
           ;------------------------------------------------

            Icg             = thresh         ; estimate of the solar continuum
            Bg[iii,jjj,*]   = contrasts
            Phig[iii,jjj,*] = 0.0; [66.d0/180.*!pi/FSR[0],-80.d0/180.*!pi/FSR[1],66.d0/180.*!pi/FSR[2]]
            
            converg = 0
            jj      = 0
           ;mu      = mu0
            
            profilef= blocker0 ; we use the averaged blocker+front window profile
            FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0
            
            WHILE (converg EQ 0) DO BEGIN
               
                line    = Icg * line0
                
                FOR j=0,nseq-1 DO BEGIN
                    
                    profileg      = 0.125d0 * profilef * (1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j])) * (1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j])) * (1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))
                  
                    Residual[j]   = Inten[iii,jjj,j] - TOTAL(line*profileg)*dlam

                    dIntendIc[j]  = TOTAL(dlinedIc*profileg)*dlam
                    dIntendPhi0[j]= TOTAL(line*(-Bg[iii,jjj,0]*FSR[0]*SIN(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef)*dlam
                    dIntendPhi1[j]= TOTAL(line*(-Bg[iii,jjj,1]*FSR[1]*SIN(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,2]*COS(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))/8.d0*profilef)*dlam
                    dIntendPhi2[j]= TOTAL(line*(-Bg[iii,jjj,2]*FSR[2]*SIN(FSR[2]*(lam+Phig[iii,jjj,2])+tuning[2,j]))*(1.d0+Bg[iii,jjj,0]*COS(FSR[0]*(lam+Phig[iii,jjj,0])+tuning[0,j]))*(1.d0+Bg[iii,jjj,1]*COS(FSR[1]*(lam+Phig[iii,jjj,1])+tuning[1,j]))/8.d0*profilef)*dlam
                ENDFOR
                

              ; Jacobian matrix
                Jac  = DBLARR(nparam,nseq)
                FOR i= 0,nseq-1 DO BEGIN
                    
                    Jac[0,i]  = dIntendIc[i]
                    Jac[1,i]  = dIntendPhi0[i]
                    Jac[2,i]  = dIntendPhi1[i]
                    Jac[3,i]  = dIntendPhi2[i]
                    
                ENDFOR

                LA_SVD,Jac,W,U,V,/DOUBLE
                
              ; regularization:
              ; filter = W*W/(W*W+mu^2.d0) ; mu is the regularization parameter
              ; Dx     = V##DIAG_MATRIX(1.d0/W*filter)##TRANSPOSE(U)##Residual
                Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
                
                temp   = TRANSPOSE(Dx)/[Icg,REFORM(Phig[iii,jjj,0]),REFORM(Phig[iii,jjj,1]),REFORM(Phig[iii,jjj,2])]
               ;err[jj]= TOTAL(temp^2.d0)/DOUBLE(nparam)
                err[jj]= MAX(ABS(temp))
                
              ; IF(FINITE(err[jj]) EQ 0 AND jj GT 0)          THEN mu = 2.5d0*mu ELSE mu = mu0
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7 )    THEN converg = 2  ; 1.d-14 
                IF(err[jj] LE 1.d-7)                          THEN converg = 1
                
                IF(converg EQ 2) THEN BEGIN
                    j               = WHERE(err EQ MIN(err))
                    Icg             = history[j[0],0]
                    Phig[iii,jjj,0] = history[j[0],1]
                    Phig[iii,jjj,1] = history[j[0],2]
                    Phig[iii,jjj,2] = history[j[0],3]
                ENDIF
                
                IF(converg NE 2) THEN BEGIN

                    Icg             = Icg            + Dx[0]
                    IF(Icg LT 0.2d0*thresh OR Icg GT 2.0d0*thresh) THEN Icg = thresh
                    Phig[iii,jjj,0] = Phig[iii,jjj,0]+ Dx[1]
                    IF(ABS(Phig[iii,jjj,0]) GT 1.25d0*!dpi/FSR[0]) THEN Phig[iii,jjj,0] = 0.;66.d0/180.*!pi/FSR[0] 
                    Phig[iii,jjj,1] = Phig[iii,jjj,1]+ Dx[2]
                    IF(ABS(Phig[iii,jjj,1]) GT 1.25d0*!dpi/FSR[1]) THEN Phig[iii,jjj,1] = 0.;-80.d0/180.*!pi/FSR[1]  
                    Phig[iii,jjj,2] = Phig[iii,jjj,2]+ Dx[3]
                    IF(ABS(Phig[iii,jjj,2]) GT 1.25d0*!dpi/FSR[2]) THEN Phig[iii,jjj,2] = 0.;66.d0/180.*!pi/FSR[2]
                    
                ENDIF
                
                history[jj,*]     = [Icg,Phig[iii,jjj,0],Phig[iii,jjj,1],Phig[iii,jjj,2]]
                lateconv[iii,jjj] = converg
                
               ;PRINT,jj,err[jj]
                jj = jj+1
                
            ENDWHILE
            Phig[iii,jjj,*]    = Phig[iii,jjj,*]*FSR[*]*180.d0/!dpi ; converted into degree units
            continuum[iii,jjj] = Icg
            PRINT,'PARAMETERS='
            FOR i=0,2 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'Ic=',Icg
            PRINT,'lateconv=',converg
            PRINT,iii,jjj

        ENDIF
        
    ENDFOR
ENDFOR

SAVE,Phig,continuum,lateconv,lam0,FILE='RESULTS/RESULTS2_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:

; COMPUTE THE AVERAGE PHASES
;--------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS2_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

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
distance[b]=0.d0

; TO SUPRESS THE STRIPES AT THE LEFT AND UPPER EDGE
distance[0:7,*]=0.0
distance[*,120:nx-1]=0.0

FOR i=0,2 DO Phig[*,*,i]=Phig[*,*,i]*distance
TVIM,distance

PRINT,'CORRECTED AVERAGED PHASES'
a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 1,COMPLEMENT=aa)
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

!P.MULTI=[0,2,2]
LOADCT,4
SET_PLOT,'ps'
DEVICE,file='yo.ps',/color,bits=24,xoffset=0,yoffset=0,xsize=21,ysize=20
titl=STRARR(3)
titl[0]='Narrow-Band Michelson'
titl[1]='Broad-Band Michelson'
titl[2]='Lyot Element E1'
FOR i=0,2 DO TVIM,Phig[*,*,i],/scale,range=[moy[i]-2.5*dis[i],moy[i]+2.5*dis[i]],tit='!17'+titl[i],stit='!17Phase (in degrees)',barwidth=0.5,pcharsize=0.75
TVIM,continuum*distance,/scale,stit='!17Continuum intensity',barwidth=0.5,pcharsize=0.75
DEVICE,/close
!P.MULTI=0


; RECONSTRUCTION OF THE DETUNE SEQUENCE
;--------------------------------------------------------------------------------------


FOR i=0,2 DO Phig[*,*,i] = Phig[*,*,i]/FSR[i]/180.d0*!dpi ; convert from degrees to radians/FSR
Bg[*,*,*]    = contrasts


; TEST OF THE CODE: WE PRODUCE AN ARTIFICIAL DETUNE SEQUENCE


FOR i=0,nseq-1 DO Inten[*,*,i] = Inten[*,*,i]*distance
a            = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0

profilef     = blocker0            ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0

Inten2  = FLTARR(nx,nx,nseq) 
FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,nx-1 DO BEGIN
        IF(distance[i,j] NE 0.0) THEN BEGIN
            line         = continuum[i,j]*line0
                        
            FOR k=0,nseq-1 DO BEGIN
                profileg      = 0.125d0 * profilef * (1.d0+Bg[i,j,0]*COS(FSR[0]*(lam+Phig[i,j,0])+tuning[0,k])) * (1.d0+Bg[i,j,1]*COS(FSR[1]*(lam+Phig[i,j,1])+tuning[1,k])) * (1.d0+Bg[i,j,2]*COS(FSR[2]*(lam+Phig[i,j,2])+tuning[2,k]))
                Inten2[i,j,k] = TOTAL(line*profileg)*dlam
            ENDFOR

        ENDIF
    ENDFOR
ENDFOR


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


SET_PLOT,'x'
!P.MULTI=0


READ,pause

END
