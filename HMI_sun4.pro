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

PRO MyHMI,X,A,F,pder

COMMON wavelength,lam0
COMMON blocker,blocker0,lam,dlam
COMMON tuning,tuning,FSR,contrasts

; A[0] = fwidth
; A[1] = Ic
; A[2], A[3], A[4] = phase
; A[5] = fdepth
A[0]  = ABS(A[0]) 
A[1]  = ABS(A[1])
A[5]  = ABS(A[5])

nseq  = N_ELEMENTS(X)
F     = FLTARR(nseq)
dpi   = 2.d0*!dpi

profilef= blocker0  ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0

templine= EXP(-(lam-lam0)^2.d0/A[0]^2.d0)
line    = A[1]-A[5]*templine

dlined1     = DBLARR(nseq)
dlined0     = DBLARR(nseq)
dlined5     = DBLARR(nseq)
dprofilegd2 = DBLARR(nseq)
dprofilegd3 = DBLARR(nseq)
dprofilegd4 = DBLARR(nseq)

FOR j=0,nseq-1 DO BEGIN
    profileg = 0.125d0 * profilef * (1.d0+contrasts*COS(FSR[0]*(lam+A[2])+tuning[0,j])) * (1.d0+contrasts*COS(FSR[1]*(lam+A[3])+tuning[1,j])) * (1.d0+contrasts*COS(FSR[2]*(lam+A[4])+tuning[2,j]))
    F[j]  = TOTAL(line*profileg)*dlam
    dprofilegd2[j] = TOTAL(0.125d0 * profilef * (-contrasts*FSR[0]*SIN(FSR[0]*(lam+A[2])+tuning[0,j])) * (1.d0+contrasts*COS(FSR[1]*(lam+A[3])+tuning[1,j])) * (1.d0+contrasts*COS(FSR[2]*(lam+A[4])+tuning[2,j]))*line)*dlam
    dprofilegd3[j] = TOTAL(0.125d0 * profilef * (-contrasts*FSR[1]*SIN(FSR[1]*(lam+A[3])+tuning[1,j])) * (1.d0+contrasts*COS(FSR[0]*(lam+A[2])+tuning[0,j])) * (1.d0+contrasts*COS(FSR[2]*(lam+A[4])+tuning[2,j]))*line)*dlam
    dprofilegd4[j] = TOTAL(0.125d0 * profilef * (-contrasts*FSR[2]*SIN(FSR[2]*(lam+A[4])+tuning[2,j])) * (1.d0+contrasts*COS(FSR[0]*(lam+A[2])+tuning[0,j])) * (1.d0+contrasts*COS(FSR[1]*(lam+A[3])+tuning[1,j]))*line)*dlam
dlined1[j] = TOTAL(profileg)*dlam
dlined0[j] = TOTAL(-A[5]*templine*(lam-lam0)^2.d0*2.d0/A[0]^3.d0*profileg)*dlam
dlined5[j] = TOTAL(-templine*profileg)*dlam    

ENDFOR

pder = [ [dlined0],[dlined1],[dprofilegd2],[dprofilegd3],[dprofilegd4],[dlined5]]



END




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

PRO HMI_sun4,draw

COMMON wavelength,lam0
COMMON blocker,blocker0,lam,dlam
COMMON tuning,tuning,FSR,contrasts

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

nseq        = 27;36
nx          = 256                 ; number of rows (we rebin the filtergrams)
ny          = 256                 ; number of columns

; MEASURED INTENSITIES
;------------------------------------------------


list='listSun061031_205938'  ; in Calmode


center=[126,129]
distance = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; READ DETUNE SEQUENCE AND GET RID OF OVERSCANS AND BAD EXPOSURES
;---------------------------------------------------------------------------------


RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
dlamdv      = lamref/2.99792458d8 ; to convert Doppler shifts into Doppler velocities
lam0        = dlamdv*lam0         ; positive velocities increase the wavelength
Inten       = imx[*,*,1:27]


;nx=1
;ny=1
;Inten=REBIN(inten,1,1,nseq)
nx=128
ny=128
Inten=REBIN(inten,nx,ny,nseq)
distance=rebin(distance,nx,ny)

;lam0=0.0018532850 ; for listSun060714_221512
;lam0=0.0024861747 ; for listSun060714_223907
;lam0=0.0032872231 ; for listSun060714_231213
;Inten[*,*,18:26] = Inten[*,*,18:26] * 0.875

; TO CORRECT FOR THE I-RIPPLE ON THE NEW NB MICHELSON
;RESTORE,'correction.ps'
;correction=correction/MEAN(correction)
;inten[0,0,*]=inten[0,0,*]/correction


; VARIABLE DEFINITIONS
;--------------------------------------------------------------------------------


Bg          = DBLARR(nx,ny,3)     ; contrasts of the tunable elements
Phig        = DBLARR(nx,ny,3)     ; relative phases
continuum   = DBLARR(nx,ny)       ; relative sun intensity
linewidth   = DBLARR(nx,nx)       ; solar linewidth
linedepth   = DBLARR(nx,nx)       ; solar linedepth
nlam        = 2250                ; number of wavelengths
dlam        = 3.6d0/1.d3          ; resolution in wavelength
lam         =(DINDGEN(nlam)-(nlam-1)/2.)*dlam ; wavelengths in Angstrom RELATIVE TO lamref
tuning      = DBLARR(3,nseq)     ; tuning sequence
lateconv    = INTARR(nx,ny)       ; to discard spurious results
contrasts   = 0.98d0 ;0.98d0
anglim      = 990.;1024;2048;980.             ; in arcseconds

; FREE SPECTRAL RANGES PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;-------------------------------------------------

FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172457d0            ; for the narrow-band Michelson
FSR[1]      = 0.344242d0            ; for the broad-band  Michelson
FSR[2]      = 0.690                 ; for E1
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


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; PROVIDED BY LOCKHEED-MARTIN
;-----------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin' ; average transmission profile
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
q            = READFITS('blocker11.fits')
blocker0     = blocker0 *INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam);centered blocker

IF draw EQ 1 THEN GOTO,draw

; LEAST-SQUARES ALGORITHM
;------------------------------------------------

SET_PLOT,'x'
WINDOW,0,RETAIN=2,xsize=900,ysize=950
!P.MULTI=[0,2,3]


; GUESS
;------------------------

FOR jjj=0,ny-1 DO BEGIN

    TVIM,Phig[*,*,0],/scale,pcharsize=1.5
    TVIM,Phig[*,*,1],/scale,pcharsize=1.5
    TVIM,Phig[*,*,2],/scale,pcharsize=1.5
    TVIM,linewidth  ,/scale,pcharsize=1.5
    TVIM,linedepth  ,/scale,pcharsize=1.5
    TVIM,continuum  ,/scale,pcharsize=1.5

    FOR iii=0,nx-1 DO BEGIN
   
        IF(distance[iii,jjj] LE anglim) THEN BEGIN

           ; GUESS VALUES
           ;------------------------------------------------

            Bg[iii,jjj,*]   = contrasts
            weight = FLTARR(nseq)+1.d0
            A      = [0.065,40000.,0.0,0.0,0.0,23000.]
            resp   = CURVEFIT(FINDGEN(nseq),REFORM(Inten[iii,jjj,*]),weight,A,FUNCTION_NAME='MyHMI',TOL=1.d-7,ITMAX=100,/DOUBLE,STATUS=stat,CHISQ=chi2)

            lateconv[iii,jjj]  = stat
            Phig[iii,jjj,*]    = A[2:4]*FSR[*]*180.d0/!dpi MOD 360; converted into degree units
            IF(A[0] LT 2.0)     THEN linewidth[iii,jjj] = A[0] ELSE linewidth[iii,jjj]=0.0
            IF(A[5] LT 100000.) THEN linedepth[iii,jjj] = A[5] ELSE linedepth[iii,jjj]=0.0
            IF(A[1] LT 100000.) THEN continuum[iii,jjj] = A[1] ELSE continuum[iii,jjj]=0.0
            PRINT,'PARAMETERS='
            FOR i=0,2 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'Ic=',A[1]
            PRINT,'fdepthg=',A[5]
            PRINT,'fwidthg=',A[0]
            PRINT,iii,jjj

        ENDIF
        
    ENDFOR
ENDFOR

SAVE,Phig,Bg,linewidth,linedepth,continuum,lam0,FILE='RESULTS/RESULTScurvefit_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

draw:

; COMPUTE THE AVERAGE PHASES
;--------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a = WHERE(lateconv EQ 0,COMPLEMENT=aa)
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

; TO SUPRESS THE STRIPES AT THE LEFT AND UPPER EDGE
distance[*,0:2]=0.0
distance[251:nx-1,*]=0.0

FOR i=0,2 DO Phig[*,*,i]=Phig[*,*,i]*distance
TVIM,distance

PRINT,'CORRECTED AVERAGED PHASES'
a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 0,COMPLEMENT=aa)
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

SET_PLOT,'ps'
IF(nx GT 1) THEN BEGIN
    !P.MULTI=[0,2,2]
    LOADCT,4
    DEVICE,file='yo.ps',/color,bits=24,xoffset=0,yoffset=0,xsize=21,ysize=20
    titl=STRARR(3)
    titl[0]='Narrow-Band Michelson'
    titl[1]='Broad-Band Michelson'
    titl[2]='Lyot Element E1'
    FOR i=0,2 DO TVIM,Phig[*,*,i],/scale,range=[moy[i]-2.5*dis[i],moy[i]+2.5*dis[i]],tit='!17'+titl[i],stit='!17Phase (in degrees)',barwidth=0.5,pcharsize=0.75
    TVIM,continuum*distance,/scale,stit='!17Continuum intensity',barwidth=0.5,pcharsize=0.75
    DEVICE,/close
    !P.MULTI=0
ENDIF


; RECONSTRUCTION OF THE DETUNE SEQUENCE
;--------------------------------------------------------------------------------------


FOR i=0,2 DO Phig[*,*,i] = Phig[*,*,i]/FSR[i]/180.d0*!dpi ; convert from degrees to radians/FSR

FOR i=0,nseq-1 DO Inten[*,*,i] = Inten[*,*,i]*distance
Inten0 = Inten
a            = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0

profilef     = blocker0            ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrasts*COS(FSR[i]*lam))/2.d0

Inten2  = FLTARR(nx,nx,nseq) 
FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,nx-1 DO BEGIN
        IF(distance[i,j] NE 0.0) THEN BEGIN
            templine     = EXP(-(lam-lam0)^2.d0/linewidth[i,j]^2.d0)
            line         = continuum[i,j] -linedepth[i,j]*templine
                        
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


;SET_PLOT,'x'
!P.MULTI=0
;WINDOW,0,RETAIN=2
device,file='yo3.ps',xsize=15,ysize=12,xoffset=0,yoffset=0
PLOT,REBIN(Inten0,1,1,nseq),xst=1,xtit='position #',ytit='intensity'
OPLOT,REBIN(Inten2,1,1,nseq),linestyle=2
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten0,1,1,nseq))-REFORM(REBIN(Inten2,1,1,nseq)))^2.0 )





READ,pause

END
