; PROGRAM TO PERFORM THE THROUGHPUT CALIBRATION TEST
; WITH THE SUN AS SOURCE
; OBJECTIVE: TO KNOW WHAT IS THE OVERALL TRANSMISSION OF THE HMI
; INSTRUMENT
;
; we write:
; ln(Ic) = ln(Ic0) - k X
; where Ic is the solar continuum intensity obtained by detune
; sequences with HMI_sun2.pro,
; Ic0 is the continuum intensity in the absence of atmospheric
; absorption/scattering, k is the extinction coefficient at
; 6173 A, and X is the airmass (returned by suninfo)
; This is the Beer-Bougher-Lambert absorption law, approximate above 500nm

;----------------------------------------------------------------------


;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------


PRO HMI_sun3

nmes       = 4 ; number of detune sequences done during the day
airmass    = DBLARR(nmes)
contin     = DBLARR(nmes)
CCDgain    = 16.3611d0   ; FOR SIDE CAMERA, FROM ANALYSIS OF FSN=1420880-1420945; to convert DNs into photoelectrons
nx         = 256.
list       = STRARR(nmes)

;list[0]    = 'listSun060222_204505'
;list[1]    = 'listSun060222_174615'
;list[2]    = 'listSun060222_225521'
;list[3]    = 'listSun060222_002209'


list[0]='listSun070905_553018'
list[1]='listSun070905_553222'
list[2]='listSun070905_548149'
list[3]='listSun070905_554777'

;year       = [2006,2006,2006,2006]
;month      = [   2,   2,   2,   2]
;day        = [  22,  22,  22,  23]
;heure      = [  20,  17,  23,   0]
;minute     = [  50,  51,   0,  27]
;seconde    = [    0,  0,   0,   0] 


year       = [2007, 2007, 2007, 2007]
month      = [   9,    9,    9,    9]
day        = [   5,    5,    4,    5]
heure      = [  16,   17,   17,   19]
minute     = [  28,   02,   46,   56]
seconde    = [  54,   04,   54,   24] 

;ratio      = MEAN([6.7/9.7,6.3/9.3,5.4/8.0]) ; ratio power inside LMSAL/power outside, measured by Jesper
;ratio      = [0.681047,0.681047,0.681047,0.681047]
ratio       = [0.551d0/1.030d0,0.694d0/1.195d0,0.694d0/1.229d0,0.703d0/1.226d0]
exposure    = [100.d0,200.d0,200.d0,200.d0] ;HMI_FSW_IMG_CMDED_EXPOSURE FOR SIDE CAMERA (CAMID=0)
Integration = MAX(exposure)/1000.d0   ; exposure time in seconds


; THE FILTER EQUIVALENT WIDTH:
nlam        = 18000               ; number of wavelengths
dlam        = 1.0d0/1.d3         ; resolution in wavelength
lam         = dlam*(DINDGEN(nlam)-(nlam-1)/2.)
FSR         = DBLARR(7)    ; FSR in Angstrom
FSR[0]      = 0.172-0.0010576d0; for the narrow-band Michelson
FSR[1]      = 0.344-0.00207683d0; for the broad-band  Michelson
FSR[2]      = 0.693+0.000483467d0  ; for E1
FSR[3]      = 1.407d0;1.405d0       ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5

lamref      = 6173.3433d0


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


profile     = blocker0
contrasts   = [0.97723184d0,0.99716162d0,0.97434975d0,0.98593564d0,0.96497714d0,0.98777667d0,0.99917262d0]
phases      = [0.d0,0.d0,0.d0,-6.138939d0,-0.17368478d0,-5.3398232d0,1.1401513d0]*!dpi/180.d0
FOR i=0,6 DO profile = profile * (1.d0+contrasts[i]*COS(2.d0*!dpi/FSR[i]*lam+phases[i]))/2.d0
;eqwidth     = TOTAL(profile)*dlam 
area        = INT_TABULATED(lam,profile) ; SURFACE AREA
eqwidth     = area/MAX(profile)          ; EQUIVALENT WIDTH IN ANGSTROM
PRINT,'HMI FILTER AREA',area
PRINT,'HMI FILTER EQUIVALENT WIDTH',eqwidth 


FOR i=0,nmes-1 DO BEGIN

    PRINT,'FILE=',list[i]


  ; RESTORE FILES CONTAINING THE CONTINUUM
  ; INTENSITY MEASURED AT DIFFERENT TIMES
  ; OF THE DAY

    RESTORE,'RESULTS/RESULTS_'+STRTRIM(list[i],1)+'_side_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
    continuum=continuum*10000.d0 ; for September 2007 data (factor in HMI_sun1.pro)
    
   ;IF(i EQ 0 OR i EQ 2) THEN TVIM,continuum,/scale,pcharsize=1.5,barwidth=0.5,tit=list[i],rmarg=2 ELSE TVIM,continuum,/scale,pcharsize=1.5,barwidth=0.5,tit=list[i],rmarg=2
    
    center   = [(nx-1.)/2.,(nx-1.)/2.] ; !!!WARNING, depends on the images
    distance=FLTARR(nx,nx)
    for i=0,nx-1 do for j=0,nx-1 do distance[i,j]=SQRT((i-center[0])^2.d0 + (j-center[1])^2.d0)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center
     a=WHERE(distance LT 900.d0,COMPLEMENT=b) ; to avoid the limb
     contin[i] = MEAN(continuum[a])/exposure[i]*MAX(exposure) ; to correct for different exposure times



                                ; COMPUTE THE EXTINCTION COEFFICIENT
                                ; OF THE ATMOSPHERE (MUST BE ON n00 TO
                                ; RUN suninfo)

    command = '/auto/home0/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year[i]),1)+'.'+STRTRIM(STRING(month[i]),1)+'.'+STRTRIM(STRING(day[i]),1)+'.'+STRTRIM(STRING(heure[i]),1)+'.'+STRTRIM(STRING(minute[i]),1)+'.'+STRTRIM(STRING(seconde[i]),1)+'| tail -10 | head -1 | cut -b19-32'
    PRINT,command
    SPAWN,command,res
    airmass[i] = DOUBLE(res)
    
ENDFOR

PRINT,'AIRMASSES=',airmass
contin = contin / ratio

; we want to determine the absolute continuum intensity
; measured by HMI, as if there were no air
; (no atmospheric extinction, ie. airmass=0)
; EXTRAPOLATE AT AIRMASS=0
contin = contin[SORT(airmass)]
airmass= airmass[SORT(airmass)]
res    = poly_fit(airmass,ALOG(contin),1,yfit=y) ; WE USE BEER-BOUGHER-LAMBERT LAW
Ic0    = EXP(res[0]) ; CONTINUUM INTENSITY IN THE ABSENCE OF ATMOSPHERIC ABSORPTION/SCATTERING (AIRMASS=0)

airmass= [0.d0,airmass]
contin = [Ic0,contin]
Ic0    = (Ic0*area) *CCDgain/ Integration ; number of photo-electrons detected on a CCD pixel per second

; photon flux per pixel per second
; computed by J. Schou (Radiometric Analysis - Photon Counts)
flux   = 7.4377396d31*5.88d-12*0.015393804d0*1.d-10*eqwidth ; number of photons per seconds per pixel in the passband, received at the entrance of the HMI instrument 

throughput = Ic0 / flux

PRINT,'TROUGHPUT=',throughput
SET_PLOT,'ps'
DEVICE,file='yo.ps',bits=24,/color,xsize=20,ysize=14,xoffset=0.5,yoffset=0
!P.MULTI=0
loadct,3
PLOT,airmass,contin *CCDgain/ Integration,psym=4,tit='!17',xtit='Airmass',ytit='Ic (photons/s/pixel/A)',charsize=1.5,thick=2;,tit='!17Throughput='+STRTRIM(STRING(throughput),1)
airmass=FINDGEN(1000)/999.*2.0
oplot,airmass,EXP(res[0]+airmass*res[1])*CCDgain/ Integration,thick=2

DEVICE,/close

READ,pause

END

