; TO CALCULATE DOPPLER VELOCITIES AND LOS MAGNETIC FIELD STRENGTHS FOR
; NOAA 11092 OF 2010/8/3 WITH I0-I5
; 0:0:0_TAI to 14:0:0_TAI
; OBS_VR varies between 1972.12  and -2608.95 m/s

FUNCTION minimize,x,y,vel
x=REFORM(x)
y=REFORM(y)
yorge=SORT(vel)
x=x[yorge]
y=y[yorge]
veli=FINDGEN(1200)/1199.*(MAX(vel)-MIN(vel))+MIN(vel)
x=INTERPOL(x,vel,veli)
y=INTERPOL(y,vel,veli)

ns=200

nx=N_ELEMENTS(x)

mini=TOTAL((x-y)^2.d0)
indexm=0

FOR i=1,ns DO BEGIN
    x2=SHIFT(x,i)
    x3=SHIFT(x,-i)
    yo=TOTAL((x2[ns:nx-ns]-y[ns:nx-ns])^2.d0)
    yo2=TOTAL((x3[ns:nx-ns]-y[ns:nx-ns])^2.d0)
    IF(yo LE mini) THEN BEGIN & mini=yo & indexm=i & ENDIF
    IF(yo2 LE mini) THEN BEGIN & mini=yo2 & indexm=-i & ENDIF
ENDFOR
PRINT,mini
PLOT,vel[yorge],SHIFT(x,indexm)
OPLOT,vel[yorge],y,col=180

RETURN,vel[nx/2+indexm]-vel[nx/2]

END

PRO lookupfit32,X,A,F,pder

COMMON filter,looktot32,looktot41,looktot50

F=INTERPOL(looktot32,X,X-A[0])*A[2]+A[1]

END

PRO lookupfit41,X,A,F,pder

COMMON filter,looktot32,looktot41,looktot50

F=INTERPOL(looktot41,X,X-A[0])*A[2]+A[1]

END

PRO lookupfit50,X,A,F,pder

COMMON filter,looktot32,looktot41,looktot50

F=INTERPOL(looktot50,X,X-A[0])*A[2]+A[1]

END

PRO analysis_NOAA11092,fitonl

COMMON filter,looktot32,looktot41,looktot50

IF(fitonl EQ 1) THEN GOTO,fitonly

;velocityLCP50=fitsio_read_image("/SUM9/D316561424/S00000/trackedv50.fits",/single)  ; NOAA 11092 in su_couvidat.v50LCP_tracked (see mtrack_script)
;velocityRCP50=fitsio_read_image("/SUM36/D316563570/S00000/trackedv50.fits",/single) ; in su_couvidat.v50RCP_tracked 

;velocityLCP50=fitsio_read_image("/SUM3/D316403443/S00000/trackedv50.fits",/single)  ; NOAA 11161 in su_couvidat.v50LCP_tracked (see mtrack_script)
;velocityRCP50=fitsio_read_image("/SUM34/D316406058/S00000/trackedv50.fits",/single) ; in su_couvidat.v50RCP_tracked 

;velocityLCP50=fitsio_read_image("/SUM34/D316301653/S00000/trackedv50.fits",/single)  ; NOAA 11243 in su_couvidat.v50LCP_tracked (see mtrack_script)
;velocityRCP50=fitsio_read_image("/SUM36/D316304575/S00000/trackedv50.fits",/single) ; in su_couvidat.v50RCP_tracked 

velocityLCP50=fitsio_read_image("/SUM6/D316586947/S00000/trackedv50.fits",/single)  ; NOAA 11306 in su_couvidat.v50LCP_tracked (see mtrack_script)
velocityRCP50=fitsio_read_image("/SUM9/D316591344/S00000/trackedv50.fits",/single) ; in su_couvidat.v50RCP_tracked

;GOTO,velocity
;-----------------------------------------------------------------------------------------

;nvel=1121
;OPENR,1,'velocitiesNOAA11092' ; file containing OBS_VR
;nvel=961
;OPENR,1,'velocitiesNOAA11161'
;nvel=1121
;OPENR,1,'velocitiesNOAA11243'
nvel=1121
OPENR,1,'velocitiesNOAA11306'
 
vel=FLTARR(nvel)
READF,1,vel
CLOSE,1

cal50=FLTARR(512,nvel)
FOR i=10,501 DO BEGIN & cal50[i,*]=TOTAL(TOTAL((velocityLCP50[i-10:i+10,185:205,*]+velocityRCP50[i-10:i+10,185:205,*])/2.,1),1)/21./21. & ENDFOR 
FOR i=0,nvel-1 DO cal50[0:9,i]=cal50[10,i]
FOR i=0,nvel-1 DO cal50[502:511,i]=cal50[501,i]


; LONGITUDE-AVERAGED LOOK-UP TABLES
caltot50=total(cal50,1)/512.0
;SAVE,vel,caltot50,FILE='NOAA11092_caltot2.bin'
;SAVE,vel,caltot50,FILE='NOAA11161_caltot2.bin'
;SAVE,vel,caltot50,FILE='NOAA11243_caltot2.bin'
SAVE,vel,caltot50,FILE='NOAA11306_caltot2.bin'


fitonly:

RESTORE,'NOAA11306_caltot2.bin'

res=POLY_FIT(vel,caltot50,3,yfit=y)
vtest=FINDGEN(821)*24.-410.*24.
d=POLY(vtest,res)

; APPLY THE LOOK-UP TABLES

velocityLCP50 = INTERPOL(vtest,d,velocityLCP50)
velocityRCP50 = INTERPOL(vtest,d,velocityRCP50)

magnetic = 1.d0/(2.d0*4.67d-5*0.000061733433d0*2.5d0*299792458.d0)
Doppler50=(velocityLCP50+velocityRCP50)/2.0
magnetic50=(velocityLCP50-velocityRCP50)*magnetic
velocityLCP50=0
velocityRCP50=0

; WRITE OUTPUT FILES
Doppler50 =FLOAT(Doppler50)
magnetic50=FLOAT(magnetic50)
writefits, 'Doppler50_NOAA11306_2.fits',Doppler50
writefits,'magnetic50_NOAA11306_2.fits',magnetic50


READ,pause

END
