; TO CALCULATE DOPPLER VELOCITIES AND LOS MAGNETIC FIELD STRENGTHS FOR
; NOAA 11243 OF 2011/7/3 AT 3 DIFFERENT ATMOSPHERIC HEIGHTS 
; 0:0:0_TAI to 12:0:0_TAI
; OBS_VR varies between -1957 m/s and  1981 m/s

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



PRO analysis_NOAA11243,fitonl

COMMON filter,looktot32,looktot41,looktot50

IF(fitonl EQ 1) THEN GOTO,fitonly

LCP0=fitsio_read_image("/SUM12/D210064977/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM9/D210122333/S00000/trackedLCP5.fits")
velocityLCP50 = float(LCP5-LCP0)/float(LCP5+LCP0)
LCP0=0.0
LCP5=0.0
RCP5=fitsio_read_image("/SUM7/D210163860/S00000/trackedRCP5.fits")
RCP0=fitsio_read_image("/SUM11/D210122533/S00000/trackedRCP0.fits")
velocityRCP50 = float(RCP5-RCP0)/float(RCP5+RCP0)
RCP0=0.0
RCP5=0.0
LCP1=fitsio_read_image("/SUM5/D210077961/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM8/D210106006/S00000/trackedLCP4.fits")
velocityLCP41 = float(LCP4-LCP1)/float(LCP4+LCP1)
LCP1=0.0
LCP4=0.0
RCP1=fitsio_read_image("/SUM15/D210135288/S00000/trackedRCP1.fits")
RCP4=fitsio_read_image("/SUM0/D210148320/S00000/trackedRCP4.fits")
velocityRCP41 = float(RCP4-RCP1)/float(RCP4+RCP1)
RCP1=0.0
RCP4=0.0
LCP2=fitsio_read_image("/SUM22/D210094364/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM11/D210105935/S00000/trackedLCP3.fits")
velocityLCP32 = float(LCP3-LCP2)/float(LCP3+LCP2)
LCP2=0.0
LCP3=0.0  
RCP2=fitsio_read_image("/SUM0/D210135384/S00000/trackedRCP2.fits")
RCP3=fitsio_read_image("/SUM3/D210148241/S00000/trackedRCP3.fits")
velocityRCP32 = float(RCP3-RCP2)/float(RCP3+RCP2)
RCP2=0.0
RCP3=0.0


OPENR,1,'velocitiesNOAA11243' ; file containing OBS_VR
vel=FLTARR(1121)
READF,1,vel
CLOSE,1

cal50=FLTARR(512,1121)
cal41=cal50
cal32=cal50
FOR i=10,501 DO BEGIN & cal50[i,*]=TOTAL(TOTAL((velocityLCP50[i-10:i+10,185:205,*]+velocityRCP50[i-10:i+10,185:205,*])/2.,1),1)/21./21.& cal41[i,*]=TOTAL(TOTAL((velocityLCP41[i-10:i+10,185:205,*]+velocityRCP41[i-10:i+10,185:205,*])/2.,1),1)/21./21. & cal32[i,*]=TOTAL(TOTAL((velocityLCP32[i-10:i+10,185:205,*]+velocityRCP32[i-10:i+10,185:205,*])/2.,1),1)/21./21. & ENDFOR
FOR i=0,960 DO cal32[0:9,i]=cal32[10,i]
FOR i=0,960 DO cal41[0:9,i]=cal41[10,i]
FOR i=0,960 DO cal50[0:9,i]=cal50[10,i]
FOR i=0,960 DO cal32[502:511,i]=cal32[501,i]
FOR i=0,960 DO cal41[502:511,i]=cal41[501,i]
FOR i=0,960 DO cal50[502:511,i]=cal50[501,i]


; LONGITUDE-AVERAGED LOOK-UP TABLES
caltot32=total(cal32,1)/512.0
caltot41=total(cal41,1)/512.0
caltot50=total(cal50,1)/512.0
SAVE,vel,caltot32,caltot41,caltot50,FILE='NOAA11243_caltot.bin'


;cal50=SMOOTH(cal50,50,/edge_truncate)
;cal41=SMOOTH(cal41,50,/edge_truncate)
;cal32=SMOOTH(cal32,50,/edge_truncate)
;
;cal502=FLTARR(512,961)
;cal412=cal50
;cal322=cal50
;FOR i=10,501 DO BEGIN & cal502[i,*]=TOTAL(TOTAL((velocityLCP50[i-10:i+10,0:20,*]+velocityRCP50[i-10:i+10,0:20,*])/2.,1),1)/21./21.& cal412[i,*]=TOTAL(TOTAL((velocityLCP41[i-10:i+10,0:20,*]+velocityRCP41[i-10:i+10,0:20,*])/2.,1),1)/21./21. & cal322[i,*]=TOTAL(TOTAL((velocityLCP32[i-10:i+10,0:20,*]+velocityRCP32[i-10:i+10,0:20,*])/2.,1),1)/21./21. & ENDFOR
;FOR i=0,960 DO cal322[0:9,i]=cal32[10,i]
;FOR i=0,960 DO cal412[0:9,i]=cal41[10,i]
;FOR i=0,960 DO cal502[0:9,i]=cal50[10,i]
;FOR i=0,960 DO cal322[502:511,i]=cal32[501,i]
;FOR i=0,960 DO cal412[502:511,i]=cal41[501,i]
;FOR i=0,960 DO cal502[502:511,i]=cal50[501,i];
;
;cal502=SMOOTH(cal502,50,/edge_truncate)
;cal412=SMOOTH(cal412,50,/edge_truncate)
;cal322=SMOOTH(cal322,50,/edge_truncate)

fitonly:

RESTORE,'NOAA11243_caltot.bin'
RESTORE,'lookup_center.bin'; look-up table calculated by at disk center [2048,2048], with calibration11=0 and phase maps of November 2010, in doppler_different_heights.pro

looktot32=INTERPOL(REFORM(ratio[0,*]),vtest,vel)
looktot41=INTERPOL(REFORM(ratio[1,*]),vtest,vel)
looktot50=INTERPOL(REFORM(ratio[2,*]),vtest,vel)


; FIT OF LOOK-UP TABLES TO MATCH THE CALTOT CURVES
!P.MULTI=[0,1,3]

parameters32= [0.d0,0.d0,1.d0]
weights   = DBLARR(N_ELEMENTS(vel))+1.d0
resp      = CURVEFIT(vel,caltot32,weights,parameters32,FUNCTION_NAME='lookupfit32',TOL=1.d-9,ITMAX=3000,/NODERIVATIVE)  

PLOT,vel,caltot32,xst=1,xrange=[-5000,5000],charsize=1.5
OPLOT,vtest,INTERPOL(REFORM(ratio[0,*]),vtest,vtest-parameters32[0])*parameters32[2]+parameters32[1],col=180

parameters41= [0.d0,0.d0,1.d0]
resp      = CURVEFIT(vel,caltot41,weights,parameters41,FUNCTION_NAME='lookupfit41',TOL=1.d-9,ITMAX=3000,/NODERIVATIVE)  

PLOT,vel,caltot41,xst=1,xrange=[-5000,5000],charsize=1.5,yrange=[-0.3,0.3],yst=1
OPLOT,vtest,INTERPOL(REFORM(ratio[1,*]),vtest,vtest-parameters41[0])*parameters32[2]+parameters41[1],col=180

parameters50= [0.d0,0.d0,1.d0]
resp      = CURVEFIT(vel,caltot50,weights,parameters50,FUNCTION_NAME='lookupfit50',TOL=1.d-9,ITMAX=3000,/NODERIVATIVE)  

PLOT,vel,caltot50,xst=1,xrange=[-5000,5000],charsize=1.5,yrange=[-0.05,0.07],yst=1
OPLOT,vtest,INTERPOL(REFORM(ratio[2,*]),vtest,vtest-parameters50[0])*parameters32[2]+parameters50[1],col=180

; APPLY THE LOOK-UP TABLES

ratio=FLOAT(ratio)

ratio[0,*]=INTERPOL(REFORM(ratio[0,*]),vtest,vtest-parameters32[0])*parameters32[2]+parameters32[1]
ratio[1,*]=INTERPOL(REFORM(ratio[1,*]),vtest,vtest-parameters41[0])*parameters41[2]+parameters41[1]
ratio[2,*]=INTERPOL(REFORM(ratio[2,*]),vtest,vtest-parameters50[0])*parameters50[2]+parameters50[1]

SAVE,ratio,vtest,file='lookup_NOAA11243.bin'


vtest=FLOAT(vtest)
maxi0=MAX(ratio[0,*])
maxi1=MAX(ratio[1,*])
maxi2=MAX(ratio[2,*])
mini0=MIN(ratio[0,*])
mini1=MIN(ratio[1,*])
mini2=MIN(ratio[2,*])
a0=WHERE(ratio[0,*] EQ maxi0)
vmax0=vtest[a0[0]]
b0=WHERE(ratio[0,*] EQ mini0,nb0)
vmin0=vtest[b0[nb0-1]]
a1=WHERE(ratio[1,*] EQ maxi1)
vmax1=vtest[a1[0]]
b1=WHERE(ratio[1,*] EQ mini1,nb1)
vmin1=vtest[b1[nb1-1]]
a2=WHERE(ratio[2,*] EQ maxi2)
vmax2=vtest[a2[0]]
b2=WHERE(ratio[2,*] EQ mini2,nb2)
vmin2=vtest[b2[nb2-1]]

; WE MODIFY ratio[2,*] SO THAT IT IS UNIQUELY DEFINED (OTHERWISE A
; GIVEN RATIO CAN MEAN SEVERAL VELOCITIES, AND THE INTERPOLATION IS
; SCREWED UP !!!!):
ratio[2,a2:N_ELEMENTS(vtest)-1]=ratio[2,a2]
ratio[2,0:b2[nb2-1]]=ratio[2,b2[nb2-1]]


a=WHERE(velocityLCP50 GT maxi2)
b=WHERE(velocityLCP50 LT mini2)
velocityLCP50 = INTERPOL(vtest,ratio[2,*],velocityLCP50)
IF(a[0] NE -1) THEN velocityLCP50[a]=vmax2
IF(b[0] NE -1) THEN velocityLCP50[b]=vmin2

a=WHERE(velocityRCP50 GT maxi2)
b=WHERE(velocityRCP50 LT mini2)
velocityRCP50 = INTERPOL(vtest,ratio[2,*],velocityRCP50)
IF(a[0] NE -1) THEN velocityRCP50[a]=vmax2
IF(b[0] NE -1) THEN velocityRCP50[b]=vmin2


; WE MODIFY ratio[1,*] SO THAT IT IS UNIQUELY DEFINED (OTHERWISE A
; GIVEN RATIO CAN MEAN SEVERAL VELOCITIES, AND THE INTERPOLATION IS
; SCREWED UP !!!!):
ratio[1,a1:N_ELEMENTS(vtest)-1]=ratio[1,a1]
ratio[1,0:b1[nb1-1]]=ratio[1,b1[nb1-1]]


a=WHERE(velocityLCP41 GT maxi1)
b=WHERE(velocityLCP41 LT mini1)
velocityLCP41 = INTERPOL(vtest,ratio[1,*],velocityLCP41)
IF(a[0] NE -1) THEN velocityLCP41[a]=vmax1
IF(b[0] NE -1) THEN velocityLCP41[b]=vmin1

a=WHERE(velocityRCP41 GT maxi1)
b=WHERE(velocityRCP41 LT mini1)
velocityRCP41 = INTERPOL(vtest,ratio[1,*],velocityRCP41)
IF(a[0] NE -1) THEN velocityRCP41[a]=vmax1
IF(b[0] NE -1) THEN velocityRCP41[b]=vmin1


; WE MODIFY ratio[0,*] SO THAT IT IS UNIQUELY DEFINED (OTHERWISE A
; GIVEN RATIO CAN MEAN SEVERAL VELOCITIES, AND THE INTERPOLATION IS
; SCREWED UP !!!!):
ratio[0,a0:N_ELEMENTS(vtest)-1]=ratio[0,a0]
ratio[0,0:b0[nb0-1]]=ratio[0,b0[nb0-1]]

a=WHERE(velocityLCP32 GT maxi0)
b=WHERE(velocityLCP32 LT mini0)
velocityLCP32 = INTERPOL(vtest,ratio[0,*],velocityLCP32)
IF(a[0] NE -1) THEN velocityLCP32[a]=vmax0
IF(b[0] NE -1) THEN velocityLCP32[b]=vmin0

a=WHERE(velocityRCP32 GT maxi0)
b=WHERE(velocityRCP32 LT mini0)
velocityRCP32 = INTERPOL(vtest,ratio[0,*],velocityRCP32)
IF(a[0] NE -1) THEN velocityRCP32[a]=vmax0
IF(b[0] NE -1) THEN velocityRCP32[b]=vmin0

magnetic = 1.d0/(2.d0*4.67d-5*0.000061733433d0*2.5d0*299792458.d0)
Doppler50=(velocityLCP50+velocityRCP50)/2.0
magnetic50=(velocityLCP50-velocityRCP50)*magnetic
velocityLCP50=0
velocityRCP50=0
Doppler41=(velocityLCP41+velocityRCP41)/2.0
magnetic41=(velocityLCP41-velocityRCP41)*magnetic
velocityLCP41=0
velocityRCP41=0
Doppler32=(velocityLCP32+velocityRCP32)/2.0
magnetic32=(velocityLCP32-velocityRCP32)*magnetic
velocityLCP32=0
velocityRCP32=0
; WRITE OUTPUT FILES
Doppler50 =FLOAT(Doppler50)
Doppler41 =FLOAT(Doppler41)
Doppler32 =FLOAT(Doppler32)
magnetic50=FLOAT(magnetic50)
magnetic41=FLOAT(magnetic41)
magnetic32=FLOAT(magnetic32)

writefits, 'Doppler50_NOAA11243.fits',Doppler50
writefits, 'Doppler41_NOAA11243.fits',Doppler41
writefits, 'Doppler32_NOAA11243.fits',Doppler32
writefits,'magnetic50_NOAA11243.fits',magnetic50
writefits,'magnetic41_NOAA11243.fits',magnetic41
writefits,'magnetic32_NOAA11243.fits',magnetic32


READ,pause

END
