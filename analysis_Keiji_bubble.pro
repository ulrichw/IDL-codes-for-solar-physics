; program to compute Doppler velocities and l.o.s. filed strength
; at 3 different heights for the datacubes related to Keiji's bubble


LCP5=fitsio_read_image("/SUM5/D153795747/S00000/trackedLCP5.fits")
LCP0=fitsio_read_image("/SUM2/D153768214/S00000/trackedLCP0.fits")
velocityLCP50 = float(LCP5-LCP0)/float(LCP5+LCP0)
LCP0=0.0
LCP5=0.0

RCP5=fitsio_read_image("/SUM11/D153840023/S00000/trackedRCP5.fits")
RCP0=fitsio_read_image("/SUM0/D153802376/S00000/trackedRCP0.fits")
velocityRCP50 = float(RCP5-RCP0)/float(RCP5+RCP0)
RCP0=0.0
RCP5=0.0

LCP4=fitsio_read_image("/SUM16/D153791091/S00000/trackedLCP4.fits")
LCP1=fitsio_read_image("/SUM22/D153773285/S00000/trackedLCP1.fits")
velocityLCP41 = float(LCP4-LCP1)/float(LCP4+LCP1)
LCP4=0.0
LCP1=0.0

RCP4=fitsio_read_image("/SUM0/D153832600/S00000/trackedRCP4.fits")
RCP1=fitsio_read_image("/SUM20/D153808415/S00000/trackedRCP1.fits")
velocityRCP41 = float(RCP4-RCP1)/float(RCP4+RCP1)
RCP4=0.0
RCP1=0.0

LCP2=fitsio_read_image("/SUM13/D153777175/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM18/D153781723/S00000/trackedLCP3.fits")
velocityLCP32 = float(LCP3-LCP2)/float(LCP3+LCP2)
LCP3=0.0
LCP2=0.0

RCP2=fitsio_read_image("/SUM7/D153814356/S00000/trackedRCP2.fits")
RCP3=fitsio_read_image("/SUM18/D153825404/S00000/trackedRCP3.fits")
velocityRCP32 = float(RCP3-RCP2)/float(RCP3+RCP2)
RCP3=0.0
RCP2=0.0

cal50=TOTAL(TOTAL((velocityLCP50+velocityRCP50)/2.,1),1)/512./256.
cal41=TOTAL(TOTAL((velocityLCP41+velocityRCP41)/2.,1),1)/512./256.
cal32=TOTAL(TOTAL((velocityLCP32+velocityRCP32)/2.,1),1)/512./256.
cal32=smooth(cal32,15,/edge_truncate)
cal41=smooth(cal41,15,/edge_truncate)
cal50=smooth(cal50,15,/edge_truncate)
OPENR,1,'velocitiesKEIJIBUBBLE'
vel=FLTARR(960)
READF,1,vel
CLOSE,1
magnetic = 1.d0/(2.d0*4.67d-5*0.000061733433d0*2.5d0*299792458.d0)

; PRODUCE LOOK-UP TABLE
;res  =poly_fit(vel,cal41,1,yfit=y)
;coeff41=res[1]
;res  =poly_fit(vel,cal50,1,yfit=y)
;coeff50=res[1]
;res  =poly_fit(vel,cal32,1,yfit=y)
;coeff32=res[1]
;Doppler50=(velocityLCP50+velocityRCP50)/2.d0/coeff50
;magnetic50=(velocityLCP50-velocityRCP50)/coeff50*magnetic
;velocityLCP50=0
;velocityRCP50=0
;Doppler41=(velocityLCP41+velocityRCP41)/2.d0/coeff41
;magnetic41=(velocityLCP41-velocityRCP41)/coeff41*magnetic
;velocityLCP41=0
;velocityRCP41=0
;Doppler32=(velocityLCP32+velocityRCP32)/2.d0/coeff32
;magnetic32=(velocityLCP32-velocityRCP32)/coeff32*magnetic
;velocityLCP32=0
;velocityRCP32=0

RESTORE,'lookup.bin' ; obtained with doppler_different_heights.pro
vtest=vtest+185.0
cal32=ratio[0,*]*1.40
cal41=ratio[1,*]*1.165
cal50=ratio[2,*]*1.15
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

a = WHERE(velocityLCP50 gt maxi2) ; we saturate the velocities
b = WHERE(velocityLCP50 lt mini2) ; we saturate the velocities
velocityLCP50 = INTERPOL(vtest[b2[nb2-1]:a2[0]],REFORM(ratio[2,b2[nb2-1]:a2[0]]),velocityLCP50)
IF(a[0] NE -1) THEN velocityLCP50[a]=vmax2
IF(b[0] NE -1) THEN velocityLCP50[b]=vmin2
a = WHERE(velocityRCP50 gt maxi2) ; we saturate the velocities
b = WHERE(velocityRCP50 lt mini2) ; we saturate the velocities
velocityRCP50 = INTERPOL(vtest[b2[nb2-1]:a2[0]],REFORM(ratio[2,b2[nb2-1]:a2[0]]),velocityRCP50)
IF(a[0] NE -1) THEN velocityRCP50[a]=vmax2
IF(b[0] NE -1) THEN velocityRCP50[b]=vmin2

a = WHERE(velocityLCP41 gt maxi1) ; we saturate the velocities
b = WHERE(velocityLCP41 lt mini1) ; we saturate the velocities
velocityLCP41 = INTERPOL(vtest[b1[nb1-1]:a1[0]],REFORM(ratio[1,b1[nb1-1]:a1[0]]),velocityLCP41)
IF(a[0] NE -1) THEN velocityLCP41[a]=vmax1
IF(b[0] NE -1) THEN velocityLCP41[b]=vmin1
a = WHERE(velocityRCP41 gt maxi1) ; we saturate the velocities
b = WHERE(velocityRCP41 lt mini1) ; we saturate the velocities
velocityRCP41 = INTERPOL(vtest[b1[nb1-1]:a1[0]],REFORM(ratio[1,b1[nb1-1]:a1[0]]),velocityRCP41)
IF(a[0] NE -1) THEN velocityRCP41[a]=vmax1
IF(b[0] NE -1) THEN velocityRCP41[b]=vmin1

a = WHERE(velocityLCP32 gt maxi0) ; we saturate the velocities
b = WHERE(velocityLCP32 lt mini0) ; we saturate the velocities
velocityLCP32 = INTERPOL(vtest[b0[nb0-1]:a0[0]],REFORM(ratio[0,b0[nb0-1]:a0[0]]),velocityLCP32)
IF(a[0] NE -1) THEN velocityLCP32[a]=vmax0
IF(b[0] NE -1) THEN velocityLCP32[b]=vmin0
a = WHERE(velocityRCP32 gt maxi0) ; we saturate the velocities
b = WHERE(velocityRCP32 lt mini0) ; we saturate the velocities
velocityRCP32 = INTERPOL(vtest[b0[nb0-1]:a0[0]],REFORM(ratio[0,b0[nb0-1]:a0[0]]),velocityRCP32)
IF(a[0] NE -1) THEN velocityRCP32[a]=vmax0
IF(b[0] NE -1) THEN velocityRCP32[b]=vmin0

Doppler50=(velocityLCP50+velocityRCP50)/2.d0
magnetic50=(velocityLCP50-velocityRCP50)*magnetic
velocityLCP50=0
velocityRCP50=0
Doppler41=(velocityLCP41+velocityRCP41)/2.d0
magnetic41=(velocityLCP41-velocityRCP41)*magnetic
velocityLCP41=0
velocityRCP41=0
Doppler32=(velocityLCP32+velocityRCP32)/2.d0
magnetic32=(velocityLCP32-velocityRCP32)*magnetic
velocityLCP32=0
velocityRCP32=0


; POLYNOMIAL CORRECTION
error=TOTAL(TOTAL(Doppler50,1),1)/512./256.-vel
res=poly_fit(vel,error,3,yfit=correction50)
error=TOTAL(TOTAL(Doppler41,1),1)/512./256.-vel
res=poly_fit(vel,error,3,yfit=correction41)
error=TOTAL(TOTAL(Doppler32,1),1)/512./256.-vel
res=poly_fit(vel,error,3,yfit=correction32)
FOR i=0,959 DO BEGIN & Doppler50[*,*,i]=Doppler50[*,*,i]-correction50[i] & Doppler41[*,*,i]=Doppler41[*,*,i]-correction41[i] & Doppler32[*,*,i]=Doppler32[*,*,i]-correction32[i] & ENDFOR

writefits, 'Doppler50_KEIJIBUBBLE.fits',Doppler50
writefits, 'Doppler41_KEIJIBUBBLE.fits',Doppler41
writefits, 'Doppler32_KEIJIBUBBLE.fits',Doppler32
writefits,'magnetic50_KEIJIBUBBLE.fits',magnetic50
writefits,'magnetic41_KEIJIBUBBLE.fits',magnetic41
writefits,'magnetic32_KEIJIBUBBLE.fits',magnetic32


READ,pause

END
