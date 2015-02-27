quietsun=1

IF(quietsun EQ 0) THEN BEGIN
; QUIET SUN
;---------------

LCP0=fitsio_read_image("/SUM9/D156478242/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM0/D156485787/S00000/trackedLCP5.fits")
velocityLCP50 = float(LCP5-LCP0)/float(LCP5+LCP0)
LCP0=0.0
LCP5=0.0
RCP5=fitsio_read_image("/SUM20/D156496127/S00000/trackedRCP5.fits")
RCP0=fitsio_read_image("/SUM14/D156504791/S00000/trackedRCP0.fits")
velocityRCP50 = float(RCP5-RCP0)/float(RCP5+RCP0)
RCP0=0.0
RCP5=0.0
LCP1=fitsio_read_image("/SUM17/D156480000/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM7/D156484565/S00000/trackedLCP4.fits")
velocityLCP41 = float(LCP4-LCP1)/float(LCP4+LCP1)
LCP1=0.0
LCP4=0.0
RCP1=fitsio_read_image("/SUM17/D156503069/S00000/trackedRCP1.fits")
RCP4=fitsio_read_image("/SUM16/D156498376/S00000/trackedRCP4.fits")
velocityRCP41 = float(RCP4-RCP1)/float(RCP4+RCP1)
RCP1=0.0
RCP4=0.0
LCP2=fitsio_read_image("/SUM12/D156481812/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM0/D156483395/S00000/trackedLCP3.fits")
velocityLCP32 = float(LCP3-LCP2)/float(LCP3+LCP2)
LCP2=0.0
LCP3=0.0  
RCP2=fitsio_read_image("/SUM15/D156501572/S00000/trackedRCP2.fits")
RCP3=fitsio_read_image("/SUM3/D156500247/S00000/trackedRCP3.fits")
velocityRCP32 = float(RCP3-RCP2)/float(RCP3+RCP2)
RCP2=0.0
RCP3=0.0

OPENR,1,'velocitiesQUIETSUN' ; file containing OBS_VR
vel=FLTARR(400)
READF,1,vel
CLOSE,1
nt=400

; THE INVERSE LOOK-UP TABLES
cal50=TOTAL(TOTAL(VELOCITYLCP50[206:306,206:306,*]+VELOCITYRCP50[206:306,206:306,*],1),1)/101./101./2.
cal41=TOTAL(TOTAL(VELOCITYLCP41[206:306,206:306,*]+VELOCITYRCP41[206:306,206:306,*],1),1)/101./101./2.
cal32=TOTAL(TOTAL(VELOCITYLCP32[206:306,206:306,*]+VELOCITYRCP32[206:306,206:306,*],1),1)/101./101./2.

ENDIF ELSE BEGIN

; QUIET SUN 2
;---------------


LCP0=fitsio_read_image("/SUM9/D156794960/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM17/D156810196/S00000/trackedLCP5.fits")
velocityLCP50 = float(LCP5-LCP0)/float(LCP5+LCP0)
LCP0=0.0
LCP5=0.0
RCP5=fitsio_read_image("/SUM18/D156827861/S00000/trackedRCP5.fits")
RCP0=fitsio_read_image("/SUM9/D156812923/S00000/trackedRCP0.fits")
velocityRCP50 = float(RCP5-RCP0)/float(RCP5+RCP0)
RCP0=0.0
RCP5=0.0
LCP1=fitsio_read_image("/SUM4/D156798543/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM13/D156806880/S00000/trackedLCP4.fits")
velocityLCP41 = float(LCP4-LCP1)/float(LCP4+LCP1)
LCP1=0.0
LCP4=0.0
RCP1=fitsio_read_image("/SUM22/D156815745/S00000/trackedRCP1.fits")
RCP4=fitsio_read_image("/SUM13/D156823785/S00000/trackedRCP4.fits")
velocityRCP41 = float(RCP4-RCP1)/float(RCP4+RCP1)
RCP1=0.0
RCP4=0.0
LCP2=fitsio_read_image("/SUM5/D156800752/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM19/D156803850/S00000/trackedLCP3.fits")
velocityLCP32 = float(LCP3-LCP2)/float(LCP3+LCP2)
LCP2=0.0
LCP3=0.0  
RCP2=fitsio_read_image("/SUM19/D156818409/S00000/trackedRCP2.fits")
RCP3=fitsio_read_image("/SUM8/D156820489/S00000/trackedRCP3.fits")
velocityRCP32 = float(RCP3-RCP2)/float(RCP3+RCP2)
RCP2=0.0
RCP3=0.0

OPENR,1,'velocitiesQUIETSUN2' ; file containing OBS_VR
vel=FLTARR(961)
READF,1,vel
CLOSE,1
nt=961

; THE INVERSE LOOK-UP TABLES
cal50=TOTAL(TOTAL(VELOCITYLCP50+VELOCITYRCP50,1),1)/128./128./2.
cal41=TOTAL(TOTAL(VELOCITYLCP41+VELOCITYRCP41,1),1)/128./128./2.
cal32=TOTAL(TOTAL(VELOCITYLCP32+VELOCITYRCP32,1),1)/128./128./2.

ENDELSE


cal502=SMOOTH(cal50,15,/EDGE_TRUNCATE)
cal412=SMOOTH(cal41,15,/EDGE_TRUNCATE)
cal322=SMOOTH(cal32,15,/EDGE_TRUNCATE)
res50=poly_fit(cal502,vel,3,yfit=y)
res41=poly_fit(cal412,vel,3,yfit=y)
res32=poly_fit(cal322,vel,3,yfit=y)
;vin=DINDGEN(1000)/999.d0*6000.d0-3000.d0
;cal50i=res50[0]+res50[1]*vin+res50[2]*vin*vin+res50[3]*vin*vin*vin
;cal41i=res41[0]+res41[1]*vin+res41[2]*vin*vin+res41[3]*vin*vin*vin
;cal32i=res32[0]+res32[1]*vin+res32[2]*vin*vin+res32[3]*vin*vin*vin

RESTORE,'lookup_center0.bin' ; with calibration11
ratio0=ratio
RESTORE,'lookup_center.bin'
ratio1=ratio
RESTORE,'temp.bin'
!P.MULTI=[0,1,3]
plot,vel,cal50,xst=1,xrange=[-3500,3000],yrange=[-0.025,0.025],yst=1
oplot,vtest,ratio1[2,*]+(cal50[521]-ratio1[2,410]),col=180
oplot,vtest,ratio0[2,*]+(cal50[521]-ratio0[2,410]),linestyle=2,thick=2
oplot,vtest,ratio[2,*]+(cal50[521]-ratio[2,410]),linestyle=2,col=180
plot,vel,cal41,xst=1,xrange=[-3500,3000],yrange=[-0.15,0.15],yst=1
oplot,vtest,ratio1[1,*]+(cal41[521]-ratio1[1,410]),col=180
oplot,vtest,ratio0[1,*]+(cal41[521]-ratio0[1,410]),linestyle=2,thick=2
oplot,vtest,ratio[1,*]+(cal41[521]-ratio[1,410]),linestyle=2,col=180
plot,vel,cal32,xst=1,xrange=[-3500,3000],yrange=[-0.25,0.25],yst=1
oplot,vtest,ratio1[0,*]+(cal32[521]-ratio1[0,410]),col=180
oplot,vtest,ratio0[0,*]+(cal32[521]-ratio0[0,410]),linestyle=2,thick=2
oplot,vtest,ratio[0,*]+(cal32[521]-ratio[0,410]),linestyle=2,col=180
$mv lookup_center.bin temp.bin


READ,PAUSE

END
