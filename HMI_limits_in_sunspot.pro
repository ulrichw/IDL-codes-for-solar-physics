; NOAA 11092
;LCP5=fitsio_read_image("/SUM15/D147734385/S00000/trackedLCP5.fits")
;LCP0=fitsio_read_image("/SUM18/D147675980/S00000/trackedLCP0.fits")
;RCP5=fitsio_read_image("/SUM21/D147829266/S00000/trackedRCP5.fits")
;RCP0=fitsio_read_image("/SUM2/D147745962/S00000/trackedRCP0.fits")
;LCP1=fitsio_read_image("/SUM9/D147683069/S00000/trackedLCP1.fits")
;LCP4=fitsio_read_image("/SUM16/D147725403/S00000/trackedLCP4.fits")
;RCP1=fitsio_read_image("/SUM14/D147761508/S00000/trackedRCP1.fits")
;RCP4=fitsio_read_image("/SUM8/D147821873/S00000/trackedRCP4.fits")
;LCP2=fitsio_read_image("/SUM13/D147695696/S00000/trackedLCP2.fits")
;LCP3=fitsio_read_image("/SUM13/D147715300/S00000/trackedLCP3.fits")
;RCP2=fitsio_read_image("/SUM16/D147801564/S00000/trackedRCP2.fits")
;RCP3=fitsio_read_image("/SUM10/D147812164/S00000/trackedRCP3.fits")

; NOAA 11161
LCP0=fitsio_read_image("/SUM20/D150904920/S00000/trackedLCP0.fits")
LCP5=fitsio_read_image("/SUM9/D150934556/S00000/trackedLCP5.fits")
velocityLCP50 = float(LCP5-LCP0)/float(LCP5+LCP0)
LCP0=0.0
LCP5=0.0
RCP5=fitsio_read_image("/SUM15/D150965082/S00000/trackedRCP5.fits")
RCP0=fitsio_read_image("/SUM2/D150944540/S00000/trackedRCP0.fits")
velocityRCP50 = float(RCP5-RCP0)/float(RCP5+RCP0)
RCP0=0.0
RCP5=0.0
LCP1=fitsio_read_image("/SUM20/D150907112/S00000/trackedLCP1.fits")
LCP4=fitsio_read_image("/SUM15/D150934493/S00000/trackedLCP4.fits")
velocityLCP41 = float(LCP4-LCP1)/float(LCP4+LCP1)
LCP1=0.0
LCP4=0.0
RCP1=fitsio_read_image("/SUM17/D150944570/S00000/trackedRCP1.fits")
RCP4=fitsio_read_image("/SUM8/D150964980/S00000/trackedRCP4.fits")
velocityRCP41 = float(RCP4-RCP1)/float(RCP4+RCP1)
RCP1=0.0
RCP4=0.0
LCP2=fitsio_read_image("/SUM22/D150925772/S00000/trackedLCP2.fits")
LCP3=fitsio_read_image("/SUM1/D150925874/S00000/trackedLCP3.fits")
velocityLCP32 = float(LCP3-LCP2)/float(LCP3+LCP2)
LCP2=0.0
LCP3=0.0  
RCP2=fitsio_read_image("/SUM19/D150955096/S00000/trackedRCP2.fits")
RCP3=fitsio_read_image("/SUM12/D150955124/S00000/trackedRCP3.fits")
velocityRCP32 = float(RCP3-RCP2)/float(RCP3+RCP2)
RCP2=0.0
RCP3=0.0

;OPENR,1,'velocitiesNOAA11092' ; file containing OBS_VR
;FOR NOAA 11092 Of Februay 19, 2011
OPENR,1,'velocitiesNOAA11161' ; file containing OBS_VR
;vel=FLTARR(1121)
vel=FLTARR(961)
READF,1,vel
CLOSE,1
; NB: vel[647] corresponds to zero velocity
RESTORE,'lookup.bin'
magnetic = 1.d0/(2.d0*4.67d-5*0.000061733433d0*2.5d0*299792458.d0)
RESTORE,'calibration32_NOAA11161.bin' ;CALIBRATION AT [0:50,0:50,*]
RESTORE,'calibration41_NOAA11161.bin'
RESTORE,'calibration50_NOAA11161.bin'

;---------------------------------------------------------------------------------------------------------------------

; IN QUIET SUN AT DIFFERENT LOCATIONS
velocity50Q0=TOTAL(TOTAL((velocityLCP50[0:50,0:50,*]+velocityRCP50[0:50,0:50,*])/2.,1),1)/51./51.
velocity41Q0=TOTAL(TOTAL((velocityLCP41[0:50,0:50,*]+velocityRCP41[0:50,0:50,*])/2.,1),1)/51./51.
velocity32Q0=TOTAL(TOTAL((velocityLCP32[0:50,0:50,*]+velocityRCP32[0:50,0:50,*])/2.,1),1)/51./51.
velocity50Q1=TOTAL(TOTAL((velocityLCP50[450:500,450:500,*]+velocityRCP50[450:500,450:500,*])/2.,1),1)/51./51.
velocity41Q1=TOTAL(TOTAL((velocityLCP41[450:500,450:500,*]+velocityRCP41[450:500,450:500,*])/2.,1),1)/51./51.
velocity32Q1=TOTAL(TOTAL((velocityLCP32[450:500,450:500,*]+velocityRCP32[450:500,450:500,*])/2.,1),1)/51./51.
LCP50Q0=TOTAL(velocityLCP50[*,0:20,*],2)/21.
RCP50Q0=TOTAL(velocityRCP50[*,0:20,*],2)/21.
LCP41Q0=TOTAL(velocityLCP41[*,0:20,*],2)/21.
RCP41Q0=TOTAL(velocityRCP41[*,0:20,*],2)/21.
LCP32Q0=TOTAL(velocityLCP32[*,0:20,*],2)/21.
RCP32Q0=TOTAL(velocityRCP32[*,0:20,*],2)/21.
LCP50Q1=TOTAL(velocityLCP50[*,185:205,*],2)/21.
RCP50Q1=TOTAL(velocityRCP50[*,185:205,*],2)/21.
LCP41Q1=TOTAL(velocityLCP41[*,185:205,*],2)/21.
RCP41Q1=TOTAL(velocityRCP41[*,185:205,*],2)/21.
LCP32Q1=TOTAL(velocityLCP32[*,185:205,*],2)/21.
RCP32Q1=TOTAL(velocityRCP32[*,185:205,*],2)/21.

velocity50Qsunspot=TOTAL(TOTAL((velocityLCP50[250:260,190:205,*]+velocityRCP50[250:260,190:205,*])/2.,1),1)/11./16.
velocity41Qsunspot=TOTAL(TOTAL((velocityLCP41[250:260,190:205,*]+velocityRCP41[250:260,190:205,*])/2.,1),1)/11./16.
velocity32Qsunspot=TOTAL(TOTAL((velocityLCP32[250:260,190:205,*]+velocityRCP32[250:260,190:205,*])/2.,1),1)/11./16.

!P.multi=[0,1,3]
PLOT,vel,SMOOTH(LCP50Q0[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP50Q0[510,*],10,/edge_truncate),col=180
OPLOT,vel,SMOOTH(LCP50Q1[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP50Q1[510,*],10,/edge_truncate),col=180
PLOT,vel,SMOOTH(LCP41Q0[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP41Q0[510,*],10,/edge_truncate),col=180
OPLOT,vel,SMOOTH(LCP41Q1[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP41Q1[510,*],10,/edge_truncate),col=180
PLOT,vel,SMOOTH(LCP32Q0[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP32Q0[510,*],10,/edge_truncate),col=180
OPLOT,vel,SMOOTH(LCP32Q1[0,*],10,/edge_truncate)
OPLOT,vel,SMOOTH(LCP32Q1[510,*],10,/edge_truncate),col=180


; VARIATION OF LOOK-UP TABLES WITH X (LONGITUDE):
; THE LOOK-UP TABLE IS SHIFTED IN VTEST, BUT THE SLOPE ALSO VARIES
!P.multi=[0,1,3]
plot,vel,TOTAL(LCP32Q0[0:15,*],1)/16.
oplot,vel,TOTAL(LCP32Q0[250:265,*],1)/16.,col=180
oplot,vel,TOTAL(LCP32Q0[496:511,*],1)/16.,col=18000
plot,vel,TOTAL(LCP41Q0[0:15,*],1)/16.
oplot,vel,TOTAL(LCP41Q0[250:265,*],1)/16.,col=180
oplot,vel,TOTAL(LCP41Q0[496:511,*],1)/16.,col=18000
plot,vel,TOTAL(LCP50Q0[0:15,*],1)/16.
oplot,vel,TOTAL(LCP50Q0[250:265,*],1)/16.,col=180
oplot,vel,TOTAL(LCP50Q0[496:511,*],1)/16.,col=18000

; ROUGHLY, THE SHIFT IN VTEST WITH X IS:
x=[7.5,37.5,57.5,77.5,107.5,137.5,157.5,177.5,207.5,227.5,257.5,277.5,307.5,357.5,407.5,457.5,477.5,503.5]
vshift=[0,-20,-55,-90,-130,-205,-225,-310,-320,-400,-490,-500,-530,-700,-720,-770,-790,-800]
;result of a third-order polynomial fit based on cal32:
res=[12.0454,-0.823745,-0.00644013,9.73507e-06]

!p.multi=[0,1,3]
PLOT,vel,SMOOTH(TOTAL(LCP41Q0[0:15,*],1)/16.,10,/edge_truncate),charsize=1.5,xrange=[-6000,6000],xst=1,yrange=[-0.25,0.25],yst=1
FOR i=0,495 DO BEGIN & dv=res[0]+res[1]*i+res[2]*i*i+res[3]*i*i*i & OPLOT,vel-dv,SMOOTH(TOTAL(LCP41Q0[i:i+15,*],1)/16.,10,/edge_truncate)
OPLOT,vtest,cal41*1.03,col=180,thick=2
PLOT,vel,SMOOTH(TOTAL(LCP50Q0[0:15,*],1)/16.,10,/edge_truncate),charsize=1.5,xrange=[-6000,6000],xst=1,yrange=[-0.25,0.25],yst=1
FOR i=0,495 DO BEGIN & dv=res[0]+res[1]*i+res[2]*i*i+res[3]*i*i*i & OPLOT,vel-dv,SMOOTH(TOTAL(LCP50Q0[i:i+15,*],1)/16.,10,/edge_truncate)
OPLOT,vtest,cal50*1.03,col=180,thick=2
PLOT,vel,SMOOTH(TOTAL(LCP32Q0[0:15,*],1)/16.,10,/edge_truncate),charsize=1.5,xrange=[-6000,6000],xst=1,yrange=[-0.25,0.25],yst=1
FOR i=0,495 DO BEGIN & dv=res[0]+res[1]*i+res[2]*i*i+res[3]*i*i*i & OPLOT,vel-dv,SMOOTH(TOTAL(LCP32Q0[i:i+15,*],1)/16.,10,/edge_truncate)
OPLOT,vtest,cal32,col=180,thick=2

;CALIBRATION AT [250:260,190:205,*] COMPARED TO CALIBRATION AT [0:50,0:50,*]


;AT ZERO VELOCITY
LCP50vel0=INTERPOL(vtest,cal50,velocityLCP50[*,*,647])
RCP50vel0=INTERPOL(vtest,cal50,velocityRCP50[*,*,647])
LCP41vel0=INTERPOL(vtest,cal41,velocityLCP41[*,*,647])
RCP41vel0=INTERPOL(vtest,cal41,velocityRCP41[*,*,647])
LCP32vel0=INTERPOL(vtest,cal32,velocityLCP32[*,*,647])
RCP32vel0=INTERPOL(vtest,cal32,velocityRCP32[*,*,647])






END
