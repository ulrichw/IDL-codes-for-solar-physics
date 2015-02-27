; paper on HMI calibration 

; temperatures:

TE1=[34.1,29.62,27.65]
TNB=[33.83,29.31,27.19]
TWB=[34.1,29.62,27.65]

TE12=FINDGEN(1000)/999.*(34.1-27.65)+27.65
TNB2=FINDGEN(1000)/999.*(33.83-27.19)+27.19
TWB2=FINDGEN(1000)/999.*(34.1-27.65)+27.65

; FSR

FSRNB      = 0.1710098d0;0.172d0-0.001057600d0 ; for the narrow-band Michelson
FSRWB      = 0.3421506d0;0.344d0-0.002076830d0 ; for the broad-band  Michelson
FSRE1      = 0.6943613d0;0.693d0+0.000483467d0 ; for E1

; phases

front_NB_OBS=-([-116.47,-79.24,-188.08]+79.24)/360.d0*FSRNB*1000.
front_WB_OBS=-([-3.79,29.20,-16.34]-29.20)/360.d0*FSRWB*1000.
front_E1_OBS=-([33.22,35.68,12.92]-35.68)/360.d0*FSRE1*1000.

side_NB_OBS=-([-118.11,-80.55,-188.85]+80.55)/360.d0*FSRNB*1000.
side_WB_OBS=-([-4.65,28.89,-16.55]-28.89)/360.d0*FSRWB*1000.
side_E1_OBS=-([33.09,35.24,12.93]-35.24)/360.d0*FSRE1*1000.

front_NB_CAL=-([-118.65,-81.91,-191.35]+81.91)/360.d0*FSRNB*1000.
front_WB_CAL=-([-5.03,29.29,-15.58]-29.29)/360.d0*FSRWB*1000.
front_E1_CAL=-([35.55,40.44,17.89]-40.44)/360.d0*FSRE1*1000.

side_NB_CAL=-([-118.39,-82.73,-191.82]+82.73)/360.d0*FSRNB*1000.
side_WB_CAL=-([-5.39,30.73,-16.02]-30.73)/360.d0*FSRWB*1000.
side_E1_CAL=-([36.50,40.85,17.83]-40.85)/360.d0*FSRE1*1000.

SET_PLOT,'PS'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color
LOADCT,3
!P.MULTI=[0,1,3]

PLOT,TNB,front_NB_OBS,psym=4,yrange=[-20.,60.],yst=1,xrange=[27.1,34.2],xst=1,charsize=1.5,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength drift (mA)'
res=POLY_FIT(TNB,front_NB_OBS,2)
OPLOT,TNB2,POLY(TNB2,res)
PRINT,res[1]+2.d0*res[2]*30.d0
OPLOT,TNB,side_NB_OBS,psym=4
res=POLY_FIT(TNB,side_NB_OBS,2)
OPLOT,TNB2,POLY(TNB2,res),linestyle=2
OPLOT,TNB,front_NB_CAL,psym=4,col=180
res=POLY_FIT(TNB,front_NB_CAL,2)
OPLOT,TNB2,POLY(TNB2,res),col=180
OPLOT,TNB,side_NB_CAL,psym=4,col=180
res=POLY_FIT(TNB,side_NB_CAL,2)
OPLOT,TNB2,POLY(TNB2,res),col=180,linestyle=2


PLOT,TWB,front_WB_OBS,psym=4,yrange=[-20.,60.],yst=1,xrange=[27.1,34.2],xst=1,charsize=1.5,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength drift (mA)'
res=POLY_FIT(TWB,front_WB_OBS,2)
PRINT,res[1]+2.d0*res[2]*30.d0
OPLOT,TWB2,POLY(TWB2,res)
OPLOT,TWB,side_WB_OBS,psym=4
res=POLY_FIT(TWB,side_WB_OBS,2)
OPLOT,TWB2,POLY(TWB2,res),linestyle=2
OPLOT,TWB,front_WB_CAL,psym=4,col=180
res=POLY_FIT(TWB,front_WB_CAL,2)
OPLOT,TWB2,POLY(TWB2,res),col=180
OPLOT,TWB,side_WB_CAL,psym=4,col=180
res=POLY_FIT(TWB,side_WB_CAL,2)
OPLOT,TWB2,POLY(TWB2,res),col=180,linestyle=2

PLOT,TE1,front_E1_OBS,psym=4,yrange=[-20.,60.],yst=1,xrange=[27.1,34.2],xst=1,charsize=1.5,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength drift (mA)'
res=POLY_FIT(TE1,front_E1_OBS,2)
PRINT,res[1]+2.d0*res[2]*30.d0
OPLOT,TE12,POLY(TE12,res)
OPLOT,TE1,side_E1_OBS,psym=4
res=POLY_FIT(TE1,side_E1_OBS,2)
OPLOT,TE12,POLY(TE12,res),linestyle=2
OPLOT,TE1,front_E1_CAL,psym=4,col=180
res=POLY_FIT(TE1,front_E1_CAL,2)
OPLOT,TE12,POLY(TE12,res),col=180
OPLOT,TE1,side_E1_CAL,psym=4,col=180
res=POLY_FIT(TE1,side_E1_CAL,2)
OPLOT,TE12,POLY(TE12,res),col=180,linestyle=2

DEVICE,/CLOSE
!P.MULTI=0

END
