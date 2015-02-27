; This program checks how the transmission profile
; of the Michelsons is affected by air pressure

PRO HMI_air

openr,1,'AirPressure.txt' ; contains the refraction index of air
data=dblarr(3,23)
readf,1,data
close,1

d2nb=8.330 ;NB Michelson vacuum leg thickness
d1nb=12.620
d2wb=4.160
d1wb=6.310
n1=1.51683

nlam = 12000
dlam = 1.d0/1.5d3       ; resolution in wavelength
lam  = (DINDGEN(nlam)-(nlam-1)/2.)*dlam+6173.d0

; in vacuum
Deltavac = cos(!dpi*2.d0*(n1*d1nb-d2nb)/lam*1.d7)^2.d0
;Deltavac = cos(!dpi*2.d0*(n1*d1nb-d2nb*1.000271576)/lam*1.d7)^2.d0 ; NB Michelson


SET_PLOT,'ps'
DEVICE,file='yo.ps',xsize=21,ysize=26,xoffset=0,yoffset=0,/color,bits=24
LOADCT,4
!p.multi=[0,3,7]
FOR i=2,22 DO BEGIN
; in air
Deltaair1 = cos(!dpi*2.d0*(n1*d1nb-d2nb*data[1,i])/lam*1.d7)^2.d0
Deltaair2 = cos(!dpi*2.d0*(n1*d1nb-d2nb*data[2,i])/lam*1.d7)^2.d0
;Deltaair = cos(!dpi*2.d0*(n1*d1nb-d2nb*1.000264086)/lam*1.d7)^2.d0 ; NB Michelson


;plot,lam,Deltavac,xrange=[6172.7,6173.3],xst=1,tit=STRING(Data[0,i])
plot ,lam,Deltavac  ,xst=1,tit=STRING(Data[0,i]),thick=2,xrange=[6172.9,6173.1]
oplot,lam,Deltaair1,linestyle=2,color=100,thick=2
oplot,lam,Deltaair2,linestyle=2,color=200,thick=2
ENDFOR

plot,lam,Deltavac,xrange=[6172.9,6173.1],xst=1,yrange=[0.5,1],yst=1

DEVICE,/CLOSE
set_plot,'x'
!p.multi=0


END
