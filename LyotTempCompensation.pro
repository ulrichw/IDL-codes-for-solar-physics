; calculate the theoretical wavelength drift for the HMI Lyot elements
; with temperature.

PRO LyotTempCompensation

; E1

d1=30.620d0 ; in mm
d2=7.056d0
n1=-0.170969d0
n2=-0.0449016d0
l0=6173.3434d-7 ; in mm
alpha1=3.25919d-6*1.d7
alpha2=1.26788d-6*1.d7
print,l0^2.d0/( (d1*n1-d2*n2)-l0*(d1*alpha1-d2*alpha2))*1.d7 ; FSR OF E1


lam=FINDGEN(10000)/9999.-0.5
lam=lam*1.d-7
lam=lam+6173.3434d-7

dt1=9.98d-6
dt2=5.23d-5
dta1=5.68d-6
dta2=4.2d-6


retard=2.d0*!dpi*(d1*n1-d2*n2)/lam
dretard=2.d0*!dpi*(d1*dt1+dta1*d1*n1-d2*dt2-dta2*n2*d2)/lam



plot,lam*1.d7,cos(retard/2.d0)^2.d0,xst=1
oplot,lam*1.d7,cos((retard+dretard)/2.d0)^2.d0,col=180

RC=1.d0/((dt1+dta1*n1)/(dt2+dta2*n2))
FSR=0.688d-7
d0=l0^2.d0/FSR/ABS(n1-l0*alpha1) ; in mm
d1b=d0/(1.d0-n2/(n1*RC))
d2b=d1b/RC
PRINT,d0,d1b,d2b

READ,pause


END
