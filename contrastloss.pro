; program to answer Jesper's question: due to the rays being nnot in
; the Michelsons, how are the phases and contrasts affected ?


n1=1.51682698726654d0
n2=1.d0
; FOR WB MICHELSON
;d1=6.310d7
;d2=4.160d7
; FOR NB MICHELSON
d1=12.620d7
d2=8.330d7
; from Zemax measurements at CF=9 for NB Michelson (PROVIDED BY JESPER
; EMAIL OF JULY 29, 2011)
extangle =3.d0*(2.599d0+0.043d0)
extangle2=3.d0*(1.321d0+1.321d0)
tilt1=1.278d0
tilt2=0.d0*1.278d0

nlam=60000
ntheta=500
T=FLTARR(nlam,ntheta) ; transmittances

lam=DINDGEN(nlam)/(nlam-1.d0)*0.25d0-0.12d0
lam=lam+6173.3433d0

theta=(DINDGEN(ntheta)/(ntheta-1.d0)*extangle-extangle/2.d0)*!dpi/180.d0
theta=theta+tilt1*!dpi/180.d0
print,min(theta)*180./!dpi,max(theta)*180./!dpi,extangle
FOR i=0,ntheta-1 DO BEGIN

    d=2.d0*(n1*d1-n2*d2)-(d1/n1^3.d0-d2/n2^3.d0)*(sin(theta[i]))^4.d0/4.d0-(d1/n1^5.d0-d2/n2^5.d0)*(sin(theta[i]))^6.d0/8.d0
    T[*,i]=0.5d0*(1.d0+cos(2.d0*!dpi*d/lam))

ENDFOR
avT=TOTAL(T^2.d0,2)/TOTAL(T,2) ; careful, it's a cone so integral= theta*2 !pi r dr
PLOT,lam,avT,xst=1;,xrange=[6173.15,6173.25],yst=1
;PRINT,(MAX(T[*,0])-MIN(T[*,0]))
PRINT,(MAX(avT)-MIN(avT))
a=WHERE(avT EQ MAX(avT))


theta=(DINDGEN(ntheta)/(ntheta-1.d0)*extangle2-extangle2/2.d0)*!dpi/180.d0
theta=theta+tilt2*!dpi/180.d0
print,min(theta)*180./!dpi,max(theta)*180./!dpi,extangle2
FOR i=0,ntheta-1 DO BEGIN

    d=2.d0*(n1*d1-n2*d2)-(d1/n1^3.d0-d2/n2^3.d0)*sin(theta[i])^4.d0/4.d0-(d1/n1^5.d0-d2/n2^5.d0)*sin(theta[i])^6.d0/8.d0
    T[*,i]=0.5d0*(1.d0+cos(2.d0*!dpi*d/lam))

ENDFOR
avT=TOTAL(T^2.d0,2)/TOTAL(T,2)
OPLOT,lam,avT,col=180,linestyle=2
;PRINT,(MAX(T[*,0])-MIN(T[*,0]))
PRINT,(MAX(avT)-MIN(avT))
b=WHERE(avT EQ MAX(avT))

PRINT,"wavelength difference:",lam[a[0]]-lam[b[0]]
PRINT,"phase difference:",(lam[a[0]]-lam[b[0]])*360.d0/0.172d0 ; for NB
PRINT,"cones",extangle,extangle2

END
