PRO limbdarkening

coeffs=[0.4430d0, 0.1390d0, 0.0410d0, 0.0125d0, 0.0019d0]

x=DINDGEN(10000)/9999.d0

costheta=SQRT(1.d0-x^2.d0)
mu=ALOG(costheta)

z=1.d0
ld=1.d0
FOR i=0,4 DO BEGIN
    z=z*mu
    ld=ld+coeffs[i]*z ; limb darkening
ENDFOR

PLOT,x,ld,xst=1,charsize=1.5,tit='!17',xtit='fractional solar radius',ytit='Relative intensity',thick=2



READ,pause

END
