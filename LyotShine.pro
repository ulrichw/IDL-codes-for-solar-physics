; this programs fits for the average phase and contrast of the Lyot
; non tunable elements using  LyotProfiles2006aug3.fits


PRO MyLyot,X,A,F,pder

A[0:3] = ABS(A[0:3])

nel = N_ELEMENTS(X)
F   = FLTARR(nel)

RESTORE,'temp.bin' ; contains wavelengths, blocker filter profile, and Iripple

;blocker0 = 0.65930921-0.0029807898*l-0.018651693*l^2.0
blocker0 = A[14] + l*A[15] + l^2.0*A[16]
bb=WHERE(blocker0 LT 0.0)
IF(bb[0] NE -1) THEN blocker0[bb]=0.0

dpi     = 2.d0*!dpi
;angles  = FINDGEN(17)*22.5 ; angle of the tuning polarizer in degrees
;angles  = angles*!dpi/180.d0

FSR     = FLTARR(4) ; in A
FSR[0]  = 1.405d0               ; for E2
FSR[1]  = 2.779d0               ; for E3
FSR[2]  = 5.682d0               ; for E4
FSR[3]  = 11.354d0              ; for E5

F       = A[0]*blocker0*(1.d0+A[1]*COS(A[9]*dpi*l/FSR[0]+A[4]))*(1.d0+A[2]*COS(A[10]*dpi*l/FSR[1]+A[5]))*(1.d0+A[3]*COS(A[11]*dpi*l/FSR[2]+A[6]))*(1.d0+A[13]*COS(A[12]*dpi*l/FSR[3]+A[7]))/32.d0+A[8] ; E2 to E5

END



PRO LyotShine

d          = READFITS('LyotProfiles2006aug3.fits')
data       = FLTARR(1280)
FOR i=0,16 DO data=data+d[*,i+1]
data       = data/MAX(data) 
;data      = data - MEAN(data[0:250]) ; dark removed
data       = REVERSE(data)

;a          = WHERE(data LT 0.d0)
;IF(a[0] NE -1) THEN data[a]       = 0.d0

;Iripple    = FLTARR(17)
;FOR i=0,16 DO Iripple[i] = TOTAL(d[*,i+1])
;Iripple    = Iripple/MAX(Iripple)

l          = FINDGEN(1280)/1279.*(6177.24d0-6168.24d0)+6168.24d0-6173.3433 ; reference wavelength
q          = READFITS('blocker11.fits')
blocker0   = INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,l)
SAVE,l,blocker0,FILE='temp.bin'

;wavelength = FLTARR(1280*17)
;FOR i=0,16 DO wavelength[i*1280:i*1280+1279] = l
;data       = FLTARR(1280*17)
;FOR i=0,16 DO data[i*1280:i*1280+1279] = REFORM(d[*,i+1])
nel        = N_ELEMENTS(data)
weight     = FLTARR(nel)+1.d0
d          = 0.0

A          = FLTARR(17)
A[0]       = MAX(data)/0.65
A[1:3]     = 1.0
A[8]       = 0.04
A[9]       = 1.0d0 ; FSRs
A[10]      = 1.0d0
A[11]      = 1.d0
A[12]      = 1.d0 ; element E2
A[13]      = 1.d0
A[14]      = 0.65930921
A[15]      = -0.0029807898
A[16]      = -0.018651693


resp = CURVEFIT(FINDGEN(nel),data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2,/NODERIVATIVE)

print,chi2,A
WINDOW,0,RETAIN=2
plot,l,data,xst=1,yrange=[0,1.05],yst=1
oplot,l,resp,color=180
oplot,l,A[14] + l*A[15] + l^2.0*A[16]

READ,pause

END
