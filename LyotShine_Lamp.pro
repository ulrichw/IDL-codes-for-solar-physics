; this programs fits for the average phase and contrast of the Lyot
; non tunable elements using  LyotProfiles2006aug3.fits


PRO MyLyot,X,A,F,pder

A[0] = ABS(A[0])
A[1] = ABS(A[1])
A[3] = ABS(A[3])
A[5] = ABS(A[5])
A[7] = ABS(A[7])
A[10]= ABS(A[10])
;A[12]= ABS(A[12])
;A[13]= ABS(A[13])

nel = N_ELEMENTS(X)
F   = FLTARR(nel)

RESTORE,'temp.bin' ; contains wavelengths, blocker filter profile, and Iripple

dpi   = 2.d0*!dpi

FSR   = FLTARR(4) ; in A
FSR[0]= 1.405d0;A[13]; 1.405d0               ; for E2
FSR[1]= 2.779d0;A[12] ;2.779d0               ; for E3
FSR[2]= 5.682d0               ; for E4
FSR[3]= 11.354d0              ; for E5

E1    = (1.d0+A[10] *COS(dpi*l/0.7039+A[11]))/2.d0  ; I try to fit for E1 if light from the source is partly polarized
;E1    = (1.d0+A[15] *COS(dpi*l/0.690))/2.d0
;E1    = 1.0 ; check or uncheck if you fit for E1
E2    = (1.d0+A[1] *COS(dpi*l/FSR[0]+A[2]))/2.d0
E3    = (1.d0+A[3] *COS(dpi*l/FSR[1]+A[4]))/2.d0
E4    = (1.d0+A[5] *COS(dpi*l/FSR[2]+A[6]))/2.d0
E5    = (1.d0+A[7] *COS(dpi*l/FSR[3]+A[8]))/2.d0

F     = A[0]*blocker0*E1*E2*E3*E4*E5+A[9]

dE1d1 = COS(dpi*l/0.7039+A[11])/2.d0
dE1d2 = -A[10]*SIN(dpi*l/0.7039+A[11])/2.d0
;dE1d1 = COS(dpi*l/0.690)/2.d0
dE2d1 = COS(dpi*l/FSR[0]+A[2])/2.d0
dE2d2 = -A[1]*SIN(dpi*l/FSR[0]+A[2])/2.d0
dE3d1 = COS(dpi*l/FSR[1]+A[4])/2.d0
dE3d2 = -A[3]*SIN(dpi*l/FSR[1]+A[4])/2.d0
dE4d1 = COS(dpi*l/FSR[2]+A[6])/2.d0
dE4d2 = -A[5]*SIN(dpi*l/FSR[2]+A[6])/2.d0
dE5d1 = COS(dpi*l/FSR[3]+A[8])/2.d0
dE5d2 = -A[7]*SIN(dpi*l/FSR[3]+A[8])/2.d0

pder  = [ [blocker0*E1*E2*E3*E4*E5] , [dE2d1*E1*E3*E4*E5*A[0]*blocker0] , [dE2d2*E1*E3*E4*E5*A[0]*blocker0] , [dE3d1*E1*E2*E4*E5*A[0]*blocker0] , [dE3d2*E1*E2*E4*E5*A[0]*blocker0] , [dE4d1*E1*E2*E3*E5*A[0]*blocker0] , [dE4d2*E1*E2*E3*E5*A[0]*blocker0] , [dE5d1*E1*E2*E3*E4*A[0]*blocker0] , [dE5d2*E1*E2*E3*E4*A[0]*blocker0]  , [FLTARR(nel)+1.d0] , [dE1d1*E2*E3*E4*E5*A[0]*blocker0] , [dE1d2*E2*E3*E4*E5*A[0]*blocker0] ];, [A[0]*blocker0*dpi*l*A[3]/FSR[1]^2.0*SIN(dpi*l/FSR[1]+A[4])/2.d0*E1*E2*E4*E5] , [A[0]*blocker0*dpi*l*A[1]/FSR[0]^2.0*SIN(dpi*l/FSR[0]+A[2])/2.d0*E1*E3*E4*E5] ]

; WE CONVOLVE BY THE "FINITE RESOLUTION WINDOW"
xx = SHIFT(DINDGEN(nel)-nel/2.d0,nel/2.d0)
xx[nel/2.d0]=-xx[nel/2.d0]
param=0.00175
kernel=1./(1.+param*ABS(xx)+param*xx^2.+param*ABS(xx)^3.)
fkernel = FFT(kernel,/DOUBLE)

F = FLOAT(FFT(FFT(F,/DOUBLE)*fkernel,/INVERSE,/DOUBLE))
FOR i=0,N_ELEMENTS(A)-1 DO pder[*,i] = FLOAT(FFT(FFT(pder[*,i],/DOUBLE)*fkernel,/INVERSE,/DOUBLE))

F = F*nel
pder = pder*nel


END



PRO LyotShine_Lamp
; we correct for the potential change in intensity
; we correct for the curvature of the spectrum due to the spectrograph
; we correct for the blocker filter shift IF BLOCKER PRESENT
; we correct for finite resolution

!P.MULTI = 0

filename = 'lyotlampsetBprofiles2006oct16.fits'

d        = READFITS(filename)
blocker  = 1

; WE GET RID OF THE LINEAR INCREASE IN THE BACKGROUND
FOR i=0,35 DO d[*,i] = d[*,i] - i*362. ; linear fit

; WE GET RID OF THE I-RIPPLE
inten    = FLTARR(36) 
FOR i=0,35 DO BEGIN
    inten[i]= TOTAL(d[*,i])
    d[*,i]  = d[*,i]/inten[i]
ENDFOR
d        = REFORM(REBIN(d,1320,1))

nelem    = N_ELEMENTS(d)
data     = d
data     = data/MAX(data) 

data     = REVERSE(data)-0.041 ; TO GET READ OF BACKGROUND DUE TO IR COMPONENT IN LAMP (?) 

l        = READFITS('ref2006oct16.fits')
spectrum = REVERSE(REFORM(l[*,1]))

l        = REVERSE(REFORM(l[*,0]))-6173.3433-0.005;-0.15   ; for lyotsunsetTprofiles2006oct16.fits

aaa      = WHERE(ABS(l) LT 10.5d0) ; WE RESTRICT THE RANGE TO FIT FOR
l        = l[aaa]
data     = data[aaa]

; We correct for the solar spectrum curvature due to the spectrograph
x        = FINDGEN(1320)
x        = x[aaa]
;data     = data/(0.935+0.000305037*x+4.03807e-08*x^2.0-1.69869e-10*x^3.0)
data     = data/MAX(data) 

IF blocker EQ 1 THEN BEGIN
    q          = READFITS('blocker11.fits')
    blocker0   = SHIFT(q[*,1],2.75/(q[1,0]-q[0,0])) ; careful, the shift depends on the focal distance
    blocker0   = INTERPOL(blocker0/100.d0,q[*,0]-6173.3433,l)
ENDIF ELSE blocker0=1.0
SAVE,l,blocker0,FILE='temp.bin'

nel        = N_ELEMENTS(data)
weight     = FLTARR(nel)+1.d0
d          = 0.0

A          = FLTARR(12)
A[0]       = MAX(data)
A[1]       = 1.0
A[3]       = 1.0
A[5]       = 1.0
A[7]       = 1.0

A[9]      = 0.08  ; dark level
;A[13]      = 1.405 ; FSR E2
;A[12]      = 2.779 ; FSR E3
A[10]      = 0.1 ; contrast E1

resp = CURVEFIT(FINDGEN(nel),data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2)

print,chi2
PRINT,'contrasts',A[1],A[3],A[5],A[7],A[10]
PRINT,'amplitude',A[0]
PRINT,'angles',A[2]*180./!pi,A[4]*180./!pi,A[6]*180./!pi,A[8]*180./!pi,A[11]*180./!pi
PRINT,'Dark Level Offset',A[9]
;PRINT,'FSRs',A[12],A[13]

SET_PLOT,'ps'
DEVICE,FILE='yo2.ps',/color,xoffset=0,yoffset=0,xsize=20,ysize=16
loadct,3
plot,l,data,xst=1,yst=1,charsize=1.5,xtit='Wavelength (A)',ytit='Normalized Intensity',thick=1,yrange=[0,1],xrange=[-5.25,4.]
;oplot,l,resp,col=180,thick=3 ; FIT RESULT
oplot,[-7.5,7.5],[0.05,0.05],linestyle=2
;oplot,[-0.5,-0.5],[0,1],linestyle=2
;oplot,[0.5,0.5],[0,1],linestyle=2
oplot,[0,0],[0,1],linestyle=2
oplot,[-7.5,7.5],[0.5,0.5],linestyle=2
;oplot,REBIN(l,330),REBIN((resp-data)/data+0.5,330),thick=2
;oplot,[-10,10],[0.5,0.5],linestyle=2

DEVICE,/CLOSE

SET_PLOT,'x'
lresp=l
data=data*blocker0
data=data/MAX(data)
SAVE,lresp,data,file='LyotLampData.bin'
READ,pause

END
