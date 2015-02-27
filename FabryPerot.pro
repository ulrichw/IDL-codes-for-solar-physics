; This program simulates a FabryPerot to understand why the blocker
; filter and front window produce some interference patterns

PRO fitting,X,A,F,pder

COMMON properties,d,n,lamref,dpi

freq  = (d*2.d0*n+lamref)/(lamref)^2.d0
F = A[0] + A[1]*COS(dpi*X*freq) + A[2]*SIN(dpi*X*freq)

pder = [ [FLTARR(N_ELEMENTS(X))+1.d0],[COS(dpi*X*freq)],[SIN(dpi*X*freq)] ]

END


PRO FabryPerot

COMMON properties,d,n,lamref,dpi

dpi         = 2.d0*!dpi
FSR         = DBLARR(7)      ; FSR in Angstrom
FSR[0]      = 0.172457d0     ; for the narrow-band Michelson
FSR[1]      = 0.344242d0     ; for the broad-band  Michelson
FSR[2]      = 0.690d0        ; for E1
FSR[3]      = 1.405d0        ; for E2
FSR[4]      = 2.779d0        ; for E3
FSR[5]      = 5.682d0        ; for E4
FSR[6]      = 11.354d0       ; for E5

n = 1.519569    ; refractive index
d = 0.005*1.d10;0.005*1.d10 ; thickness in A
T = 0.95d0      ; transmissivity
R = 0.05d0      ; reflectivity
nx=256
ny=256
Inten = FLTARR(nx,ny)
;theta = FLTARR(nx,ny)
;center = [-100,360]
;dtheta = 0.0085*!pi/180.;0.006*!pi/180.
;FOR i=0,nx-1 DO FOR j=0,ny-1 DO theta[i,j] =
;SQRT((i-center[0])^2.d0+(j-center[1])^2.d0)*dtheta
dista =DIST(2.*nx,2.*ny)*0.023*1.d10/FLOAT(nx)
dista = dista[0:nx-1,0:ny-1]
theta= 0.0115*!pi/180.d0
dista = dista*SIN(theta)

mask = SHIFT(DIST(nx,ny),nx/2,ny/2)*0.5d0*4096.d0/nx
a    = WHERE(mask LE 990.,COMPLEMENT=b)
mask[a] = 1.0
mask[b] = 0.0


nlam = 2250                ; number of wavelengths
dlam = 3.6d0/1.d3          ; resolution in wavelength
lamref= 6173.3433d0
lam  =(DINDGEN(nlam)-(nlam-1)/2.)*dlam+lamref ; wavelengths in Angstrom
q    = READFITS('blocker11.fits')
profile     = INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,lam)
FOR i=0,6 DO profile     = profile*(1.d0+COS(dpi*lam/FSR[i]))/2.d0


FOR i=0,nx-1 DO BEGIN
    PRINT,i
    FOR j=0,ny-1 DO BEGIN
        delta = 2.d0*!dpi/lam*(2.d0*(d+dista[i,j])*n);*COS(theta[i,j])
        Inten0 = T^2.d0/(1.d0+R^2.d0-2.d0*R*COS(delta)) ; intensity of outgoing light
        Inten[i,j] = TOTAL(Inten0*profile)
    ENDFOR
ENDFOR


;weight = FLTARR(nlam)+1.d0
;A = [1.0,0.,0.]
;resp      = CURVEFIT(lam,I,weight,A,FUNCTION_NAME='fitting',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2)

;WINDOW,0,RETAIN=2
;!P.MULTI=[0,1,2]
;PLOT,lam,I,xst=1,yst=1
;oplot,lam,resp,col=180
;PLOT,lam,(resp-I)/I,xst=1,yst=1
;PLOT,theta,inten
;TVIM,Inten*mask

SET_PLOT,'PS'
DEVICE,FILE='yo.ps',BITS=24,xoffset=0,yoffset=0,xsize=20,ysize=20,/COLOR
LOADCT,3
Inten=Inten/MEAN(Inten)
TVIM,Inten*mask,pcharsize=1.5,xtit='Pixels',ytit='Pixels',stit='Amplitude',/scale,range=[0.99*MIN(Inten),MAX(Inten)]
DEVICE,/close
SET_PLOT,'X'


READ,pause

END
