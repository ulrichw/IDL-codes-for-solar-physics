; This program shows what is the output intensity when we tune NB and
; WB with the tuning polarizer instead of the two tuning waveplates

PRO doublesinewave,X,A,F,pder
COMMON wavelength,lam
FSR0 = 0.172d0
FSR1 = 0.344d0
B    = [0.98,0.97]

COS0 = COS(2.d0*!dpi*lam/FSR0+A[0]+2.*X)
COS1 = COS(2.d0*!dpi*lam/FSR1+A[1]-2.*X)
SIN0 =-SIN(2.d0*!dpi*lam/FSR0+A[0]+2.*X)
SIN1 =-SIN(2.d0*!dpi*lam/FSR1+A[1]-2.*X)

F = ALOG((1.d0+B[0]*COS0)*(1.d0+B[1]*COS1)*A[2]/4.d0)

pder = [ [B[0]*SIN0*(1.d0+B[1]*COS1)*A[2]/4.d0],[SIN1*B[1]*(1.d0+B[0]*COS0)*A[2]/4.d0],[(1.d0+B[0]*COS0)*(1.d0+B[1]*COS1/4.d0)] ]



END


PRO doublesinewave2,X,A,F,pder

FSR0 = 0.172d0
FSR1 = 0.344d0

lam  = 0.0

COS0 = COS(2.d0*!dpi*lam/FSR0+2.*X)
SIN0 = SIN(2.d0*!dpi*lam/FSR0+2.*X)
COS1 = COS(2.d0*!dpi*lam/FSR1-2.*X)
SIN1 = SIN(2.d0*!dpi*lam/FSR1-2.*X)

F = (1.d0+A[0]*COS0-A[1]*SIN0)*(1.d0+A[2]*COS1-A[3]*SIN1)*A[4]

pder = [ [COS0*(1.d0+A[2]*COS1-A[3]*SIN1)*A[4]],[-SIN0*(1.d0+A[2]*COS1-A[3]*SIN1)*A[4]],[(1.d0+A[0]*COS0-A[1]*SIN0)*COS1*A[4]],[-SIN1*(1.d0+A[0]*COS0-A[1]*SIN0)*A[4]],[(1.d0+A[0]*COS0-A[1]*SIN0)*(1.d0+A[2]*COS1-A[3]*SIN1)] ]


END


PRO HMI_tuningpolarizer

COMMON wavelength,lam

WINDOW,0,retain=2
FSR  = [0.172d0,0.344d0]

nangle = 24
;nlam = 2000
;dlam = 3.6d0/1.0d3      ; resolution in wavelength
;lam  = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
lam  =  0.008d0
profile = 1454.

B    = [.98,.97] ; contrasts
Phi  = [-91.,151.456]*!dpi/180.d0 ; phases

INB  = FLTARR(nangle)
IWB  = INB

angle= FINDGEN(nangle)/(nangle-1.)*2.*!dpi
sign = [1.,-1.]

FOR i=0,nangle-1 DO BEGIN

    INB[i]  = (1.d0+B[0]*COS(2.d0*!dpi*lam/FSR[0]+Phi[0]+sign[0]*2.*angle[i]))/2.
    IWB[i]  = (1.d0+B[1]*COS(2.d0*!dpi*lam/FSR[1]+Phi[1]+sign[1]*2.*angle[i]))/2.

ENDFOR
data = ALOG(INB*IWB*profile)
weights = FLTARR(nangle)+1.0

;A = [Phi[0]*0.8,Phi[1]*1.25,MAX(data)]
A = [0.,0.,ALOG(MAX(data))]
PRINT,'GUESS=',A[0]*180./!pi MOD 360.,A[1]*180./!pi MOD 360.,A[2]
;A=[1.,0.,1.,0.,1./4.]

resp      = CURVEFIT(angle,data,weights,A,FUNCTION_NAME='doublesinewave',TOL=1.d-9,ITMAX=15000,/DOUBLE,CHISQ=chi2,STATUS=stat,/NODERIVATIVE)
;resp      = CURVEFIT(angle,data,weights,A,FUNCTION_NAME='doublesinewave2',TOL=1.d-7,ITMAX=15000,/DOUBLE,CHISQ=chi2,STATUS=stat)
PRINT,'RESULT=',A[0]*180./!pi MOD 360.,A[1]*180./!pi MOD 360.,A[2]

doublesinewave,angle,A,F,pder

PLOT,angle,data,xst=1
OPLOT,angle,F,LINESTYLE=2,THICK=2
PRINT,'STATUS=',stat
READ,pause

END
