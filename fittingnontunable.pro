; to fit non-tunasble profile using data of September 2007


PRO MyLyot,X,A,F,pder

COMMON blocker,blocker0,B

A[0]  = ABS(A[0])
;A[6]  = ABS(A[6])
;A[7]  = ABS(A[7])

nel   = N_ELEMENTS(X)
F     = FLTARR(nel)

dpi   = 2.d0*!dpi

FSR   = FLTARR(4) ; in A
FSR[0]= 1.405d0; A[6]; 1.405d0               ; for E2
FSR[1]= 2.779d0; A[7] ;2.779d0               ; for E3
FSR[2]= 5.682d0               ; for E4
FSR[3]= 11.354d0              ; for E5

E2    = (1.d0+B[0] *COS(dpi*X/FSR[0]+A[1]))/2.d0
E3    = (1.d0+B[1] *COS(dpi*X/FSR[1]+A[2]))/2.d0
E4    = (1.d0+B[2] *COS(dpi*X/FSR[2]+A[3]))/2.d0
E5    = (1.d0+B[3] *COS(dpi*X/FSR[3]+A[4]))/2.d0

F     = A[0]*blocker0*E2*E3*E4*E5+A[5]

dE2d2 = -B[0]*SIN(dpi*X/FSR[0]+A[0])/2.d0
dE3d2 = -B[1]*SIN(dpi*X/FSR[1]+A[1])/2.d0
dE4d2 = -B[2]*SIN(dpi*X/FSR[2]+A[2])/2.d0
dE5d2 = -B[3]*SIN(dpi*X/FSR[3]+A[3])/2.d0

pder  = [ [blocker0*E2*E3*E4*E5] , [dE2d2*E3*E4*E5*A[0]*blocker0] , [dE3d2*E2*E4*E5*A[0]*blocker0] , [dE4d2*E2*E3*E5*A[0]*blocker0] , [dE5d2*E2*E3*E4*A[0]*blocker0] , [FLTARR(nel)+1.d0]]; , [A[0]*blocker0*dpi*X*B[0]/FSR[0]^2.0*SIN(dpi*X/FSR[0]+A[1])/2.d0*E3*E4*E5] ,[A[0]*blocker0*dpi*X*B[1]/FSR[1]^2.0*SIN(dpi*X/FSR[1]+A[2])/2.d0*E2*E4*E5] ]


END

PRO fittingnontunable

COMMON blocker,blocker0,B

dpi        = 2.d0*!dpi
lamref     = 6173.3433d0
RESTORE,'temp.bin' ; must include wavelength and intensity

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin'
blocker0   = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
q          = READFITS('blocker11.fits')
blocker0   = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam) ; I center the profile ;+3.61328

nel        = N_ELEMENTS(data)
weight     = FLTARR(nel)+1.d0

A          = FLTARR(6)
A[0]       = MAX(data)
A[1]       = 0.d0
A[2]       = 0.d0
A[3]       = 0.d0
A[4]       = 0.d0
A[5]       = 0.d0
;A[6]       = 1.405d0
;A[7]       = 2.779d0

B          = FLTARR(4)+0.98d0


resp = CURVEFIT(lam,data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-12,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2)
PRINT,'PHASES:',A[1:4]*180./!dpi
;PRINT,'FSRs:',A[6:7]
PRINT,'AMPLITUDE:',A[0]

FSR   = FLTARR(4) ; in A
FSR[0]= 1.405d0 ; A[6]               ; for E2
FSR[1]= 2.779d0 ; A[7]              ; for E3
FSR[2]= 5.682d0               ; for E4
FSR[3]= 11.354d0   
E2    = (1.d0+B[0] *COS(dpi*lam/FSR[0]+A[1]))/2.d0
E3    = (1.d0+B[1] *COS(dpi*lam/FSR[1]+A[2]))/2.d0
E4    = (1.d0+B[2] *COS(dpi*lam/FSR[2]+A[3]))/2.d0
E5    = (1.d0+B[3] *COS(dpi*lam/FSR[3]+A[4]))/2.d0

F     = A[0]*blocker0*E2*E3*E4*E5+A[5]
SET_PLOT,'x'
!P.MULTI=0
PLOT,lam,data,psym=4
OPLOT,lam,F,col=180,psym=4
READ,pause


END
