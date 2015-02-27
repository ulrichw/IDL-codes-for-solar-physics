; Jones calculus for the HMI Lyot elements
;---------------------------------------------------------------------------


PRO sinewave,X,A,F,pder

A[0]=ABS(A[0])

F = A[0] + (A[1]*cos(X)+A[2]*sin(X))^2.d0

pder= [[FLTARR(N_ELEMENTS(X))+1.d0],[2.d0*A[1]*cos(X)^2.d0+2.d0*A[2]*cos(X)*sin(X)],[2.d0*A[2]*sin(X)^2.d0+2.d0*A[1]*cos(X)*sin(X)] ]


END

PRO JonesLyot

nlam= 30000 ; NUMBER OF WAVELENGTHS
lam = (FINDGEN(nlam)/(nlam-1.)*150.-75.);+0.086d0
FSR = 0.690d0 ; FOR E1, in Angstroms

theta = 0.0*!dpi/180.d0 ; tilt of the entrance polarizer (NORMALLY AT 0 DEGREE)
;Rpol=1.0
;Jpol=0.0
;Pol= [[cos(theta), sin(theta)],[-sin(theta),cos(theta)]]##[ [SQRT(Rpol),0.],[0.,SQRT(Jpol)] ]##[[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]
;JV = [1.,0.]  ; INPUT JONES VECTOR

JV = [cos(theta),sin(theta)] ; INPUT JONES VECTOR

err    = 0.d0*!dpi/180.d0;tilt
theta  = !dpi/4.d0+err ;angle of first WF block
ROTpi4 = [ [cos(theta),-sin(theta)],[ sin(theta),cos(theta)] ]
ROTpi4i= [ [cos(theta), sin(theta)],[-sin(theta),cos(theta)] ] ; for -theta

theta2 =  -!dpi/4.d0 ;angle of second WF block
ROTpi42=  [ [cos(theta2),-sin(theta2)],[ sin(theta2),cos(theta2)] ]
ROTpi42i= [ [cos(theta2), sin(theta2)],[-sin(theta2),cos(theta2)] ]

beta = 0.0*!dpi/180.   ; tilt error
RET14= COMPLEX(1,1)/SQRT(2.)*[[1.,0.],[0.,-COMPLEX(0.,1)]]  ; JONES MATRIX FOR THE EXIT 1/4-WAVE PLATE
RET14= [ [cos(beta), sin(beta)],[-sin(beta),cos(beta)] ] ##RET14##[ [cos(beta),-sin(beta)],[ sin(beta),cos(beta)] ]

delt=0.0*!dpi/180.d0 ;retardation error
gam =2.0*!dpi/180.d0 ;tilt error
RET12=cos(!dpi/2.d0+delt)*[[1.,0.],[0.,1.]]+COMPLEX(0,1)*sin(!dpi/2.d0+delt)*[[1.,0.],[0.,-1.]]##[[cos(2.d0*gam),-sin(2.d0*gam)],[sin(2.d0*gam),cos(2.d0*gam)]] ; JONES MATRIX OF HALF-WAVE PLATE

A0=1.0
B0=1.0


nphi  = 300
Inten = FLTARR(nlam,nphi)
psi=2.d0*!dpi/(nphi-1.)*FINDGEN(nphi)

FOR j=0l,nphi-1l DO BEGIN

    PRINT,j
    phi= psi[j]            ; tuning angle
;TUNING DONE BY POLARIZER
    PJ = [ [cos(phi)^2.d0,-cos(phi)*sin(phi)],[-cos(phi)*sin(phi),sin(phi)^2.d0] ] ;JONES MATRIX OF TUNING POLARIZER
    
;TUNING DONE BY 1/2-WAVE PLATE
    delt=0.0*!dpi/180.d0        ;retardation error
    gam = phi ; tuning angle
    HJ=cos(!dpi/2.d0+delt)*[[1.,0.],[0.,1.]]+COMPLEX(0,1)*sin(!dpi/2.d0+delt)*[[1.,0.],[0.,-1.]]##[[cos(2.d0*gam),-sin(2.d0*gam)],[sin(2.d0*gam),cos(2.d0*gam)]] ; JONES MATRIX OF TUNING HALF-WAVE PLATE
    phi2= 0.d0 ;angle of fixed polarizer after 1/2-wave plate
    Pol = [ [cos(phi2)^2.d0,-cos(phi2)*sin(phi2)],[-cos(phi2)*sin(phi2),sin(phi2)^2.d0] ]

    FOR i=0l,nlam-1l DO BEGIN
        JU1 = [ [SQRT(A0)*exp(complex(0,1)*(!dpi*lam[i]/FSR/2.)),0.0],[0.0,SQRT(B0)*exp(-complex(0,1)*(!dpi*lam[i]/FSR/2.))] ] ; JONES MATRIX FOR THE BEAM PROPAGATION, 1st WF BLOCK
        JU1 = ROTpi4i##JU1##ROTpi4
        JU2 = [ [SQRT(A0)*exp(complex(0,1)*(!dpi*lam[i]/FSR/2.)),0.0],[0.0,SQRT(B0)*exp(-complex(0,1)*(!dpi*lam[i]/FSR/2.))] ] ; JONES MATRIX FOR THE BEAM PROPAGATION, 2nd WF BLOCK
        JU2 = ROTpi42i##JU2##ROTpi42


       ;TJ =  PJ##RET14##JU2##RET12##JU1; TOTAL JONES MATRIX (TUNING WITH POLARIZER)
        TJ = Pol##HJ##RET14##JU2##RET12##JU1 ; TOTAL JONES MATRIX (TUNING WITH 1/2-WAVE PLATE)
        JVout = TJ##JV
        JVout = REFORM(ABS(JVout))
        Inten[i,j]=JVout[0]^2.d0+JVout[1]^2.d0 ; transmittance
    ENDFOR

ENDFOR

inten2=fltarr(nphi)
FOR i=0l,nphi-1l DO inten2[i]=INT_TABULATED(lam,inten[*,i])

inten2=inten2/mean(inten2)
plot,psi,inten2,xst=1,yst=1,thick=3

weights=FLTARR(nphi)+1.d0
A=[1.d0,-0.1,0.1]
resp      = CURVEFIT(psi,inten2,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2)
PRINT,A
funct=A[0] + (A[1]*cos(psi)+A[2]*sin(psi))^2.d0
OPLOT,psi,funct,col=180,linestyle=2,thick=3
PRINT,'I-RIPPLE=',(MAX(inten2)-MIN(inten2))/MEAN(inten2)

; maxi-mini is supposed to be constant
mini=fltarr(300)
maxi=mini
for i=0,299 do mini[i]=MIN(inten[*,i])
for i=0,299 do maxi[i]=MAX(inten[*,i])
temp=maxi-mini

READ,pause

END
