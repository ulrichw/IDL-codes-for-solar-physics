; JONES CALCULUS FOR THE HMI MICHELSONS
;-----------------------------------------------------------------------------------------------------------------------------------



PRO sinewave2,X,A,F,pder

A[0]=ABS(A[0])

F = A[0] + (A[1]*cos(X+A[3])+A[2]*sin(X+A[3]))^2.d0

pder= [[FLTARR(N_ELEMENTS(X))+1.d0],[2.d0*A[1]*cos(X+A[3])^2.d0+2.d0*A[2]*cos(X+A[3])*sin(X+A[3])],[2.d0*A[2]*sin(X+A[3])^2.d0+2.d0*A[1]*cos(X+A[3])*sin(X+A[3])], [-A[1]^2.d0*2.d0*cos(X+A[3])*sin(X+A[3])+A[2]^2.d0*2.d0*sin(X+A[3])*cos(X+A[3])-2.d0*A[1]*A[2]*sin(X+A[3])^2.d0+2.d0*A[1]*A[2]*cos(X+A[3])^2.d0] ]


END


PRO sinewave,X,A,F,pder

A[0]=ABS(A[0])

F = A[0] + (A[1]*cos(X)+A[2]*sin(X))^2.d0

pder= [[FLTARR(N_ELEMENTS(X))+1.d0],[2.d0*A[1]*cos(X)^2.d0+2.d0*A[2]*cos(X)*sin(X)],[2.d0*A[2]*sin(X)^2.d0+2.d0*A[1]*cos(X)*sin(X)] ]


END

PRO JonesMichelson

nlam= 30000 ; NUMBER OF WAVELENGTHS
lam = (FINDGEN(nlam)/(nlam-1.)*150.-75.);+0.086d0
FSR = 0.172d0 ; FOR NB MICHELSON, in Angstroms

Rcs = 1.0
Rcp = 0.0
JR  = [[sqrt(Rcs),0.0],[0.0,sqrt(Rcp)]] ; JONES MATRIX FOR THE BEAM-SPLITTER IN REFLECTION

Tcs = 0.0
Tcp = 1.0
JT  = [[sqrt(Tcs),0.0],[0.0,sqrt(Tcp)]] ; JONES MATRIX FOR THE BEAM-SPLITTER IN TRANSMISSION

Rm  = 1.0
JPm = [[sqrt(Rm),0.0],[0.0,sqrt(Rm)]]   ; JONES MATRIX FOR THE PLANE MIRROR IN VACUUM ARM

RmS = 1.0
JPmS= [[sqrt(RmS),0.0],[0.0,sqrt(RmS)]]   ; JONES MATRIX FOR THE PLANE MIRROR IN SOLID ARM


epsE = 0.0*!dpi/180.d0 ; retardation error
betaE= 0.0*!dpi/180.   ; tilt error
Retard025E = [ [exp(COMPLEX(0,1)*(!dpi/4.d0+epsE))*cos(!dpi/4.d0+betaE)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsE))*sin(!dpi/4.d0+betaE)^2.d0,cos(!dpi/4.d0+betaE)*sin(!dpi/4.d0+betaE)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsE))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsE)))],[cos(!dpi/4.d0+betaE)*sin(!dpi/4.d0+betaE)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsE))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsE))),exp(COMPLEX(0,1)*(!dpi/4.d0+epsE))*sin(!dpi/4.d0+betaE)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsE))*cos(!dpi/4.d0+betaE)^2.d0] ]   ; JONES MATRIX FOR THE EXIT 1/4-WAVE PLATE

epsV = 0.0
betaV= 0.0*!dpi/180.
Retard025V = [ [exp(COMPLEX(0,1)*(!dpi/4.d0+epsV))*cos(!dpi/4.d0+betaV)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsV))*sin(!dpi/4.d0+betaV)^2.d0,cos(!dpi/4.d0+betaV)*sin(!dpi/4.d0+betaV)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsV))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsV)))],[cos(!dpi/4.d0+betaV)*sin(!dpi/4.d0+betaV)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsV))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsV))),exp(COMPLEX(0,1)*(!dpi/4.d0+epsV))*sin(!dpi/4.d0+betaV)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsV))*cos(!dpi/4.d0+betaV)^2.d0] ]   ; JONES MATRIX FOR THE VACUUM 1/4-WAVE PLATE

epsS = 0.0
betaS= 0.0*!dpi/180.
Retard025S = [ [exp(COMPLEX(0,1)*(!dpi/4.d0+epsS))*cos(!dpi/4.d0+betaS)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsS))*sin(!dpi/4.d0+betaS)^2.d0,cos(!dpi/4.d0+betaS)*sin(!dpi/4.d0+betaS)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsS))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsS)))],[cos(!dpi/4.d0+betaS)*sin(!dpi/4.d0+betaS)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsS))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsS))),exp(COMPLEX(0,1)*(!dpi/4.d0+epsS))*sin(!dpi/4.d0+betaS)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsS))*cos(!dpi/4.d0+betaS)^2.d0] ]   ; JONES MATRIX FOR THE SOLID 1/4-WAVE PLATE

;JV    = 1.d0/SQRT(2.)*[1.,1.]  ; INPUT JONES VECTOR

theta = 10.*!dpi/180.d0 ;tilt of the entrance polarizer
;PolR= 1.0
;PolW= 0.0
;Pol = [[cos(!dpi/4.+theta),-sin(!dpi/4.+theta)],[sin(!dpi/4.+theta),cos(!dpi/4.+theta)]]##[ [sqrt(PolR),0.],[0.,sqrt(PolW)] ]##[[cos(!dpi/4.+theta),sin(!dpi/4.+theta)],[-sin(!dpi/4.+theta),cos(!dpi/4.+theta)]]  ; JONES MATRIX OF INPUT POLARIZER

JV = [cos(!dpi/4.+theta),sin(!dpi/4.+theta)]

epsBS = 0.0*!dpi/180.
betaBS= 0.0*!dpi/180.
Retard025BS = [ [exp(COMPLEX(0,1)*(!dpi/4.d0+epsBS))*cos(!dpi/4.d0+betaBS)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsBS))*sin(!dpi/4.d0+betaBS)^2.d0,cos(!dpi/4.d0+betaBS)*sin(!dpi/4.d0+betaBS)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsBS))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsBS)))],[cos(!dpi/4.d0+betaBS)*sin(!dpi/4.d0+betaBS)*(exp(COMPLEX(0,1)*(!dpi/4.d0+epsBS))-exp(-COMPLEX(0,1)*(!dpi/4.d0+epsBS))),exp(COMPLEX(0,1)*(!dpi/4.d0+epsBS))*sin(!dpi/4.d0+betaBS)^2.d0+exp(-COMPLEX(0,1)*(!dpi/4.d0+epsBS))*cos(!dpi/4.d0+betaBS)^2.d0] ]   ; JONES MATRIX FOR THE BEAM-SPLITTER 1/4-WAVE PLATE

nphi  = 300
Inten = FLTARR(nlam,nphi)

FOR j=0l,nphi-1l DO BEGIN

    phi= 2.d0*!dpi/(nphi-1.)*j                   ; tuning angle

;TUNING DONE BY POLARIZER
   ;PJ = [ [cos(phi+!pi/4.)^2.d0,cos(phi+!pi/4.)*sin(phi+!pi/4.)],[cos(phi+!pi/4.)*sin(phi+!pi/4.),sin(phi+!pi/4.)^2.d0] ] ;JONES MATRIX OF TUNING POLARIZER (!pi/4 ADDED BECAUSE WHAT MATTERS IS THE ANGLE WITH THE 1/4-WAVE PLATE)
    PJ = [ [cos(phi)^2.d0,cos(phi)*sin(phi)],[cos(phi)*sin(phi),sin(phi)^2.d0] ]

;TUNING DONE BY 1/2-WAVE PLATE
    delt=10.0*!dpi/180.d0        ;retardation error
    gam = phi ; tuning angle
    HJ=cos(!dpi/2.d0+delt)*[[1.,0.],[0.,1.]]+COMPLEX(0,1)*sin(!dpi/2.d0+delt)*[[1.,0.],[0.,-1.]]##[[cos(2.d0*gam),sin(2.d0*gam)],[-sin(2.d0*gam),cos(2.d0*gam)]] ; JONES MATRIX OF TUNING HALF-WAVE PLATE
    phi2= 0.d0 ;angle of fixed polarizer after 1/2-wave plate
    Pol = [ [cos(phi2)^2.d0,cos(phi2)*sin(phi2)],[cos(phi2)*sin(phi2),sin(phi2)^2.d0] ]


    FOR i=0l,nlam-1l DO BEGIN
        JU = [ [exp(complex(0,1)*(!dpi*lam[i]/FSR)),0.0],[0.0,exp(complex(0,1)*(!dpi*lam[i]/FSR))] ] ; JONES MATRIX FOR THE BEAM PROPAGATION
        TJ = (Retard025E##JT##Retard025V##JPm##Retard025V##JR) + (Retard025E##JR##JU##Retard025S##JPmS##Retard025S##JU##JT)
       ;JVout = Retard025BS##(PJ##TJ)##JV ; TUNING BY POLARIZER
        JVout = Retard025BS##(Pol##HJ##TJ)##JV ; TUNING BY 1/2-WAVE PLATE
        JVout = REFORM(ABS(JVout))
        Inten[i,j]=JVout[0]^2.d0+JVout[1]^2.d0 ; transmittance
    ENDFOR

ENDFOR

inten2=fltarr(nphi)
FOR i=0l,nphi-1l DO inten2[i]=INT_TABULATED(lam,inten[*,i])

inten2=inten2/mean(inten2)
psi=2.d0*!dpi/(nphi-1.)*FINDGEN(nphi)
plot,psi,inten2,xst=1,yst=1,thick=3

weights=FLTARR(nphi)+1.d0
;A=[1.d0,-0.1,0.1,-139.d0/2.d0/180.d*!dpi]
A=[1.d0,-0.1,0.1]
resp      = CURVEFIT(psi,inten2,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2)
PRINT,A
;funct=A[0] + (A[1]*cos(psi+A[3])+A[2]*sin(psi+A[3]))^2.d0
funct=A[0] + (A[1]*cos(psi)+A[2]*sin(psi))^2.d0
OPLOT,psi,funct,col=180,linestyle=2,thick=3
PRINT,'I-RIPPLE=',(MAX(inten2)-MIN(inten2))/MEAN(inten2)


; maxi-mini is supposed to be constant
mini=fltarr(300)
maxi=mini
for i=0,299 do mini[i]=MIN(inten[*,i])
for i=0,299 do maxi[i]=MAX(inten[*,i])
temp=maxi-mini

read,pause


; TO MEASURE THE ERROR WE MAKE BY ASSUMING THAT \int NT T dl/(\int NT dl) = 0.5
; NON-TUNABLE TRANSMISSION PROFILE
;---------------------------------------------------------------------------------------------------------------------------------------
FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172-0.0010576d0; for the narrow-band Michelson
FSR[1]      = 0.344-0.00207683d0; for the broad-band  Michelson
FSR[2]      = 0.693+0.000483467d0; for E1
FSR[3]      = 1.407d0       ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5
FSR         = 2.d0*!dpi/FSR


nlam        = 2880*5*2
dlam        = .5d-10 ; in mm
l           = (DINDGEN(nlam)-(nlam-1)/2.)*dlam+0.00061733433d0 ; wavelength in mm
lam         = l*1.d7 ; wavelength in Angstroms
lam0        = 6173.3433d0              ; HMI Fe line in Angstrom
lam         = lam-lam0
transmission= FLTARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength  = REFORM(transmission[0,*])
transmission= REFORM(transmission[1,*])
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0,lam+lam0)
q           = READFITS('blocker11.fits')
blocker     = blocker0*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6d0,lam+lam0)
contrast    = [0.969,0.997,0.974d0,0.986,0.965,0.988,0.999]; contrast of E1 is adjusted to match the I-ripple effect as well as possible
phase       = [10.,-20.,80.,-6.13,-0.174,-5.34,1.14]*!dpi/180.d0

profilef     = blocker ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrast[i]*COS(FSR[i]*lam+phase[i]))/2.d0 ; NON-TUNABLE TRANSMISSION PROFILE



; TUNING SEQUENCE IN WHICH E1 AND WB ARE FIXED, BUT NB IS TUNED


profilef2=profilef*(1.d0+contrast[2]*COS(FSR[2]*lam+phase[2]))/2.0*(1.d0+contrast[1]*COS(FSR[1]*lam+phase[1]))/2.0

nphi=500
phi= 2.d0*!dpi/(nphi-1.)*FINDGEN(nphi)
integ=FLTARR(nphi)
TOTALI=TOTAL(profilef2)

FOR i=0,nphi-1 DO integ[i]=TOTAL(profilef2*(1.d0+contrast[0]*COS(FSR[0]*lam+phase[0]+phi[i]))/2.0)/TOTALI
PRINT,(MAX(integ)-MIN(integ))/MEAN(integ)


; TUNING SEQUENCE IN WHICH NB AND WB ARE FIXED, BUT E1 IS TUNED


profilef2=profilef*(1.d0+contrast[0]*COS(FSR[0]*lam+phase[0]))/2.0*(1.d0+contrast[1]*COS(FSR[1]*lam+phase[1]))/2.0

nphi=500
phi= 2.d0*!dpi/(nphi-1.)*FINDGEN(nphi)
integ=FLTARR(nphi)
TOTALI=TOTAL(profilef2)

FOR i=0,nphi-1 DO integ[i]=TOTAL(profilef2*(1.d0+contrast[2]*COS(FSR[2]*lam+phase[2]+phi[i]))/2.0)/TOTALI
PRINT,(MAX(integ)-MIN(integ))/MEAN(integ)



END
