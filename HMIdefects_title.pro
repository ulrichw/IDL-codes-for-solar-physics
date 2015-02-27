; use equations of Alan Title to simulate behaviour of Lyot elements
; and understand the origin of an I-ripple

PRO sinewave,X,A,F,pder

T   = ABS(A[2])
A[0]= ABS(A[0]) ; we want a positive amplitude

F = A[0]*SIN(X/T+A[1])+1.d0

pder= [[SIN(X/T+A[1])],[A[0]*COS(X/T+A[1])],[-A[0]/T^2.0*X*COS(X/T+A[1])]]

END


PRO sinewave2,X,A,F,pder

A[0]=ABS(A[0])
A[1]=ABS(A[1])
A[2]=ABS(A[2])

F = A[0] + (A[1]*cos(X)+A[2]*sin(X))^2.d0

pder= [[FLTARR(N_ELEMENTS(X))+1.d0],[2.d0*A[1]*cos(X)^2.d0+2.d0*A[2]*cos(X)*sin(X)],[2.d0*A[2]*sin(X)^2.d0+2.d0*A[1]*cos(X)*sin(X)] ]


END



;-----------------------------------------------------------------------------------------------------------------------------------------
;
;
; MAIN PROGRAM
;
;
;-----------------------------------------------------------------------------------------------------------------------------------------


PRO HMIdefects_title



;-----------------------------------------------------------------------------------------------------------------------------------------
;
; TRANSMISSION PROFILE OF E1 INCLUDING 1/2-WAVE AND 1/4-WAVE PLATES ERRORS
;
;-----------------------------------------------------------------------------------------------------------------------------------------

; WARNING: TUNING WITH A POLARIZER INSTEAD OF 1/2-WAVE PLATE

; FOR THE WF LYOT ELEMENT E1

nlam   = 2880*5*2
dlam   = .5d-10 ; in mm
l      = (DINDGEN(nlam)-(nlam-1)/2.)*dlam+0.00061733433d0 ; wavelength in mm

FSR    = (0.693+0.000483467d0)*1.d-7 ; FSR in mm of the Lyot element
d1     = 15.310*2.d0             ; calcite thickness
d2     = 3.528 *2.d0             ; ADP thickness
Dn2    = -0.0449016d0            ; birefringence of ADP at 25 C and 6173 Angstrom
Dn1    = -(6173.3433d-7)^2.d0/FSR/d1  + d2 * Dn2/d1; birefringence of calcite from expected FSR
d      = !dpi * d1 * Dn1         ; retardance of Calcite + ADP block
d      = d/l

nang   = 50    ; NUMBER OF TUNING ANGLES
ang    = FINDGEN(nang)/(nang-1.d0)*2.d0*!dpi ; TUNING ANGLE


; ERRORS IN E1
beta   = 0.d0*!dpi/180.d0 ; ERROR OF 1/4-WAVE PLATE IN TILT
eps    = 0.d0*!dpi/4.d0 ; ERROR OF 1/4-WAVE PLATE IN RETARDATION

gam    = -2.1d0*!dpi/180.d0 ; ERROR OF 1/2-WAVE PLATE IN TILT
delt   = 0.0d0*!dpi/180.d0 ; ERROR OF 1/2-WAVE PLATE IN RETARDATION


; NON-TUNABLE TRANSMISSION PROFILE
FSR         = DBLARR(7)             ; FSR in Angstrom
FSR[0]      = 0.172-0.0010576d0; for the narrow-band Michelson
FSR[1]      = 0.344-0.00207683d0; for the broad-band  Michelson
FSR[2]      = 0.693+0.000483467d0; for E1
FSR[3]      = 1.407d0       ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5
inttune     = 5./2                 ; Set number of tuning positions over wnarrow. (1/spacing)
dtune       = FSR[0]/inttune      ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
FSR         = 2.d0*!dpi/FSR

lam          = l*1.d7 ; wavelength in Angstroms
lam0         = 6173.3433d0              ; HMI Fe line in Angstrom
lam          = lam-lam0
transmission = FLTARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0,lam+lam0)
q            = READFITS('blocker11.fits')
blocker      = blocker0*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854,lam+lam0) ; +2.6 A
phase        = [0.0d0,0.d0,0.0d0,-6.09d0,-0.011d0,-4.8d0,1.6d0]*!dpi/180.d0 
contrast     = [0.978,0.997,0.974,0.985,0.965,0.988,1.0]  

profilef     = blocker ; we use the averaged blocker+front window profile
FOR i=3,6 DO profilef = profilef * (1.d0+contrast[i]*COS(FSR[i]*lam+phase[i]))/2.d0 ; NON-TUNABLE TRANSMISSION PROFILE

ntune  = 6                    ; Set number of tuning positions (NUMBER OF HMI FILTERS)
tune   = DBLARR(ntune)
tune   =(-(ntune-1)/2.+dindgen(ntune))*dtune 
dlamdv = lam0/3d8
dvdlam = 1/dlamdv
dvtune = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s

; PROFILE OF E1

Tr14   = FLTARR(nlam,nang)
Tr12   = Tr14
int14  = FLTARR(nang) ; integral values of transmission profile
int12  = int14


FOR i=0,nang-1 DO BEGIN

; SECOND-ORDER EXPANSION OF TRANSMISSION INCLUDING ERRORS IN THE
; 1/4-WAVE PLATE
    Tr14[*,i]  = (1.d0-2.d0*beta^2.d0-eps^2.d0)*cos(ang[i]-d)^2.d0-beta*sin(2.d0*ang[i]-d)-beta*eps*2.d0*sin(2.d0*ang[i])*cos(d)+eps^2.d0*cos(ang[i]+d)^2.d0+2.d0*beta^2.d0*(cos(d)^2.d0*sin(ang[i])^2.d0+cos(ang[i])^2.d0*sin(d)^2.d0)
    int14[i]   = TOTAL(Tr14[*,i]*profilef)*dlam

; SECOND-ORDER EXPANSION OF TRANSMISSION INCLUDING ERRORS IN THE
; 1/2-WAVE PLATE
    Tr12[*,i]  = cos(2.d0*gam)^2.d0*cos(delt)^2.d0*cos(ang[i]-d)^2.d0+(sin(delt)*cos(ang[i])+sin(2.d0*gam)*cos(delt)*sin(ang[i]))^2.d0
    int12[i]   = TOTAL(Tr12[*,i]*blocker)*dlam

ENDFOR

!P.MULTI=[0,1,3]
PLOT,l*1.d7,Tr14[*,0],xrange=[6170,6176.5],xst=1,thick=2
FOR i=1,nang-1 DO OPLOT,l*1.d7,Tr14[*,i],col=i*20.,thick=2
PLOT,l*1.d7,Tr12[*,0],xrange=[6170,6176.5],xst=1,thick=2
FOR i=1,nang-1 DO OPLOT,l*1.d7,Tr12[*,i],col=i*20.,thick=2
OPLOT,l*1.d7,blocker
int12=int12/MEAN(int12)
int14=int14/MEAN(int14)
PLOT,ang,int12,xst=1,thick=2
OPLOT,ang,int14,col=180,thick=2

;------------------------------------------------------------------------------------------------------------------------
; WE FIT THE I-RIPPLE WITH A SINE WAVE
;------------------------------------------------------------------------------------------------------------------------
weights = FLTARR(nang)+1.0

;A = [0.02,-!dpi/2.,.5]
A=[1.d0,0.1,0.1]

;resp      = CURVEFIT(ang,int12,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-10,ITMAX=5000,/DOUBLE,CHISQ=chi2)
resp      = CURVEFIT(ang,int12,weights,A,FUNCTION_NAME='sinewave2',TOL=1.d-13,ITMAX=5000,/DOUBLE,CHISQ=chi2)

PRINT,A
ang2=FINDGEN(10000)/9999.*2.d0*!dpi
;funct=A[0]*SIN(ang2/A[2]+A[1])+1.d0
funct=A[0] + (A[1]*cos(ang2)+A[2]*sin(ang2))^2.d0
PLOT,ang,int12,xst=1,yrange=[MIN(funct),MAX(funct)],yst=1,charsize=1.5
OPLOT,ang2,funct,col=180,linestyle=2
; I-RIPPLE VALUES
PRINT,'I-RIPPLE (from fit) =',(max(funct)-MIN(funct))/MEAN(funct)
PRINT,'I-RIPPLE (measured) =',(MAX(int12)-MIN(int12))/MEAN(int12)
REAd,pause


;------------------------------------------------------------------------------------------------------------------------
;
; WE DETERMINE THE ERROR ON THE DOPPLER VELOCITY PRODUCED BY USING THE
; WRONG MODEL FOR E1 (NO I-RIPPLE ACCOUNTED FOR)
;
;------------------------------------------------------------------------------------------------------------------------


;------------------------------------------------------------------------------------------------
; FOR THE DOPPLER VELOCITY DETERMINATION
; WE GENERATE ntest DOPPLER SHIFTED LINE PROFILES
;------------------------------------------------------------------------------------------------

dlamdv = lam0/3d8
OPENR,1,'Ulrich_Fe_0.txt'
data = FLTARR(2,98)
READF,1,data
lamp = REFORM(data[0,*])
nlamp= N_ELEMENTS(lamp)
line = REFORM(data[1,*])
close,1
ntest         = 501 ; Set test velocities: we will use ntest different velocities, resulting
                                ; in ntest Doppler shifted line profiles
dvtest        = 100.; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1)/2)
lines         = DBLARR(nlam,ntest)
dlinesdv      = DBLARR(nlam,ntest)
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
q2            = [MIN(lamp)-10,lamp,MAX(lamp)+10]
FOR i = 0,ntest-1 DO lines(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)

; AFTER JULY 2007, FRONT WINDOW S/N 1 WAS REPLACED BY S/N 3

filters = DBLARR(nlam,ntune) ; contain the filter profiles (Blocker+Lyot+Michelsons)
filters2=filters

; REFERENCE TRANSMISSION PROFILES FOR THE 6 HMI FILTERS: CONTRASTS=1, PHASES=0
;------------------------------------------------------------------------------

FOR itune = 0,ntune-1 DO filters[*,itune] =profilef*(1.d0+contrast[2]*COS(FSR[2]*(lam+tune[itune])+phase[2]))/2.0 ;E1 MODEL, EXCLUDING I-RIPPLE
FOR itune = 0,ntune-1 DO FOR i  = 0,1 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(FSR[i]*(lam+tune[itune])+phase[i]))/2 ; MICHELSONS

FOR itune = 0,ntune-1 DO filters2[*,itune]=profilef*(cos(2.d0*gam)^2.d0*cos(delt)^2.d0*(1.d0+contrast[2]*COS(FSR[2]*(lam+tune[itune])+phase[2]))/2.0+(sin(delt)*cos(FSR[2]*tune[itune]/2.d0)+sin(2.d0*gam)*cos(delt)*sin(FSR[2]*tune[itune]/2.d0))^2.d0) ;E1 MODEL, EXCLUDING I-RIPPLE
FOR itune = 0,ntune-1 DO FOR i  = 0,1 DO filters2[*,itune] = filters2[*,itune]*(1.d0+contrast[i]*COS(FSR[i]*(lam+tune[itune])+phase[i]))/2 ; MICHELSONS



; REFERENCE DOPPLER VELOCITY ERROR COMPUTED FROM A MDI-LIKE ALGORITHM
; WITH 6 FILTERS. WE APPLY A FACTOR SQRT(2) TO SIMULATE FOR THE 
; LCP AND RCP MEASUREMENTS
; WE ALSO PROVIDE THE RESULTS FOR 2 AVERAGED MDI-LIKE ALGORITHMS,
; USING THE PHASES OF THE 1ST AND 2ND FOURIER COMPONENTS OF THE LINE
;--------------------------------------------------------------------
pv1       = dvtune*inttune*2.0
pv2       = dvtune*inttune
x         = 2.0*!dpi*(-(ntune-1)/2.+DINDGEN(ntune))/ntune
well      = 2.d5;1.75d5

inten     = TRANSPOSE(filters)#lines
inten2    = TRANSPOSE(filters2)#lines
texp      = well/MAX(inten) 

; LOOK-UP TABLES OBTAINED FROM THE MODEL OF E1 WE USE (NO I-RIPPLE
; ACCOUNTED FOR)
c1        = REFORM(COS(x)#inten)  ; coefficient a1
s1        = REFORM(SIN(x)#inten)  ; coefficient b1 (for the sine)
c2        = REFORM(COS(2.d0*x)#inten)  ; coefficient a2
s2        = REFORM(SIN(2.d0*x)#inten)  ; coefficient b2 (for the sine)
phi1      = ATAN(-s1,-c1)
phi2      = ATAN(-s2,-c2)
vel1      = phi1*pv1/2.d0/!dpi
vel2      = phi2*pv2/2.d0/!dpi
vel1a     = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; The look-up table will be (vtest,vel1) or (vtest,vel1a)
vel2a     = (vel2-vtest+10.5d0*pv2) MOD pv2-pv2/2.d0+vtest

; LOOK-UP TABLES OBTAINED FROM THE ACTUAL E1 (WITH I-RIPPLE)
c1        = REFORM(COS(x)#inten2)  ; coefficient a1
s1        = REFORM(SIN(x)#inten2)  ; coefficient b1 (for the sine)
c2        = REFORM(COS(2.d0*x)#inten2)  ; coefficient a2
s2        = REFORM(SIN(2.d0*x)#inten2)  ; coefficient b2 (for the sine)
phi1      = ATAN(-s1,-c1)
phi2      = ATAN(-s2,-c2)
vel1      = phi1*pv1/2.d0/!dpi
vel2      = phi2*pv2/2.d0/!dpi
vel1b     = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; The look-up table will be (vtest,vel1) or (vtest,vel1a)
vel2b     = (vel2-vtest+10.5d0*pv2) MOD pv2-pv2/2.d0+vtest

nt=10000
           ;WE COMPUTE THE ERRORS: SYSTEMATIC + RMS VARIATION DUE TO PHOTON NOISE
            sigmat1 = FLTARR(ntest)
            sigmat=sigmat1
            systematicerror= FLTARR(ntest)
            systematicerror2=systematicerror
            iseed   = 1000l
            randomge = FLTARR(nt,ntune)

            FOR i= 0,ntune-1 DO randomge[*,i] = RANDOMN(iseed,nt)

            FOR itest= 0,ntest-1 DO BEGIN
 
                vt   = vtest[itest]
                i0   = inten2[*,itest]*texp

                c1   = TOTAL(COS(x)*i0)
                s1   = TOTAL(SIN(x)*i0)
                c2   = TOTAL(COS(2*x)*i0)
                s2   = TOTAL(SIN(2*x)*i0)
                phi1t = ATAN(-s1,-c1)
                phi2t = ATAN(-s2,-c2)
                vel1t = phi1t*pv1/2.0/!dpi
                vel2t = phi2t*pv2/2/!dpi
                vel1at= vel1t
                vel2at= (vel2t-vel1t+10.5*pv2) MOD pv2-pv2/2+vel1t
                v1tSYSa = INTERPOL(vtest,vel1a,vel1at) ; measured velocity using 1st look-up table (one from model with no I-ripple)
                v2tSYSa = INTERPOL(vtest,vel2a,vel2at)
                v1tSYSb = INTERPOL(vtest,vel1b,vel1at) ; measured velocity using 2nd look-up table (one from real model with I-ripple) 
                v2tSYSb = INTERPOL(vtest,vel2b,vel2at)

                systematicerror[itest] = ((v1tSYSa-vt)+(v2tSYSa-vt))/2.d0
                systematicerror2[itest]= ((v1tSYSb-vt)+(v2tSYSb-vt))/2.d0

                intent=DBLARR(ntune,nt)
                FOR i= 0,ntune-1 DO intent(i,*)=i0(i)+SQRT(i0(i))*randomge[*,i] ; SQRT(i[0]) is photon noise
                c1   = REFORM(COS(x)#intent)
                s1   = REFORM(SIN(x)#intent)
                c2   = REFORM(COS(2*x)#intent)
                s2   = REFORM(SIN(2*x)#intent)
                phi1t = ATAN(-s1,-c1)
                phi2t = ATAN(-s2,-c2)
                vel1t = phi1t*pv1/2.0/!dpi
                vel2t = phi2t*pv2/2/!dpi
                vel1at= vel1t
                vel2at= (vel2t-vel1t+10.5*pv2) MOD pv2-pv2/2+vel1t
                v1t = INTERPOL(vtest,vel1a,vel1at) ; measured velocity
                v2t = INTERPOL(vtest,vel2a,vel2at)
                sigmat1[itest] = SQRT(REBIN( ( ((v1t-v1tSYSa)+(v2t-v2tSYSa))/2.d0 )^2.0,1))/SQRT(2.d0) ; dispersion on the measured velocity for what we think are the good filters. I DIVIDE BY SQRT(2.d0) BECAUSE I DO LC AND RCP SEPARATELY, AND THEN AVERAGE THEIR VALUE
                v1t = INTERPOL(vtest,vel1b,vel1at) ; measured velocity
                v2t = INTERPOL(vtest,vel2b,vel2at)
                sigmat [itest] = SQRT(REBIN( ( ((v1t-v1tSYSb)+(v2t-v2tSYSb))/2.d0 )^2.0,1))/SQRT(2.d0) ; dispersion on the measured velocity from the actual filters

            ENDFOR 

!P.MULTI=0
PLOT,vtest,vel1a,xrange=[-8500,8500],xst=1,charsize=1.5
OPLOT,vtest,vel1b

SET_PLOT,'PS'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=28,/color,bits=24
LOADCT,3
!P.MULTI=[0,1,2]
plot,vtest,systematicerror,xst=1,xrange=[-6500,6500],thick=3,tit='!17 MDI-like algorithm, 6 filters, at disk center',xtit='Input Velocity (m/s)',ytit='Systematic Error (m/s)',charsize=1.5
a=WHERE(vtest ge -6500 and vtest le 6500)
res=poly_fit(vtest[a],systematicerror[a],1,yfit=y)
oplot,vtest[a],y,col=180
plot,vtest,sigmat1,xst=1,xrange=[-6500,6500],thick=3,tit='!17',xtit='Input Velocity (m/s)',ytit='Standard Deviation (m/s)',charsize=1.5
;OPLOT,vtest,sigmat,col=180
DEVICE,/CLOSE

SET_PLOT,'X'
READ,pause


END
