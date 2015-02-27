;------------------------------------------------------------------------
;
; THIS PROGRAM RETURNS THE DOPPLER VELOCITY ERROR DUE TO PHOTON NOISE
; IN m/s, FOR ntest2 INPUT VELOCITIES, FOR THE HMI FILTER SYSTEM
; (5 LYOT ELEMENTS + 2 MICHELSONS)
; THE DATA USED ARE THE RELATIVE PHASE MAPS OF THE LYOT AND MICHELSON
; NO CONTRAST MAP IS ACCOUNTED FOR
; ELEMENTS,  PROVIDED BY RICHARD SHINE (LMSAL)
;
; AUTHOR: S. COUVIDAT, BASED ON A CODE BY J. SCHOU
;
;------------------------------------------------------------------------


FUNCTION i_f, C, fac
; C is REAL (no imaginary part)
; C must be periodic (C[0]=C[N-1])
; C must have an even number of points
; the sampling rate must be uniform
; and initial and interpolated domain must cover THE SAME RANGE (NT=N0T0)                
                                                                
nt0 = DOUBLE(N_ELEMENTS(C))
nt  = nt0*fac
C_i = DBLARR(nt)
G   = DCOMPLEXARR(nt)
                                                                                
C2              = FFT(C,/DOUBLE)
G[0:nt0/2]      = C2[0:nt0/2]
G[nt0/2+1:nt/2] = DCOMPLEX(0,0)
G[nt/2+1:nt-1]  = CONJ(REVERSE(G[1:nt/2-1])) ; function is real
C_i=FFT(G,/INVERSE,/DOUBLE)
RETURN,DOUBLE(C_i)
                                                                                
END

PRO lineprofile,X,A,F,pder

A[0] = ABS(A[0])
A[2] = ABS(A[2])

F = 1.d0-A[0]*exp(-(X-A[1])^2.d0/A[2]^2.d0)

pder = [ [-exp(-(X-A[1])^2.d0/A[2]^2.d0)], [-A[0]*exp(-(X-A[1])^2.d0/A[2]^2.d0)*(-2.d0/A[2]^2.d0*(A[1]-X))], [-A[0]*exp(-(X-A[1])^2.d0/A[2]^2.d0)*(2.d0*(A[1]-X)^2.d0/A[2]^3.d0)]]

END


PRO HMI2

algo = 0 ;algo=0 for MDI-like algorithm, algo=1 for least-squares fit

; PARAMETERS OF THE SIMULATION
;-----------------------------
nx     = 200 ; the relative phase maps are rebined into nx*nx maps
ntest2 = 5   ; number of input velocities for which we want the error
loc    = FLTARR(ntest2) ; test input velocities
loc[0] = 37 ;-6.5 km/s
loc[1] = 38 ;-6 km/s
loc[2] = 50 ; 0
loc[3] = 62 ; 6 km/s
loc[4] = 63 ; 6.5 km/s
nlam   = 1550 ; number of points that sample the HMI filter transmission profiles
dlam   = 4/1d3; resolution on these profiles
lam    = (DINDGEN(nlam)-(nlam-1)/2.)*dlam

; COMMON VARIABLES
;-----------------------------
erreur = FLTARR(nx,nx,2*(ntest2+1)) ; will contain the errors on the input velocities in m/s
                                ; and the zero-velocity error
; Iron line property (from Graham et al., see HMI web page) 
lam0   = 6173.d0              ; HMI Fe line in Angstrom
fdepth = 0.623/0.53           ; depth of the line according to Yang's web page
fwidth = 0.10/0.12            ; width of the line

well   = 1.25d5               ; Set full well depth for the CCDs (or actually the maximum exposure level) in electrons
ntune  = 6                    ; Set number of tuning positions (NUMBER OF HMI FILTERS)
inttune= 5./2                 ; Set number of tuning positions over wnarrow. (1/spacing)
                              ; WARNING: Don't mess with this unless you know what you are doing.
FSRmich= 0.189D0              ; Free Spectral Range of the narrowest Michelson for MDI at the wavelength 6768 A 
wnarrow= FSRmich*lam0/6768.0  ; FSR of the narrow-band Michelson

; PARAMETERS OF THE LYOT AND MICHELSON ELEMENTS
;----------------------------------------------
lyotw0 = 8*wnarrow   
lyotw  = lyotw0*[0.5,1,2,4,8] ; elements E1, E2, E3, E4, & E5
nlyot  = N_ELEMENTS(lyotw)
mlist  = [1,2]               
nmich  = n_elements(mlist)
wmich  = mlist*wnarrow        ; elements M2 (NB Michelson), & M1 (WB Michelson)

; TUNING POSITIONS
;-----------------
dtune  = wnarrow/inttune      ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
tune   = DBLARR(ntune)
tune   =(-(ntune-1)/2.+dindgen(ntune))*dtune 

; SIZE OF THE CONE OF RAYS OVER WHICH THE LIGHT IS INTEGRATED
; WHEN IT REACHES A CCD PIXEL 
;------------------------------------------------------------
largeur    = DBLARR(5)          ; the sizes are provided by Jesper Schou
largeur[0] = 2.01d0             ; radius (in mm) of the beam on the E1 Lyot element aperture
largeur[1] = 0.85d0             ; E2
largeur[2] = 1.16d0             ; E3
largeur[3] = 1.39d0             ; E4
largeur[4] = 1.59d0             ; E5
largeurM   = DBLARR(2)          ; idem for the Michelsons
largeurM[1]= 4.63               ; Wide Band Michelson (at the reflection off the solid leg)
largeurM[0]= 8.69               ; Narrow Band Michelson (at the reflection off the solid leg)

; LYOT ELEMENTS
;------------------------------------------------------------
e1 = MRDFITS('element1phase.fits')*!pi/180.0
a  = WHERE(e1[*,400] NE 0.0)
b  = WHERE(e1[400,*] NE 0.0)
e1 = e1[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]
e2 = MRDFITS('element2phase.fits')*!pi/180.0
a  = WHERE(e2[*,400] NE 0.0)
b  = WHERE(e2[400,*] NE 0.0)
e2 = e2[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]
e3 = MRDFITS('element3phase.fits')*!pi/180.0
a  = WHERE(e3[*,400] NE 0.0)
b  = WHERE(e3[400,*] NE 0.0)
e3 = e3[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]
e4 = MRDFITS('element4phase.fits')*!pi/180.0
a  = WHERE(e4[*,400] NE 0.0)
b  = WHERE(e4[400,*] NE 0.0)
e4 = e4[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]
e5 = MRDFITS('element5phase.fits')*!pi/180.0
a  = WHERE(e5[*,400] NE 0.0)
b  = WHERE(e5[400,*] NE 0.0)
e5 = e5[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]


e1[*,*] = 5.d0*!dpi/180.d0;1.d-10 ; !!!!! WARNING !!!!!!
e5[*,*] = 5.d0*!dpi/180.d0;1.d-10
e2[*,*] = 5.d0*!dpi/180.d0;1.d-10
e3[*,*] = 5.d0*!dpi/180.d0;1.d-10
e4[*,*] = 5.d0*!dpi/180.d0;1.d-10
;contrast=[.9,.9,.9,.9,.9,.9,.9]
contrast=[1.,1.,1.,1.,1.,1.,1.]

; MICHELSON ELEMENTS
;--------------------------------
;M2 = MRDFITS('michN00_14jul2005.fits')*!pi/180.0 ;NARROW-BAND MICHELSON (M2) at 23 degrees Celsius
M2 = MRDFITS('michN00_18jul2005.fits')*!pi/180.0 ;NARROW-BAND MICHELSON (M2) at 30 degrees Celsius
a  = WHERE(M2[*,430] NE 0.0)
b  = WHERE(M2[430,*] NE 0.0)
M2 = M2[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]
;M1 = MRDFITS('michB00_19jul2005.fits')*!pi/180.0 ;BROAD-BAND MICHELSON (M1) at 23 degrees Celsius
M1 = MRDFITS('michBE00_20jul2005.fits')*!pi/180.0 ;BROAD-BAND MICHELSON (M1) at 30 degrees Celsius
; WARNING: orientations of the figures may vary between 23 and 30
; degrees Celsius
a  = WHERE(M1[*,430] NE 0.0)
b  = WHERE(M1[430,*] NE 0.0)
M1 = M1[a[0]:a[N_ELEMENTS(a)-1],b[0]:b[N_ELEMENTS(b)-1]]


M2[*,*]=5.d0*!dpi/180.d0;1.d-10 ; !!! WARNING !!!
M1[*,*]=5.d0*!dpi/180.d0;1.d-10


; WE REBIN THE RELATIVE-PHASE MAPS
;---------------------------------
e        = DBLARR(5,nx,nx)
e[0,*,*] = CONGRID(e1,nx,nx)
e1       = 0.0
e[1,*,*] = CONGRID(e2,nx,nx)
e2       = 0.0
e[2,*,*] = CONGRID(e3,nx,nx)
e3       = 0.0
e[3,*,*] = CONGRID(e4,nx,nx)
e4       = 0.0
e[4,*,*] = CONGRID(e5,nx,nx)
e5       = 0.0
dx       = 30.d0/DOUBLE(nx) ; we assume that all the relative-phase maps have a diameter of 30 mm
ny       = DBLARR(5)
FOR i=0,4 DO ny[i]    = CEIL(largeur[i]/dx*2.d0)

sizeM1   = 30.06 ; we assume the relative-phase disc has a diameter of 30.06 mm
fac1     = 1.20  ; actual number: 1.19772; ratio between the size of the solar image in M1 and in the Lyot:
                 ; 25.2/20.8
sizeM2   = 29.58 ; we assume the relative-phase disc has a diameter of 29.58 mm
fac2     = 2.06  ; actual number: 2.06557; ratio between the size of the solar image in M2 and in the Lyot:
                 ; 25.2/12.2

M1       = CONGRID(M1,nx*fac1,nx*fac1)
dx1      = sizeM1*fac1/DOUBLE(N_ELEMENTS(M1[0,*]))
nyM1     = CEIL(largeurM[1]*fac1/dx1*2.d0)
M2       = CONGRID(M2,nx*fac2,nx*fac2)
dx2      = sizeM2*fac2/DOUBLE(N_ELEMENTS(M2[0,*]))
nyM2     = CEIL(largeurM[0]*fac2/dx2*2.d0)

maxi     = ny[0]
maxlargeur = largeur[0]

; Fe LINE PARAMETERS
;-------------------

GOTO,skipo
; THEORETICAL PROFILE FROM TED'S APPROXIMATE EQUATION
npar  = 3    ; Number of line parameters
nlamp = 1223 ; Set wavelength points for line profile
lamp  = -0.3D0+DINDGEN(nlamp)*0.6D/(nlamp-1)

t      = -2.7726*(lamp/0.125/fwidth)^2 ; Line profile from Ted's approximate equation
t      = -25. + (t+25.)*(t GT -25.)
line   = 1-0.6*exp(t)*fdepth
skipo:

GOTO,skipo2
; LINE PROFILE FROM THE Kitt Peak NSO Atlas
RESTORE,'atlas_profile.sav'
line  = i_f(pr1[218-20:218+20],32)
line  = line/MAX(line)
nlamp = N_ELEMENTS(line)
line2 = FLTARR(nlamp+1) ; to have an odd number
line2[0:nlamp-1]=line
line2[nlamp]=line2[nlamp-1]
line  = line2
line2 = 0.0
nlamp = nlamp+1
a     = WHERE(line EQ MIN(line))
line  = SHIFT(line,nlamp/2-a[0])
mini2 = MIN(lb1[218-20])
maxi2 = MAX(lb1[218+20])
mili2 = (mini2+maxi2)/2.0
lamp  = (mini2-mili2)+DINDGEN(nlamp)*(maxi2-mini2)/(nlamp-1)
lb1   = 0.0
pr1   = 0.0
skipo2:

; LINE PROFILE FROM ROGER ULRICH'S WEBSITE
OPENR,1,'Ulrich_Fe_0.txt'
;OPENR,1,'Ulrich_Fe_45.txt' ;center-to-limb angle = 45 degrees
;OPENR,1,'Ulrich_Fe_60.txt'
;data = FLTARR(2,100) ; for angle=0
data = FLTARR(2,98)  ; for angle=45 or 60
READF,1,data
lamp = REFORM(data[0,*])
nlamp= N_ELEMENTS(lamp)
line = REFORM(data[1,*])
close,1


; Convert wavelength to velocity
;--------------------------------------------
dlamdv = lam0/3d8
dvdlam = 1/dlamdv
vel    = lam*dvdlam   ; DOPPLER shift in m/s
dvel   = dlam*dvdlam
dvtune = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune  = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s


;------------------------------------------------------------------------------------------------
; BLOCKER FILTER PARAMETERS & PROFILE
;------------------------------------------------------------------------------------------------

blocker = DBLARR(nlam)                 ; contains the smooth blocker profile
blocker(WHERE(ABS(lam) LE lyotw0)) = 1 ; Simple boxcar blocker to arbitrary null in Lyot
bfac    = 1.6
q       = readx('smooth.txt')
blocker = INTERPOL(q(1,*)/100,bfac*(q(0,*)-6170),lam,/spline)

;------------------------------------------------------------------------------------------------
; FOR THE DOPPLER VELOCITY DETERMINATION
; WE GENERATE ntest DOPPLER SHIFTED LINE PROFILES
;------------------------------------------------------------------------------------------------

ntest         = 101 ; Set test velocities: we will use ntest different velocities, resulting
                                ; in ntest Doppler shifted line profiles
dvtest        = 500.; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1)/2)
lines         = DBLARR(nlam,ntest)
dlinesdv      = DBLARR(nlam,ntest)
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
;q2            = [MIN(lamp)-10,lamp,MAX(lamp)+10] ; Following two lines take care of continuum outside of where profile is given
q2            = [MIN(lamp)-4,lamp,MAX(lamp)+4] ; !!!!MODIFICATION
FOR i = 0,ntest-1 DO lines(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)

;------------------------------------------------------------------------------------------------
; LYOT FILTER PROFILE
;------------------------------------------------------------------------------------------------

clyot   = 2*!dpi/lyotw
cmich   = 2*!dpi/wmich
lyot    = blocker ; Set Lyot+Blocker filter profile. No Doubled elements.
filters = DBLARR(nlam,ntune) ; contain the filter profiles (Blocker+Lyot+Michelsons)

; REFERENCE TRANSMISSION PROFILES FOR THE 6 HMI FILTERS: CONTRASTS=1, PHASES=0
;------------------------------------------------------------------------------
FOR i     = 1,nlyot-1 DO lyot = lyot*(1.0+COS(clyot[i]*lam))/2.0 ; reference transmission profile
FOR itune = 0,ntune-1 DO filters[*,itune] =lyot*(1.0+COS(clyot[0]*(lam+tune[itune])))/2.0
FOR itune = 0,ntune-1 DO FOR i  = 0,nmich-1 DO filters[*,itune] = filters[*,itune]*(1+COS(cmich[i]*(lam+tune[itune])))/2



IF (algo EQ 0) THEN BEGIN
; REFERENCE DOPPLER VELOCITY ERROR COMPUTED FROM A MDI-LIKE ALGORITHM
; WITH 6 FILTERS. WE APPLY A FACTOR SQRT(2) TO SIMULATE FOR THE 
; LCP AND RCP MEASUREMENTS
; WE ALSO PROVIDE THE RESULTS FOR 2 AVERAGED MDI-LIKE ALGORITHMS,
; USING THE PHASES OF THE 1ST AND 2ND FOURIER COMPONENTS OF THE LINE
;--------------------------------------------------------------------
inten     = TRANSPOSE(filters)#lines
texp      = well/MAX(inten)
x         = 2.0*!dpi*(-(ntune-1)/2.+DINDGEN(ntune))/ntune
c1        = REFORM(COS(x)#inten)  ; coefficient a1
s1        = REFORM(SIN(x)#inten)  ; coefficient b1 (for the sine)
c2        = REFORM(COS(2.d0*x)#inten)  ; coefficient a2
s2        = REFORM(SIN(2.d0*x)#inten)  ; coefficient b2 (for the sine)
pv1       = dvtune*inttune*2.0
pv2       = dvtune*inttune
phi1      = ATAN(-s1,-c1)
phi2      = ATAN(-s2,-c2)
vel1      = phi1*pv1/2.d0/!dpi
vel2      = phi2*pv2/2.d0/!dpi
vel1a     = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; The look-up table will be (vtest,vel1) or (vtest,vel1a)
vel2a     = (vel2-vtest+10.5d0*pv2) MOD pv2-pv2/2.d0+vtest
read,pause
nt        = 8000 ; number of random intensities
randomge  = FLTARR(nt,ntune)
iseed     = 1000l
FOR i     = 0,ntune-1 DO randomge[*,i] = RANDOMN(iseed,nt)
FOR k = 0,ntest2-1 DO BEGIN
    itest = loc[k]
    vt    = vtest(itest)
    i0    = inten(*,itest)*texp
    intent= DBLARR(ntune,nt)
    FOR ii= 0,ntune-1 DO intent(ii,*)=i0(ii)+SQRT(i0(ii))*randomge[*,ii]
    c1    = REFORM(COS(x)#intent)
    s1    = REFORM(SIN(x)#intent)
    c2    = REFORM(COS(2.d0*x)#intent)
    s2    = REFORM(SIN(2.d0*x)#intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel2t = phi2t*pv2/2.d0/!dpi
    vel1at= (vel1t-vt+10.5d0*pv1) MOD pv1-pv1/2.d0+vt
    vel2at= (vel2t-vt+10.5d0*pv2) MOD pv2-pv2/2.d0+vt
    v1tb  = INTERPOL(vtest,vel1a,vel1at) ; measured velocity
    v2tb  = INTERPOL(vtest,vel2a,vel2at)
;   v1tb  = INTERPOL(vtest[20:80],vel1[20:80],vel1t) ; measured velocity
    sigmat1b = SQRT(REBIN((v1tb-vtest[itest])^2.0,1)) ; error on the measured velocity   
    sigmat2b = SQRT(REBIN((v2tb-vtest[itest])^2.0,1))
    PRINT,'reference error at (MDI-like)',vt,sigmat1b/SQRT(2.d0)
;   PRINT,'reference error at (weighted MDI-like)',vt,1.d0/SQRT((1.d0/sigmat1b^2.d0+1.d0/sigmat2b^2.d0))/SQRT(2.d0) 
    PRINT,'reference error at (averaged MDI-like)',vt,SQRT(REBIN(((v1tb+v2tb)/2.d0-vtest[itest])^2.0,1))/SQRT(2.d0)
ENDFOR
ENDIF ELSE BEGIN

; USING A LEAST-SQUARE ALGORITHM
inten     = TRANSPOSE(filters)#lines
texp      = well/MAX(inten)
weight    = DBLARR(6)+1.d0
vel1      = DBLARR(ntest)

FOR k = 0,ntest-1 DO BEGIN 
A       = [0.75,0.0,0.05]
res     = CURVEFIT(tune,inten[*,k]/MAX(inten[*,k]),weight,A,FUNCTION_NAME='lineprofile',itmax=1000,TOL=1.d-6,/DOUBLE,CHISQ=chi1)
B       = [0.75,-0.123,0.05]
res     = CURVEFIT(tune,inten[*,k]/MAX(inten[*,k]),weight,B,FUNCTION_NAME='lineprofile',itmax=1000,TOL=1.d-6,/DOUBLE,CHISQ=chi2)
C       = [0.75, 0.123,0.05]
res     = CURVEFIT(tune,inten[*,k]/MAX(inten[*,k]),weight,C,FUNCTION_NAME='lineprofile',itmax=1000,TOL=1.d-6,/DOUBLE,CHISQ=chi3)
IF(chi1 LE chi2 AND chi1 LE chi3) THEN vel1[k] = A[1]*dvdlam
IF(chi2 LE chi1 AND chi2 LE chi3) THEN vel1[k] = B[1]*dvdlam
IF(chi3 LE chi1 AND chi3 LE chi2) THEN vel1[k] = C[1]*dvdlam

ENDFOR
plot,vtest,vel1-vtest,xrange=[-7500,7500],xst=1,yst=1

nt        = 8000 ; number of random intensities
randomge  = FLTARR(nt,ntune)
iseed     = 1000l
FOR i     = 0,ntune-1 DO randomge[*,i] = RANDOMN(iseed,nt)
FOR k = 0,ntest2-1 DO BEGIN
    itest = loc[k]
    vt    = vtest(itest)
    i0    = inten(*,itest)*texp
    intent= DBLARR(ntune,nt)
    vel2  = DBLARR(nt)
    FOR ii= 0,ntune-1 DO  intent(ii,*)=i0(ii)+SQRT(i0(ii))*randomge[*,ii]
    FOR ii= 0,nt-1 DO BEGIN
        A = [0.75,0.0,0.05]
        res = CURVEFIT(tune,intent[*,k]/MAX(intent[*,k]),weight,A,FUNCTION_NAME='lineprofile',itmax=400,/DOUBLE)
        vel2[ii]=A[1]*dvdlam
        ;vel2[ii]=INTERPOL(vtest[35:65],vel1[35:65],vel2[ii])
    ENDFOR
    PRINT,'reference error at (MDI-like)',vt,SQRT(REBIN((vel2-vtest[itest])^2.0,1))/SQRT(2.d0)
ENDFOR
ENDELSE

discf    = DBLARR(7,maxi,maxi)
discc    = DBLARR(2,maxi,maxi)
disc     = SHIFT(DIST(maxi,maxi)*dx,maxi/2,maxi/2)
aa       = WHERE(disc le maxlargeur,COMPLEMENT=bb)
disc[aa] = 1.d0
na       = FLOAT(N_ELEMENTS(aa))
disc[bb] = 0.d0

FOR i=0,nx-1 DO BEGIN
    FOR j=0,nx-1 DO BEGIN
        
        IF(SQRT( (i-nx/2.)^2.0+(j-nx/2.)^2.0 )*dx LE 12.5) THEN BEGIN ; 12.6 mm is the radius of the image on a CCD

        PRINT,i,j
        ; NON-TUNABLE AND TUNABLE ELEMENTS OF THE LYOT
        ;---------------------------------------------
        FOR k=0,4 DO BEGIN
            IF(ny[k] MOD 2 EQ 0) THEN discp = e[k,i-ny[k]/2:i+ny[k]/2-1,j-ny[k]/2:j+ny[k]/2-1] ELSE discp = e[k,i-ny[k]/2:i+ny[k]/2,j-ny[k]/2:j+ny[k]/2]
            discf[k,*,*] = CONGRID(discp,1,maxi,maxi)
            discf[k,*,*] = disc[*,*]*discf[k,*,*]
        ENDFOR
        ; MICHELSONS
        ;---------------------------------------------
        i2 = sizeM1*fac1/dx1/2.0-30.0/dx1/2.0+i*dx/dx1
        j2 = sizeM1*fac1/dx1/2.0-30.0/dx1/2.0+j*dx/dx1
        IF(nyM1 MOD 2 EQ 0) THEN discp = M1[i2-nyM1/2:i2+nyM1/2-1,j2-nyM1/2:j2+nyM1/2-1] ELSE discp = M1[i2-nyM1/2:i2+nyM1/2,j2-nyM1/2:j2+nyM1/2]
        discM1  = CONGRID(discp,maxi,maxi)
        discM1  = disc[*,*]*discM1[*,*]
        i3 = sizeM2*fac2/dx2/2.0-30.0/dx2/2.0+i*dx/dx2
        j3 = sizeM2*fac2/dx2/2.0-30.0/dx2/2.0+j*dx/dx2
        IF(nyM2 MOD 2 EQ 0) THEN discp = M2[i3-nyM2/2:i3+nyM2/2-1,j3-nyM2/2:j3+nyM2/2-1] ELSE discp = M2[i3-nyM2/2:i3+nyM2/2,j3-nyM2/2:j3+nyM2/2]
        discM2  = CONGRID(discp,maxi,maxi)
        discM2  = disc[*,*]*discM2[*,*]
        
        temp  = DBLARR(nlam,ntune)
        lyot2 = DBLARR(nlam,ntune)
        FOR a=0,ny[0]-1 DO BEGIN
            FOR b=0,ny[0]-1 DO BEGIN
                IF(discf[0,a,b] NE 0.0) THEN BEGIN
                lyot2[*,0]= blocker[*]
                FOR k     = 1,4 DO lyot2[*,0] = lyot2[*,0]*(1.0+contrast[k+2]*COS(clyot[k]*lam[*]+discf[k,a,b]))/2.0
                FOR itune = 1,ntune-1 DO lyot2[*,itune] = lyot2[*,0]
                FOR itune = 0,ntune-1 DO lyot2[*,itune] = lyot2[*,itune]*(1.0+contrast[2]*COS(clyot[0]*(lam[*]+tune[itune])+discf[0,a,b]))/2.0*(1.0+contrast[0]*COS(cmich[0]*(lam[*]+tune[itune])+discM2[a,b]))/2*(1.0+contrast[1]*COS(cmich[1]*(lam[*]+tune[itune])+discM1[a,b]))/2
                temp[*,*] = temp[*,*] + lyot2[*,*]
                ENDIF
            ENDFOR
        ENDFOR
        nonnul= DOUBLE(N_ELEMENTS(WHERE(discf[0,*,*]) NE 0.0))
        temp  = temp/nonnul
   
        ; VELOCITY ERROR FOR DIFFERENT INPUT VELOCITIES
        ;----------------------------------------------
        inten2  = TRANSPOSE(temp)#lines
        c1      = REFORM(COS(x)#inten2) ; coefficient a1
        s1      = REFORM(SIN(x)#inten2) ; coefficient b1 (for the sine)
        c2      = REFORM(COS(2.d0*x)#inten2) ; coefficient a2
        s2      = REFORM(SIN(2.d0*x)#inten2) ; coefficient b2 (for the sine)
        phi1    = ATAN(-s1,-c1)
        phi2    = ATAN(-s2,-c2)
        vel1b   = phi1*pv1/2.d0/!dpi ; velocity measurement with an MDI-like algorithm
        vel2b   = phi2*pv2/2.d0/!dpi
        vel1bb  = (vel1b-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; corrected look-up tables
        vel2bb  = (vel2b-vtest+10.5d0*pv2) MOD pv2-pv2/2.d0+vtest

        FOR k=0,ntest2-1 DO BEGIN
            itest= loc[k]
            vt   = vtest[itest]
            i0   = inten2[*,itest]*texp
            intent=DBLARR(ntune,nt)
            FOR ii= 0,ntune-1 DO intent[ii,*]=i0[ii]+SQRT(i0[ii])*randomge[*,ii]
            c1   = REFORM(COS(x)#intent)
            s1   = REFORM(SIN(x)#intent)
            c2   = REFORM(COS(2.d0*x)#intent)
            s2   = REFORM(SIN(2.d0*x)#intent)
            phi1t = ATAN(-s1,-c1)
            phi2t = ATAN(-s2,-c2)
            vel1t = phi1t*pv1/2.d0/!dpi
            vel2t = phi2t*pv2/2.d0/!dpi
            vel1at= (vel1t-vt+10.5d0*pv1) MOD pv1-pv1/2.d0+vt
            vel2at= (vel2t-vt+10.5d0*pv2) MOD pv2-pv2/2.d0+vt
            v1tb = INTERPOL(vtest,vel1bb,vel1at) ; measured velocity obtained with the corrected look-up table
            v1uncorrected = INTERPOL(vtest,vel1a,vel1at) ; !!! WARNING !!!measured velocity obtained with the uncorrected look-up table
;           v1tb = INTERPOL(vtest[20:80],vel1b[20:80],vel1t) ; measured velocity
            v2tb = INTERPOL(vtest,vel2bb,vel2at) ; measured velocity
            v2uncorrected = INTERPOL(vtest,vel2a,vel2at) ;!!! WARNING !!!
            sigmat1b = SQRT(REBIN((v1tb-vtest[itest])^2.0,1)) ; error on the measured velocity   
            ;sigmat2b = SQRT(REBIN((v2tb-vtest[itest])^2.0,1)) ; error on the measured velocity   
            erreur[i,j,k]=sigmat1b/SQRT(2.0);/sigmat1a[itest]
            erreur[i,j,k+ntest2+1]=SQRT(REBIN(((v1tb+v2tb)/2.0-vtest[itest])^2.0,1))/SQRT(2.d0) ;1.d0/SQRT( (1.d0/sigmat1b^2.d0)+(1.d0/sigmat2b^2.d0) )/SQRT(2.0);for the weighted MDI-like algorithm
            IF(itest EQ 50) THEN print,i,j,MEAN((v1uncorrected+v2uncorrected)/2.d0),SIGMA((v1uncorrected+v2uncorrected)/2.d0)/SQRT(2.d0) ; !!! WARNING !!!
        ENDFOR

        ; ZERO-VELOCITY ERROR
        ;---------------------
        itest= 50
        c1   = REFORM(COS(x)#inten2)
        s1   = REFORM(SIN(x)#inten2)
        c2   = REFORM(COS(2.d0*x)#inten2)
        s2   = REFORM(SIN(2.d0*x)#inten2)  
        phi1t = ATAN(-s1,-c1)
        phi2t = ATAN(-s2,-c2)
        vel1t = phi1t*pv1/2.d0/!dpi
        vel2t = phi2t*pv2/2.d0/!dpi
        erreur[i,j,ntest2]=vel1t[itest] ;/sigmat1a[itest]
        erreur[i,j,2*(ntest2+1)-1]=(vel1t[itest]+vel2t[itest])/2.d0 ;/sigmat1a[itest]

        PRINT,TRANSPOSE(erreur[i,j,*])

    ENDIF ELSE erreur[i,j,*]=0.0

    ENDFOR
ENDFOR

;SAVE,erreur,FILE='lyot_NSOAtlas.bin'
;SAVE,erreur,FILE='lyot_Roger_60.bin'
;SAVE,erreur,FILE='lyot_Roger_0_30deg.bin'
set_plot,'ps'
device,file='yo.ps',bits=24
a=WHERE(erreur EQ 0.0,COMPLEMENT=b)
erreur[a]=-100.0
FOR i=0,2*(ntest2+1)-1 DO BEGIN
mini=MIN(ABS(erreur[*,*,i]))
IF(i EQ ntest2 OR i EQ 2*(ntest2+1)-1) THEN mini=MIN(erreur[*,*,i])
maxi=MAX(erreur[*,*,i])
tvim,erreur[*,*,i],/scale,tit='!17',xtit='x (mm)',ytit='y (mm)',xrange=[0,dx*200],yrange=[0,dx*200],stit='Error (m/s)',range=[mini,maxi],pcharsize=1.5,/rct
ENDFOR
device,/close
END
