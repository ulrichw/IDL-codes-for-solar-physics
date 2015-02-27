; a lighter version of HMI.pro to compute errors on the Doppler
; velocity (systematic errors + photon noise errors) resulting from a
; wrong calibration
; computes error for:
; ---MDI-like algorithm with 6 positions, using the first 2 Fourier
; coefficients, and the RCP+LCP filtergrams
; ---Least-Squares Fit with 4 parameters, fitting RCP and LCP separately
;
; WARNING: THE SIGN CONVENTION FOR THE VELOCITY IS REVERSED FROM WHAT
; IS USED IN HMI_observables.c
; POSITIVE VELOCITIES ARE BLUESHIFT, I.E. VELOCITIES TOWARD THE OBSERVER
;
; FUNCTION THAT COMPUTES THE LIMB DARKENING
; FROM H. NECKEL (2005)
;------------------------------------------------------------------------

FUNCTION limb,long

mu0  = COS(long*!dpi/180.d0)
mu   = [1.d0,mu0,mu0^2.d0,mu0^3.d0,mu0^4.d0,mu0^5.d0]
lamref=0.61733433d0

A    = DBLARR(6)
A[0] =  0.75267d0-0.265577d0/lamref
A[1] =  0.93874d0+0.265577d0/lamref-0.004095d0/lamref^5.d0
A[2] = -1.89287d0+0.012582d0/lamref^5.d0
A[3] =  2.42234d0-0.017117d0/lamref^5.d0
A[4] = -1.71150d0+0.011977d0/lamref^5.d0
A[5] =  0.49062d0-0.003347d0/lamref^5.d0

ratio= TOTAL(A*mu); ratio intensity at long over intensity at the center


RETURN,ratio

END



PRO solarline,X,A,F,pder

COMMON filtre,filters2,lam2,dlamdv

A[0]   = ABS(A[0])
A[1]   = ABS(A[1])
A[3]   = ABS(A[3])

F1     = A[0]*(1.d0-A[1]*EXP(-(lam2+A[2]*dlamdv)^2.d0/A[3]^2.d0))

F      = filters2#F1
Fd     = filters2#EXP(-(lam2+A[2]*dlamdv)^2.d0/A[3]^2.d0)
Fdd    = filters2#((lam2+A[2]*dlamdv)*EXP(-(lam2+A[2]*dlamdv)^2.d0/A[3]^2.d0))
Fddd   = filters2#((lam2+A[2]*dlamdv)^2.d0*EXP(-(lam2+A[2]*dlamdv)^2.d0/A[3]^2.d0))


pder   = [ [F/A[0]],[-A[0]*Fd],[A[0]*A[1]*2.d0*dlamdv/A[3]^2.d0*Fdd],[-A[0]*A[1]*2.d0/A[3]^3.d0*Fddd] ]

END

;-----------------------------------------------------------------
;
; MAIN PROGRAM
;
;-----------------------------------------------------------------


PRO HMIcalibrationerror

COMMON filtre,filters2,lam2,dlamdv

lam0    = 6173.3433d0
well    = 163000.d0;200000.d0;150000.d0 ; CCD full well

ntune   = 6      ; Set number of tuning positions (NUMBER OF FILTERS)
ntunep  = 11     ; All positions for the purpose of plotting.
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)
dvdb    = 1./2.13d12*lam0*2.5d0*299792458.d0  ; conversion from magnetic field to Doppler velocity

; PARAMETERS
;-----------------------------------------------------------------------

long       = 0. ; longitude (0=solar center, 90=limb)
field      = 0. ; magnetic field in Gauss
; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
phase      = [0.0d0,0.d0,0.0d0,-6.60914,-0.154691,-8.61043,-5.15019]*!dpi/180.d0 
contrast   = [0.976811,0.999908,0.976513,0.986981,0.965252,0.989889,0.999402]  
; PHASES AND CONTRASTS WE USE TO OBTAIN THE LOOK-UP TABLE
phase0     = [0.0d0,0.d0,0.0d0,-6.60914,-0.154691,-8.61043,-5.15019]*!dpi/180.d0 
contrast0  = [0.976811,0.999908,0.976513,0.986981,0.965252,0.989889,0.999402]  

;-----------------------------------------------------------------------

;wmich      = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]          ; FSRs of the tunable elements
wmich      = [0.169,0.337d0,0.695d0]
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; theoretical FSRs non-tunable elements
lyotwseb   = [1.407d0,2.779d0,5.682d0,11.354d0] ; real FSRs


dtune  = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune


nlamp  = 30061 ; Set wavelength points for line profile
lamp   = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

ROGER  = 1
; SOLAR LINE
IF(ROGER EQ 1) THEN BEGIN
    OPENR,1,'Ulrich_Fe_0.txt'
   ;OPENR,1,'Ulrich_Fe_45.txt'
   ;OPENR,1,'Ulrich_Fe_60.txt'
    roger  = FLTARR(2,98)
    READF,1,roger
    CLOSE,1
    nroger = 98
    rlam   = REFORM(roger(0,*))
    rp     = REFORM(roger(1,*))
    rp[97] = 1.d0
    rlam   = [-10.d0,rlam,10.d0]
    rp     = [1.d0,rp,1.d0]
    line   = INTERPOL(rp,rlam,lamp)
ENDIF ELSE BEGIN
; SOLAR LINE FROM AN ATLAS (DON'T REMEMBER WHICH ONE)
    data   = FLTARR(30002)
    wavelength = FLTARR(15001)
    spectrum   = FLTARR(15001)
    OPENR,1,'spectre_6160_30.dat'
    READF,1,data
    CLOSE,1
    FOR i=0,15000 DO BEGIN
        wavelength[i]=data[i*2]
        spectrum[i]=data[i*2+1]
    ENDFOR
    wavelength= wavelength-6173.3433d0
    spectrum  = spectrum/MAX(spectrum)
    line   = INTERPOL(spectrum,wavelength,lamp)
    a      = WHERE(lamp LE -7.8) ; to get rid of the line that is half-cut
    line[a]= 1.d0
ENDELSE
; MODEL OF LINE FOR THE LEAST-SQUARES FIT
width  = 0.098d0/2.d0/SQRT(ALOG(2.d0));0.102d0/2.d0/SQRT(ALOG(2.d0)) ; w in A^2 (part added by seb)
depth  = 0.66d0;0.62d0
line2  = 1.d0-depth*EXP(-(lamp)^2.d0/width^2.d0)
dlinedd= -EXP(-lamp^2.d0/width^2.d0)
dlinedi= line2
dlinedw=-depth*2.d0*lamp^2.d0/width^3.d0*EXP(-lamp^2.d0/width^2.d0)

; GRID WE WANT IN WAVELENGTH
nlam          = 90000.;38500.
dlam          = 1.d0/3.5d3;1.d0/1.75d3
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = lam0/299792458.d0
dvdlam        = 1.d0/dlamdv
vel           = lam*dvdlam   ; DOPPLER shift in cm/s
dvel          = dlam*dvdlam
dvtune        = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune         = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; BLOCKER FILTER
RESTORE,'frontwindow3.bin' ; front window
blocker       = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
q             = READFITS('blocker11.fits')
blockerseb    = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6-6173.3433d0,lam) ; real blocker
blocker       = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6-6173.3433d0,lam) ; blocker used for look-up tables

; INPUT VELOCITIES
ntest         = 501 ; Set test velocities
dvtest        = 27.75d0;dlam/dlamdv;100.d0 ; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
lines         = DBLARR(nlam,ntest)
lines2        = lines
lines3        = lines
;dlinesdv      = DBLARR(nlam,ntest) 
;dlinesdd      = DBLARR(nlam,ntest) ; Get derivative with respect to depth
;dlinesdi      = DBLARR(nlam,ntest) ; Get derivative with respect to overall intensity
;dlinesdw      = DBLARR(nlam,ntest) ; Get derivative with respect to linewidth (added by seb)
q2            = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles FROM LINE
FOR i = 0,ntest-1 DO lines[*,i]    = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
q2            = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles FROM LINE 
FOR i = 0,ntest-1 DO lines2[*,i]   = INTERPOL(q1,q2,lam+(vtest[i]+field*dvdb)*dlamdv) ;LCP
FOR i = 0,ntest-1 DO lines3[*,i]   = INTERPOL(q1,q2,lam+(vtest[i]-field*dvdb)*dlamdv) ;RCP

;read,pause

;dv            = 100.d0                     ; Get derivative with respect to v
;FOR i = 0,ntest-1 DO dlinesdv[*,i] =(INTERPOL(q1,q2,lam+(vtest[i]+dv)*dlamdv)-INTERPOL(q1,q2,lam+(vtest[i]-dv)*dlamdv))/2.d0/dv
;q1            = [MAX(dlinedd),dlinedd,MAX(dlinedd)]
;for i = 0,ntest-1 DO dlinesdd[*,i] = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
;q1            = [MAX(dlinedi),dlinedi,MAX(dlinedi)]
;for i = 0,ntest-1 DO dlinesdi[*,i] = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
;q1            = [0.0,dlinedw,0.0]
;for i = 0,ntest-1 DO dlinesdw[*,i] = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)

; NON-TUNABLE PROFILE
lyot          = blocker
lyotseb       = blockerseb
FOR i = 0,3 DO BEGIN
    lyot      = lyot   *(1.d0+contrast0[i+3]*COS(2.d0*!dpi/lyotw[i]   *lam+phase0[i+3]))/2.d0
    lyotseb   = lyotseb*(1.d0+contrast [i+3]*COS(2.d0*!dpi/lyotwseb[i]*lam+phase [i+3]))/2.d0
ENDFOR

; TUNABLE PROFILE
cmich     = 2.d0*!dpi/wmich
filterseb = DBLARR(nlam,ntune)
filters   = DBLARR(nlam,ntune)
FOR itune = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast0[i]*COS(cmich[i]*(lam+tune[i,itune])+phase0[i]))/2.d0
ENDFOR
FOR itune = 0,ntune-1 DO BEGIN
    filterseb[*,itune] = lyotseb
    FOR i  = 0,2 DO BEGIN 
        filterseb[*,itune] = filterseb[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0 ; The Michelsons
    ENDFOR
ENDFOR


; INTENSITIES
inten      = TRANSPOSE(filters)  #lines ; intensities from what we believe are the filters
inten      = inten
;inten2     = TRANSPOSE(filters)  #lines2
inten2     = TRANSPOSE(filterseb)#lines2 ; intensities from actual filters LCP
inten3     = TRANSPOSE(filterseb)#lines3 ; intensities from actual filters RCP
texp       = well/MAX(inten)            ; "exposure" time CALCULATED BEFORE WE APPLY THE LIMB DARKENING FUNCTION


; LIMB DARKENING (AFTER CALCULATION OF TEXP!)
inten      = inten    *limb(long)
inten2     = inten2   *limb(long)
inten3     = inten3   *limb(long)
filters    = filters  *limb(long)
filterseb  = filterseb*limb(long)
well       = well*limb(long) ; CAUSE WE WANT TO KNOW THE ACTUAL INTENSITY AT A GIVEN LONGITUDE


nt         = 7500;200.
iseed      = 1000l
iseed2     = 1002l
randomge   = FLTARR(nt,ntune)
randomge2  = randomge
FOR i      = 0,ntune-1 DO BEGIN
    randomge[*,i]  = RANDOMN(iseed,nt) ;FOR PHOTON NOISE ON LCP
    randomge2[*,i] = RANDOMN(iseed2,nt);FOR PHOTON NOISE ON RCP
ENDFOR

; LINEARIZED INVERSE WITH FITTING:
; LEAST-SQUARES FITTING (ACTUAL FIT) 
;------------------------------------------------------------------------------------------------

;funcname = 'solarline'
;weights  = FLTARR(ntune)+1.d0
;A0       = [1.d0,0.62015688,0.0d,0.062]
;vel1a    = FLTARR(ntest,nt)
;
;aaaa     = WHERE(vtest GE -4000. AND vtest LE 4000.,na)
;
;selection= FINDGEN(10000)*2.+27500
;filters2 = TRANSPOSE(filters[selection,*])
;lam2     = lam[selection]
;
;FOR itest= aaaa[0],aaaa[na-1] DO BEGIN
;    vt   = vtest[itest]
;    i0   = inten[*,itest]*texp
;    intent=DBLARR(ntune,nt)
;    FOR i= 0,ntune-1 DO intent[i,*]=i0[i]+SQRT(i0[i])*randomge[*,i] ; SQRT(i[0]) is photon noise
;    FOR ii=0,nt-1 DO BEGIN
;        A= A0
;        res=CURVEFIT(FINDGEN(6),intent[*,ii]/well,weights,A,FUNCTION_NAME=funcname,TOL=1.d-11,ITMAX=4000,STATUS=stat2)
;        IF(stat2 EQ 1) THEN BEGIN
;            A=A0
;            A[2]=-6500.
;            res=CURVEFIT(FINDGEN(6),intent[*,ii]/well,weights,A,FUNCTION_NAME=funcname,TOL=1.d-11,ITMAX=4000,STATUS=stat2)
;        ENDIF
;        IF(stat2 EQ 1) THEN BEGIN
;            A=A0
;            A[2]=+6500.
;            res=CURVEFIT(FINDGEN(6),intent[*,ii]/well,weights,A,FUNCTION_NAME=funcname,TOL=1.d-11,ITMAX=4000,STATUS=stat2)
;        ENDIF
;        vel1a[itest,ii]= A[2]
;   ENDFOR
;   PRINT,itest,MEAN(vel1a[itest,*]),vt,SIGMA(vel1a[itest,*])
;   PLOT,vel1a[itest,*]
;ENDFOR


;npar     = 4    ; Number of line parameters (3 for Jesper)  
;dintendv = TRANSPOSE(filters)#dlinesdv ; from what we believe are the filters
;dintendd = TRANSPOSE(filters)#dlinesdd
;dintendi = TRANSPOSE(filters)#dlinesdi 
;dintendw = TRANSPOSE(filters)#dlinesdw

GOTO,mdi
; LINEARIZED INVERSE WITH FITTING: 
; LEAST-SQUARES FIT (JUST THEORETICAL COMPUTATION OF HESSIAN MATRIX)
;------------------------------------------------------------------------------------------------
sigma1   = DBLARR(ntest)
a        = DBLARR(ntune,npar)
FOR itest= 0,ntest-1 DO BEGIN
    i0     = texp*inten[*,itest]   ; Calculate intensities
    si0    = 1.d0/SQRT(i0)         ; Calculate inverse sigmas
    a[*,0] = texp*dintendv[*,itest]*si0 ; Calculate derivative with respect to velocity
    a[*,1] = texp*dintendi[*,itest]*si0
    a[*,2] = texp*dintendd[*,itest]*si0
    a[*,3] = texp*dintendw[*,itest]*si0
    cov    = INVERT(TRANSPOSE(a)#a)
    sigma1[itest] = SQRT(cov[0,0]) ; Calculate least squares estimate of noise. (This is just a weighted mean.)
ENDFOR

dintendv2= TRANSPOSE(filterseb)#dlinesdv ; from actual filters
dintendd2= TRANSPOSE(filterseb)#dlinesdd
dintendi2= TRANSPOSE(filterseb)#dlinesdi
dintendw2= TRANSPOSE(filterseb)#dlinesdw 
sigma1b  = DBLARR(ntest)
a        = DBLARR(ntune,npar)
FOR itest= 0,ntest-1 DO BEGIN
    i0     = texp*inten2[*,itest] ; Calculate intensities
    si0    = 1.d0/SQRT(i0)         ; Calculate inverse sigmas
    a[*,0] = texp*dintendv2[*,itest]*si0 ; Calculate derivative with respect to velocity
    a[*,1] = texp*dintendi2[*,itest]*si0
    a[*,2] = texp*dintendd2[*,itest]*si0
    a[*,3] = texp*dintendw2[*,itest]*si0
    cov    = INVERT(TRANSPOSE(a)#a)
    sigma1b[itest] = SQRT(cov[0,0]) ; Calculate least squares estimate of noise. (This is just a weighted mean.)
ENDFOR        


SET_PLOT,'ps'
DEVICE,file='yo2.ps',xoffset=0,yoffset=0,xsize=25,ysize=25,/color
LOADCT,3
!P.MULTI=[0,2,2]

a=WHERE(ABS(vtest) EQ MIN(ABS(vtest)))
PLOT,lam,lyotseb/max(lyotseb),xst=1,xtit='Wavelength (A)',ytit='Transmittance',charsize=1.15,xrange=[-5.5,5.5]
OPLOT,lam,lyot/max(lyot),linestyle=2
OPLOT,[0,0],[0,1]
OPLOT,lam,lines[*,a],thick=2
OPLOT,lam,blockerseb,col=180
OPLOT,lam,blocker,col=180,linestyle=2

PLOT,lam,lines[*,a],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.15,tit='!17',xtit='!7k!17 (A)',ytit='Transmittance' ,yst=1,yrange=[0,1]
loadct,33
for i=0,ntune-1 do OPLOT,lam,filterseb[*,i]/MAX(filterseb),color=i*(254./ntune)+1
  
loadct,0
PLOT ,vtest,sigma1b/SQRT(2.0),xrange=[-8000,8000],yrange=[0,45],charsize=1.15,xst=1,yst=1,tit='!17',xtit='Input velocity (m/s)',ytit='Error due to photon noise (m/s)'
OPLOT,vtest,sigma1/SQRT(2.0),linestyle=2
OPLOT,[-6500,-6500],[0,200],linestyle=3
OPLOT,[ 6500, 6500],[0,200],linestyle=3


DEVICE,/close
SET_PLOT,'x'
!P.MULTI=0

; LEAST-SQUARE FIT WITH
; SIMULTANEOUS FIT OF LCP AND RCP
;------------------------------------------------------------------------------------------------

    vt1      = DBLARR(ntest,ntest) ; simultaneous fit of LCP and RCP assuming same
    vt2      = DBLARR(ntest,ntest) ; intensity and depth.
    sigmas   = DBLARR(npar+1,ntest,ntest)
    a        = DBLARR(2*ntune,npar+1)
    FOR it1  = 0,ntest-1 DO BEGIN
        FOR it2= 0,ntest-1 DO BEGIN
            vt1[it1,it2] = vtest[it1]
            vt2[it1,it2] = vtest[it2]
            i0           = texp*[inten[*,it1],inten[*,it2]] ;   Calculate intensities
            si0          = 1.d0/SQRT(i0) ;   Calculate inverse sigmas
            a[*,0]       = texp*[dintendv[*,it1],0*dintendv[*,it2]]*si0 ;derivative / velocity
            a[*,1]       = texp*[0*dintendv[*,it1],dintendv[*,it2]]*si0
            a[*,2]       = texp*[dintendi[*,it1],dintendi[*,it2]]*si0
            a[*,3]       = texp*[dintendd[*,it1],dintendd[*,it2]]*si0
            a[*,4]       = texp*[dintendw[*,it1],dintendw[*,it2]]*si0
            cov          = invert(transpose(a)#a)
            FOR i = 0,npar DO sigmas[i,it1,it2] = SQRT(cov[i,i]) ;   Calculate least squares estimate of noise. (This is just a weighted mean.)
        ENDFOR
    ENDFOR

    vta     = [vt1+vt2]/2.d0
    vts     = [vt2-vt1]/2.d0
    vt1r    = REFORM(vt1,ntest^2.d0)
    vt2r    = REFORM(vt2,ntest^2.d0)
    vtar    = REFORM(vta,ntest^2.d0)
    vtsr    = REFORM(vts,ntest^2.d0)
    sigmasr = REFORM(sigmas,npar+1,ntest^2.d0)
    sigmavr = SQRT(sigmasr[0,*]^2.d0+sigmasr[1,*]^2.d0)/2.d0
    sigmav  = REFORM(sigmavr,ntest,ntest)

vi0 = FINDGEN(ntest)
; This is just the diagonal elements. ii=lindgen(ntest) & sig0=sigmav(ii,ii)
sig0= INTERPOLATE(sigmav,vi0,vi0)
vix = vi0-dvdb*1000.d0/dvtest ; in presence of a 1000 Gauss field
viy = vi0+dvdb*1000.d0/dvtest
sig1000=INTERPOLATE(sigmav,vix,viy)
vix = vi0-dvdb*2000.d0/dvtest
viy = vi0+dvdb*2000.d0/dvtest
sig2000=INTERPOLATE(sigmav,vix,viy)
vix = vi0-dvdb*3000.d0/dvtest
viy = vi0+dvdb*3000.d0/dvtest
sig3000=INTERPOLATE(sigmav,vix,viy)
vix = vi0-dvdb*4000.d0/dvtest
viy = vi0+dvdb*4000.d0/dvtest
sig4000=INTERPOLATE(sigmav,vix,viy)

READ,PAUSE
mdi:

; WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

period     = 5.d0*dtune
Kfourier   = dtune/period*2.d0

; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd
; COEFFICIENT) FOR WHAT WE THINK ARE THE CORRECT FILTER PROFILES
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune
c1         = REFORM(COS(x)#inten)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)#inten)  ; coefficient b1 (for the sine)
linedepths = sqrt(c1^2.d0+s1^2.d0)
c2         = REFORM(COS(2.d0*x)#inten)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)#inten)  ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0   ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi    ; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi    ; velocity measurement
vel1a      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; the lookup tables for the filters we think are correct: actual velocities vs. measured ones
vel2a      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1

; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd
; COEFFICIENT) FOR THE ACTUAL FILTER PROFILES (VALID ONLY WHEN B=0)
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune
c1         = REFORM(COS(x)#inten2)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)#inten2)  ; coefficient b1 (for the sine)
linedepths = sqrt(c1^2.d0+s1^2.d0)
c2         = REFORM(COS(2.d0*x)#inten2)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)#inten2)  ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0   ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi    ; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi    ; velocity measurement
vel1b      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; the lookup tables for the filters we think are correct: actual velocities vs. measured ones
vel2b      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1

READ,pause


; COMPUTE DOPPLER VELOCITY ERROR DUE TO PHOTON NOISE AND SYSTEMATICS

sigmata  = FLTARR(ntest)
sigmataB = FLTARR(ntest)
sigmata2  = FLTARR(ntest)
sigmataB2 = FLTARR(ntest)
systematicerror1a= FLTARR(ntest)
velocity = FLTARR(ntest)
velocity2 = FLTARR(ntest)
systematicerror1a2= FLTARR(ntest)
linewidth= FLTARR(ntest)
linedepth= FLTARR(ntest)
Ic       = FLTARR(ntest)
meanL    = FLTARR(ntest)
meanR    = meanL
sigmaIc  = FLTARR(ntest)
sigmawidth = FLTARR(ntest)
sigmadepth  = FLTARR(ntest)
FOR itest= 0,ntest-1 DO BEGIN
    vt   = vtest[itest]

    ; LCP
    ;-----------------------------------------------------------------------------------
    i0   = inten2[*,itest]*texp
    intent=DBLARR(ntune,nt)
    FOR i= 0,ntune-1 DO intent[i,*]=i0[i]+SQRT(i0[i])*randomge[*,i] ; SQRT(i[0]) is photon noise
    c1   = REFORM(COS(x)#intent)
    s1   = REFORM(SIN(x)#intent)
    c2   = REFORM(COS(2.d0*x)#intent)
    s2   = REFORM(SIN(2.d0*x)#intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel2t = phi2t*pv2/2.d0/!dpi
    vel1at= vel1t
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1tL = INTERPOL(vtest,vel1a,vel1at) ; velocity obtained from look-up table
    v2tL = INTERPOL(vtest,vel2a,vel2at) ; velocity obtained from look-up table
    ;covt12a         = MEAN((v1t-MEAN(v1t))*(v2t-MEAN(v2t)))
    ;Weigth1 = (sigmat2a[itest]-covt12a)/(sigmat1a[itest]+sigmat2a[itest]-2.d0*covt12a)
    ;Weigth2 = 1.d0-Weigth1
    ;sigmata[itest] = SQRT(REBIN(((v1t-MEAN(v1t))*Weigth1+(v2t-MEAN(v2t))*Weigth2)^2.0,1))


    c1=c1*Kfourier
    s1=s1*Kfourier
    c2=c2*Kfourier
    s2=s2*Kfourier
   ;We compute the linewidth (in Angstroms)
    tempL        = period/!dpi*sqrt(1.0/6.0*alog((c1*c1+s1*s1)/(c2*c2+s2*s2))) ;
    widthgL      = tempL*2.d0*sqrt(alog(2.d0)) ; //we want the FWHM not the sigma of the Gaussian	      
   ;We compute the linedepth
    temp2L       = period/2.d0*sqrt(c1*c1+s1*s1)/sqrt(!dpi)/tempL*exp(!dpi*!dpi*tempL*tempL/period/period) ;
    IdgL         = temp2L/WELL
   ;We compute the continuum intensity
    temp3L       = (v1tL+v2tL)/2.0/dvdlam
    FOR i= 0,ntest-1 DO BEGIN
        meanL[i]=MEAN(intent[*,i])
        FOR j=0,5 DO meanL[i] = meanL[i]+temp2L[i]/6.d0*exp(-(tune[0,j]-temp3L[i])*(tune[0,j]-temp3L[i])/tempL[i]/tempL[i]) ; 
    ENDFOR
    meanL        = meanL/WELL

    ; RCP
    ;-----------------------------------------------------------------------------------
    i0   = inten3[*,itest]*texp
    intent=DBLARR(ntune,nt)
    FOR i= 0,ntune-1 DO intent[i,*]=i0[i]+SQRT(i0[i])*randomge2[*,i] ; SQRT(i[0]) is photon noise
    c1   = REFORM(COS(x)#intent)
    s1   = REFORM(SIN(x)#intent)
    c2   = REFORM(COS(2.d0*x)#intent)
    s2   = REFORM(SIN(2.d0*x)#intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel2t = phi2t*pv2/2.d0/!dpi
    vel1at= vel1t
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1tR = INTERPOL(vtest,vel1a,vel1at) ; velocity obtained from look-up table
    v2tR = INTERPOL(vtest,vel2a,vel2at) ; velocity obtained from look-up table

    c1=c1*Kfourier
    s1=s1*Kfourier
    c2=c2*Kfourier
    s2=s2*Kfourier
   ;We compute the linewidth (in Angstroms)
    tempR        = period/!dpi*sqrt(1.0/6.0*alog((c1*c1+s1*s1)/(c2*c2+s2*s2))) ;
    widthgR      = tempR*2.d0*sqrt(alog(2.d0)) ; //we want the FWHM not the sigma of the Gaussian	      
   ;We compute the linedepth
    temp2R       = period/2.d0*sqrt(c1*c1+s1*s1)/sqrt(!dpi)/tempR*exp(!dpi*!dpi*tempR*tempR/period/period) ;
    IdgR         = temp2R/WELL
   ;We compute the continuum intensity
    temp3R       = (v1tR+v2tR)/2.0/dvdlam
    FOR i= 0,ntest-1 DO BEGIN
        meanR[i]=MEAN(intent[*,i])
        FOR j=0,5 DO meanR[i] = meanR[i]+temp2R[i]/6.d0*exp(-(tune[0,j]-temp3R[i])*(tune[0,j]-temp3R[i])/tempR[i]/tempR[i]) ; 
    ENDFOR
    meanR        = meanR/WELL

   ;compute observables
   ;--------------------------------------------------------------------------

    velocity[itest]=MEAN(v1tL+v2tL+v1tR+v2tR)/4.d0
    sigmata[itest] =SQRT(REBIN(((v1tL+v2tL+v1tR+v2tR)/4.d0-velocity[itest])^2.d0,1)) ; STD DEVIATION ON DOPPLER VELOCITY
    systematicerror1a[itest] = velocity[itest]-vtest[itest]                          ; SYSTEMATIC ERROR ON DOPPLER VELOCITY
    sigmataB[itest] = sigmata[itest]*2.d0*0.2314d0                                   ; STD DEVIATION ON L.O.S. FIELD STRENGTH


    velocity2[itest]=MEAN(v1tL+v1tR)/2.d0                                  ; WE ONLY USE THE 1ST FOURIER COEFFICIENTS
    sigmata2[itest]  =SQRT(REBIN(((v1tL+v1tR)/2.d0-velocity[itest])^2.d0,1)) ; STD DEVIATION ON DOPPLER VELOCITY
    systematicerror1a2[itest] = velocity[itest]-vtest[itest]               ; SYSTEMATIC ERROR ON DOPPLER VELOCITY
    sigmataB2[itest] = sigmata2[itest]*2.d0*0.2314d0                       ; STD DEVIATION ON L.O.S. FIELD STRENGTH


    linewidth[itest]=MEAN((widthgL+widthgR)/2.d0)
    sigmawidth[itest]=SQRT(REBIN(((widthgL+widthgR)/2.d0-linewidth[itest])^2.d0,1))  ; STD DEVIATION OF LINEWIDTH

    linedepth[itest]=MEAN((IdgL+IdgR)/2.d0)
    sigmadepth[itest]=SQRT(REBIN(((IdgL+IdgR)/2.d0-linedepth[itest])^2.d0,1))        ; STD DEVIATION OF LINEWIDTH

    Ic[itest]    = MEAN((meanL+meanR)/2.d0)
    sigmaIc[itest]=SQRT(REBIN(((meanL+meanR)/2.d0-Ic[itest])^2.d0,1))                ; STD DEVIATION OF CONTINUUM
 
ENDFOR

SET_PLOT,'PS'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=28,/color
!P.MULTI=[0,1,2]
LOADCT,3
PLOT,vtest,systematicerror1a,xst=1,xrange=[-6500,6500],thick=3,tit='!17 MDI-like algorithm, 6 filters, at disk center',xtit='Input Velocity (m/s)',ytit='Systematic Error (m/s)',charsize=1.5
PLOT,vtest,sigmata2,xst=1,xrange=[-6500,6500],thick=3,tit='!17',xtit='Input Velocity (m/s)',ytit='Standard Deviation (m/s)',charsize=1.5
DEVICE,/CLOSE

DEVICE,FILE='yo2.ps',xoffset=0,yoffset=0,xsize=20,ysize=28,/color
LOADCT,3
;PLOT,vtest,sigmata,xst=1,xrange=[0,6500],thick=3,tit='!17',xtit='Input Velocity (m/s)',ytit='Standard Deviation (m/s)',charsize=1.5
;PLOT,vtest,sigmata*2.d0
;PLOT,vtest/dvdb,sigmata*2.d0*0.2314d0,xst=1,xrange=[-6500,6500]/dvdb,thick=3,tit='!17 MDI-like algorithm, 6 filters, v=0',xtit='l.o.s. Magnetic Field (G)',ytit='Standard Deviation (m/s)',charsize=1.5
PLOT,vtest,Ic,xst=1,xrange=[-6500,6500],thick=3,tit='!17',xtit='Input Velocity (m/s)',ytit='Systematic Error on Ic',charsize=1.5,yst=1
PLOT,vtest,sigmaIc*1000.,xst=1,xrange=[-6500,6500],thick=3,tit='!17',xtit='Input Velocity (m/s)',ytit='Standard Deviation of Ic (*10!U-3!N)',charsize=1.5,yst=1
DEVICE,/CLOSE

!P.MULTI=0
SET_PLOT,'X'


READ,pause

;c1      = REFORM(COS(x)#inten2) ; coefficient a1
;s1      = REFORM(SIN(x)#inten2) ; coefficient b1 (for the sine)
;c2      = REFORM(COS(2.d0*x)#inten2) ; coefficient a1
;s2      = REFORM(SIN(2.d0*x)#inten2) ; coefficient b1 (for the sine)
;phi1    = ATAN(-s1,-c1)
;phi2    = ATAN(-s2,-c2)
;vel1b   = phi1*pv1/2.d0/!dpi
;vel2b   = phi2*pv2/2.d0/!dpi
;vel1bb  = (vel1b-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; correct look-up table
;vel2bb  = (vel2b-vel1b+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1b ; correct look-up table
;sigmat1b = FLTARR(ntest)
;sigmat2b = FLTARR(ntest)
;sigmatb  = FLTARR(ntest)
;systematicerror1b= FLTARR(ntest)
;FOR itest= 0,ntest-1 DO BEGIN
;    vt   = vtest[itest]
;    i0   = inten2[*,itest]*texp
;    intent=DBLARR(ntune,nt)
;    FOR i= 0,ntune-1 DO intent[i,*]=i0[i]+SQRT(i0[i])*randomge[*,i] 
;    c1   = REFORM(COS(x)#intent)
;    s1   = REFORM(SIN(x)#intent)
;    c2   = REFORM(COS(2.d0*x)#intent)
;    s2   = REFORM(SIN(2.d0*x)#intent)
;    phi1t= ATAN(-s1,-c1)
;    phi2t= ATAN(-s2,-c2)
;    vel1t= phi1t*pv1/2.d0/!dpi
;    vel2t= phi2t*pv2/2.d0/!dpi
;    vel1at= (vel1t-vt+10.5d0*pv1) MOD pv1-pv1/2.d0+vt
;    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
;    v1tb = INTERPOL(vtest,vel1a,vel1at,/SPLINE)
;    v2tb = INTERPOL(vtest,vel2a,vel2at,/SPLINE)
;    sigmat1b[itest] = SQRT(REBIN((v1tb-MEAN(v1tb))^2.d0,1)) ; photon noise
;    sigmat2b[itest] = SQRT(REBIN((v2tb-MEAN(v2tb))^2.d0,1)) ; photon noise
;    covt12b         = MEAN((v1tb-MEAN(v1tb))*(v2tb-MEAN(v2tb)))
;    Weigth1 = (sigmat2b[itest]-covt12b)/(sigmat1b[itest]+sigmat2b[itest]-2.d0*covt12b)
;    Weigth2 = 1.d0-Weigth1
;    sigmatb[itest] = SQRT(REBIN(((v1tb-MEAN(v1tb))*Weigth1+(v2tb-MEAN(v2tb))*Weigth2)^2.0,1))
;    systematicerror1b[itest] = MEAN(v1tb-vtest[itest]) ; systematic error
;ENDFOR

;-------------------------------------------------------------------------------------

; PLOT
;WINDOW,0,RETAIN=2,xsize=600,ysize=900
;SET_PLOT,'ps'
;DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=25,ysize=25,/color
;LOADCT,3
;!P.MULTI=[0,2,2]
;a  = WHERE(vtest ge -3200. and vtest le 3200.)
;bb = WHERE(vtest ge -6500. and vtest le 6500.)
;res= POLY_FIT(vtest[a],systematicerror1b[a],2,yfit=y)
;y  = res[0]+res[1]*vtest+res[2]*vtest^2.d0
;PLOT,vtest,systematicerror1b,xrange=[-6500,6500],xst=1,xtit='Input velocity (m/s)',ytit='Error (m/s)',charsize=1.15,yrange=[MIN([y[bb],systematicerror1b[bb]]),MAX([y[bb],systematicerror1b[bb]])]
;OPLOT,vtest,systematicerror1a,linestyle=2
;OPLOT,vtest,y,col=180
;
;a=WHERE(ABS(vtest) EQ MIN(ABS(vtest)))
;PLOT,lam,lyotseb/max(lyotseb),xst=1,xtit='Wavelength (A)',ytit='Transmittance',charsize=1.15,xrange=[-5.5,5.5]
;OPLOT,lam,lyot/max(lyot),linestyle=2
;OPLOT,[0,0],[0,1]
;OPLOT,lam,lines[*,a],thick=2
;OPLOT,lam,blockerseb,col=180
;OPLOT,lam,blocker,col=180,linestyle=2
;
;PLOT,lam,lines[*,a],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.15,tit='!17',xtit='!7k!17 (A)',ytit='Transmittance' ,yst=1,yrange=[0,1]
;loadct,33
;for i=0,ntune-1 do OPLOT,lam,filterseb[*,i]/MAX(filterseb),color=i*(254./ntune)+1
;  
;loadct,0
;PLOT ,vtest,sigmatb/SQRT(2.0),xrange=[-6500,6500],yrange=[0,45],charsize=1.15,xst=1,yst=1,tit='!17',xtit='Input velocity (m/s)',ytit='Error due to photon noise (m/s)'
;OPLOT,vtest,sigmata/SQRT(2.0),linestyle=2
;
;; FWHM CALCULATION
;lyotseb = lyotseb/max(lyotseb)
;lyot    = lyot/max(lyot)
;a = where(lyotseb GE 0.9975*0.5d0 AND lyotseb LE 1.0025*0.5d0)
;;PRINT,lam[a]
;b=WHERE(lam[a] LE 0.,complement=c)
;FWHM=MEAN(lam[a[c]])-MEAN(lam[a[b]])
;a=WHERE(lyotseb EQ 1.0)
;PRINT,'ACTUAL FILTER',FWHM,lam[a]
;a = where(lyot GE 0.9975*0.5d0 AND lyot LE 1.0025*0.5d0)
;;PRINT,lam[a]
;b=WHERE(lam[a] LE 0.,complement=c)
;FWHM2=MEAN(lam[a[c]])-MEAN(lam[a[b]])
;a=WHERE(lyot EQ 1.0)
;PRINT,'THEORETICAL FILTER USED FOR LOOK-UP TABLES',FWHM2,lam[a]
;DEVICE,/close
;SET_PLOT,'x'
;!P.MULTI=0
;read,pause
END


; FOR DRAWING THE RESULTS WITH DIFFERENT LONGITUDES

;restore,'sigma0.bin'
;s0=sigmatb
;restore,'sigma10.bin'
;s10=sigmatb
;restore,'sigma20.bin'
;s20=sigmatb
;restore,'sigma30.bin'
;s30=sigmatb
;restore,'sigma40.bin'
;s40=sigmatb
;restore,'sigma50.bin'
;s50=sigmatb
;restore,'sigma60.bin'
;s60=sigmatb
;restore,'sigma70.bin'
;s70=sigmatb
;restore,'sigma80.bin'
;s80=sigmatb
;restore,'sigma90.bin'
;s90=sigmatb
;restore,'sigmat0.bin'
;s0=sigmatb
;restore,'sigmat45.bin'
;s45=sigmatb
;restore,'sigmat60.bin'
;s60=sigmatb
;set_plot,'ps'
;device,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color
;loadct,4
;ntest         = 501 ; Set test velocities
;dvtest        = 100.d0 ; in m/s
;vtest         = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
;plot,vtest,s60/sqrt(2.),xst=1,xrange=[-6500,6500],charsize=1.5,tit='!17',xtit='Velocity (m/s)',ytit='Error due to photon noise (m/s)',yrange=[0,25],yst=1
;oplot,vtest,s0/sqrt(2.),color=240
;;;oplot,vtest,s10/sqrt(2.),color=210
;;;oplot,vtest,s20/sqrt(2.),color=180
;;;oplot,vtest,s30/sqrt(2.),color=150
;;;oplot,vtest,s40/sqrt(2.),color=120
;oplot,vtest,s45/sqrt(2.),color=120
;;;oplot,vtest,s50/sqrt(2.),color=90
;;;oplot,vtest,s60/sqrt(2.),color=60
;;;oplot,vtest,s70/sqrt(2.),color=30
;;;oplot,vtest,s80/sqrt(2.),color=10
;;;legend,['0','10','20','30','40','50','60','70','80','90'],color=[240,210,180,150,120,90,60,30,10,0],linestyle=[0,0,0,0,0,0,0,0,0,0],/bottom
;legend,['0','45','60'],color=[240,120,0],linestyle=[0,0,0],/bottom
;;
;device,/close
