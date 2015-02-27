; program to test how we can reconstruct the look-up tables for a
; MDI-like algorithm
; by using the orbital velocity of SDO

PRO sinus,X,A,F,pder

Nx=N_ELEMENTS(X)

F=A[0]*SIN(2.d0*!dpi*A[1]*X+A[2])+A[3]
;F=A[0]*TAN(2.d0*!dpi*A[1]*X+A[2])+A[3]



pder=[ [SIN(2.d0*!dpi*A[1]*X+A[2])],[2.d0*!dpi*A[0]*X*COS(2.d0*!dpi*A[1]*X+A[2])],[A[0]*COS(2.d0*!dpi*A[1]*X+A[2])],[FLTARR(Nx)+1.d0] ]
;pder=[ [TAN(2.d0*!dpi*A[1]*X+A[2])],[2.d0*!dpi*A[0]*X/COS(2.d0*!dpi*A[1]*X+A[2])^2.d0],[A[0]/COS(2.d0*!dpi*A[1]*X+A[2])^2.d0],[FLTARR(Nx)+1.d0] ]

END


PRO tangent,X,A,F,pder

Nx=N_ELEMENTS(X)
F=A[0]*ATAN(A[1]*X+A[2])+A[3]

pder=[ [ATAN(A[1]*X+A[2])],[A[0]*X/(1.d0+(A[1]*X+A[2])^2.d0)],[A[0]/(1.d0+(A[1]*X+A[2])^2.d0)],[FLTARR(Nx)+1.d0] ]

END

PRO HMI_orbitalvel

lam0    = 6173.3433d0
ntune   = 6           ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0   ; Set number of tuning positions over wnarrow. (1/spacing)
ntime   = 600.*3.        ; duration of observation in minutes
nframe  = ntune*2     ; total number of filtergrams taken per observable (each filter is used twice
                      ; because of LCP and RCP) 
cadence = 45.d0
dt      = cadence/nframe ; time interval between 2 filtergrams, or integration time (?)
dtime   = 60.d0/dt
nt      = ntime*dtime        ; total number of filtergrams
nt2     = nt/nframe   ; total number of observables
time    = FINDGEN(nt)*dt

; Parameters for the solar line
;------------------------------------------------------------------------

nlam       = 18000.;21500.
dlam       = 1.d0/1.75d3
lam        = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

; Fe I LINE PROFILE FROM ROGER ULRICH'S WEBSITE
;-------------------------------------------------------------------------

;OPENR,1,'Ulrich_Fe_0.txt'
;roger  = FLTARR(2,98)
;READF,1,roger
;CLOSE,1
;rlam   = REFORM(roger[0,*])
;rp     = REFORM(roger[1,*])
;rlam   = [-10.d0,rlam,-10.d0]
;rp     = [1.d0,rp,1.d0]
;line   = INTERPOL(rp,rlam,lam)
;line   = line/INTERPOL([line[0],line[nlam-1]],[lam[0],lam[nlam-1]],lam)

; MODEL OF SOLAR LINE
;-------------------------------------------------------------------------
width  = 0.102d0/2.d0/SQRT(ALOG(2.d0)) ; w in A^2 (part added by seb)
depth  = 0.62d0
line  = 1.d0-depth*EXP(-lam^2.d0/width^2.d0)
dlinedv= depth*2.d0*lam/width^2.d0*EXP(-lam^2.d0/width^2.d0)

; REFERENCE LOOK-UP TABLES WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

phase      = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
contrast   = [1.,1.,1.,1.,1.,1.,1.]
wmich      = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]  
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0]

dtune      = wmich[0]/inttune  ; Set tuning position interval
tune       = DBLARR(3,ntune)
FOR  i     = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

dlamdv     = lam0/299792458.d0
dvdlam     = 1.d0/dlamdv
vel        = lam*dvdlam   ; DOPPLER shift in cm/s
dvel       = dlam*dvdlam
dvtune     = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune      = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; theoretical blocker filter
blocker    = EXP(-lam^2.d0/2.d0/(8.d0/2.354820d0)^2.d0)

; NON-TUNABLE PROFILE
lyot       = blocker
FOR i = 0,3 DO lyot = lyot*(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0

; TUNABLE PROFILE
cmich      = 2.d0*!dpi/wmich
filters    = DBLARR(nlam,ntune)
FOR itune  = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i  = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
ENDFOR

; the filters are in this order:
;--------------------------------------------------------
;filters[*,0] is F0 centered on +171 mA
;filters[*,1] is F1 centered on +103 mA
;filters[*,2] is F2 centered on +034 mA
;filters[*,3] is F3 centered on -034 mA
;filters[*,4] is F4 centered on -103 mA
;filters[*,5] is F5 centered on -171 mA
;--------------------------------------------------------

ntest      = 501 ; Set test velocities
dvtest     = dlam/dlamdv;100.d0 ; in m/s
vtest      = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
lines      = DBLARR(nlam,ntest)
FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+(vtest[i])*dlamdv)
inten      = TRANSPOSE(filters)#lines

x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune
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
vel1a      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; the lookup tables
vel2a      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1

; SECOND SET OF LOOK-UP TABLES WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

; CHANGE THE SOLAR LINE (TO SIMULATE DIFFERENCE BETWEEN OBSMODE AND CALMODE)
;width  = 1.2d0*0.102d0/2.d0/SQRT(ALOG(2.d0)) ; w in A^2 (part added by seb)
;depth  = 0.8d0*0.62d0
;line  = 1.d0-depth*EXP(-lam^2.d0/width^2.d0)
;dlinedv= depth*2.d0*lam/width^2.d0*EXP(-lam^2.d0/width^2.d0)

; PHASE AND CONTRAST FOR THE SECON D SET OF LOOKUP TABLES
phase      = [6.0d0,-6.d0,6.d0,6.0700685d0,3.4029859d0,8.2042572d0,20.1685274d0]*!dpi/180.d0 
contrast   = [0.98,0.98,0.97,0.952,0.964,0.987,1.0]  
; FRONT WINDOW AND BLOCKER FILTER TRANSMISSION PROFILES
RESTORE,'frontwindow.bin' ; front window
blocker    = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
q          = READFITS('blocker11.fits')
blocker    = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) 

; TUNABLE PROFILE
cmich      = 2.d0*!dpi/wmich
filters    = DBLARR(nlam,ntune)
FOR itune  = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i  = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
ENDFOR

; we simulate the intensities actually measured by HMI (IF IN OBSMODE,
; WE ADD RANDOM GAUSSIAN NOISE TO SIMULATE THE IMPACT OF P-MODES + GRANULATION
seed=1l
noiselevel = 226.d0 ; 226=rms variation of Doppler velocities on solar disk due to p modes and granulation
noiselevel = 50.d0
RESTORE,'velocity_pmodes.bin' ; from quiet Sun data at [256,256], I just added the same data 3 times in a row.

ntest      = 501 ; Set test velocities
dvtest     = dlam/dlamdv;100.d0 ; in m/s
vtest      = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
lines      = DBLARR(nlam,ntest)
;FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+(vtest[i]+RANDOMN(seed)*noiselevel)*dlamdv)
FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+(vtest[i]+velocity_pmodes[i])*dlamdv)
inten      = TRANSPOSE(filters)#lines ; actual intensities obtained by HMI

; I apply the MDI-like algorithm on the intensities measured by HMI
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune
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
vel1b      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; velocity returned by the MDI-like algorithm
vel2b      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1  ; velocity returned by the MDI-like algorithm


; WHAT THE ACTUAL LOOK-UP TABLES SHOULD LOOK LIKE IN ABSENCE OF P
; MODES AND GRANULATION

FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+vtest[i]*dlamdv) ; 226=rms variation of Doppler velocities on solar disk due to p modes and granulation
inten      = TRANSPOSE(filters)#lines ; actual intensities obtained by HMI

; I apply the MDI-like algorithm on the intensities measured by HMI
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune
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
vel1bb      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; velocity returned by the MDI-like algorithm
vel2bb      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1  ; velocity returned by the MDI-like algorithm



WINDOW,0,RETAIN=2,xsize=2000,YSIZE=1000
!p.multi=[0,4,2]

; FIT BY POLYNOMIALS
;-----------------------------------

read,pause

; LOW-PASS FILTER TO TRY TO REMOVE THE P-MODE NOISE ON THE LOOK-UP TABLE
filter = exp(-vtest^2.d0/1000000.d0)
filter = shift(filter,-250)
vel1b=fft(fft(vel1b)*fft(filter/TOTAL(filter)),/inverse)*501

moishe=WHERE(abs(vtest) le 3000.,nmoishe)

ntest=1440.;6 ; number of cotunes: 1440 means 24h of cotunes at a cadence of 60 seconds
; IN OBSMODE: SMOOTH vel1b FIRST
;vel1b=SMOOTH(vel1b,8)
;vel2b=SMOOTH(vel2b,8)
selection=ROUND(findgen(ntest)/(ntest-1.d0)*(nmoishe-1.d0))+moishe[0]

res1=POLY_FIT(vtest[selection],vel1b[selection]-vel1a[selection],1,yerror=error1)
res2=POLY_FIT(vtest[selection],vel1b[selection]-vel1a[selection],2,yerror=error2)
error2=10000000000.d0; so that it's always a linear fit
IF(error1 LE error2) THEN vel1b_rec=res1[0]+vtest*res1[1] ELSE vel1b_rec=res2[0]+vtest*res2[1]+vtest^2.d0*res2[2]


res1b=POLY_FIT(vtest[selection],vel1bb[selection]-vel1a[selection],1,yerror=error1)
res2b=POLY_FIT(vtest[selection],vel1bb[selection]-vel1a[selection],2,yerror=error2)
error2=10000000000.d0; so that it's always a linear fit
IF(error1 LE error2) THEN vel1bb_rec=res1b[0]+vtest*res1b[1] ELSE vel1bb_rec=res2b[0]+vtest*res2b[1]+vtest^2.d0*res2b[2]

plot,vtest,vel1b,xst=1,yst=1,charsize=1.5
oplot,vtest,vel1b_rec+vel1a,col=180
oplot,vtest,vel1bb_rec+vel1a,col=180,linestyle=2,thick=2
oplot,vtest[selection],vel1b_rec[selection]+vel1a[selection],psym=4

plot,vtest,vel1bb-(vel1b_rec+vel1a),xst=1,yst=1,charsize=1.5
oplot,vtest,vel1bb-(vel1bb_rec+vel1a),col=180

res1=POLY_FIT(vtest[selection],vel2b[selection]-vel2a[selection],1,yerror=error1)
res2=POLY_FIT(vtest[selection],vel2b[selection]-vel2a[selection],2,yerror=error2)
error2=10000000000.d0; so that it's always a linear fit
IF(error1 LE error2) THEN vel2b_rec=res1[0]+vtest*res1[1] ELSE vel2b_rec=res2[0]+vtest*res2[1]+vtest^2.d0*res2[2] 

res1b=POLY_FIT(vtest[selection],vel2bb[selection]-vel2a[selection],1,yerror=error1)
res2b=POLY_FIT(vtest[selection],vel2bb[selection]-vel2a[selection],2,yerror=error2)
error2=10000000000.d0; so that it's always a linear fit
IF(error1 LE error2) THEN vel2bb_rec=res1b[0]+vtest*res1b[1] ELSE vel2bb_rec=res2b[0]+vtest*res2b[1]+vtest^2.d0*res2b[2] 


plot,vtest,vel2b,xst=1,yst=1,charsize=1.5
oplot,vtest,vel2b_rec+vel2a,col=180
oplot,vtest,vel2bb_rec+vel2a,col=180,linestyle=2,thick=2

plot,vtest,vel2bb-(vel2b_rec+vel2a),xst=1,yst=1,charsize=1.5
oplot,vtest,vel2bb-(vel2bb_rec+vel2a),col=180


; FIT BY SINE
;-----------------------------------

vmax=10000.d0
vtest=vtest/vmax
vel1b=vel1b/vmax
vel2b=vel2b/vmax

A=[1.,1.d0/40000.,0.,0.]

weightsp=FLTARR(ntest)+1.d0

resp=CURVEFIT(vtest[selection],vel1b[selection],weightsp[selection],A,FUNCTION_NAME='sinus',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat1)

PLOT,vtest*vmax,vel1b*vmax,xrange=[-6500,6500],xst=1,charsize=1.5
OPLOT,vtest*vmax,(A[0]*SIN(2.d0*!dpi*A[1]*vtest+A[2])+A[3])*vmax,col=180
Plot,vtest*vmax,(vel1b-(A[0]*SIN(2.d0*!dpi*A[1]*vtest+A[2])+A[3]))*vmax,xrange=[-6500,6500],xst=1,yst=1,charsize=1.5

A=[1.,1.d0/40000.,0.,0.]

resp=CURVEFIT(vtest[selection],vel2b[selection],weightsp[selection],A,FUNCTION_NAME='sinus',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat1)

PLOT,vtest*vmax,vel2b*vmax,xrange=[-6500,6500],xst=1,charsize=1.5
OPLOT,vtest*vmax,(A[0]*SIN(2.d0*!dpi*A[1]*vtest+A[2])+A[3])*vmax,col=180
Plot,vtest*vmax,(vel2b-(A[0]*SIN(2.d0*!dpi*A[1]*vtest+A[2])+A[3]))*vmax,xrange=[-6500,6500],xst=1,yst=1,charsize=1.5
;OPLOT,vtest*vmax,(A[0]*ATAN(A[1]*vtest+A[2])+A[3])*vmax,col=180
;Plot,vtest*vmax,(vel2b-(A[0]*ATAN(A[1]*vtest+A[2])+A[3]))*vmax,xrange=[-6500,6500],xst=1,yst=1,charsize=1.5




READ,pause


;order1=7
;res=POLY_FIT(vtest[a],vel1a[a],order1)
;vel1a_rec=res[0]+vtest*res[1]
;FOR i=2,order1 DO vel1a_rec=vel1a_rec+vtest^i*res[i]

;order2=15
;res=POLY_FIT(vtest[a],vel2a[a],order2)
;vel2a_rec=res[0]+vtest*res[1]
;FOR i=2,order2 DO vel2a_rec=vel2a_rec+vtest^i*res[i]

;t=FINDGEN(15)/14.*na+a[0]
;Y2=SPL_INIT(vtest[t],vel2a[t])
;res2=SPL_INTERP(vtest[t],vel2a[t],Y2,vtest)

;plot,vtest,vel2a-res2,xst=1




END
