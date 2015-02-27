; THIS PROGRAM IS DESIGNED TO CHECK WHICH FILTERGRAM ORDER (FRAMELIST)
; PROVIDES THE LEAST NOISE ON THE OBSERVABLES
; THE NOISE IS ASSUMED TO COME FROM 2 MAIN SOURCES:
; 1)  FROM ONE FILTERGRAM TO THE OTHER THE SOLAR SURFACE HAS ROTATED A
; LITTLE BIT
; 2) THE FILTERGRAMS NEED TO BE INTERPOLATED IN TIME
; 
; FOR THE TEMPORAL INTERPOLATION, I USE accel.pro PROVIDED BY JESPER SCHOU
; AND WHICH DOES A WIENER INTERPOLATION
; FOR THE SOLAR ROTATION I USE ~richard/idlwork/hmi/level1/rotation.pro
; PROVIDED BY RICHARD WACHTER
;
; THE DIFFERENCE WITH HMI_framelist.pro IS THAT HERE I USE SINES
; INSTEAD OF MDI DATA
;-------------------------------------------------------------------------

; SHIFT IN THE FOURIER DOMAIN
FUNCTION decalage,f,lag,dx

; FUNCTION f MUST BE PERIODICAL (f[0]=f[N-1])
; f MUST HAVE AN EVEN NUMBER OF POINTS

N=N_ELEMENTS(f)
k=DINDGEN(N)/DOUBLE(N)/dx
k[N/2+1:N-1]=-REVERSE(k[1:N/2-1])
f2=FFT((FFT(f,/DOUBLE)*exp(-2.0*!dpi*DCOMPLEX(0,1)*k*lag)),/INVERSE,/DOUBLE)

RETURN,DOUBLE(f2)

END


;-------------------------------------------------------------------------
;
;
; MAIN PROGRAM
;
;-------------------------------------------------------------------------


PRO HMI_framelist2

lam0    = 6173.3433d0
ntune   = 6           ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0   ; Set number of tuning positions over wnarrow. (1/spacing)
dvdb    = 1./2.13d12*lam0*2.5d0*299792458.d0  ; conversion from magnetic field to Doppler velocity
ntime   = 720.        ; duration of observation in minutes
nframe  = ntune*2     ; total number of filtergrams taken per observable (each filter is used twice
                      ; because of LCP and RCP) 
cadence = 45.d0
dt      = cadence/nframe ; time interval between 2 filtergrams, or integration time (?)
LANDE   = 1.0/(2.0*4.67e-5*0.000061733433*2.5*299792458.0) ; Lande factor
dtime   = 60.d0/dt
nt      = ntime*dtime        ; total number of filtergrams
nt2     = nt/nframe   ; total number of observables
time    = FINDGEN(nt)*dt

; Parameters for the solar line
;------------------------------------------------------------------------

nlam       = 21500.
dlam       = 1.d0/1.75d3
;nlam       = 45000.
;dlam       = 1.d0/4.0d3
lam        = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam


; MODEL OF SOLAR LINE
;-------------------------------------------------------------------------
width  = 0.102d0/2.d0/SQRT(ALOG(2.d0)) ; w in A^2 (part added by seb)
depth  = 0.62d0
line  = 1.d0-depth*EXP(-lam^2.d0/width^2.d0)

; HMI FILTER TRANSMISSION PROFILES
;--------------------------------------------------------------------------
 
;phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
phase      = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
;contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
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


; FRONT WINDOW AND BLOCKER FILTER TRANSMISSION PROFILES
;RESTORE,'frontwindow.bin' ; front window
;blocker    = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
;q          = READFITS('blocker11.fits')
;blocker    = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) 

; theoretical blocker filter
blocker    = EXP(-lam^2.d0/2.d0/(8.d0/2.354820d0)^2.d0)

; NON-TUNABLE PROFILE
lyot       = blocker
FOR i = 0,3 DO lyot = lyot*(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0

; TUNABLE PROFILE
cmich      = 2.d0*!dpi/wmich
filters    = DBLARR(nlam,ntune*2)
FOR itune  = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i  = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
    filters[*,itune+ntune]=filters[*,itune]
ENDFOR



; WE DEFINE THE FRAMELIST
;------------------------------------------------------------------------------------------

; the filters are in this order:
;--------------------------------------------------------
;filters[*,0] is F0 centered on +171 mA
;filters[*,1] is F1 centered on +103 mA
;filters[*,2] is F2 centered on +034 mA
;filters[*,3] is F3 centered on -034 mA
;filters[*,4] is F4 centered on -103 mA
;filters[*,5] is F5 centered on -171 mA
;--------------------------------------------------------

;tframelist='F0L,F0R,F1R,F1L,F2L,F2R,F3R,F3L,F4L,F4R,F5R,F5L'
;framelist=[0,6,7,1,2,8,9,3,4,10,11,5]

;tframelist='F2L,F2R,F3R,F3L,F4L,F4R,F1R,F1L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,4,10,7,1,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F1L,F1R,F4R,F4L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,1,7,10,4,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F5L,F5R,F0R,F0L,F1L,F1R,F4R,F4L'
;framelist=[2,8,9,3,5,11,6,0,1,7,10,4]

;tframelist='F0L,F0R,F5R,F5L,F1L,F1R,F4R,F4L,F2L,F2R,F3R,F3L'
;framelist=[0,6,11,5,1,7,10,4,2,8,9,3]

;tframelist='F0L,F0R,F2R,F2L,F4L,F4R,F1R,F1L,F3L,F3R,F5R,F5L'
;framelist=[0,6,8,2,4,10,7,1,3,9,11,5]

;tframelist='F2L,F2R,F3R,F3L,F1L,F1R,F4R,F4L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,1,7,10,4,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F0L,F0R,F5R,F5L,F1L,F1R,F4R,F4L'
;framelist=[2,8,9,3,0,6,11,5,1,7,10,4]

;tframelist='F3L,F3R,F2R,F2L,F5L,F5R,F0R,F0L,F4L,F4R,F1R,F1L'
;framelist=[3,9,8,2,5,11,6,0,4,10,7,1]

;tframelist='F3L,F3R,F2L,F2R,F5L,F5R,F0L,F0R,F4L,F4R,F1L,F1R'
;framelist=[3,9,2,8,5,11,0,6,4,10,1,7]

;tframelist='F3L,F3R,F2R,F2L,F5L,F5R,F0R,F0L,F4L,F4R,F1R,F1L'
;framelist=[3,9,8,2,5,11,6,0,4,10,7,1]

;tframelist='F2L,F2R,F4R,F4L,F1L,F1R,F3R,F3L,F0L,F0R,F5R,F5L'
;framelist=[2,8,10,4,1,7,9,3,0,6,11,5]

;tframelist='F3L,F3R,F4R,F4L,F0L,F0R,F5R,F5L,F1L,F1R,F2R,F2L'
;framelist=[3,9,10,4,0,6,11,5,1,7,8,2]

tframelist='F3L,F3R,F1L,F1R,F5L,F5R,F0L,F0R,F4L,F4R,F2L,F2R'
framelist=[3,9,1,7,5,11,0,6,4,10,2,8]
;framelist=[3,9,7,1,5,11,6,0,4,10,8,2]


; LOOK-UP TABLES WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

ntest      = 501 ; Set test velocities
dvtest     = dlam/dlamdv;100.d0 ; in m/s
vtest      = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
lines      = DBLARR(nlam,ntest)
FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+(vtest[i])*dlamdv)
inten      = TRANSPOSE(filters[*,0:ntune-1])#lines

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


; DOPPLER VELOCITY DATA
;----------------------------------------------------------------------------

nu = 0.0035 ; temporal frequency in Hz
A  = 1200.d0; m/s

; p-modes (5-minute oscillations)
;phase=!dpi/4.3256d0
phase= 0.d0
velocity = A*SIN(2.d0*!dpi*nu*time+phase)

; SDO orbital velocity ; between -3.1 to +3.1 km/s in 12 hours (I
; assume a circular orbit)
phase=0.d0;!pi;!dpi/3.d0
vSDO=0.d0;3079.89d0*COS(2.d0*!dpi*time/24.d0/3600.d0+phase)
;vSDO=FINDGEN(nt)/(nt-1.d0)*13000.20-6500.d0

; solar rotation
offset= 0.d0 ; disk center
;offset= 2000 ; m/s at the limb

velocity=velocity+vSDO+offset

; INTERPOLATION COEFFICIENTS
; interpolation coefficients were obtained from accel.pro
RESTORE,'interpolation_coefficients.bin'
dtint= tint[1]-tint[0]
nint = N_ELEMENTS(tint)


inten    = DBLARR(nt2,nframe)
inten2   = inten
FOR t=0,nt2-1 DO BEGIN          ; observable index    
    FOR j=0,nframe-1 DO BEGIN
        velo=-velocity[t*nframe+j]*dlamdv
        line2 = decalage(line,velo,dlam)
        inten[t,j] = TOTAL(filters[*,framelist[j]]*line2)
    ENDFOR
ENDFOR

v1LCP= DBLARR(nt2)
v2LCP= v1LCP
v1RCP= v1LCP
v2RCP= v1LCP

; WE LOOP OVER THE TIME AND SPACE COORDINATES
;------------------------------------------------------------------------------


FOR t=3,nt2-3 DO BEGIN ; we start at 2 because with the interpolation algorithm we need 2 before and 3 after

    PRINT,t

; WE CORRECT THE FILTERGRAMS FOR THE FACT THAT THEY ARE NOT TAKEN AT
; THE EXACT SAME TIME AND THAT BECAUSE OF SOLAR ROTATION THEY ARE NOT
; TAKEN AT THE EXACT SAME LOCATION
;------------------------------------------------------------------------------------

; WE INTERPOLATE IN TIME
; the first frame (filtergram) is taken at relative time t=0
; the second frame is taken at t=dt
; the third at t=2*dt, and so on....


    FOR j=0,nframe-1 DO BEGIN
         dtime=j*dt
         dtime=ROUND(dtime/dtint)
         ; WE RE-INTERPOLATE ONTO THE FIRST FRAME, i.e. AT TIMES t=0,45,90,... seconds
          inten2[t,j] = coeffjx[0,nint-1-dtime]*inten[t-3,j]+coeffjx[1,nint-1-dtime]*inten[t-2,j]+coeffjx[2,nint-1-dtime]*inten[t-1,j]+coeffjx[3,nint-1-dtime]*inten[t,j]+coeffjx[4,nint-1-dtime]*inten[t+1,j]+coeffjx[5,nint-1-dtime]*inten[t+2,j]
         ;inten2[t,j]=INTERPOL(inten[*,j],FINDGEN(nt2)*cadence+dt*j,t*cadence,/SPLINE)
         ; WE RE-INTERPOLATE ONTO THE SIXTH FRAME, i.e. AT TIMES t=22.5,67.5,112.5,... seconds
    ENDFOR

; MDI-LIKE ALGORITHM TO OBTAIN THE DOPPLER VELOCITY
;-----------------------------------------------------------------------------------

    temp=SORT(framelist)

    ; LCP
    intent=REFORM(inten2[t,temp[0:nframe/2-1]]) ; depends on framelist definition
    c1   = TOTAL(COS(x)*intent)
    s1   = TOTAL(SIN(x)*intent)
    c2   = TOTAL(COS(2.d0*x)*intent)
    s2   = TOTAL(SIN(2.d0*x)*intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel1at=vel1t
    vel2t = phi2t*pv2/2.d0/!dpi
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1LCP[t] = INTERPOL(vtest,vel1a,vel1at,/SPLINE) ; velocity obtained from look-up table
    v2LCP[t] = INTERPOL(vtest,vel2a,vel2at,/SPLINE) ; velocity obtained from look-up table

    ; RCP
    intent=REFORM(inten2[t,temp[nframe/2:nframe-1]]) ; depends on framelist definition
    c1   = TOTAL(COS(x)*intent)
    s1   = TOTAL(SIN(x)*intent)
    c2   = TOTAL(COS(2.d0*x)*intent)
    s2   = TOTAL(SIN(2.d0*x)*intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel1at=vel1t
    vel2t = phi2t*pv2/2.d0/!dpi
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1RCP[t] = INTERPOL(vtest,vel1a,vel1at,/SPLINE) ; velocity obtained from look-up table
    v2RCP[t] = INTERPOL(vtest,vel2a,vel2at,/SPLINE) ; velocity obtained from look-up table


    
ENDFOR
PRINT,'DONE'

velocity=velocity[INDGEN(nt2)*nframe] ; re-interpolation at times 0, 45, 90...

velocity = velocity[3:nt2-3]

time=time[INDGEN(nt2)*nframe]
time=2.d0*!dpi*nu*time[3:nt2-3]

estvel = (v1LCP+v1RCP+v2LCP+v2RCP)/4.d0 ; estimated velocity
estvel = estvel[3:nt2-3]

SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xoffset=1,yoffset=0,xsize=20,ysize=26,/color,bits=24
LOADCT,3
!P.MULTI=[0,2,2]

PLOT,velocity,estvel-velocity,psym=1,xst=1,tit='!17'+tframelist,xtit='Velocity (m/s)',ytit='Error (m/s)',charsize=1.,yrange=[-15,15],yst=1

PLOT,FINDGEN(nt2-5)*cadence/60.d0,velocity,tit='!17',xtit='Time (min)',ytit='Velocity (m/s)',charsize=1.,yst=1,xst=1


;t=sort(velocity)
;PLOT,velocity[t],estvel[t]-velocity[t],xst=1,tit='!17'+tframelist,xtit='Velocity (m/s)',ytit='Error (m/s)',charsize=1.5,yrange=[-15,15],yst=1

res=HISTOGRAM(estvel-velocity,BINSIZE=.25,LOCATIONS=xh,MAX=25,MIN=-25)
;histo,estvel-velocity,-50,50,1.
PLOT,xh,res,charsize=1.,tit='!17'+STRING(MIN(estvel-velocity))+STRING(MAX(estvel-velocity)),xtit='Error (m/s)',ytit='Histogram',xst=1

;histo,(estvel-velocity)/velocity,-0.04,0.04,.001
;res=HISTOGRAM((estvel-velocity)/velocity,BINSIZE=.00025,LOCATIONS=xh,MAX=0.04,MIN=-0.04)
;PLOT,xh,res,charsize=1.,tit='!17'+STRING(SIGMA(estvel-velocity))+STRING(MEAN(estvel-velocity)),xtit='Relative error (m/s)',ytit='Histogram',xst=1

plot,time MOD (2.d0*!dpi),estvel-velocity,xst=1,psym=4,tit='!17',xtit='Phase (rad)',ytit='Error (m/s)',charsize=1.,yrange=[-15,15],yst=1



DEVICE,/CLOSE

READ,pause

END
