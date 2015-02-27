; simulates polarization

PRO temp6,angle,integral,contrast,thru

FSR      = [172.*8.,172.*16.,172.*32.,172.*64.,172.*4.]/1000. ; in A for the Lyot filter
; in the order E2, E3, E4, E5, and E1
nl       = 8001
l        = FINDGEN(nl)/FLOAT(nl-1)*16.*FSR[4]-8.*FSR[4]; wavelengths in A

;nangle   = 100 ; even number to have the zero !!!!
;angle    = FINDGEN(nangle)/FLOAT(nangle)*2.*!pi+!pi/2.
nangle   = 3
angle    = [0.,120.,240.]*!pi/180.
integral = FLTARR(nangle)
contrast = FLTARR(nl)
thru     = integral
a0       = WHERE(angle EQ 0.0 OR angle EQ 2.*!pi)
IF(a0[0] EQ -1) THEN STOP
Ep       = COMPLEXARR(nangle,nl)
Eo12     = Ep
Ee12     = Ep
Ep2      = Ep


angle0   =   0.*!pi/180. ; angle of the entrance polarizer of E1. should be zero
angle12  =   11.75*!pi/180.;11.75*!pi/180. ; angle of E1 1/2-wave plate, should be zero or 180.
retshift =   1. ; variation in the retardance of the E1 1/2-wave plate


; WE ADD THE BLOCKER FILTER + FRONT WINDOW
;--------------------------------------------------------------------------------------------

A            = FLTARR(nl)+1.0
RESTORE,'frontwindow.bin'
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433,l);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,l);,/LSQUADRATIC)
A            = A*SQRT(blocker0) ; because blocker0 is for the intensity


; COMPLETE LYOT FILTER
;----------------------------------------------------------------------------------------------

; Lyot entrance 1/4-wave plate (we assume unpolarized light arrives)
Eo14 =  A*EXP(+COMPLEX(0,1)*!pi/4.)
Ee14 =  A*EXP(-COMPLEX(0,1)*!pi/4.)

FOR j=0,4 DO BEGIN

  ; entrance polarizer
    IF j EQ 0 THEN A    = 1./SQRT(2.)*Eo14+1./SQRT(2.)*Ee14 
    IF j EQ 4 THEN A    = -COS(angle0)*Ee14 + SIN(angle0)*Eo14
    IF j EQ 1 OR j EQ 2 OR j EQ 3 THEN A    = -Ee14

  ; first calcite bloc
    IF j EQ 4 Then angleA = angle0 ELSE angleA = 0.0
    Ef   = COS(!pi/4.+angleA) * A   * EXP(+COMPLEX(0,1)*!pi*l/FSR[j]/2.) ; fast axis of calcite
    Es   = COS(!pi/4.-angleA) * A   * EXP(-COMPLEX(0,1)*!pi*l/FSR[j]/2.) ; slow axis of calcite

  ; half-wave plate

    IF j EQ 4 Then angleA =  angle12 ELSE angleA = 0.0
    IF j EQ 4 Then ret = retshift ELSE ret = 1.0
    Eo   = ( COS(!pi/4.+angleA) * Ef + SIN(!pi/4.+angleA) * Es) * EXP(+COMPLEX(0,1)*!pi/2.*ret)
    Ee   = (-SIN(!pi/4.+angleA) * Ef + COS(!pi/4.+angleA) * Es) * EXP(-COMPLEX(0,1)*!pi/2.*ret)

  ; second calcite bloc

    Ef   = ( COS(!pi/4.-angleA) * Eo + SIN(!pi/4.-angleA) * Ee) * EXP(+COMPLEX(0,1)*!pi*l/FSR[j]/2.)
    Es   = (-SIN(!pi/4.-angleA) * Eo + COS(!pi/4.-angleA) * Ee) * EXP(-COMPLEX(0,1)*!pi*l/FSR[j]/2.)

  ; quarter-wave plate
    IF j EQ 3 THEN ret2=1.0 ELSE ret2=1.0
    Eo14 = 1./SQRT(2.) * ( Ef + Es) * EXP(+COMPLEX(0,1)*!pi/4.*ret2) ; achromatic 1/4 wave plate
    Ee14 = 1./SQRT(2.) * (-Ef + Es) * EXP(-COMPLEX(0,1)*!pi/4.*ret2)

ENDFOR

; if polarizer
FOR i=0,nangle-1 DO Ep[i,*]   = SIN(angle[i])*Eo14 - COS(angle[i])*Ee14

; if 1/2-wave plate
FOR i=0,nangle-1 DO BEGIN
Eo12[i,*] = ( SIN(angle[i])*Eo14 - COS(angle[i])*Ee14 ) * EXP(+COMPLEX(0,1)*!pi/2.)
Ee12[i,*] = ( COS(angle[i])*Eo14 + SIN(angle[i])*Ee14 ) * EXP(-COMPLEX(0,1)*!pi/2.)
Ep2[i,*]  = ( COS(angle[i])*Eo12[i,*] - SIN(angle[i])*Ee12[i,*] ) ; waveplate followed by polarizer
ENDFOR

FOR i=0,nangle-1 DO BEGIN
    integral[i]=TOTAL(ABS(Ep[i,*])^2.0)*(l[1]-l[0])
;    contrast[i]=MAX(ABS(Ep[i,*])^2.0)/MEAN(ABS(Ep[i,*])^2.0)-1.0
    thru[i]    =MAX(ABS(Ep[i,*])^2.0)
ENDFOR

; TO REPRODUCE LASER DATA
FOR i=0,nl-1 DO BEGIN
    lam0     = 6173.3433+l[i]
    mat      = DBLARR(2,2)
    vec      = DBLARR(1,2)
    mat[0,0] =  COS(2.0*!dpi/FSR[4]*lam0)
    mat[1,0] = -SIN(2.0*!dpi/FSR[4]*lam0)
    mat[0,1] =  COS(2.0*!dpi/FSR[4]*lam0 + 2.0*!dpi/3.d0)
    mat[1,1] = -SIN(2.0*!dpi/FSR[4]*lam0 + 2.0*!dpi/3.d0)
    mat[*,*] =  LA_INVERT(REFORM(mat[*,*]),/DOUBLE)
    
    tot      = ABS(Ep[0,i])^2.0 + ABS(Ep[1,i])^2.0 + ABS(Ep[2,i])^2.0
    vec[0,0] = ABS(Ep[0,i])^2.0*3.d0/tot-1.d0
    vec[0,1] = ABS(Ep[1,i])^2.0*3.d0/tot-1.d0
    vec      = REFORM(mat[*,*])##vec
    contrast[i]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)
ENDFOR


; WE PLOT THE LYOT TRANSMITTANCE
;--------------------------------------------------------------------------------------------

;WINDOW,0,RETAIN=2,xsize=700,ysize=900
LOADCT,4
SET_PLOT,'ps'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color
!p.multi=[0,1,2]

PLOT,l,ABS(Ep[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[0,FSR[4]],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
oplot,[0,0],[0,1]
oplot,l,blocker0,linestyle=2,thick=3
FOR i=0,nangle-1 DO OPLOT,l,ABS(Ep[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2

PLOT,l,ABS(Ep2[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[0,FSR[4]],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
oplot,[0,0],[0,1]
FOR i=0,nangle-1 DO OPLOT,l,ABS(Ep2[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2


; WE ADD THE WB MICHELSON
;--------------------------------------------------------------------------------------------

FSRM         = 0.172*2.

retardquartv = !dpi / 4.d0 ; plate in the vacuum leg
retardquarts = !dpi / 4.d0 ; plate in the solid leg
retardquartM = !dpi / 4.d0 ; output quarter-wave plate
alpha        = 1.d0        ; reflectivity of the vacuum leg mirror
beta         = 1.d0        ; reflectivity of the solid leg mirror
Rs           = 1.d0
Ts           = 0.d0
Rp           = 0.d0
Tp           = 1.d0
angle0       = !dpi/4.d0  ; orientation of the input polarizer
angle1       = !dpi/4.d0  ; orientation of the vacuum leg waveplate
angle2       = !dpi/4.d0  ; orientation of the solid leg waveplate
angle3       = !dpi/4.d0  ; orientation of the output waveplate
Tpara0       = 1.d0
Tperp0       = 0.d0
angleout     = 0.0       ; for tuning: angle relative to the fast axis of the output 1/4-waveplate 

Trtemp       = FLTARR(nangle,nl)

FOR iii=0,nangle-1 DO BEGIN

    A = REFORM(Ep2[iii,*])
  
; after the input polarizer
        Epara0 = COS(angle0)*A*Tpara0-SIN(angle0)*A*Tperp0
        Eperp0 = SIN(angle0)*A*Tpara0+COS(angle0)*A*Tperp0

; after beamsplitter in vacuum leg
        Eperp  = Rs * Eperp0
        Epara  = Rp * Epara0
; after two passages through the 1/4-wave plate (THE !dpi PHASE SHIFT
; IS DUE TO THE MIRRORS)
        E0v    = (Epara*COS(angle1) + Eperp*SIN(angle1)) * exp(DCOMPLEX(0,1)*  (2.d0*retardquartv+!dpi))* exp(DCOMPLEX(0,1) *!pi*l/FSRM)*alpha
        Eev    =(-Epara*SIN(angle1) + Eperp*COS(angle1)) * exp(DCOMPLEX(0,1)* (-2.d0*retardquartv+!dpi))* exp(DCOMPLEX(0,1) *!pi*l/FSRM)*alpha
; after second passage through beamsplitter
        Eperpv = (E0v*SIN(angle1)   + Eev*COS(angle1))*Ts
        Eparav = (E0v*COS(angle1)   - Eev*SIN(angle1))*Tp
        
;after beamsplitter in solid leg
        Eperp  = Ts * Eperp0
        Epara  = Tp * Epara0
; after two passages through the 1/4-wave plate
        E0s    = (Epara*COS(angle2) + Eperp*SIN(angle2)) * exp(DCOMPLEX(0,1)*  (2.d0*retardquarts+!dpi))* exp(-DCOMPLEX(0,1)*!pi*l/FSRM)*beta
        Ees    =(-Epara*SIN(angle2) + Eperp*COS(angle2)) * exp(DCOMPLEX(0,1)* (-2.d0*retardquarts+!dpi))* exp(-DCOMPLEX(0,1)*!pi*l/FSRM)*beta
; after second passage through beamsplitter
        Eperps = (E0s*SIN(angle2)   + Ees*COS(angle2))*Rs
        Eparas = (E0s*COS(angle2)   - Ees*SIN(angle2))*Rp

; after the passage through the exit quarter-wave plate
        E0    = ((Eperps+Eperpv)*SIN(angle3) + (Eparas+Eparav)*COS(angle3)) * exp( DCOMPLEX(0,1) * retardquartM)
        Ee    =(-(Eparas+Eparav)*SIN(angle3) + (Eperps+Eperpv)*COS(angle3)) * exp(-DCOMPLEX(0,1) * retardquartM)
        
; output polarizer oriented at angleout degrees
;Trtemp = ABS(E0)^2.d0
        Trtemp[iii,*] = ABS(E0*COS(angleout)+Ee*SIN(angleout))^2.d0
        
    ENDFOR

; WE ADD THE NB MICHELSON CONSIDERED PERFECT
;-------------------------------------------------------------------------------------------

FOR i=0,nangle-1 DO Trtemp[i,*]=Trtemp[i,*]*COS(!pi*l/0.172)^2.0


; WE PLOT THE COMPLETE FILTER TRANSMITTANCE
;--------------------------------------------------------------------------------

;PLOT,l,Trtemp[a0[0],*],xst=1,thick=3,yrange=[0,1],yst=1,tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
;FOR i=0,nangle-1 DO OPLOT,l,Trtemp[i,*],color=i*256./FLOAT(nangle),thick=2


DEVICE,/close
SET_PLOT,'X'
!P.MULTI=0

END
