;----------------------------------------------------------------------------------------------
;
; simulates the HMI Lyot filter
;
; WITH THE CORRECT MODEL FOR THE LYOT ELEMENTS (MODEL OF 2005, NOT 2004)
;----------------------------------------------------------------------------------------------


PRO E1Model,angle,integral,contrast,thru,angle0,angle12,retshift,angle14,retshift2,angleD,angleE

; PROBLEMS SIMULATED
angle0   =  angle0 *!pi/180.   ; angle of the entrance polarizer of E1. should be zero (vertical direction)
angle12  =  angle12*!pi/180. ; angle of E1 1/2-wave plate, should be zero or 180 (vertical direction).
;retshift =   1.              ; error in the retardance of the E1 1/2-wave plate
angle14  =  angle14*!pi/180. ; E1 1/4-wave plate angle compared to 90 degrees. Should be close to 0
;retshift2=  1.-2.*0.190986   ; E1 1/4-wave plate retardance error
angleD   =  angleD *!pi/180.   ; misalignment of the 1st calcite+ADP bloc of E1
angleE   =  angleE *!pi/180.   ; misalignment of the 2nd calcite+ADP bloc of E1




; theoretical FSRs
;FSR      = [172.*8.,172.*16.,172.*32.,172.*64.,172.*4.]/1000. ; in A for the Lyot filter

; measured FSRs
FSR      = [1.405d0, 2.779d0, 5.682d0, 11.354d0, 0.7039d0]    ; in the order E2, E3, E4, E5, and E1

nl       = 8001
l        = FINDGEN(nl)/FLOAT(nl-1)*16.*FSR[4]-8.*FSR[4]; wavelengths in A

;nangle   = 8 ; even number to have the zero !!!!
;angle    = FINDGEN(nangle)/FLOAT(nangle-1)*2.*!pi ; angle of the polarizer relative to the 1/4-wave plate fast axis
nangle   = 5
angle    = -[0.,!pi/8.,!pi/4.,3.*!pi/8.,!pi/2.]
integral = FLTARR(nangle)
;contrast = FLTARR(nl)
contrast = integral
thru     = integral
a0       = WHERE(ABS(angle) LT !pi/180.)
IF(a0[0] EQ -1) THEN STOP
Ep       = COMPLEXARR(nangle,nl)
Epp      = Ep
Eo12     = Ep
Ee12     = Ep
Ep2      = Ep



; WE ADD THE BLOCKER FILTER + FRONT WINDOW
;--------------------------------------------------------------------------------------------

q            = READFITS('blocker11.fits')
blocker0     = FLTARR(nl)+1.d0
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,l);,/LSQUADRATIC)
A            = SQRT(blocker0) ; because blocker0 is for the intensity

; COMPLETE LYOT FILTER
;----------------------------------------------------------------------------------------------

; Lyot entrance 1/4-wave plate (we assume unpolarized light arrives)
Eo14 =  A*EXP(+COMPLEX(0,1)*!pi/4.)
Ee14 =  A*EXP(-COMPLEX(0,1)*!pi/4.)

FOR j=4,4 DO BEGIN ; start from 4 if E1 alone

  ; entrance polarizer
    IF j EQ 0 THEN A    = 1./SQRT(2.)*Eo14+1./SQRT(2.)*Ee14 
    IF j EQ 4 THEN A    = FLTARR(nl)+1.0; -COS(angle0)*Ee14 + SIN(angle0)*Eo14
    IF j EQ 1 OR j EQ 2 OR j EQ 3 THEN A    = -Ee14

  ; first calcite+ADP bloc
    IF j EQ 4 Then angleA = angle0 ELSE angleA = 0.0
    Ef   = COS(angleA-(angleD-!pi/4.)) * A   * EXP(+COMPLEX(0,1)*!pi*l/FSR[j]/2.) ; fast axis of calcite
    Es   = SIN(angleA-(angleD-!pi/4.)) * A   * EXP(-COMPLEX(0,1)*!pi*l/FSR[j]/2.) ; slow axis of calcite

  ; half-wave plate

    IF j EQ 4 Then angleA =  angle12 ELSE angleA = 0.0
    IF j EQ 4 Then ret = retshift ELSE ret = 1.0
    Eo   = ( COS(angleA-(angleD-!pi/4.)) * Ef + SIN(angleA-(angleD-!pi/4.)) * Es) * EXP(+COMPLEX(0,1)*!pi/2.*ret)
    Ee   = (-SIN(angleA-(angleD-!pi/4.)) * Ef + COS(angleA-(angleD-!pi/4.)) * Es) * EXP(-COMPLEX(0,1)*!pi/2.*ret)

  ; second calcite+ADP bloc
    IF j EQ 4 Then angleB =  angle12 ELSE angleB = 0.0
    Ef   = ( COS(!pi/4.+angleE-angleB) * Eo + SIN(!pi/4.+angleE-angleB) * Ee) * EXP(+COMPLEX(0,1)*!pi*l/FSR[j]/2.)
    Es   = (-SIN(!pi/4.+angleE-angleB) * Eo + COS(!pi/4.+angleE-angleB) * Ee) * EXP(-COMPLEX(0,1)*!pi*l/FSR[j]/2.)

  ; quarter-wave plate
    IF j EQ 4 THEN BEGIN
        ret2=retshift2 
        angleC = angle14
    ENDIF ELSE BEGIN
        ret2=1.0
        angleC = 0.0
    ENDELSE
    Eo14 = ( COS(!pi/4.+angleC-angleE) * Ef + SIN(!pi/4.+angleC-angleE) * Es) * EXP(+COMPLEX(0,1)*!pi/4.*ret2) ; achromatic 1/4 wave plate
    Ee14 = (-SIN(!pi/4.+angleC-angleE) * Ef + COS(!pi/4.+angleC-angleE) * Es) * EXP(-COMPLEX(0,1)*!pi/4.*ret2)

ENDFOR

; if polarizer (0 degrees is vertical position)
transmit=.9d0 ; transmissivity of the polarizer in the right polarization state.
FOR i=0,nangle-1 DO BEGIN
    Ep[i,*]   = transmit*(COS(!pi/2.-angle[i]+angle14)*Eo14 - SIN(!pi/2.-angle[i]+angle14)*Ee14)
    Epp[i,*]  = (1.d0-transmit)*(SIN(!pi/2.-angle[i]+angle14)*Eo14 + COS(!pi/2.-angle[i]+angle14)*Ee14)
ENDFOR

; if 1/2-wave plate
lastret=.15d0 ; to simulate a problem with the waveplate
FOR i=0,nangle-1 DO BEGIN
Eo12[i,*] = ( SIN(angle[i])*Eo14 - COS(angle[i])*Ee14 ) * EXP(+COMPLEX(0,1)*(!pi/2.+lastret))
Ee12[i,*] = ( COS(angle[i])*Eo14 + SIN(angle[i])*Ee14 ) * EXP(-COMPLEX(0,1)*(!pi/2.+lastret))
Ep2[i,*]  = ( COS(angle[i])*Eo12[i,*] - SIN(angle[i])*Ee12[i,*] ) ; waveplate followed by a perfect polarizer
ENDFOR


; FOR TUNING WITH POLARIZER
FOR i=0,nangle-1 DO BEGIN
    integral[i]=TOTAL(ABS(Ep[i,*])^2.0+ABS(Epp[i,*])^2.0)*(l[1]-l[0])
    contrast[i]=MAX(ABS(Ep[i,*])^2.0+ABS(Epp[i,*])^2.0)/MEAN(ABS(Ep[i,*])^2.0+ABS(Epp[i,*])^2.0)-1.0 ;NB: has no meaning if we consider all the Lyot stages
    thru[i]    =MAX(ABS(Ep[i,*])^2.0+ABS(Epp[i,*])^2.0)
ENDFOR

; FOR TUNING WITH WAVEPLATE
FOR i=0,nangle-1 DO BEGIN
    integral[i]=TOTAL(ABS(Ep2[i,*])^2.0)*(l[1]-l[0])
    contrast[i]=MAX(ABS(Ep2[i,*])^2.0)/MEAN(ABS(Ep2[i,*])^2.0)-1.0 ;NB: has no meaning if we consider all the Lyot stages
    thru[i]    =MAX(ABS(Ep2[i,*])^2.0)
ENDFOR



; TO REPRODUCE LASER DATA
;FOR i=0,nl-1 DO BEGIN
;    lam0     = 6173.3433+l[i]
;    mat      = DBLARR(2,2)
;    vec      = DBLARR(1,2)
;    mat[0,0] =  COS(2.0*!dpi/FSR[4]*lam0)
;    mat[1,0] = -SIN(2.0*!dpi/FSR[4]*lam0)
;    mat[0,1] =  COS(2.0*!dpi/FSR[4]*lam0 + 2.0*!dpi/3.d0)
;    mat[1,1] = -SIN(2.0*!dpi/FSR[4]*lam0 + 2.0*!dpi/3.d0)
;    mat[*,*] =  LA_INVERT(REFORM(mat[*,*]),/DOUBLE)
;    
;    tot      = ABS(Ep[0,i])^2.0 + ABS(Ep[1,i])^2.0 + ABS(Ep[2,i])^2.0
;    vec[0,0] = ABS(Ep[0,i])^2.0*3.d0/tot-1.d0
;    vec[0,1] = ABS(Ep[1,i])^2.0*3.d0/tot-1.d0
;    vec      = REFORM(mat[*,*])##vec
;    contrast[i]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)
;ENDFOR


; WE PLOT THE LYOT TRANSMITTANCE
;--------------------------------------------------------------------------------------------

;WINDOW,0,RETAIN=2,xsize=700,ysize=900
LOADCT,4
SET_PLOT,'ps'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color
!p.multi=[0,1,2]

PLOT,l,ABS(Ep[a0[0],*])^2.d0+ABS(Epp[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[-1.*FSR[4],1.*FSR[4]],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
oplot,[0,0],[0,1]
contras=1.0;0.9395d0
OPLOT,l,transmit^2.d0*contras*COS(!dpi/FSR[4]*l)^2.d0+(1.d0-transmit)^2.d0*contras*SIN(!dpi/FSR[4]*l)^2.d0+transmit^2.d0*(1.d0-contras)/2.d0,linestyle=2,thick=3
;oplot,l,blocker0,linestyle=2,thick=3
FOR i=1,nangle-1 DO BEGIN
    OPLOT,l,ABS(Ep[i,*])^2.d0+ABS(Epp[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2
    OPLOT,l,transmit^2.d0*contras*COS(!dpi/FSR[4]*l+angle[i])^2.d0+(1.d0-transmit)^2.d0*contras*SIN(!dpi/FSR[4]*l+angle[i])^2.d0+transmit^2.d0*(1.d0-contras)/2.d0,linestyle=2,thick=3,color=i*256./FLOAT(nangle)
ENDFOR

PLOT,l,ABS(Ep2[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[-1.*FSR[4],1.*FSR[4]],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
oplot,[0,0],[0,1]
OPLOT,l,(1.d0-lastret^2.d0/2.d0)^2.d0*contras*COS(!dpi/FSR[4]*l)^2.d0+lastret^2.d0*contras*COS(!dpi/FSR[4]*l)^2.d0+(1.d0-contras)/2.d0,linestyle=2,thick=3
FOR i=1,nangle-1 DO BEGIN
    OPLOT,l,ABS(Ep2[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2
    OPLOT,l,(1.d0-lastret^2.d0/2.d0)^2.d0*contras*COS(!dpi/FSR[4]*l+2.d0*angle[i])^2.d0+lastret^2.d0*contras*COS(!dpi/FSR[4]*l)^2.d0+(1.d0-contras)/2.d0,linestyle=2,thick=3,color=i*256./FLOAT(nangle)
ENDFOR

DEVICE,/close
SET_PLOT,'X'
!P.MULTI=0

GOTO,endo

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




endo:

END
