; simulates polarization

PRO temp7,angle,integral

FSR  = 0.688        ; in A
nl   = 8001
l    = FINDGEN(nl)/FLOAT(nl-1)*16.*FSR-8.*FSR; wavelengths in A

nangle = 100 ; even number to have the zero !!!!
angle  = FINDGEN(nangle)/FLOAT(nangle)*2.*!pi+!pi/2.*0.9
integral = FLTARR(nangle)
a0   = WHERE(ABS(angle) LE 2.*!pi/FLOAT(nangle) OR ABS(angle-2.*!pi) LE 2.*!pi/FLOAT(nangle))
IF(a0[0] EQ -1) THEN STOP
Ep   = COMPLEXARR(nangle,nl)
;Eo14 = Ep
;Ee14 = Ep
;Ef   = Ep
;Es   = Ep
;Eo   = Ep
;Ee   = Ep

angle14 = (-!pi/4.)-4.1*!pi/180. ; must be close to -!pi/4
ret     = 1.0

; entrance light
A   = 1.0

; quarter-wave plate
Eo14 = EXP(+COMPLEX(0,1)*!pi/4.)*A ; achromatic 1/4 wave plate
Ee14 = EXP(-COMPLEX(0,1)*!pi/4.)*A
; second calcite bloc
Ef   = 1./SQRT(2.0)*(Eo14 - Ee14) * EXP(+COMPLEX(0,1)*!pi*l/FSR/2.)
Es   = 1./SQRT(2.0)*(Eo14 + Ee14) * EXP(-COMPLEX(0,1)*!pi*l/FSR/2.)
; half-wave plate
Eo   = 1./SQRT(2.) * ( Ef - Es)     * EXP(+COMPLEX(0,1)*!pi/2.)
Ee   = 1./SQRT(2.) * ( Ef + Es)     * EXP(-COMPLEX(0,1)*!pi/2.)
; first calcite bloc
Ef   = 1./SQRT(2.) * (Eo - Ee)      * EXP(+COMPLEX(0,1)*!pi*l/FSR/2.) ; fast axis of calcite
Es   = 1./SQRT(2.) * (Eo + Ee)      * EXP(-COMPLEX(0,1)*!pi*l/FSR/2.) ; slow axis of calcite
; entrance polarizer
Ep2  = 1./SQRT(2.) * ( Ef + Es)

    
; Lyot entrance 1/4-wave plate
Eo14 =  COS(angle14)*Ep2*EXP(+COMPLEX(0,1)*!pi/4.*ret)
Ee14 =  SIN(angle14)*Ep2*EXP(-COMPLEX(0,1)*!pi/4.*ret)

; tuning polarizer
FOR i=0,nangle-1 DO Ep[i,*]   = COS(angle[i]+angle14) * Eo14 + SIN(angle[i]+angle14) * Ee14

FOR i=0,nangle-1 DO integral[i] = TOTAL(ABS(Ep[i,*])^2.0)*(l[1]-l[0])
PRINT,''

; BLOCKER FILTER + FRONT WINDOW
;--------------------------------------------------------------------------------------------

RESTORE,'frontwindow.bin'
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433,l);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]-6169.8,l);,/LSQUADRATIC)

; COMPLETE LYOT FILTER
;--------------------------------------------------------------------------------------------

FSRe=[172.*8.,172.*16.,172.*32.,172.*64.]/1000.
T = FLTARR(nl)+1.0
FOR i=0,3 DO T=T*COS(!pi*l/FSRe[i])^2.0
T = T *blocker0

; WE PLOT THE E1 TRANSMITTANCE
;--------------------------------------------------------------------------------------------

;WINDOW,0,RETAIN=2,xsize=700,ysize=900
LOADCT,4
SET_PLOT,'ps'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color
!p.multi=[0,1,2]

PLOT,l,ABS(Ep[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[0,FSR],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
oplot,[0,0],[0,1]
FOR i=0,nangle-1 DO OPLOT,l,ABS(Ep[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2

;PLOT,l,ABS(Ep2[a0[0],*])^2.0,xst=1,thick=3,yrange=[0,1],yst=1,xrange=[0,FSR],tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
;oplot,[0,0],[0,1]
;FOR i=0,nangle-1 DO OPLOT,l,ABS(Ep2[i,*])^2.0,color=i*256./FLOAT(nangle),thick=2


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

    A = REFORM(Ep[iii,*])
  
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

; WE add the NB Michelson CONSIDERD PERFECT
;-------------------------------------------------------------------------------

FOR i=0,nangle-1 DO Trtemp[i,*]=Trtemp[i,*]*COS(!pi*l/0.172)^2.0*T


; WE PLOT THE COMPLETE FILTER TRANSMITTANCE
;--------------------------------------------------------------------------------

PLOT,l,Trtemp[a0[0],*],xst=1,thick=3,yrange=[0,1],yst=1,tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
FOR i=0,nangle-1 DO OPLOT,l,Trtemp[i,*],color=i*256./FLOAT(nangle),thick=2

;PLOT,l,ABS(Ep2[a0[0],*])^2.0*T,xst=1,thick=3,yrange=[0,1],yst=1,tit='!17',xtitle='!7k!17 (mA)',ytit='Transmittance',charsize=1.5
;FOR i=0,nangle-1 DO OPLOT,l,ABS(Ep2[i,*])^2.0*T,color=i*256./FLOAT(nangle),thick=2

;FOR i=0,nangle-1 DO PRINT,'Integral =',TOTAL(ABS(Ep2[i,*])^2.0*T)*(l[1]-l[0])

DEVICE,/close
SET_PLOT,'X'
!P.MULTI=0

END
