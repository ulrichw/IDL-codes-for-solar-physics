; THIS PROGRAM TESTS THE IMPACT OF:
; ---imperfect polarizers and half-wave plates in the Lyot elements
; ---imperfect polarizers, beamsplitters, and waveplates in the
; Michelson
; NB: WAVEPLATES ARE CONSIDERED ACHROMATIC
; the quarter-wave plate of the Lyot stages is considered perfect
; (retardance and orientation)
; the calcite and ADP/KDP blocks are considered perfectly oriented
; the polarizing beamsplitter is considerer perfectly oriented

; the Michelsons can be tuned thru the variable "anglout"

; NB: THE MODEL FOR THE LYOT ELEMENTS IS NOT CORRECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PRO HMIdefects,angleout,returnedt,returnedc
; angleout is used to tune the Michelsons
; returned=FLTARR(2) is the quantities returned, throughput and
; contrast for the NB

set_plot,'ps'
device,file='yo.ps',xsize=22,ysize=26,xoffset=-0.75,yoffset=0.5
!p.multi =[0,3,3]

; LYOT FILTER
;------------------------------------------------------------------------

; LYOT element
nlam   = 2880*5*2
dlam   = .5d-10 ; in mm
l      = (DINDGEN(nlam)-(nlam-1)/2.)*dlam+0.0006173d0 ; wavelength in mm

FSR    = 0.189d0*6173.d0/6768.d0*4.d0*1.d-7 ; FSR in mm of the narrowest Lyot element

retardcalcite1 = DBLARR(5)
retardcalcite2 = DBLARR(5)
retardhalf1    = DBLARR(5)
retardhalf2    = DBLARR(5)
retardADPKDP1  = DBLARR(5)
retardADPKDP2  = DBLARR(5)
retardquart    = DBLARR(5)
angle0         = DBLARR(5)
angle1         = DBLARR(5)
angle2         = DBLARR(5)
Tpara          = DBLARR(5)
Tperp          = DBLARR(5)


; element E1
;-----------


d1                = 15.310d0*2.d0          ; calcite thickness
d2                = 3.528*2.d0             ; ADP thickness
Dn2               = -0.0449016d0           ; birefringence of ADP at 25 C and 6173 Angstrom
Dn1               = -((6173.d-7)^2.d0/FSR - d2 * Dn2)/d1 ; birefringence of calcite derived from the expected FSR
retardcalcite1[0] = !dpi * (d1/2.d0) * Dn1 ; retardance of the 1st calcite block
retardcalcite2[0] = retardcalcite1[0]      ; retardance of the 2nd calcite block
retardADPKDP1[0]  = !dpi *  (d2/2.d0) *Dn2 ; retardance of the 1st ADP block
retardADPKDP2[0]  = retardADPKDP1[0]       ; retardance of the 2nd ADP block

retardhalf1[0]    = !dpi / 2.d0            ; retardance of the 1st half-wave plate
retardhalf2[0]    = !dpi / 2.d0            ; retardance of the 2nd half-wave plate
retardquart[0]    = !dpi / 4.d0            ; retardance of the output quarter-wave plate
angle0[0]         = 0.d0                   ; orientation of the input polarizer (varies from 0 to !pi/4)
angle1[0]         = 0.d0                   ; orientation of the fast axis of the 1st waveplate (varies from 0 to !pi/4)
angle2[0]         = 0.d0                   ; orientation of the fast axis of the second waveplate (varies from 0 to !pi/4)
Tpara[0]          = 1.d0                   ; transmitivity in the direction parallel to the input polarizer axis
Tperp[0]          = 0.d0                  ; transmitivity in the direction perpendicular to the input polarizer axis

; element E2
;-----------

retardcalcite1[1] = retardcalcite1[0]/2.d0
retardcalcite2[1] = retardcalcite1[1]
retardADPKDP1[1]  = retardADPKDP1[0]/2.d0
retardADPKDP2[1]  = retardADPKDP1[1]

retardhalf1[1]    = !dpi / 2.d0
retardhalf2[1]    = !dpi / 2.d0
retardquart[1]    = !dpi / 4.d0
angle0[1]         = 0.d0
angle1[1]         = 0.d0
angle2[1]         = 0.d0
Tpara[1]          = 1.d0
Tperp[1]          = 0.d0

; element E3
;-----------

retardcalcite1[2] = retardcalcite1[1]/2.d0
retardcalcite2[2] = retardcalcite1[2]
retardADPKDP1[2]  = retardADPKDP1[1]/2.d0
retardADPKDP2[2]  = retardADPKDP1[2]

retardhalf1[2]    = !dpi / 2.d0
retardhalf2[2]    = !dpi / 2.d0
retardquart[2]    = !dpi / 4.d0
angle0[2]         = 0.d0
angle1[2]         = 0.d0
angle2[2]         = 0.d0
Tpara[2]          = 1.d0
Tperp[2]          = 0.d0

; element E4
;-----------

FSR               = FSR*8.d0
d1                = 2.331d0*2.d0
d2                = 2.25*2.d0
Dn2               = (d1 * Dn1 + (6173.d-7)^2.d0/FSR)/d2
retardcalcite1[3] = !dpi * (d1/2.d0) * Dn1
retardcalcite2[3] = retardcalcite1[3]
retardADPKDP1[3]  = !dpi *  (d2/2.d0) *Dn2
retardADPKDP2[3]  = retardADPKDP1[3]

retardhalf1[3]    = !dpi / 2.d0
retardhalf2[3]    = !dpi / 2.d0
retardquart[3]    = !dpi / 4.d0
angle0[3]         = 0.d0
angle1[3]         = 0.d0
angle2[3]         = 0.d0
Tpara[3]          = 1.d0
Tperp[3]          = 0.d0

; element E5
;-----------

retardcalcite1[4] = retardcalcite1[3]/2.d0
retardcalcite2[4] = retardcalcite1[4]
retardADPKDP1[4]  = retardADPKDP1[3]/2.d0
retardADPKDP2[4]  = retardADPKDP1[4]

retardhalf1[4]    = !dpi / 2.d0
retardhalf2[4]    = !dpi / 2.d0
retardquart[4]    = !dpi / 4.d0
angle0[4]         = 0.d0
angle1[4]         = 0.d0
angle2[4]         = 0.d0
Tpara[4]          = 1.d0
Tperp[4]          = 0.d0

;E0= DCOMPLEXARR(nlam)
;Ee= E0
Tr= 1.d0

FOR k=0,4 DO BEGIN

A = 1.d0 ; amplitude of incident wave after the input polarizer

CASE k OF
0: i = 1  ; element E2 comes first
1: i = 2  ; E3 
2: i = 3  ; E4
3: i = 4  ; E5
4: i = 0  ; E1
ENDCASE

; after the input polarizer
;IF(i NE 0) THEN A = -SIN(angle0[i])*E0-COS(angle0[i])*Ee
E0=A*COS(!dpi/4.d0-angle0[i])*Tpara[i]-A*SIN(!dpi/4.d0-angle0[i])*Tperp[i]
Ee=A*SIN(!dpi/4.d0-angle0[i])*Tpara[i]+A*COS(!dpi/4.d0-angle0[i])*Tperp[i]

; after the first calcite block
E0=E0*exp(DCOMPLEX(0,1)*(-retardcalcite1[i])/l)
Ee=Ee*exp(DCOMPLEX(0,1)*  retardcalcite1[i] /l)

; after the first 1/2-wave plate
Etemp = E0
E0 = (COS(!dpi/4.d0-angle1[i])*E0 + SIN(!dpi/4.d0-angle1[i])*Ee)   * exp(DCOMPLEX(0,1)*  retardhalf1[i])
Ee = (COS(!dpi/4.d0-angle1[i])*Ee - SIN(!dpi/4.d0-angle1[i])*Etemp)* exp(DCOMPLEX(0,1)*(-retardhalf1[i]))

; after the second calcite block
Etemp = E0
E0 = (COS(!dpi/4.d0+angle1[i])*E0 + COS(!dpi/4.d0-angle1[i])*Ee)   * exp(DCOMPLEX(0,1)*(-retardcalcite2[i])/l)
Ee = (COS(!dpi/4.d0+angle1[i])*Ee - COS(!dpi/4.d0-angle1[i])*Etemp)* exp(DCOMPLEX(0,1)*  retardcalcite2[i]/l)

; after the first ADP/KDP block
Etemp = E0
E0 = -Ee*exp(DCOMPLEX(0,1)*(-retardADPKDP1[i])/l)
Ee = Etemp*exp(DCOMPLEX(0,1)*retardADPKDP1[i] /l)

; after the second 1/2-wave plate
Etemp = E0
E0 = (COS(!dpi/4.d0-angle2[i])*E0 + SIN(!dpi/4.d0-angle2[i])*Ee)   * exp(DCOMPLEX(0,1)*  retardhalf2[i])
Ee = (COS(!dpi/4.d0-angle2[i])*Ee - SIN(!dpi/4.d0-angle2[i])*Etemp)* exp(DCOMPLEX(0,1)*(-retardhalf2[i]))

; after the second ADP/KDP block
Etemp = E0
E0 = (COS(!dpi/4.d0+angle2[i])*E0 + COS(!dpi/4.d0-angle2[i])*Ee)   * exp(DCOMPLEX(0,1)*(-retardADPKDP2[i])/l)
Ee = (COS(!dpi/4.d0+angle2[i])*Ee - COS(!dpi/4.d0-angle2[i])*Etemp)* exp(DCOMPLEX(0,1)*  retardADPKDP2[i]/l)

; after the 1/4-wave plate
Etemp = E0
E0 = (1.d0/SQRT(2.d0)*E0 + 1.d0/SQRT(2.d0)*Ee)   * exp(DCOMPLEX(0,1)*  retardquart[i])
Ee = (1.d0/SQRT(2.d0)*Ee - 1.d0/SQRT(2.d0)*Etemp)* exp(DCOMPLEX(0,1)*(-retardquart[i]))

; output polarizer oriented at 0 degrees
;Trtemp = ABS(COS(angle0[i])*Ee+SIN(angle0[i]*E0))^2.d0
Trtemp = ABS(Ee)^2.d0
IF i EQ 4 THEN Trtemp = ABS(E0)^2.d0

Tr = Tr * Trtemp
plot,l*1.d7,Trtemp,xst=1;,xrange=[6170,6176.5],tit='!17 E'+STRTRIM(STRING(i+1),1),xtit='Wavelength (Angstrom)',ytit='Transmittance',charsize=1.5
aa=WHERE(ABS(l*1.d7 - 6173.2) le .6)
maxi=MAX(Trtemp[aa])
centr=l[WHERE(ABS(l*1.d7 - 6173.2) le .6 AND Trtemp EQ maxi)]*1.d7
;PRINT,centr
ENDFOR

;Tr = ABS(E0/SQRT(2.d0)+Ee/SQRT(2.d0))^2.d0 ; transmittance
plot,l*1.d7,Tr,xst=1,xrange=[6170.5,6176],tit='!17 Lyot transmittance',xtit='Wavelength (Angstrom)',ytit='Transmittance',charsize=1.5


; MICHELSONS
;--------------------------------------------------------------------------

retardquartv = DBLARR(2)
retardquarts = DBLARR(2)
retardM      = DBLARR(2)
retardquartM = DBLARR(2)
alpha        = DBLARR(2)
beta         = DBLARR(2)
R            = beta
T            = beta
Rs           = beta
Ts           = beta
Rv           = beta
Tv           = beta
angle0       = beta
angle1       = beta
angle2       = beta
angle3       = beta
Tpara0       = beta
Tperp0       = beta
;angleout     = beta
FSRM         = beta

; Michelson M1 (narrow Michelson)
;--------------------------------

FSRM[0]         = 0.172d-7 ; expected FSR of the narrow Michelson
Dm              = (6173.d-7)^2.d0/FSRM[0]
retardM[0]      = !dpi * Dm   ; phase difference due to the path difference in the Michelson

retardquartv[0] = !dpi / 4.d0 ; plate in the vacuum leg
retardquarts[0] = !dpi / 4.d0 ; plate in the solid leg
retardquartM[0] = !dpi / 4.d0 ; output quarter-wave plate
alpha[0]        = 1.d0        ; reflectivity of the vacuum leg mirror
beta[0]         = 1.d0        ; reflectivity of the solid leg mirror
R[0]            = .9d0        ; reflectivity of the beamsplitter on the entrance side
T[0]            = .9d0        ; transmittivity of the beamsplitter on the entrance side
Rs[0]           = .9d0        ; reflectivity of the beamsplitter IN THE SOLID LEG COMING FROM THE MIRROR
Ts[0]           = .9d0        ; transmittivity of the beamsplitter IN THE SOLID LEG COMING FROM THE MIRROR
angle0[0]       = !dpi/4.d0   ; orientation of the input polarizer
angle1[0]       = !dpi/4.d0   ; orientation of the vacuum leg waveplate
angle2[0]       = !dpi/4.d0   ; orientation of the solid leg waveplate
angle3[0]       = !dpi/4.d0   ; orientation of the output waveplate
Tpara0[0]       = 1.d0        ; transmittivity of the input polarizer in the direction parallel to the optic axis
Tperp0[0]       = 0.d0
;angleout[0]     = !pi/6.    ; for tuning: angle relative to the fast axis of the output 1/4-waveplate ?

; Michelson M2 (wide Michelson)
;------------------------------

FSRM[1]         = FSRM[0]*2.d0
Dm              = (6173.d-7)^2.d0/FSRM[1]
retardM[1]      = !dpi * Dm   ; phase difference due to the path difference in the Michelson

retardquartv[1] = !dpi / 4.d0 ; plate in the vacuum leg
retardquarts[1] = !dpi / 4.d0 ; plate in the solid leg
retardquartM[1] = !dpi / 4.d0 ; output quarter-wave plate
alpha[1]        = 1.d0        ; reflectivity of the vacuum leg mirror
beta[1]         = 1.d0        ; reflectivity of the solid leg mirror
R[1]            = 1.d0
T[1]            = 1.d0
Rs[1]           = 1.d0        ; reflectivity of the beamsplitter IN THE SOLID LEG COMING FROM THE MIRROR
Ts[1]           = 1.d0        ; transmittivity of the beamsplitter IN THE SOLID LEG COMING FROM THE MIRROR
angle0[1]       = !dpi/4.d0  ; orientation of the input polarizer
angle1[1]       = !dpi/4.d0  ; orientation of the vacuum leg waveplate
angle2[1]       = !dpi/4.d0  ; orientation of the solid leg waveplate
angle3[1]       = !dpi/4.d0  ; orientation of the output waveplate
Tpara0[1]       = 1.d0
Tperp0[1]       = 0.d0
;angleout[1]     = 0.d0       ; for tuning: angle relative to the fast axis of the output 1/4-waveplate ?

A = 1.d0

FOR k=0,1 DO BEGIN

IF k EQ 0 THEN i=1 ELSE i=0

;    IF i NE 0 THEN BEGIN
;        A = COS(angle0[i]-angle3[i-1]) * E0 + COS(!dpi/2.d0 + angle0[i] - angle3[i-1])
;    ENDIF
; after the input polarizer
Epara0 = COS(angle0[i])*A*Tpara0[i]+SIN(angle0[i])*A*Tperp0[i]
Eperp0 = SIN(angle0[i])*A*Tpara0[i]-COS(angle0[i])*A*Tperp0[i]

; after beamsplitter in vacuum leg
Eperp = R[i]         * Eperp0
Epara = (1.d0-T[i])  * Epara0
; after two passages through the 1/4-wave plate (THE !dpi PHASE SHIFT
; IS DUE TO THE MIRRORS)
E0v   = (Epara*COS(angle1[i]) + Eperp*SIN(angle1[i])) * exp(DCOMPLEX(0,1)*  (2.d0*retardquartv[i]+!dpi))* exp(DCOMPLEX(0,1) *retardM[i]/l)*alpha[i]
Eev   = (Epara*SIN(angle1[i]) - Eperp*COS(angle1[i])) * exp(DCOMPLEX(0,1)* (-2.d0*retardquartv[i]+!dpi))* exp(DCOMPLEX(0,1) *retardM[i]/l)*alpha[i]

;after beamsplitter in solid leg
Eperp = (1.d0-R[i])   * Eperp0
Epara = T[i]          * Epara0
; after two passages through the 1/4-wave plate
E0s   = (Epara*COS(angle2[i]) + Eperp*SIN(angle2[i])) * exp(DCOMPLEX(0,1)*  (2.d0*retardquarts[i]+!dpi))* exp(-DCOMPLEX(0,1)*retardM[i]/l)*beta[i]
Ees   = (Epara*SIN(angle2[i]) - Eperp*COS(angle2[i])) * exp(DCOMPLEX(0,1)* (-2.d0*retardquarts[i]+!dpi))* exp(-DCOMPLEX(0,1)*retardM[i]/l)*beta[i]

;if i EQ 0 THEN read,pause
; after the second passage through the beamsplitter
Eperp = (E0s*SIN(angle2[i])-Ees*COS(angle2[i]))*Rs[i]+(1.-Ts[i])*(E0s*COS(angle2[i])+Ees*SIN(angle2[i]))
Epara = (E0v*COS(angle1[i])+Eev*SIN(angle1[i]))*T[i] +(1.-R[i]) *(E0v*SIN(angle1[i])-Eev*COS(angle1[i]))

; after the passage through the exit quarter-wave plate
E0    = (Eperp*SIN(angle3[i]) + Epara*COS(angle3[i])) * exp( DCOMPLEX(0,1) * retardquartM[i])
Ee    = (Epara*SIN(angle3[i]) - Eperp*COS(angle3[i])) * exp(-DCOMPLEX(0,1) * retardquartM[i])

; output polarizer oriented at angleout degrees
;Trtemp = ABS(E0)^2.d0
Trtemp = ABS(E0*COS(angleout[i])+Ee*SIN(angleout[i]))^2.d0

aa=WHERE(ABS(l - 6173.27*1.d-7) le 2.*FSRM[i])
contrast = MAX(Trtemp[aa])/MEAN(Trtemp[aa])-1.0
plot,l*1.d7,Trtemp,xst=1,xrange=[6173.27-2.*FSRM[i]*1.d7,6173.27+2.*FSRM[i]*1.d7],tit='!17 M'+STRTRIM(STRING(i+1),1),charsize=1.5,xtit='Wavelength (Angstrom)',ytit='Transmittance'
;IF(i EQ 0) THEN PRINT,angleout[0],angleout[1],'  Throughput=',MAX(Trtemp[aa]),'  contrast=',contrast ; print for NB
IF(i EQ 0) THEN BEGIN
    PRINT,angleout[0],MAX(Trtemp[aa]),contrast ; print for NB
    returnedt=MAX(Trtemp[aa])
    returnedc=contrast                    ;quantities returned by the code
ENDIF

Tr = Tr * Trtemp

ENDFOR

a=WHERE(ABS(l*1.d7 - 6173.2627) le .6)
maxi=MAX(Tr[a])
mini=MIN(Tr[a])
a=WHERE(Tr GE maxi/2.d0 AND ABS(l*1.d7 - 6173.2627) le .6)
fwhm=FLOAT(N_ELEMENTS(a))*(l[1]-l[0])*1.d7
centr=l[WHERE(ABS(l*1.d7 - 6173.2627) le 0.6 AND Tr EQ maxi)]*1.d7
;PRINT,'minimum=',mini

plot,l*1.d7,Tr,xst=1,xrange=[6172.8,6173.72],tit='!17A='+STRTRIM(STRING(maxi,format='(f7.3)'),1)+', !7dk!17='+STRTRIM(STRING(fwhm,format='(f8.4)'),1)+', !7k!17='+STRTRIM(STRING(centr,format='(f8.3)'),1),xtit='Wavelength (Angstrom)',ytit='Transmittance',charsize=1.5

device,/close
set_plot,'x'
!P.MULTI=0
;WINDOW,0,RETAIN=2

;READ,pause






END
