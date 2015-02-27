; FUNCTION THAT DRAWS THE NON-UNIFORMITY MAPS
;---------------------------------------------------------------------------

FUNCTION drawmaps
set_plot,'ps'
loadct,4
FOR element=2,6 DO BEGIN
device,file='yo.ps',bits=24,xsize=20,ysize=26,xoffset=0,yoffset=0.5,/color
!p.multi=[0,2,3]
IF(element GE 2) THEN title='E'+STRTRIM(STRING(element-1),1)+' v=-6000' ; Lyot
IF(element LT 2) THEN title='M'+STRTRIM(STRING(element+1),1)+' v=-6000' ; Michelson
FILENAME='erreurv-6000.0000'+STRTRIM(STRING(element),1)+'.bin'
RESTORE,filename
a=WHERE(erreur NE 0.0)
tvim,erreur,range=[MIN(erreur[a]),MAX(erreur[a])],/rct,/scale,tit='!17'+title,xtit='x',ytit='y',pcharsize=1.5,stit='!17Error ratio'
histo,erreur[a],MIN(erreur[a]),MAX(erreur[a]),0.001
xyouts,0.7,0.925,'Mean='+STRTRIM(STRING(MEAN(erreur[a])),1),charsize=1.5,/NORMAL
FILENAME='erreurv0.0000000'+STRTRIM(STRING(element),1)+'.bin'
RESTORE,filename
a=WHERE(erreur NE 0.0)
tvim,erreur,range=[MIN(erreur[a]),MAX(erreur[a])],/rct,/scale,tit='!17v=0',xtit='x',ytit='y',pcharsize=1.5,stit='!17Error ratio'
histo,erreur[a],MIN(erreur[a]),MAX(erreur[a]),0.001
xyouts,0.7,0.6,'Mean='+STRTRIM(STRING(MEAN(erreur[a])),1),charsize=1.5,/NORMAL
FILENAME='erreurv6000.0000'+STRTRIM(STRING(element),1)+'.bin'
RESTORE,filename
a=WHERE(erreur NE 0.0)
tvim,erreur,range=[MIN(erreur[a]),MAX(erreur[a])],/rct,/scale,tit='!17v=6000',xtit='x',ytit='y',pcharsize=1.5,stit='!17Error ratio'
histo,erreur[a],MIN(erreur[a]),MAX(erreur[a]),0.001
xyouts,0.7,0.275,'Mean='+STRTRIM(STRING(MEAN(erreur[a])),1),charsize=1.5,/NORMAL
device,/close
READ,pause
ENDFOR
set_plot,'x'
RETURN,'done'
END

; FUNCTION THAT COMPUTES THE LIMB DARKENING
; FROM H. NECKEL (2005)
;------------------------------------------------------------------------

FUNCTION limb,long

mu0 = long;COS(long*!dpi/180.d0)
mu  = [1.d0,mu0,mu0^2.d0,mu0^3.d0,mu0^4.d0,mu0^5.d0]
lamref=0.61733433d0

A = DBLARR(6)
A[0] = 0.75267d0-0.265577d0/lamref
A[1] = 0.93874d0+0.265577d0/lamref-0.004095d0/lamref^5.d0
A[2] = -1.89287d0+0.012582d0/lamref^5.d0
A[3] = 2.42234d0-0.017117d0/lamref^5.d0
A[4] = -1.71150d0+0.011977d0/lamref^5.d0
A[5] = 0.49062d0-0.003347d0/lamref^5.d0

ratio = TOTAL(A*mu); ratio intensity at the limb over intensity at the center


RETURN,ratio

END

;------------------------------------------------------------------------
;
; Program to set filter profiles and calculate expected noise performance.
; Author: Jesper Schou (based on sim.pro)
; Modified by S. Couvidat
;
;------------------------------------------------------------------------

PRO HMI,phasesdonnees,contrastdonnees,erreur0,erreur1,erreur2,long


; PARAMETERS
;------------------------------------------------------------------------

iline=1         ; Set line to look at (0=6768, 1=6173)

IF (iline EQ 0) THEN BEGIN ; Ni line, MDI case
  lam0  = 6767.68d0   ; Wavelength in Angstrom
  geff  = 1.426    ; Effective Lande factor
  fdepth= 1.0      ; line depth relative to Ni line
  fwidth= 1.0      ; line width relative to Ni line
  title = 'Ni6768'
ENDIF
IF (iline EQ 1) THEN BEGIN ; Fe line, HMI case
  lam0  = 6173.3433d0
  geff  = 2.500
  fdepth= 0.623/0.53;According to Yang's web page
  fwidth= 0.10/0.12
  title = 'Fe6173'
ENDIF

dvdb    = 1.35*(lam0/6767.68)*(geff/1.426) ; Conversion factor from field to velocity
well    = 1.50d5 ;1.25e5                        ; Set full well depth (or actually the maximum exposure level) in electrons
ntune   = 6                              ; Set number of tuning positions (NUMBER OF FILTERS)
ntunep  = 11                             ; All positions for the purpose of plotting.
inttune = 5./2                           ; Set number of tuning positions over wnarrow. (1/spacing)
                                         ; Don't mess with this unless you know what you are doing.


;parameters added by S. Couvidat
darkening  = 0     ; limb darkening ? 1=yes, 0=no
;long       = 75.0   ; angular distance from disk center
noangle    = 1     ; to test the impact of the incidence angle (0 or 1)
uniformity = 1     ; to test the impact of the spatial non-uniformity on the Lyot element aperture
                   ; =0 no test of spatial uniformity
                   ; =1 test 1 value
                   ; =2 test all the values for a specific map E*.bin or for all the maps
                   ; =3 give the zero-velocity error for a specific map E*.bin
IF(noangle EQ 0 AND uniformity NE 0) THEN STOP
phase      = [0.0d0,0.d0,0.d0,-4.7d0,0.08d0,-1.45d0,-5.0d0]*!dpi/180.d0 ; one per filter element
contrast   = [0.98,0.99,0.98,0.93,0.96,0.97,0.93]  
; only the first 3 elements of each array are used if the actual Lyot
; profile from HMI is used (see below)
phase0     = [0.0d0,0.d0,0.d0,-6.7d0,0.08d0,-1.45d0,-5.0d0]*!dpi/180.d0
contrast0  = [0.98,0.99,0.98,0.93,0.96,0.97,0.93]
recommence = 0
abscisse   = 0
ordonnee   = 0
nx         = 1



; DIFFERENT CASE
;---------------------------------------------------------------------------------------------------------------

IF (uniformity EQ 2) THEN BEGIN ; we produce maps of the velocity error
    element    = 3    ; element to which the phase and contrast apply  0=NB Michelson, 1=WB Michelson, 2=E1, 3=E2, 4=E3, 5=E4, 6=E5 7=ALL THE ELEMENTS
    IF(element NE 7) THEN BEGIN
        IF(element GE 2) THEN RESTORE,'E'+STRTRIM(STRING(element-2),1)+'.bin' ; Lyot
        IF(element LT 2) THEN RESTORE,'/scr109a/couvidat/HMI_RESULTS/M'+STRTRIM(STRING(element),1)+'.bin' ; Michelson
        nx = N_ELEMENTS(resultat[0,0,*])
    ENDIF ELSE BEGIN
        resultati   = FLTARR(6,2,1,1) ; for compatibility with the E*.bin files
        nx = 800
        FOR i=0,1 DO BEGIN
            RESTORE,'/scr109a/couvidat/HMI_RESULTS/M'+STRTRIM(STRING(element),1)+'.bin'
            resultati[i,0,*,*]=resultat[0,*,*]
            resultati[i,1,*,*]=resultat[1,*,*]
        ENDFOR
        FOR i=2,6 DO BEGIN
            RESTORE,'/scr109a/couvidat/HMI_RESULTS/E'+STRTRIM(STRING(element),1)+'.bin'
            resultati[i,0,*,*]=resultat[0,*,*]
            resultati[i,1,*,*]=resultat[1,*,*]
        ENDFOR
    ENDELSE
    loc        = 62;38             ;50 corresponds to zero input velocity
    recommence = 1
    erreur     = FLTARR(nx,nx) ; will contain the error ratio at a specific velocity
ENDIF

IF (uniformity EQ 3) THEN BEGIN ; we produce maps of the velocity error
    element    = 1    ; element to which the phase and contrast apply 0=NB Michelson, 1=WB Michelson, 2=E1, 3=E2, 4=E3, 5=E4, 6=E5
    IF(element GE 2) THEN RESTORE,'~/E'+STRTRIM(STRING(element-2),1)+'.bin' ; Lyot
    IF(element LT 2) THEN RESTORE,'~/M'+STRTRIM(STRING(element),1)+'.bin' ; Michelson
    nx         = N_ELEMENTS(resultat[0,0,*])
    loc        = 50
    recommence = 1
    erreur     = FLTARR(nx,nx) ; will contain the error ratio at a specific velocity
ENDIF

;IF (uniformity EQ 1) THEN BEGIN
;    phase[0]   =  phasesdonnees[0]*!pi/180.; relative phase
;    phase[1]   =  phasesdonnees[1]*!pi/180.; relative phase
;    phase[2]   =  phasesdonnees[2]*!pi/180.; relative phase
;   ;phase[3]   =  3.0*!pi/180.; relative phase
;   ;phase[4]   =  -4.58d0*!pi/180.; relative phase
;   ;phase[5]   =   2.07d0*!pi/180.; relative phase
;   ;phase[6]   =  15.22d0*!pi/180.; relative phase
;    contrast[0] = contrastdonnees[0]
;    contrast[1] = contrastdonnees[1]
;    contrast[2] = contrastdonnees[2]
;   ;contrast[3] = 0.94d0
;   ;contrast[4] = 0.94d0
;   ;contrast[5] = 0.95d0
;   ;contrast[6] = 0.93d0
;ENDIF


;------------------------------------------------------------------------------------------------
; MICHELSON PARAMETERS
;------------------------------------------------------------------------------------------------

FSRmich = 0.189D0                       ; Free Spectral Range of the narrowest Michelson for MDI at the wavelength 6768 A 
wnarrow = FSRmich*lam0/6767.68           ; Set period of narrowest Michelson
mlist  = [1,2,4]                        ; Set widths of tunable elements in units of wnarrow
                                        ; [1,2] is present MDI. [1,2,4] is with a third tunable element
                                        ; the last element is the Lyot tunable element
nmich  = n_elements(mlist)
;wmich  = mlist*wnarrow
wmich  = [0.172d0,0.344d0,0.693d0]; actual FSR (January 2007)

dtune  = wmich[0]/inttune                ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
tune   = DBLARR(nmich,ntune)

testtuning= -1                          ; Set to 1 to test the effect of off-tuning, to -1 to use Jesper original code
IF (noangle EQ 0) THEN testtuning=0
phi     = FLTARR(nmich,ntune)           ; angles of the three half-wave plates in degrees
IF(testtuning EQ 0) THEN BEGIN
    phi[*,0]= [18,9,4.5]
    phi[*,1]= [54,27,13.5]
    phi[*,2]= [90,45,22.5]
    phi[*,3]= [-18,-9,-4.5]
    phi[*,4]= [-54,-27,-13.5]
    phi[*,5]= [-90,-45,-22.5]
    FOR i = 0,nmich-1 DO tune[i,*] = (phi[i,*]*!dpi/180.)*2./!dpi*wmich[i]
ENDIF
IF(testtuning EQ 1) THEN  BEGIN
    phi[*,1]= [18.,9,4.5]
    phi[*,2]= [18.25,9.4,4.5]            ; 0.0269444 degress = 97 arcseconds
    phi[*,3]= [18,9,4.725]
    phi[*,0]= [18,9,4.5]
    phi[*,4]= [18.9,8.55,4.5]
    phi[*,5]= [18.9,8.55,4.725]
    FOR i = 0,nmich-1 DO tune[i,*] = (phi[i,*]*!dpi/180.)*2./!dpi*wmich[i]
ENDIF
IF(testtuning EQ -1) THEN FOR  i = 0,nmich-1 DO tune(i,*)=(-(ntune-1)/2.+dindgen(ntune))*dtune ;TUNING POSITIONS FOR JESPER SCHOU

tunep  = (-(ntunep-1)/2.+dindgen(ntunep))*dtune

;------------------------------------------------------------------------------------------------
; LYOT PARAMETERS
;------------------------------------------------------------------------------------------------

lyotw0 = 8*wnarrow                      ; Set period of narrowest Lyot element
IF (iline EQ 0) THEN BEGIN
    nlyot = 4
    lyotw = DBLARR(nlyot)
    lyotw = lyotw0*[1,1,2,2,4,8]  ; MDI
ENDIF
IF (iline EQ 1) THEN BEGIN
   ;lyotw = lyotw0*[1,2,4,8] ; HMI
    lyotw = [1.405d0,2.779d0,5.682d0,11.354d0] ; actual FSRs, January 2007
    nlyot = N_ELEMENTS(lyotw)
ENDIF

lyottune  = DBLARR(nlyot)
fixexp    = 1 ; Identical exposures ?

;------------------------------------------------------------------------------------------------
; LINE PARAMETERS
;------------------------------------------------------------------------------------------------

npar  = 4    ; Number of line parameters (3 for Jesper)
nlamp = 1223 ; Set wavelength points for line profile
lamp  = -0.3D0+DINDGEN(nlamp)*0.6D/(nlamp-1)

IF 0 THEN BEGIN
    restore,'/scr305/schou/iquv.sav' ; Line profile from Yang's file
    profile = iquv0000(*,0:nlamp-1,*)
    line    = REFORM(profile(0,*,0)) ; contain the line profile
ENDIF

IF 1 THEN BEGIN

    OPENR,1,'Ulrich_Fe_0.txt'
    roger  = FLTARR(2,98)
    READF,1,roger
    CLOSE,1
    nroger = 98
    rlam   = REFORM(roger(0,*))
    rp     = REFORM(roger(1,*))
    rp1    = rp/INTERPOL([rp(0),rp(nroger-1)],[rlam(0),rlam(nroger-1)],rlam)
    line   = INTERPOL(rp,rlam,lamp)
    line   = line/INTERPOL([line(0),line(nlamp-1)],[lamp(0),lamp(nlamp-1)],lamp)

    ; for the least-squares fit
    ; I assume the fit the line by the
    ; function I0*(1-d*exp(-lamp^2/w))
    ; with line=(1-d*exp(-lamp^2/w))
    dlinedd= -(1-line) ; just 1-line for Jesper
    dlinedd= dlinedd/ABS(MIN(dlinedd))
    dlinedi= line
    width  = 0.102d0/2.d0/SQRT(ALOG(2.d0))   ; w in A^2 (part added by seb)
    dlinedw= -(1-line)*lamp^2.d0*2.d0/width^3.d0

    lineseb = line
    dlinesebdd= -(1-line) 
    dlinesebdi= line
    dlinesebdw=-(1-line)*lamp^2.d0*2.d0/width^3.d0
    ;RESTORE,'interpolatedsolarline.bin' ; inter/extrapolated line with FUNCTION line in HMI_erreur.pro
    ;lineseb   = lineint
    ;dlinesebdd= -(1-lineseb)
    ;dlinesebdd= dlinesebdd/ABS(MIN(dlinesebdd))
    ;dlinesebdi= lineseb
    ;dlinesebdw= -(1-lineseb)*lamp^2.d0*2.d0/width^3.d0
ENDIF

IF 0 THEN BEGIN
    t      = -2.7726*(lamp/0.125/fwidth)^2 ; Line profile from Ted's approximate equation
    t      = -25. + (t+25.)*(t GT -25.)
    line   = (1-0.5615*exp(t)*fdepth) ;0.6
    dlinedd= exp(t)*fdepth
    dlinedi= line
ENDIF

; Set wavelength mesh for filter calculations
;--------------------------------------------
IF(uniformity EQ 0) THEN nlam = 2880*5*2 ELSE nlam = 6200; 720*5*2 THAT IS NOT ENOUGH !!! WE NEED TO GO FAR FROM THE TARGET WAVELENGTH TAKE nlam=105000 INSTEAD, AND dlam=1./3.5d3
dlam   = 1/1d3
lam    = (DINDGEN(nlam)-(nlam-1)/2.)*dlam

; Convert wavelength to velocity
;--------------------------------------------
dlamdv = lam0/3d8
dvdlam = 1/dlamdv
vel    = lam*dvdlam   ; DOPPLER shift in cm/s
dvel   = dlam*dvdlam
dvtune = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune  = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; Number of test filter profiles to test sensitivities
; ONLY THE BLOCKER FILTER PROFILE IS ALTERED AND WE TEST HOW IT
; AFFECTS THE TOTAL PROFILE
;-----------------------------------------------------
ntestx = 1       ; if ntestx=1, no test
iseed  = 1000l
tunex  = DBLARR(nmich,ntune,ntestx)
FOR  i = 0,ntune-1 DO FOR j = 0,nmich-1 DO tunex(j,i,*) = tune(j,i)

;------------------------------------------------------------------------------------------------
; BLOCKER FILTER PARAMETERS & PROFILE
; NB:NOT USED IF WE USE THE ACTUAL LYOT PROFILE (SEE BELOW)
;------------------------------------------------------------------------------------------------

;blocker = DBLARR(nlam)                 ; contains the blocker profile
;blocker(WHERE(ABS(lam) LE lyotw0)) = 1 ; Simple boxcar blocker to arbitrary null in Lyot
;bfac    = 1.6
;q       = readx('smooth.txt')
;blocker = INTERPOL(q(1,*)/100,bfac*(q(0,*)-6170),lam,/spline)

RESTORE,'frontwindow.bin' ; front windo
blocker  = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
q        = READFITS('blocker11.fits')
blockerseb= blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam)
blocker  = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) ; uncentered blocker (correction of blocker11.fits because of f8 instead of collimated light)
blockerseb  = blockerseb/MAX(blockerseb) ; to normalize
blocker  = blocker/MAX(blocker) ; to normalize

;------------------------------------------------------------------------------------------------
; FOR THE DOPPLER VELOCITY DETERMINATION
; WE GENERATE ntest DOPPLER SHIFTED LINE PROFILES
;------------------------------------------------------------------------------------------------

ntest         = 101 ; Set test velocities: we will use ntest different velocities, resulting in ntest Doppler shifted line profiles
dvtest        = 500.
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1)/2)
lines         = DBLARR(nlam,ntest)
dlinesdv      = DBLARR(nlam,ntest)
q1            = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
q2            = [MIN(lamp)-10,lamp,MAX(lamp)+10] ; Following two lines take care of continuum outside of where profile is given
FOR i = 0,ntest-1 DO lines(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dv            = 100D0           ; Get derivative with respect to v
FOR i = 0,ntest-1 DO dlinesdv(*,i) = (INTERPOL(q1,q2,lam+(vtest(i)+dv)*dlamdv)-INTERPOL(q1,q2,lam+(vtest(i)-dv)*dlamdv))/2/dv
dlinesdd = DBLARR(nlam,ntest)   ; Get derivative with respect to depth
q1            = [MAX(dlinedd),dlinedd,MAX(dlinedd)]
for i = 0,ntest-1 DO dlinesdd(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dlinesdi = DBLARR(nlam,ntest) ; Get derivative with respect to overall intensity
q1            = [MAX(dlinedi),dlinedi,MAX(dlinedi)]
for i = 0,ntest-1 DO dlinesdi(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dlinesdw = DBLARR(nlam,ntest) ; Get derivative with respect to linewidth (added by seb)
q1            = [0.0,dlinedw,0.0]
for i = 0,ntest-1 DO dlinesdw(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)


linesseb         = DBLARR(nlam,ntest)
dlinessebdv      = DBLARR(nlam,ntest)
q1               = [MAX(lineseb),lineseb,MAX(lineseb)] ; Calculate Doppler shifted line profiles
q2               = [MIN(lamp)-10,lamp,MAX(lamp)+10] ; Following two lines take care of continuum outside of where profile is given
FOR i = 0,ntest-1 DO linesseb(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dv               = 100D0           ; Get derivative with respect to v
FOR i = 0,ntest-1 DO dlinessebdv(*,i) = (INTERPOL(q1,q2,lam+(vtest(i)+dv)*dlamdv)-INTERPOL(q1,q2,lam+(vtest(i)-dv)*dlamdv))/2/dv
dlinessebdd      = DBLARR(nlam,ntest)   ; Get derivative with respect to depth
q1               = [MAX(dlinesebdd),dlinesebdd,MAX(dlinesebdd)]
for i = 0,ntest-1 DO dlinessebdd(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dlinessebdi      = DBLARR(nlam,ntest) ; Get derivative with respect to overall intensity
q1               = [MAX(dlinesebdi),dlinesebdi,MAX(dlinesebdi)]
for i = 0,ntest-1 DO dlinessebdi(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)
dlinessebdw = DBLARR(nlam,ntest) ; Get derivative with respect to linewidth (added by seb)
q1            = [0.0,dlinesebdw,0.0]
for i = 0,ntest-1 DO dlinessebdw(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)


WHILE (recommence GE 0) DO BEGIN
  
    IF(uniformity EQ 2 OR uniformity EQ 3) THEN BEGIN
        phase[element]    = resultat[1,abscisse,ordonnee]
        contrast[element] = resultat[0,abscisse,ordonnee]
        IF(phase[element] EQ 0.0 AND contrast[element] EQ 0.0 AND abscisse+ordonnee NE 0) THEN GOTO,unif3 
    ENDIF

;------------------------------------------------------------------------------------------------
; LYOT FILTER PROFILE
;------------------------------------------------------------------------------------------------

    clyot   = 2*!dpi/lyotw
    lyot    = blocker ; Set Lyot+Blocker filter profile. No Doubled elements.
    lyotseb = blockerseb

    IF(noangle EQ 0) THEN  BEGIN
        theta   = 0.d0          ; incidence angle in degrees
        TARGET  = 0.690d0       ; target FWHM of the E2 element
        T       = 30.0d0        ; Temperature in Celsius degrees

        T       = T-30.d0
        theta   = theta*!dpi/180.d0
        l       = (lam+lam0)*1.d-4
    
                                ;SELLMEIER'S EQUATIONS
        nO1          = SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0) ; Calcite
        ne1          = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)
        nO2          = DBLARR(nlyot,nlam)
        ne2          = nO2
        FOR i=0,1 DO BEGIN      ; ADP
            nO2[i,*] = SQRT(2.3041d0+0.0111d0/(l^2.0d0-0.0133d0)+15.1086d0*l^2.0d0/(l^2.0d0-400.d0))
            ne2[i,*] = SQRT(2.1643d0+0.0097d0/(l^2.0d0-0.0129d0)+5.8057d0*l^2.0d0/(l^2.0d0-400.d0))
        ENDFOR
        FOR i=2,nlyot-1 DO BEGIN ; KDP
            nO2[i,*] = SQRT(2.2576d0+0.0101d0/(l^2.0d0-0.0142d0)+1.7623d0*l^2.0d0/(l^2.0d0-57.8984d0))
            ne2[i,*] = SQRT(2.1295d0+0.0097d0/(l^2.0d0-0.0014d0)+0.7580d0*l^2.0d0/(l^2.0d0-127.0535d0))
        ENDFOR
        
                                ;THERMO-OPTICAL COEFFICIENTS AT 6173 A
        dno1dt=-12.975d-6 ;-10.54d-6            ; per Celsius degree The value has been corrected to obtain thermal compensation
        dne1dt= 0.d0           ;-0.56d-6             ; corrected value
        dneadt=-2.040863d-6
        dnoadt=-5.4321289d-5
        dnekdt=-3.0100346d-5
        dnokdt=-4.43d-5        ;-4.4419765d-5        ; corrected value
                                ;THERMAL EXPANSION COEFFICIENTS
        alp1=5.68d-6
        alpa=4.2d-6
        alpk=44.d-6
        
                                ;FINE-TUNING OF THE LYOT ELEMENT THICKNESS AT T=30 CELSIUS DEGREES AND 6173 A
        l            = lam0*1.d-4
        dc           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        da           = SQRT(2.1643d0+0.0097d0/(l^2.0d0-0.0129d0)+5.8057d0*l^2.0d0/(l^2.0d0-400.d0))-SQRT(2.3041d0+0.0111d0/(l^2.0d0-0.0133d0)+15.1086d0*l^2.0d0/(l^2.0d0-400.d0))
        dk           = SQRT(2.1295d0+0.0097d0/(l^2.0d0-0.0014d0)+0.7580d0*l^2.0d0/(l^2.0d0-127.0535d0))-SQRT(2.2576d0+0.0101d0/(l^2.0d0-0.0142d0)+1.7623d0*l^2.0d0/(l^2.0d0-57.8984d0))
        
        thickness1    = DBLARR(4) ; calcite block thicknesses
        thickness2    = DBLARR(4) ; ADP or KDP block thicknesses
        RCADP         = (dneadt-dnoadt+da*alpa)/(dne1dt-dno1dt+dc*alp1) ; compensation ratios Calcite/ADP 
        RCKDP         = (dnekdt-dnokdt+dk*alpk)/(dne1dt-dno1dt+dc*alp1) ; compensation ratios Calcite/KDP
        
                                ; thickness of the first calcite block (E2)
        l             = (lam0+TARGET/2.d0)*1.d-4 ; target FWHM=TARGET
        dc2           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        L0            = ABS(0.25d0/(dc/lam0-dc2/(lam0+TARGET/2.d0))) ; from the exact definition of the FWHM
        L0            = ROUND(L0/ABS(lam0/dc))*ABS(lam0/dc) ; from the condition that dc L0=k lam0 where k is an integer
        thickness1[0] = L0/(1.d0-da/(dc*RCADP))
        
                                ; thickness of the second calcite block (E3)
        l             = (lam0+TARGET)*1.d-4 ; target FWHM=TARGET*2
        dc2           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        L0            = ABS(0.25d0/(dc/lam0-dc2/(lam0+TARGET))) ; from the exact definition of the FWHM
        L0            = ROUND(L0/ABS(lam0/dc))*ABS(lam0/dc) ; from the condition that dc L0=k lam0 where k is an integer
        thickness1[1] = L0/(1.d0-da/(dc*RCADP))
        
                                ; thickness of the third calcite block (E3)
        l             = (lam0+TARGET*2.d0)*1.d-4
        dc2           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        L0            = ABS(0.25d0/(dc/lam0-dc2/(lam0+TARGET*2.d0))) ; from the exact definition of the FWHM
        L0            = ROUND(L0/ABS(lam0/dc))*ABS(lam0/dc) ; from the condition that dc L0=k lam0 where k is an integer
        thickness1[2] = L0/(1.d0-dk/(dc*RCKDP))
        
                                ; thickness of the last calcite block (E4)
        l             = (lam0+TARGET*4.d0)*1.d-4
        dc2           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        L0            = ABS(0.25d0/(dc/lam0-dc2/(lam0+TARGET*4.d0))) ; from the exact definition of the FWHM
        L0            = ROUND(L0/ABS(lam0/dc))*ABS(lam0/dc) ; from the condition that dc L0=k lam0 where k is an integer
        thickness1[3] = L0/(1.d0-dk/(dc*RCKDP))
        
                                ; TUNABLE LYOT ELEMENT
                                ; thickness of the calcite block for the tunable element (E1)
        l             = (lam0+TARGET/4.d0)*1.d-4 ; target FWHM=TARGET/2.0
        dc2           = SQRT(2.18438d0+0.0087309d0/(l^2.0d0-0.01018d0)-0.0024411d0*l^2.0d0)-SQRT(2.69705d0+0.0192064d0/(l^2.0d0-0.01820d0)-0.0151624d0*l^2.0d0)
        L0            = ABS(0.25d0/(dc/lam0-dc2/(lam0+TARGET/4.d0))) ; from the exact definition of the FWHM
        L0            = ROUND(L0/ABS(lam0/dc))*ABS(lam0/dc) ; from the condition that dc L0=k lam0 where k is an integer
        thicknessc    = L0/(1.d0-da/(dc*RCADP))
        
                                ; thicknesses of the ADP and KDP blocks
        thickness2[0] = thickness1[0]/RCADP
        thickness2[1] = thickness1[1]/RCADP
        thickness2[2] = thickness1[2]/RCKDP
        thickness2[3] = thickness1[3]/RCKDP
        thicknessa    = thicknessc   /RCADP
        PRINT,'Crystal thicknesses at T=30 Celsius degrees, in Angstrom'
        PRINT,thickness1
        PRINT,thickness2
        PRINT,'Crystal thickness for the tunable Lyot element',thicknessc,thicknessa
        
        
        IF(T GT 0) THEN BEGIN
                                ;THERMAL EXPANSION OF THE CRYSTALS
            thickness1      = thickness1*(1.d0+T*alp1)
            thickness2[0:1] = thickness2[0:1]*(1.d0+T*alpa)
            thickness2[2:3] = thickness2[2:3]*(1.d0+T*alpk)
            hicknessc       = thicknessc*(1.d0+T*alp1)
            thicknessa      = thicknessa*(1.d0+T*alpa)  
                                ;BIREFRINGENCE CHANGE DUE TO TEMPERATURE VARIATION
            nO1=nO1+T*dno1dt
            ne1=ne1+T*dne1dt
            nO2[0:1,*]=nO2[0:1,*]+T*dnoadt
            ne2[0:1,*]=ne2[0:1,*]+T*dneadt
            nO2[2:3,*]=nO2[2:3,*]+T*dnokdt
            ne2[2:3,*]=ne2[2:3,*]+T*dnekdt
        ENDIF
        
        retardance=DBLARR(nlam)
        FOR   i = 0,nlyot-1 DO BEGIN
            retardance[*]=2.d0*!dpi*thickness1[i]/(lam+lam0)*(ne1[*]-nO1[*])*(1.0d0-theta^2.d0/4.d0/nO1[*]*(1.d0/nO1[*]-1.d0/ne1[*]))-2.d0*!dpi*thickness2[i]/(lam+lam0)*(ne2[i,*]-nO2[i,*])*(1.0d0-theta^2.d0/4.d0/nO2[i,*]*(1.d0/nO2[i,*]-1.d0/ne2[i,*]))
            lyot[*]=lyot[*]*(1.d0+COS(retardance[*]))/2.d0
        ENDFOR


    ENDIF ELSE BEGIN            ; Jesper initial code
        FOR i = 0,nlyot-1 DO begin
            lyot = lyot*(1.0+contrast0[i+3]*COS(clyot[i]*(lam+lyottune(i))+phase0[i+3]))/2.0
            IF(uniformity NE 0) THEN lyotseb= lyotseb*(1.0+contrast[i+3]*COS(clyot[i]*(lam+lyottune(i))+phase[i+3]))/2.0
        ENDFOR

    ENDELSE

;--------------------------------------------------------------------
; ACTUAL (FROM HMI)
;--------------------------------------------------------------------

; LYOT PROFILE OF E2-E5 OBTAINED BY DICK SHINE WITH LAMP LIGHT
; THIS PART COMPLETELY CANCELS ANYTHING DONE BEFORE REGARDING THE LYOT
; PROFILE !!!
; INCLUDES A FRONT WINDOW SAMPLE BUT NO BLOCKER FILTER
    ;RESTORE,'Lyot2.bin'
    ;lyot = INTERPOL(data,lresp,lam)
    ;lyot = lyot/MAX(lyot)*blocker
    ;a=where(lyot lt 0.0)
    ;lyot[a]=0.0
    ;RESTORE,'Lyot.bin'
    ;lyotseb = INTERPOL(f,l,lam)
    ;lyotseb = lyotseb/MAX(lyotseb);*blocker
    ;read,pause

;------------------------------------------------------------------------------------------------
; MICHELSON FILTER PROFILE
;------------------------------------------------------------------------------------------------

    IF (noangle EQ 0) THEN BEGIN

        filters   = DBLARR(nlam,ntune) ; contain the filter profiles (Blocker+Lyot+Michelsons)
        FOR itune=0,ntune-1 DO filters[*,itune]=lyot
        delta=DBLARR(nmich,nlam) ; contains the path difference
        
                                ;TUNABLE LYOT ELEMENT
        delta[0,*]=thicknessc*(ne1[*]-nO1[*])*(1.0d0-theta^2.d0/4.d0/nO1[*]*(1.d0/nO1[*]-1.d0/ne1[*]))-thicknessa*(ne2[0,*]-nO2[0,*])*(1.0d0-theta^2.d0/4.d0/nO2[0,*]*(1.d0/nO2[0,*]-1.d0/ne2[0,*]))
        
        
                                ;MICHELSON ELEMENTS
        
        d1=[6.310d7,12.620d7]   ; estimated width of the solid leg
        d2=[4.160d7,8.330d7]    ; estimated width of the vacuum leg
        n1=DBLARR(nlam)+1.51682698726654 ; at 6173 A, obtained by the condition of optical coincidence: d1/n1=d2/n2
        n2=DBLARR(nlam)+1.0     ; vacuum
        dn1dl=-3.65402d-6 ; from the LightMAchinery website, in A^{-1}
        dn2dl=0.d0
        
                                ;l=TARGET/8.d0               ; target FWHM=TARGET/4.0
                                ;d1[0]=-1.d0/8.d0/( (n1[0]+dn1dl*l)/(l+lam0)+n2[0]^2.d0/n1[0]/lam0-n1[0]/lam0-n2[0]*(n2[0]+dn2dl*l)/n1[0]/(l+lam0) ) ; obtained from the definitionof the FWHM + condition of optical coincidence + assuming target wavelength=lam0
        d1[0]=lam0*17434.d0/(n1[0]-n2[0]^2.0/n1[0])/2.d0 ; k=17434 comes from: delta/lam0= k * 2 !pi
                                ;l=TARGET/16.d0
                                ;d1[1]=-1.d0/8.d0/(n1[0]+dn1dl*l)/(l+lam0)+n2[0]^2.d0/n1[0]/lam0-n1[0]/lam0-n2[0]*(n2[0]+dn2dl*l)/n1[0]/(l+lam0) )    
        d1[1]=lam0*34868.d0/(n1[0]-n2[0]^2.0/n1[0])/2.d0
        d2[0]=n2[0]*d1[0]/n1[0]
        d2[1]=n2[0]*d1[1]/n1[0]
        
        ;n2=DBLARR(nlam)+1.0002763d0 ; !!! WARNING !!! dry air at sea level and 600 nm
        n1=n1+dn1dl*lam
        n2=n2+dn2dl*lam
        PRINT,'thickness of the Michelson',d1,d2
        
        delta[1,*]=2.0d0*(n1*d1[0]-n2*d2[0])-(d1[0]/n1-d2[0]/n2)*sin(theta)^2.d0-(d1[0]/n1^3.0d0-d2[0]/n2^3.0d0)*sin(theta)^4.0d0/4.0d0-(d1[0]/n1^5.0d0-d2[0]/n2^5.0d0)*sin(theta)^6.0d0/8.0d0
        delta[2,*]=2.0d0*(n1*d1[1]-n2*d2[1])-(d1[1]/n1-d2[1]/n2)*sin(theta)^2.d0-(d1[1]/n1^3.0d0-d2[1]/n2^3.0d0)*sin(theta)^4.0d0/4.0d0-(d1[1]/n1^5.0d0-d2[1]/n2^5.0d0)*sin(theta)^6.0d0/8.0d0
        
        tune[*,0]=0.d0
        mdi=fltarr(2,nlam)
        plot,lam+lam0,lyot,xrange=[6172,6174]
        FOR itune = 0,ntune-1 DO BEGIN
            filters[*,itune]=filters[*,itune]*(1+COS(2.d0*!dpi*delta[0,*]/(lam0+lam[*]+tune[0,itune])))/2.d0     
            IF(itune EQ 0) THEN oplot,lam+lam0,(1+COS(2.d0*!dpi*delta[0,*]/(lam0+lam[*]+tune[0,itune])))/2.d0,linestyle=3
            FOR i  = 1,nmich-1 DO BEGIN filters[*,itune]=filters[*,itune]*(1.0d0+COS(2.0d0*!dpi*delta[i,*]/(lam0+lam[*]+tune[i,itune])))/2.d0
                IF(itune EQ 0) THEN mdi[i-1,*]=(1+COS(2.d0*!dpi*delta[i,*]/(lam0+lam[*]+tune[i,itune])))/2.d0 ;,linestyle=3+i
            ENDFOR
        ENDFOR
        
    ENDIF ELSE BEGIN            ; original code by Jesper
        

        IF(uniformity EQ 0) THEN BEGIN
            cmich     = 2*!dpi/wmich
            filters   = DBLARR(nlam,ntune) ; contain the filter profiles (Blocker+Lyot+Michelsons)
            FOR itune = 0,ntune-1 DO BEGIN
                filters(*,itune) = lyot
                FOR i  = 0,nmich-1 DO filters(*,itune) = filters(*,itune)*(1+contrast0[i]*COS(cmich(i)*(lam+tune(i,itune))+phase0[i]))/2
            ENDFOR

        ENDIF ELSE BEGIN
            cmich     = 2*!dpi/wmich
            filterseb = DBLARR(nlam,ntune)
            filters   = DBLARR(nlam,ntune)
            FOR itune = 0,ntune-1 DO BEGIN
                filters(*,itune) = lyot
                FOR i  = 0,nmich-1 DO filters(*,itune) = filters(*,itune)*(1+contrast0[i]*COS(cmich(i)*(lam+tune(i,itune))+phase0[i]))/2
            ENDFOR
            FOR itune = 0,ntune-1 DO BEGIN
                filterseb(*,itune) = lyotseb
                FOR i  = 0,nmich-1 DO BEGIN 
                    filterseb(*,itune) = filterseb(*,itune)*(1+contrast[i]*COS(cmich(i)*(lam+tune(i,itune))+phase[i]))/2 ; The Michelsons
                ENDFOR
            ENDFOR
        ENDELSE
        
    ENDELSE
 
; WE ADD LIMB DARKENING
IF darkening EQ 1 THEN filterseb=filterseb*limb(long)

   
;------------------------------------------------------------------------------------------------
; VARIABLE EXPOSURE
; THE GOAL OF VARIABLE EXPOSURE IS TO INCREASE/DECREASE THE EXPOSURE
; TIME AS A FUNCTION OF A TUNING POSITION, SO THAT THE MAXIMUM
; TRANSMISSION IS ABOUT THE SAME AT ALL POSITIONS
;------------------------------------------------------------------------------------------------

    IF (fixexp EQ 0) THEN BEGIN ; Adjust filters if variable exposures are allowed
        trans = DBLARR(ntune)   ; Find total filter transmissions
        FOR itune = 0,ntune-1 DO trans(itune) = TOTAL(filters(*,itune))
        trans = trans/max(trans) ; Keep max transmission at roughly 1
        FOR itune = 0,ntune-1 DO filters(*,itune)    = filters(*,itune)/trans(itune)
    ENDIF

;------------------------------------------------------------------------------------------------
; CALCULATION OF NOISE ON THE VELOCITY DETERMINATION
;------------------------------------------------------------------------------------------------

; Calculate sigma for a single velocity determination from linearized inverse without fitting for depth or intensity. Calculate intensities

    IF(abscisse EQ 0 AND ordonnee EQ 0) THEN BEGIN
        inten         = TRANSPOSE(filters)#lines    ; filters contain the perfect filter (contrast=1, phase=0)
        dintendv      = TRANSPOSE(filters)#dlinesdv ; And derivatives
        texp          = well/MAX(inten) ; "exposure" time: MAX(inten) is the continuum texp is the number of photons received by unit intensity
    ENDIF

    IF (uniformity EQ 2 OR uniformity EQ 3) THEN GOTO,unif2
; LINEARIZED INVERSE WITH NO FITTING: sigma
;------------------------------------------------------------------------------------------------
    
    sigma         = DBLARR(ntest)
    FOR itest     = 0,ntest-1 DO BEGIN
        i0        = texp*inten(*,itest) ; Calculate intensities in number of photons
        sens      = texp*dintendv(*,itest) ; Calculate derivative with respect to velocity
        sigma(itest)=1.0/SQRT(TOTAL(sens^2/i0)) ; Calculate least squares estimate of noise. (This is just a weighted mean.)
    ENDFOR
    
; LINEARIZED INVERSE WITH FITTING: sigma1
;------------------------------------------------------------------------------------------------
; we do a least-squares fitting and compute the Hessian matrix
    
    dintendd = TRANSPOSE(filters)#dlinesdd ; Calculate sigma(1) for a single velocity determination from linearized
    dintendi = TRANSPOSE(filters)#dlinesdi ; inverse also fitting for depth and intensity.
    dintendw = TRANSPOSE(filters)#dlinesdw
    sigma1   = DBLARR(ntest)
    a        = DBLARR(ntune,npar)
    FOR itest= 0,ntest-1 DO BEGIN
        i0     = texp*inten(*,itest) ; Calculate intensities
        si0    = 1/SQRT(i0)     ; Calculate inverse sigmas
        a(*,0) = texp*dintendv(*,itest)*si0 ; Calculate derivative with respect to velocity
        a(*,1) = texp*dintendi(*,itest)*si0
        a(*,2) = texp*dintendd(*,itest)*si0
        a[*,3] = texp*dintendw[*,itest]*si0
        cov    = INVERT(TRANSPOSE(a)#a)
        sigma1(itest) = SQRT(cov(0,0)) ; Calculate least squares estimate of noise. (This is just a weighted mean.)
    ENDFOR
    
    IF(uniformity EQ 1) THEN BEGIN
        inten2   = TRANSPOSE(filterseb)#linesseb    ; filterseb contains the actual filters (contrast NE 1, phase NE 0)
        dintendv2= TRANSPOSE(filterseb)#dlinessebdv
        dintendd2= TRANSPOSE(filterseb)#dlinessebdd ; Calculate sigma(1) for a single velocity determination from linearized
        dintendi2= TRANSPOSE(filterseb)#dlinessebdi ; inverse also fitting for depth and intensity.
        dintendw2= TRANSPOSE(filterseb)#dlinessebdw 
        sigma1b  = DBLARR(ntest)
        a        = DBLARR(ntune,npar)
        FOR itest= 0,ntest-1 DO BEGIN
            i0     = texp*inten2(*,itest) ; Calculate intensities
            si0    = 1/SQRT(i0) ; Calculate inverse sigmas
            a(*,0) = texp*dintendv2(*,itest)*si0 ; Calculate derivative with respect to velocity
            a(*,1) = texp*dintendi2(*,itest)*si0
            a(*,2) = texp*dintendd2(*,itest)*si0
            a[*,3] = texp*dintendw2[*,itest]*si0
            cov    = INVERT(TRANSPOSE(a)#a)
            sigma1b(itest) = SQRT(cov(0,0)) ; Calculate least squares estimate of noise. (This is just a weighted mean.)
        ENDFOR        

        erreur0 = sigma1b[37]/SQRT(2.0) ; erreur a -6.5 km/s
        erreur1 = sigma1b[50]/SQRT(2.0) ; erreur a  0 km/s
        erreur2 = sigma1b[63]/SQRT(2.0) ; erreur a  6.5 km/s


;GOTO,endo


        IF (uniformity EQ 1) THEN BEGIN
            set_plot,'ps'
            device,file='yo2.ps',xsize=21,ysize=16,xoffset=0.2,yoffset=7.0,/color,bits=24
            !p.multi=[0,2,2]
            loadct,0
            bozo=STRTRIM(STRING(phase[0]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[1]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[2]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[3]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[4]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[5]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase[6]*180./!dpi,format='(F6.2)'),1)       ;Lyot prof.'
            bozo2=STRTRIM(STRING(contrast[0],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[1],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[2],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[3],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[4],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[5],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast[6],format='(F6.3)'),1)        ;Lyot prof.'
            bozo3=STRTRIM(STRING(phase0[0]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[1]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[2]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[3]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[4]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[5]*180./!dpi,format='(F6.2)'),1)+';'+STRTRIM(STRING(phase0[6]*180./!dpi,format='(F6.2)'),1)       ;Lyot prof.'
            bozo4=STRTRIM(STRING(contrast0[0],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[1],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[2],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[3],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[4],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[5],format='(F6.3)'),1)+';'+STRTRIM(STRING(contrast0[6],format='(F6.3)'),1)        ;Lyot prof.'

            plot,vtest,sigma1b/sigma1,xrange=[-1.e4,1.e4],yrange=[0,2],charsize=1.,xst=1,yst=1,tit='!17Phase='+STRTRIM(bozo3,1),xtit='input velocity (m/s)',ytit='Error ratio'
            plot,vtest,sigma1b/SQRT(2.0),xrange=[-6500,6500],yrange=[0,40],charsize=1.,xst=1,yst=1,tit='!17'+STRTRIM(bozo,1),xtit='input velocity (m/s)',ytit='!7r!17 (m/s)'
            oplot,vtest,sigma1/SQRT(2.0),linestyle=2
            plot,lam,lines[*,50],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.,tit='!17Contrast='+STRTRIM(bozo4,1),xtit='!7k!17 (A)',ytit='transmission'
            loadct,33
            for i=0,ntune-1 do oplot,lam,filters[*,i],color=i*(254./ntune)+1
            loadct,0
            plot,lam,lines[*,50],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.,tit='!17'+STRTRIM(bozo2,1),xtit='!7k!17 (A)',ytit='transmission'
            loadct,33 
            for i=0,ntune-1 do oplot,lam,filterseb[*,i],color=i*(254./ntune)+1
            loadct,0
            device,/close
            !p.multi=0
            set_plot,'x'
        ENDIF
    ENDIF

    
; SIMULTANEOUS FIT OF LCP AND RCP: sigmav
;------------------------------------------------------------------------------------------------

    vt1      = DBLARR(ntest,ntest) ; Calculate sigma(v) for simultaneous fit of LCP and RCP assuming same
    vt2      = DBLARR(ntest,ntest) ; intensity and depth.
    sigmas   = DBLARR(npar+1,ntest,ntest)
    a        = DBLARR(2*ntune,npar+1)
    FOR it1  = 0,ntest-1 DO BEGIN
        FOR it2= 0,ntest-1 DO BEGIN
            vt1(it1,it2) = vtest(it1)
            vt2(it1,it2) = vtest(it2)
            i0           = texp*[inten(*,it1),inten(*,it2)] ;   Calculate intensities
            si0          = 1/SQRT(i0) ;   Calculate inverse sigmas
            a(*,0)       = texp*[dintendv(*,it1),0*dintendv(*,it2)]*si0 ;   Calculate derivative with respect to velocity
            a(*,1)       = texp*[0*dintendv(*,it1),dintendv(*,it2)]*si0
            a(*,2)       = texp*[dintendi(*,it1),dintendi(*,it2)]*si0
            a(*,3)       = texp*[dintendd(*,it1),dintendd(*,it2)]*si0
            a[*,4]       = texp*[dintendw(*,it1),dintendw(*,it2)]*si0
            cov          = invert(transpose(a)#a)
            FOR i = 0,npar DO sigmas(i,it1,it2) = SQRT(cov(i,i)) ;   Calculate least squares estimate of noise. (This is just a weighted mean.)
        ENDFOR
    ENDFOR

    vta     = (vt1+vt2)/2
    vts     = (vt2-vt1)/2
    vt1r    = REFORM(vt1,ntest^2)
    vt2r    = REFORM(vt2,ntest^2)
    vtar    = REFORM(vta,ntest^2)
    vtsr    = REFORM(vts,ntest^2)
    sigmasr = REFORM(sigmas,npar+1,ntest^2)
    sigmavr = SQRT(sigmasr(0,*)^2+sigmasr(1,*)^2)/2
    sigmav  = REFORM(sigmavr,ntest,ntest)

;------------------------------------------------------------------------
; MDI-LIKE VELOCITY MEASUREMENT TECHNIQUE WITH 6 SYMMETRICAL POINTS 
;------------------------------------------------------------------------
    unif2:
    IF(abscisse EQ 0 AND ordonnee EQ 0) THEN BEGIN
        x       = 2.0*!dpi*(-(ntune-1)/2.+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune
       ;c0      = REFORM(COS(0*x)#inten) ; coefficient a0 (for the cosine) of the Fourier series expansion
        c1      = REFORM(COS(x)#inten) ; coefficient a1
        s1      = REFORM(SIN(x)#inten) ; coefficient b1 (for the sine)
       ;c2      = REFORM(COS(2*x)#inten) ; coefficient a2
       ;s2      = REFORM(SIN(2*x)#inten) ; coefficient b2
        pv1     = dvtune*inttune*2.0 ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
       ;pv2     = dvtune*inttune ; half of the dynamic range
        phi1    = ATAN(-s1,-c1)
       ;phi2    = ATAN(-s2,-c2)
        vel1    = phi1*pv1/2.0/!dpi ; velocity measurement with an MDI-like algorithm
       ;vel2    = phi2*pv2/2/!dpi
        vel1a   = (vel1-vtest+10.5*pv1) MOD pv1-pv1/2.0+vtest ; the lookup tables for the filters we think are correct
       ;vel2a   = (vel2-vtest+10.5*pv2) MOD pv2-pv2/2+vtest
    ENDIF    
; thanks to the previous paragraph, we have the relationship between
; velocity measurements given by the MDI-like algorithm and actual
; velocity when the HMI filters are perfect (contrast=1, phase=0)
        
    IF(uniformity EQ 3) THEN BEGIN
        itest = loc
        ;vt    = vtest(itest)
        inten2= TRANSPOSE(filterseb)#linesseb
        ;i0    = inten2(*,itest)*texp
        ;nt    = 7500          ;10000
        ;intent= DBLARR(ntune,nt)
        ;IF(abscisse EQ 0 AND ordonnee EQ 0) THEN BEGIN
        ;    randomge = FLTARR(nt,ntune)
        ;    iseed = 1000l
        ;    FOR i = 0,ntune-1 DO randomge[*,i] = RANDOMN(iseed,nt)
        ;ENDIF
        ;FOR i = 0,ntune-1 DO intent(i,*)=i0(i)+SQRT(i0(i))*randomge[*,i]
        c1    = REFORM(COS(x)#inten2)
        s1    = REFORM(SIN(x)#inten2)
        phi1t = ATAN(-s1,-c1)
        vel1t = phi1t*pv1/2.0/!dpi
        ;vel1at= (vel1t-vt+10.5*pv1) MOD pv1-pv1/2.0+vt
        ;v1t   = INTERPOL(vtest,vel1a,vel1at) ; measured velocity
        ;sigmat1a = SQRT(REBIN((v1t-vt)^2.0,1)) ; error on the measured velocity
        erreur[abscisse,ordonnee]=vel1t[itest]-vel1[itest];sigmat1a/SQRT(2.)
    ENDIF

  ; IMPACT OF THE NON-UNIFORMITY OF THE LYOT ELEMENTS ON THE DOPPLER
  ; VELOCITY DETERMINATION
    IF(uniformity EQ 1 OR uniformity EQ 2) THEN BEGIN
        IF(abscisse EQ 0 AND ordonnee EQ 0) THEN BEGIN
            nt      = 7500      ;10000
            sigmat1 = FLTARR(ntest)
           ;sigmat2 = FLTARR(ntest)
            sigmat1a= FLTARR(ntest)
            systematicerror1a= FLTARR(ntest)
           ;sigmat2a= FLTARR(ntest)
            iseed   = 1000l
            randomge = FLTARR(nt,ntune)
            FOR i= 0,ntune-1 DO randomge[*,i] = RANDOMN(iseed,nt)
            FOR itest= 0,ntest-1 DO BEGIN
                vt   = vtest(itest)
                i0   = inten(*,itest)*texp
                intent=DBLARR(ntune,nt)
                FOR i= 0,ntune-1 DO intent(i,*)=i0(i)+SQRT(i0(i))*randomge[*,i] ; SQRT(i[0]) is photon noise
               ;FOR i= 0,ntune-1 DO intent(i,*)=RANDOMN(iseed,nt,POISSON=i0[i])
               ;c0   = REFORM(COS(0*x)#intent)
                c1   = REFORM(COS(x)#intent)
                s1   = REFORM(SIN(x)#intent)
               ;c2   = REFORM(COS(2*x)#intent)
               ;s2   = REFORM(SIN(2*x)#intent)
                phi1t = ATAN(-s1,-c1)
               ;phi2t = ATAN(-s2,-c2)
                vel1t = phi1t*pv1/2.0/!dpi
               ;vel2t = phi2t*pv2/2/!dpi
                vel1at= (vel1t-vt+10.5*pv1) MOD pv1-pv1/2.0+vt
               ;vel2at= (vel2t-vt+10.5*pv2) MOD pv2-pv2/2+vt
                sigmat1(itest) = SQRT(REBIN((vel1at-vel1a(itest))^2.0,1))
               ;sigmat2(itest) = SQRT(REBIN((vel2at-vel2a(itest))^2,1))
                v1t = INTERPOL(vtest,vel1a,vel1at,/QUADRATIC) ; measured velocity
               ;v2t = INTERPOL(vtest,vel2a,vel2at)
                sigmat1a(itest) = SQRT(REBIN((v1t-vtest(itest))^2.0,1)) ; error on the measured velocity for what we think are the good filters
                systematicerror1a[itest] = MEAN(v1t-vtest[itest])
               ;sigmat2a(itest) = SQRT(REBIN((v2t-vtest(itest))^2,1))
            ENDFOR 
        ENDIF
            
        inten2  = TRANSPOSE(filterseb)#lines
        c1      = REFORM(COS(x)#inten2) ; coefficient a1
        s1      = REFORM(SIN(x)#inten2) ; coefficient b1 (for the sine)
        phi1    = ATAN(-s1,-c1)
        vel1b   = phi1*pv1/2.0/!dpi  
        vel1bb  = (vel1b-vtest+10.5*pv1) MOD pv1-pv1/2.0+vtest ; look-up table for actual filters (contrast NE 1, phase NE 0)
       ;inten2  = TRANSPOSE(filterseb)#lines

        IF (uniformity EQ 1) THEN BEGIN
            sigmat1 = FLTARR(ntest)
            sigmat1b= FLTARR(ntest)
            systematicerror1b= FLTARR(ntest)
            FOR itest= 0,ntest-1 DO BEGIN
                vt   = vtest(itest)
                i0   = inten2(*,itest)*texp
                intent=DBLARR(ntune,nt)
                FOR i= 0,ntune-1 DO intent(i,*)=i0(i)+SQRT(i0(i))*randomge[*,i] ;!!!WARNING: PHOTON NOISE CANCELED HERE
               ;FOR i= 0,ntune-1 DO intent(i,*)=RANDOMN(iseed,nt,POISSON=i0[i])
                c1   = REFORM(COS(x)#intent)
                s1   = REFORM(SIN(x)#intent)
                phi1t = ATAN(-s1,-c1)
                vel1t = phi1t*pv1/2.0/!dpi
                vel1at= (vel1t-vt+10.5*pv1) MOD pv1-pv1/2.0+vt
                sigmat1(itest) = SQRT(REBIN((vel1at-vel1a(itest))^2.0,1)) ;vel1bb
                v1tb = INTERPOL(vtest,vel1a,vel1at,/QUADRATIC) ; velocity obtained with the wrong table (vel1a instead of vel1bb)
                sigmat1b(itest) = SQRT(REBIN((v1tb-vtest(itest))^2.0,1)) ; error on the measured velocity for actual filters
                systematicerror1b[itest] = MEAN(v1tb-vtest[itest])
            ENDFOR
        ENDIF ELSE BEGIN
            itest= loc
            vt   = vtest(itest)
            i0   = inten2(*,itest)*texp
            intent=DBLARR(ntune,nt)
            FOR i= 0,ntune-1 DO intent(i,*)=i0(i)+SQRT(i0(i))*randomge[*,i]
           ;FOR i= 0,ntune-1 DO intent(i,*)=RANDOMN(iseed,nt,POISSON=i0[i])
            c1   = REFORM(COS(x)#intent)
            s1   = REFORM(SIN(x)#intent)
            phi1t = ATAN(-s1,-c1)
            vel1t = phi1t*pv1/2.0/!dpi
            vel1at= (vel1t-vt+10.5*pv1) MOD pv1-pv1/2.0+vt
            v1tb = INTERPOL(vtest,vel1bb,vel1at) ; measured velocity
            sigmat1b = SQRT(REBIN((v1tb-vtest(itest))^2.0,1)) ; error on the measured velocity   
            erreur[abscisse,ordonnee]=sigmat1b/sigmat1a[itest]
        ENDELSE
    ENDIF
        
    IF( uniformity EQ 1) THEN BEGIN
        set_plot,'ps'
        device,file='yo.ps',xsize=18,ysize=26,xoffset=1,yoffset=0.5,/color,bits=24
        !p.multi=[0,2,3]
        loadct,0
        plot,vtest,vel1,xrange=[-1.e4,1.e4],charsize=1.5,xst=1,yst=1,tit='!17Phase='+STRTRIM(bozo3,1),xtit='input velocity (m/s)',ytit='output velocity (m/s)'
        oplot,vtest,vel1b,linestyle=2
        ;plot,vtest,vel1b,xrange=[-1.e4,1.e4],charsize=1.5,xst=1,yst=1,tit='!17Phase='+STRTRIM(bozo,1),xtit='input velocity (m/s)',ytit='output velocity (m/s)'
        plot,vtest,systematicerror1b,xrange=[-6500,6500],charsize=1.5,xst=1,yst=1,tit='!17Phase='+STRTRIM(bozo,1),xtit='input velocity (m/s)',ytit='Systematic error on velocity (m/s)'
        oplot,vtest,systematicerror1a,linestyle=2
        a=where(vtest ge -6500 and vtest le 6500)
        res3=poly_fit(vtest[a],systematicerror1b[a],4,yfit=y3)
        loadct,3
        oplot,vtest[a],y3,col=180
;        plot,vtest,sigmat1b/sigmat1a,xrange=[-1.e4,1.e4],yrange=[0,2],charsize=1.5,xst=1,yst=1,tit='!17',xtit='input velocity (m/s)',ytit='Error ratio'
        ;plot,vtest,vel1,xrange=[-500,500],charsize=1.5,xst=1,yst=1,tit='!17Relative phase=0',xtit='input velocity (m/s)',ytit='output velocity (m/s)'
        plot,vtest,sigmat1b/sigmat1a,xrange=[-6500,6500],charsize=1.5,xst=1,yst=1,tit='!17',xtit='input velocity (m/s)',ytit='Erreur Ratio'
        plot,vtest,sigmat1b/SQRT(2.0),xrange=[-6500,6500],yrange=[0,45],charsize=1.5,xst=1,yst=1,tit='!17',xtit='input velocity (m/s)',ytit='Error due to photon noise (m/s)'
        oplot,vtest,sigmat1a/SQRT(2.0),linestyle=2
        plot,lam,lines[*,50],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.5,tit='!17B='+STRTRIM(bozo4,1),xtit='!7k!17 (A)',ytit='transmission',yst=1,yrange=[0,1]
        oplot,lam,blocker,col=180
        loadct,33
        for i=0,ntune-1 do oplot,lam,filters[*,i],color=i*(254./ntune)+1
        loadct,0
        plot,lam,lines[*,50],xrange=[-0.6,0.6],thick=2,xst=1,charsize=1.5,tit='!17B='+STRTRIM(bozo2,1),xtit='!7k!17 (A)',ytit='transmission' ,yst=1,yrange=[0,1]
        oplot,lam,blocker,col=180
        loadct,33
        for i=0,ntune-1 do oplot,lam,filterseb[*,i],color=i*(254./ntune)+1
        loadct,0
        device,/close
    ENDIF
        

;read,pause

    !p.multi=0
    set_plot,'x'
    
    unif3:
    IF abscisse EQ nx-1 AND ordonnee EQ nx-1 THEN recommence=0
    IF ordonnee LT nx-1 THEN abscisse=abscisse ELSE BEGIN 
        abscisse=abscisse+1 
        PRINT,abscisse          ;,ordonnee,erreur[abscisse,ordonnee]
    ENDELSE
    IF ordonnee LT nx-1 THEN ordonnee=ordonnee+1 ELSE ordonnee=0
    IF recommence EQ 0 THEN recommence = -1
    
ENDWHILE
    
IF(uniformity EQ 2 OR uniformity EQ 3) THEN BEGIN
    IF(element GE 2) THEN title='E'+STRTRIM(STRING(element-1),1)+' v='+STRTRIM(STRING(vtest[loc]),1) ; Lyot
    IF(element LT 2) THEN title='M'+STRTRIM(STRING(element+1),1)+' v='+STRTRIM(STRING(vtest[loc]),1) ; Michelson
    set_plot,'ps'
    device,file='yo.ps',bits=24
    erreur[0,0]=0.0
    a=WHERE(erreur NE 0.0)
    IF(uniformity EQ 2) THEN tvim,erreur,range=[MIN(erreur[a]),MAX(erreur[a])],/rct,/scale,tit='!17'+title,xtit='x',ytit='y',pcharsize=1.5,stit='!17Error ratio'
    IF(uniformity EQ 3) THEN tvim,erreur,range=[MIN(erreur[a]),MAX(erreur[a])],/rct,/scale,tit='!17'+title,xtit='x',ytit='y',pcharsize=1.5,stit='Offset (m/s)';'!7r!17 (m/s)'
    device,/close
    set_plot,'x'
    IF(uniformity EQ 2) THEN SAVE,erreur,FILE='erreurv_new'+STRTRIM(STRING(vtest[loc]),1)+STRTRIM(STRING(element),1)+'.bin'
    IF(uniformity EQ 3) THEN SAVE,erreur,FILE='erreurzero2_'+STRTRIM(STRING(element),1)+'.bin'
ENDIF

endo:


END


