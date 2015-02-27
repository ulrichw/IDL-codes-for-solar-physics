;------------------------------------------------------
;
; PROGRAM TO COMPUTE SOME OF THE OBSERVABLES OF HMI:
; DOPPLER VELOCITY
; L.O.S. MAGNETIC-FIELD STRENGTH
; LINEWIDTH
; LINEDEPTH
; SOLAR CONTINUUM INTENSITY
;
; BASED ON AN ANSI C CODE
; AUTHOR: S. COUVIDAT BASED ON A CODE BY J. SCHOU
;
;------------------------------------------------------

; The following function is only used to test the main program Observables.pro
FUNCTION testing

lam0    = 6173.3433d0

ntune   = 6      ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)

; PARAMETERS
;-----------------------------------------------------------------------
; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
;phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
;contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
phase      = FLTARR(7)
contrast   = FLTARR(7)+1.0

wmich      = [0.172d0,0.344d0,0.693d0]          ; FSRs of the tunable elements
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements

dtune  = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

nlamp  = 30061 ; Set wavelength points for line profile
lamp   = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

; SOLAR LINE
RESTORE,'/scr20/richard/hmi/line_profile_library.sav'
rlam  = lambda-6173.3433d0
inclination = 0
field = 0
rp  = REFORM(rcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables

;OPENR,1,'Ulrich_Fe_0.txt'
;roger  = FLTARR(2,98)
;READF,1,roger
;CLOSE,1
;nroger = 98
;rlam   = REFORM(roger(0,*))
;rp     = REFORM(roger(1,*))
 
;rp = 1.d0-0.624*EXP(-rlam^2.d0/0.0675^2.d0)
   

; GRID WE WANT IN WAVELENGTH
nlam          = 38500.
dlam          = 1.d0/1.75d3
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = 2.059205672212074294d-5
dvdlam        = 1.d0/dlamdv
vel           = lam*dvdlam   ; DOPPLER shift in cm/s
dvel          = dlam*dvdlam
dvtune        = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune         = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; BLOCKER FILTER
RESTORE,'frontwindow.bin' ; front window
blocker       = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
q             = READFITS('blocker11.fits')
blocker       = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) ; blocker used for look-up tables

;blocker       = 0.57471215*exp(-lam^2.d0/4.5d0^2.d0) ;perfect blocker

; INPUT VELOCITIES
ntest         = 501 ; Set test velocities
dvtest        = dlam/dlamdv ; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)

; NON-TUNABLE PROFILE
lyot          = blocker
FOR i = 0,3 DO BEGIN
    lyot      = lyot   *(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0
ENDFOR

; TUNABLE PROFILE
cmich     = 2.d0*!dpi/wmich
filters   = DBLARR(nlam,ntune)
FOR itune = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
ENDFOR


line= INTERPOL(rp,rlam,lamp)
a   = WHERE(lamp ge MAX(rlam))
line[a] = 1.d0
a   = WHERE(lamp le MIN(rlam))
line[a] = 1.d0
PLOT,lamp,line,xrange=[-1,1],xst=1

lines   = DBLARR(nlam,ntest)
q2      = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1      = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
FOR i   = 0,ntest-1 DO lines[*,i]    = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
; CONVENTION: VELOCITIES > 0 FOR BLUESHIFTS

; INTENSITIES
inten      = TRANSPOSE(filters)#lines

; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd COEFFICIENT)
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune
c1         = REFORM(COS(x)#inten)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)#inten)  ; coefficient b1 (for the sine)
c2         = REFORM(COS(2.d0*x)#inten)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)#inten)  ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0   ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi; velocity measurement
vel1a      = ((vel1-vtest+10.5d0*pv1) MOD pv1)-pv1/2.d0+vtest ; look-up table for 1st Fourier coefficient
vel2a      = ((vel2-vel1a+10.5d0*pv2) MOD pv2)-pv2/2.d0+vel1a ; look-up table for 2nd Fourier coefficient

SAVE,vtest,vel1a,vel2a,FILE='lookup.bin'

inclination = 0
field   = 26
rp      = REFORM(rcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables
line    = INTERPOL(rp,rlam,lamp)
a       = WHERE(lamp ge MAX(rlam))
line[a] = 1.d0
a       = WHERE(lamp le MIN(rlam))
line[a] = 1.d0
PLOT,lamp,line,xrange=[-1,1],xst=1
lines   = DBLARR(nlam,ntest)
q2      = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1      = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
FOR i   = 0,ntest-1 DO lines[*,i]    = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
intenR  = TRANSPOSE(filters)#lines

rp      = REFORM(lcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables
line    = INTERPOL(rp,rlam,lamp)
a       = WHERE(lamp ge MAX(rlam))
line[a] = 1.d0
a       = WHERE(lamp le MIN(rlam))
line[a] = 1.d0
OPLOT,lamp,line,col=180
lines   = DBLARR(nlam,ntest)
q2      = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1      = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
FOR i   = 0,ntest-1 DO lines[*,i]    = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
intenL  = TRANSPOSE(filters)#lines

SAVE,intenR,intenL,file='intensity.bin'

RETURN,0

END



PRO observables

nx      = 1            ; number of rows in the filtergrams
ny      = 501            ; number of columns
lam0    = 6173.3433d0     ; Fe I central wavelength
N       = 6               ; number of tuning positions (NUMBER OF HMI FILTERS)
FSRNB   = 0.172d0         ; FSR of Narrow-Band Michelson, in Angstroms (MIGHT CHANGE WHEN ACTUAL VALUE IS MEASURED AFTER LAUNCH)

inttune = (N-1.)/2.       ; number of tuning intervals
dtune   = FSRNB/inttune   ; wavelength separation between each tuning position
dv      = 299792458.0/lam0; conversion factor from wavelength to velocity
dvtune  = dtune*dv
tune    = [-2.5,-1.5,-0.5,0.5,1.5,2.5] ; tuning positions in the range [-T/2,+T/2]
                                       ; where T is the period (we assume
                                       ; solar line is periodic because we
                                       ; compute the Fourier coefficients.
                                       ; we take filter positions that are
                                       ; symmetrical around 0
magnetic= 1.d0/(2.d0*4.67d-5*0.000061733433d0*2.5d0*299792458.d0) ; conversion factor from velocity to magnetic field
                                                                  ; Lande factor=2.5 for Fe I line
angle   = tune*2.d0*!dpi/double(N)
cosi    = COS(angle)
sini    = SIN(angle)
cosi2   = COS(2.d0*angle)
sini2   = SIN(2.d0*angle)
pv1     = dvtune*inttune*2.d0 ; conversion factor from phase of the 1st Fourier coefficients to velocity
pv2     = pv1/2.d0            ; conversion factor from phase of the 2nd Fourier coefficients to velocity
period  = inttune*2.d0*dtune


; READ FILTERGRAMS
; the exact sequence that will be used by HMI has not been decided yet
; here I assume that we will have 6 LCP filtergrams and 6 RCP filtergrams
;-----------------------------------------------------------------------

r       = testing()
RESTORE,'intensity.bin'
intensityL0=intenL[0,*]
intensityL1=intenL[1,*]
intensityL2=intenL[2,*]
intensityL3=intenL[3,*]
intensityL4=intenL[4,*]
intensityL5=intenL[5,*]
intensityR0=intenR[0,*]
intensityR1=intenR[1,*]
intensityR2=intenR[2,*]
intensityR3=intenR[3,*]
intensityR4=intenR[4,*]
intensityR5=intenR[5,*]


; LOOK-UP TABLES:
;----------------------------------------------------------------
; here we read the 2 look-up tables: 1 for the 1st Fourier
; coefficient, 1 for the 2nd Fourier coefficient.
; By look-up "table" we mean the Doppler velocity obtained by the MDI
; algorithm as a function of the actual Doppler velocity
; HERE, FOR SIMPLICITY, WE USE ONLY 1 PAIR OF TABLES FOR THE ENTIRE 4096*4096 FILTERGRAMS
; BUT ONCE SDO IS LAUNCHED WE WILL HAVE, MAYBE, 1 PAIR OF TABLES FOR
; EACH HMI PIXEL

RESTORE,'lookup.bin'  ; idl binary file contains the lookup tables
                      ; vtest is the actual (input) velocity,
                      ; vel1a is the output velocity
                      ; obtained from the phase of the 1st Fourier coefficient
                      ; vel2a is the output velocity
                      ; obtained from the phase of the 2nd Fourier coefficient
                      ; AFTER LAUNCH THE FORMAT OF THE
                      ; TABLES WILL PROBABLY BE .FITS

; DEFINITION OF THE ARRAYS CONTAINING THE RESULT:
lam0g = INTARR(nx,ny) ; will contain the Doppler velocity
B0g   = INTARR(nx,ny) ; will contain the l.o.s. magnetic field strength
widthg= INTARR(nx,ny) ; will contain the linewidth
Idg   = INTARR(nx,ny) ; will contain the linedepth
I0g   = INTARR(nx,ny) ; will contain the continuum intensity

; LOOP OVER THE 4096*4096 PIXELS
FOR iii=0,nx-1 DO BEGIN
    FOR jjj=0,ny-1 DO BEGIN

        ; CALCULATION OF THE 1st FOURIER COEFFICIENTS
        f1LCPc  = cosi[0]*intensityL0[iii,jjj]+cosi[1]*intensityL1[iii,jjj]+cosi[2]*intensityL2[iii,jjj]+cosi[3]*intensityL3[iii,jjj]+cosi[4]*intensityL4[iii,jjj]+cosi[5]*intensityL5[iii,jjj] ; 1st Fourier coefficient for cosine and LCP filtergrams 
        f1RCPc  = cosi[0]*intensityR0[iii,jjj]+cosi[1]*intensityR1[iii,jjj]+cosi[2]*intensityR2[iii,jjj]+cosi[3]*intensityR3[iii,jjj]+cosi[4]*intensityR4[iii,jjj]+cosi[5]*intensityR5[iii,jjj] ; 1st Fourier coefficient for cosine and RCP filtergrams 
        
        f1LCPs  = sini[0]*intensityL0[iii,jjj]+sini[1]*intensityL1[iii,jjj]+sini[2]*intensityL2[iii,jjj]+sini[3]*intensityL3[iii,jjj]+sini[4]*intensityL4[iii,jjj]+sini[5]*intensityL5[iii,jjj] ; 1st Fourier coefficient for sine and LCP filtergrams 
        f1RCPs  = sini[0]*intensityR0[iii,jjj]+sini[1]*intensityR1[iii,jjj]+sini[2]*intensityR2[iii,jjj]+sini[3]*intensityR3[iii,jjj]+sini[4]*intensityR4[iii,jjj]+sini[5]*intensityR5[iii,jjj] ; 1st Fourier coefficient for sine and RCP filtergrams 
        
        vLCP    = ATAN(-f1LCPs,-f1LCPc)*pv1/2.d0/!dpi   ; -f1LCPs and -f1LCPc because the first Fourier coefficients are -2/T*...
        vRCP    = ATAN(-f1RCPs,-f1RCPc)*pv1/2.d0/!dpi

        ; CALCULATION OF THE 2nd FOURIER COEFFICIENTS
        f2LCPc  = cosi2[0]*intensityL0[iii,jjj]+cosi2[1]*intensityL1[iii,jjj]+cosi2[2]*intensityL2[iii,jjj]+cosi2[3]*intensityL3[iii,jjj]+cosi2[4]*intensityL4[iii,jjj]+cosi2[5]*intensityL5[iii,jjj] ; 2nd Fourier coefficient for cosine and LCP filtergrams 
        f2RCPc  = cosi2[0]*intensityR0[iii,jjj]+cosi2[1]*intensityR1[iii,jjj]+cosi2[2]*intensityR2[iii,jjj]+cosi2[3]*intensityR3[iii,jjj]+cosi2[4]*intensityR4[iii,jjj]+cosi2[5]*intensityR5[iii,jjj] ; 2nd Fourier coefficient for cosine and RCP filtergrams 
        
        f2LCPs  = sini2[0]*intensityL0[iii,jjj]+sini2[1]*intensityL1[iii,jjj]+sini2[2]*intensityL2[iii,jjj]+sini2[3]*intensityL3[iii,jjj]+sini2[4]*intensityL4[iii,jjj]+sini2[5]*intensityL5[iii,jjj] ; 2nd Fourier coefficient for sine and LCP filtergrams 
        f2RCPs  = sini2[0]*intensityR0[iii,jjj]+sini2[1]*intensityR1[iii,jjj]+sini2[2]*intensityR2[iii,jjj]+sini2[3]*intensityR3[iii,jjj]+sini2[4]*intensityR4[iii,jjj]+sini2[5]*intensityR5[iii,jjj] ; 2nd Fourier coefficient for sine and RCP filtergrams 
        
        v2LCP   = ATAN(-f2LCPs,-f2LCPc)*pv2/2.d0/!dpi   ; -f2LCPs and -f2LCPc because the second Fourier coefficients are -2/T*...
        v2RCP   = ATAN(-f2RCPs,-f2RCPc)*pv2/2.d0/!dpi 
        v2LCP   = ((v2LCP-vLCP+10.5d0*pv2) MOD pv2)-pv2/2.d0+vLCP ; we use the uncorrected velocity, i.e. phase, of the 1st Fourier coefficient to correct for the estimate of v2LCP and v2RCP, because the range of velocities obtained with the second Fourier coefficient is half the range of the first Fourier coefficient
        v2RCP   = ((v2RCP-vRCP+10.5d0*pv2) MOD pv2)-pv2/2.d0+vRCP


        ; We interpolate in the look-up tables  to retrieve the actual velocities
	vLCP   = INTERPOL(vtest,vel1a,vLCP)
	vRCP   = INTERPOL(vtest,vel1a,vRCP)
	v2LCP  = INTERPOL(vtest,vel2a,v2LCP)
	v2RCP  = INTERPOL(vtest,vel2a,v2RCP)

        ;We compute the Doppler velocity (result is rounded because it must be an integer to save space)
        temp0           = (vLCP+vRCP+v2LCP+v2RCP)/4.d0
        lam0g[iii,jjj]  = ROUND(temp0) ; simple average between 1st and 2nd Fourier coeff. Need weights?

	;We compute the l.o.s. magnetic field
        B0g[iii,jjj]    = ROUND((vRCP-vLCP+v2RCP-v2LCP)/2.d0*magnetic)



       ;THE FOLLOWING IS AN ATTEMPT AT QUICKLY ESTIMATING THE LINEWIDTH, LINEDEPTH, AND CONTINUUM INTENSITY
       ;THIS IS MOST LIKELY NOT HOW WE WILL ESTIMATE THESE PARAMETERS ONCE SDO IS LAUNCHED

        f1LCPc = -f1LCPc*dtune/period*2.d0
        f1RCPc = -f1RCPc*dtune/period*2.d0
        f1LCPs = -f1LCPs*dtune/period*2.d0
        f1RCPs = -f1RCPs*dtune/period*2.d0
        f2LCPc = -f2LCPc*dtune/period*2.d0
        f2RCPc = -f2RCPc*dtune/period*2.d0
        f2LCPs = -f2LCPs*dtune/period*2.d0
        f2RCPs = -f2RCPs*dtune/period*2.d0
              

	;We compute the linewidth (in Angstroms) and multiply it by 10000 because the type is SHORT INT. 
        ; Based on analytical expression of the linewidth for a Gaussian profile
        ; here linewidth is the \sigma value of the Gaussian profile: I=I0g-Idg*exp(-(lam-lam0g)^2/sigma^2)
        temp            = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(f1LCPc^2.d0+f1LCPs^2.d0)/SQRT(f2LCPc^2.d0+f2LCPs^2.d0)))
        tempbis         = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(f1RCPc^2.d0+f1RCPs^2.d0)/SQRT(f2RCPc^2.d0+f2RCPs^2.d0)))
	widthg[iii,jjj] = ROUND(10000.d0*(temp+tempbis)/2.d0)
	      
        ;We compute the linedepth. Based on analytical expression of the linedepth for a Gaussian profile
        temp2           = period/2.d0*SQRT(f1LCPc^2.d0+f1LCPs^2.d0)/SQRT(!dpi)/temp   *EXP(!dpi^2.d0*temp^2.d0   /period^2.d0)
        temp2bis        = period/2.d0*SQRT(f1RCPc^2.d0+f1RCPs^2.d0)/SQRT(!dpi)/tempbis*EXP(!dpi^2.d0*tempbis^2.d0/period^2.d0)
        Idg[iii,jjj]    = ROUND((temp2+temp2bis)/2.d0)

        ;We compute the continuum intensity. Based on analytical expression of the continuum for a Gaussian profile
        temp3           = (vLCP+v2LCP)/2.d0/dv
        temp3bis        = (vRCP+v2RCP)/2.d0/dv
        meanL           = (intensityL0[iii,jjj]+intensityL1[iii,jjj]+intensityL2[iii,jjj]+intensityL3[iii,jjj]+intensityL4[iii,jjj]+intensityL5[iii,jjj])/6.d0
        meanR           = (intensityR0[iii,jjj]+intensityR1[iii,jjj]+intensityR2[iii,jjj]+intensityR3[iii,jjj]+intensityR4[iii,jjj]+intensityR5[iii,jjj])/6.d0
        meanL           = meanL+temp2/6.d0*(EXP(-(tune[0]*dtune-temp3)^2.d0/temp^2.d0)+EXP(-(tune[1]*dtune-temp3)^2.d0/temp^2.d0)+EXP(-(tune[2]*dtune-temp3)^2.d0/temp^2.d0)+EXP(-(tune[3]*dtune-temp3)^2.d0/temp^2.d0)+EXP(-(tune[4]*dtune-temp3)^2.d0/temp^2.d0)+EXP(-(tune[5]*dtune-temp3)^2.d0/temp^2.d0))
        meanR           = meanR+temp2bis/6.d0*(EXP(-(tune[0]*dtune-temp3bis)^2.d0/tempbis^2.d0)+EXP(-(tune[1]*dtune-temp3bis)^2.d0/tempbis^2.d0)+EXP(-(tune[2]*dtune-temp3bis)^2.d0/tempbis^2.d0)+EXP(-(tune[3]*dtune-temp3bis)^2.d0/tempbis^2.d0)+EXP(-(tune[4]*dtune-temp3bis)^2.d0/tempbis^2.d0)+EXP(-(tune[5]*dtune-temp3bis)^2.d0/tempbis^2.d0))
        I0g[iii,jjj]    = ROUND((meanL+meanR)/2.d0)


    ENDFOR
ENDFOR
plot,vtest,lam0g,xrange=[-6500,6500],xst=1

; we save the results
SAVE,lam0g,B0g,widthg,Idg,I0g,FILE='result.bin'

END
