; a lighter version of HMI.pro to compute errors on the Doppler
; velocity (systematic errors + photon noise errors)
; FOR LINE PROFILES OF RICHARD WACHTER (different fields
; and inclinations)
; TO TEST IMPACT OF MAGNETIC FIELD ON DOPPLER VELOCITY ERROR
; computes error for:
; ---MDI-like algorithm with 6 positions, using the first 2 Fourier
; coefficients, and the RCP+LCP filtergrams
; NB: ALGORITHM HAS BEEN MODIFIED: WE NORMALIZE BY FILTERS' INTEGRAL VALUE
; TO IMPROVE THE LINEWIDTH, LINEDEPTH, AND CONTINUUM INTENSITY MEASUREMENTS
;
;-----------------------------------------------------------------


PRO gaussianline,X,A,F,pder

A[0]   = ABS(A[0])
A[1]   = ABS(A[1])
A[3]   = ABS(A[3])

F1     = A[0]*(1.d0-A[1]*EXP(-(X+A[2])^2.d0/A[3]^2.d0))

F      = F1
Fd     = EXP(-(X+A[2])^2.d0/A[3]^2.d0)
Fdd    = ((X+A[2])*EXP(-(X+A[2])^2.d0/A[3]^2.d0))
Fddd   = ((X+A[2])^2.d0*EXP(-(X+A[2])^2.d0/A[3]^2.d0))


pder   = [ [F/A[0]],[-A[0]*Fd],[A[0]*A[1]*2.d0/A[3]^2.d0*Fdd],[-A[0]*A[1]*2.d0/A[3]^3.d0*Fddd] ]

END

;-----------------------------------------------------------------
;
; MAIN PROGRAM
;
;-----------------------------------------------------------------


PRO HMIrichard

SET_PLOT,'x'
!P.MULTI=0
WINDOW,0,RETAIN=2

lam0    = 6173.34d0  ; not 6173.3433 otherwise Richard's lines are not centered at [0,0] !!!
ntune   = 6          ; Set number of tuning positions (NUMBER OF FILTERS)
ntunep  = 11         ; All positions for the purpose of plotting.
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)
dvdb    = 1./2.13d12*lam0*2.5d0*299792458.d0  ; conversion from magnetic field to Doppler velocity


; PARAMETERS
;-----------------------------------------------------------------------
; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
;phase      = [25.0d0,-20.d0,20.d0,8.0700685d0,0.4029859d0,4.2042572d0,8.1685274d0]*!dpi/180.d0 
contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
;phase     = FLTARR(7)
;contrast  = FLTARR(7)+1.0

 wmich     = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]          ; FSRs of the tunable elements
 lyotw     = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements
;wmich     = [0.172d0,0.344d0,0.688d0] ; ideal FSRs
;lyotw     = [wmich[2]*2.d0,wmich[2]*4.d0,wmich[2]*8.d0,wmich[2]*16.d0]

dtune      = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
period     = inttune*2.d0*dtune
tune       = DBLARR(3,ntune)
FOR  i     = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune


nlamp      = 30061 ; Set wavelength points for line profile
lamp       = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

; GRID WE WANT IN WAVELENGTH
nlam          = 6500.;38500./7.
dlam          = 0.004;1.d0/1.75d3*7.5
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = lam0/299792458.d0
dvdlam        = 1.d0/dlamdv
vel           = lam*dvdlam   ; DOPPLER shift in cm/s
dvel          = dlam*dvdlam
dvtune        = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune         = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; BLOCKER FILTER
RESTORE,'frontwindow.bin' ; front window
blocker       = INTERPOL(transmission/100.d0,wavelength*10.d0-lam0,lam)
q             = READFITS('blocker11.fits')
blocker       = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lam0,lam) ; blocker used for look-up tables

;blocker       = 0.57471215*exp(-lam^2.d0/4.5d0^2.d0) ;perfect blocker

; INPUT VELOCITIES
ntest         = 901;2801;501 ; Set test velocities
dvtest        = 21.;7.5;100.d0 ; in m/s
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

; NORM OF FILTERS
norm = DBLARR(ntune)
FOR i=0,ntune-1 DO norm[i]=INT_TABULATED(lam,filters[*,i],/DOUBLE)
;FOR i=0,ntune-1 DO filters[*,i]=filters[*,i]/norm[i]

; MAGNETIC FIELD STRENGTH AND INCLINATION USED TO PRODUCE LOOK-UP TABLES
field      = 0
inclination= 0
PRINT,inclination,field

ROGER = 0 
RESTORE,'/scr20/richard/hmi/line_profile_library.sav'
rlam       = lambda-lam0

IF ROGER NE 1 THEN BEGIN       
    rp  = REFORM(rcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables
    line= INTERPOL(rp,rlam,lam)
    a   = WHERE(lam ge MAX(rlam))
    line[a] = 1.d0
    a   = WHERE(lam le MIN(rlam))
    line[a] = 1.d0
ENDIF ELSE BEGIN
    OPENR,1,'Ulrich_Fe_0.txt'
;   OPENR,1,'Ulrich_Fe_45.txt'
;   OPENR,1,'Ulrich_Fe_60.txt'
    roger  = FLTARR(2,98)
    READF,1,roger
    CLOSE,1
    nroger = 98
    rlam   = REFORM(roger[0,*])
    rp     = REFORM(roger[1,*])
    a      = WHERE(rlam GE 0.4)
    rp[a]  = 1.d0
    a      = WHERE(rlam LE -0.4)
    rp[a]  = 1.d0
    line   = INTERPOL(rp,rlam,lam)
    ;aaa    = [1.d0,0.6,0.0,0.062]
    ;weights= FLTARR(nlam)+1.0
    ;res    = CURVEFIT(lam,line,weights,aaa,FUNCTION_NAME='gaussianline')
ENDELSE

PLOT,lam,line,xrange=[-1,1],xst=1
lines   = DBLARR(nlam,ntest)
q2      = [MIN(lam)-10.d0,lam,MAX(lam)+10.d0]
q1      = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles
FOR i   = 0,ntest-1 DO lines[*,i]    = INTERPOL(q1,q2,lam+(vtest[i])*dlamdv)
; CONVENTION: VELOCITIES > 0 FOR BLUESHIFTS

; INTENSITIES
PRINT,'COMPUTATION OF INTENSITIES'
;inten = FLTARR(ntune,ntest)
;FOR i=0,ntune-1 DO FOR j=0,ntest-1 DO inten[i,j] = INT_TABULATED(lam,filters[*,i]*lines[*,j])
inten = (TRANSPOSE(filters)  #lines)*dlam
PRINT,'DONE'
rlam       = lambda-lam0

; WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd COEFFICIENT)
FOR i=0,ntune-1 DO inten[i,*]=inten[i,*]/norm[i]
x          = 2.d0*!dpi*(-(ntune-1.d0)/2.d0+DINDGEN(ntune))/ntune
c1         = REFORM(COS(x)#inten) ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)#inten) ; coefficient b1 (for the sine)
c2         = REFORM(COS(2.d0*x)#inten) ; coefficient a2
s2         = REFORM(SIN(2.d0*x)#inten) ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0 ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi ; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi ; velocity measurement
vel1a      = ((vel1-vtest+10.5d0*pv1) MOD pv1)-pv1/2.d0+vtest ; look-up table for 1st Fourier coefficient
vel2a      = ((vel2-vel1 +10.5d0*pv2) MOD pv2)-pv2/2.d0+vel1 ; look-up table for 2nd Fourier coefficient
read,pause
        
        
;GOTO,endo
velocity   = FLTARR(ntest)
velocity2  = FLTARR(ntest)
magnetic   = FLTARR(ntest)
magnetic2  = FLTARR(ntest)
widthL     = FLTARR(ntest)
widthR     = FLTARR(ntest)
depthL     = widthL
depthR     = widthR
continuumL = widthL
continuumR = widthR
        
richard2   = 0 ; MAGNETIC FIELD STRENGTH
richard1   = 0  ; MAGNETIC FIELD INCLINATION
contin     = 1.d0
        
FOR i=0,ntest-1 DO BEGIN

    ; LCP
    rp         = REFORM(lcp0[*,richard1,richard2])*contin
    line2      = INTERPOL(rp,rlam,lam+vtest[i]*dlamdv)
    a          = WHERE(lam ge MAX(rlam))
    line2[a]   = contin
    a          = WHERE(lam le MIN(rlam))
    line2[a]   = contin
    inten2     = (filters##line2)*dlam
    FOR ii=0,ntune-1 DO inten2[0,ii]=inten2[0,ii]/norm[ii]
    c1         = REFORM(COS(x)##inten2)      ; coefficient a1 (for the cosine)
    s1         = REFORM(SIN(x)##inten2)      ; coefficient b1 (for the sine)
    c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
    s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
    phi1       = ATAN(-s1,-c1)
    phi2       = ATAN(-s2,-c2)
    vel1x      = phi1*pv1/2.d0/!dpi          ; velocity measurement
    vel2x      = phi2*pv2/2.d0/!dpi          ; velocity measurement
    res1       = INTERPOL(vtest,vel1a,vel1x)
    vel2x      = (vel2x-vel1x+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1x
    res3       = INTERPOL(vtest,vel2a,vel2x)
        
    c1         = -c1*dtune/period*2.d0
    s1         = -s1*dtune/period*2.d0
    c2         = -c2*dtune/period*2.d0
    s2         = -s2*dtune/period*2.d0
    widthL[i]  = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(c1^2.d0+s1^2.d0)/SQRT(c2^2.d0+s2^2.d0)))
    depthL[i]  = period/2.d0/SQRT(!dpi)*SQRT(c1^2.d0+s1^2.d0)/widthL[i]*EXP(!dpi^2.d0*widthL[i]^2.d0/period^2.d0)
    velo       = (res1+res3)/2.d0
    continuumL[i]= TOTAL(inten2)/6.d0+depthL[i]/6.d0*(EXP(-(tune[0,0]-velo*dlamdv)^2.d0/widthL[i]^2.d0)+EXP(-(tune[0,1]-velo*dlamdv)^2.d0/widthL[i]^2.d0)+EXP(-(tune[0,2]-velo*dlamdv)^2.d0/widthL[i]^2.d0)+EXP(-(tune[0,3]-velo*dlamdv)^2.d0/widthL[i]^2.d0)+EXP(-(tune[0,4]-velo*dlamdv)^2.d0/widthL[i]^2.d0)+EXP(-(tune[0,5]-velo*dlamdv)^2.d0/widthL[i]^2.d0))
    widthL[i]  = widthL[i]*2.d0*SQRT(ALOG(2.d0)) ; FWHM instead of \sigma


    ; RCP
    rp         = REFORM(rcp0[*,richard1,richard2])*contin
    line2      = INTERPOL(rp,rlam,lam+vtest[i]*dlamdv)
    a          = WHERE(lam ge MAX(rlam))
    line2[a]   = contin
    a          = WHERE(lam le MIN(rlam))
    line2[a]   = contin
    inten2     = (filters##line2)*dlam
    FOR ii=0,ntune-1 DO inten2[0,ii]=inten2[0,ii]/norm[ii]
    c1         = REFORM(COS(x)##inten2)      ; coefficient a1 (for the cosine)
    s1         = REFORM(SIN(x)##inten2)      ; coefficient b1 (for the sine)
    c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
    s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
    phi1       = ATAN(-s1,-c1)
    phi2       = ATAN(-s2,-c2)
    vel1y      = phi1*pv1/2.d0/!dpi          ; velocity measurement
    vel2y      = phi2*pv2/2.d0/!dpi          ; velocity measurement
    res2       = INTERPOL(vtest,vel1a,vel1y)
    vel2y      = (vel2y-vel1y+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1y
    res4       = INTERPOL(vtest,vel2a,vel2y)

    c1         = -c1*dtune/period*2.d0
    s1         = -s1*dtune/period*2.d0
    c2         = -c2*dtune/period*2.d0
    s2         = -s2*dtune/period*2.d0
    widthR[i]  = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(c1^2.d0+s1^2.d0)/SQRT(c2^2.d0+s2^2.d0)))
    depthR[i]  = period/2.d0/SQRT(!dpi)*SQRT(c1^2.d0+s1^2.d0)/widthR[i]*EXP(!dpi^2.d0*widthR[i]^2.d0/period^2.d0)
    velo       = (res2+res4)/2.d0
    continuumR[i]= TOTAL(inten2)/6.d0+depthR[i]/6.d0*(EXP(-(tune[0,0]-velo*dlamdv)^2.d0/widthR[i]^2.d0)+EXP(-(tune[0,1]-velo*dlamdv)^2.d0/widthR[i]^2.d0)+EXP(-(tune[0,2]-velo*dlamdv)^2.d0/widthR[i]^2.d0)+EXP(-(tune[0,3]-velo*dlamdv)^2.d0/widthR[i]^2.d0)+EXP(-(tune[0,4]-velo*dlamdv)^2.d0/widthR[i]^2.d0)+EXP(-(tune[0,5]-velo*dlamdv)^2.d0/widthR[i]^2.d0))
    widthR[i]  = widthR[i]*2.d0*SQRT(ALOG(2.d0)) ; FWHM instead of \sigma

    ; VELOCITY + L.O.S. MAGNETIC FIELD STRENGTH
    velocity [i] = (res1+res2)/2. ;((res1+res2)+(res3+res4))/4.
    magnetic [i] = (res2-res1)/dvdb/2.d0 ;((res2-res1)+(res4-res3))/4./dvdb
    velocity2[i] = (res3+res4)/2. ;((res1+res2)+(res3+res4))/4.
    magnetic2[i] = (res4-res3)/dvdb/2.d0 ;((res2-res1)+(res4-res3))/4./dvdb

ENDFOR

a = WHERE(vtest GE -6500 AND vtest LE 6500)
width=(widthL+widthR)/2.d0
depth=(depthL+depthR)/2.d0
joseph=(1.d0+MIN(line2))/2.d0
mini  = WHERE(line2 EQ MIN(line2))
b=WHERE(lam LE lam[mini[0]] AND ABS(line2-joseph)/joseph LT 0.007*contin )
c=WHERE(lam GT lam[mini[0]] AND ABS(line2-joseph)/joseph LT 0.007*contin )
FWHM=lam[c[0]]-lam[b[0]]
!p.multi=0
SET_PLOT,'PS'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color
LOADCT,4
PLOT,vtest,(velocity+velocity2)/2.d0-vtest,xrange=[-6500,6500],xst=1,tit='!17',charsize=1.5,xtit='Input velocity (m/s)',ytit='Error (m/s)'
PLOT,vtest,(magnetic+magnetic2)/2.d0-richard2*100.*COS(richard1*10./180.*!dpi),xrange=[-6500,6500],xst=1,tit='!17',charsize=1.5,xtit='Input velocity (m/s)',ytit='Error on l.o.s. field strength (G)'
PLOT,vtest,width,xrange=[-6500,6500],xst=1,tit='!17',charsize=1.5,xtit='Input velocity (m/s)',ytit='Linewidth (FWHM) in A'
XYOUTS,-5000,0.02,'Actual linewidth='+STRTRIM(STRING(FWHM),1),charsize=1.5
OPLOT,[-6500,6500],[MEAN(width[a]),MEAN(width[a])],linestyle=2
XYOUTS,-5000,0.01,'Peak-to-peak variation='+STRTRIM(STRING((MAX(width[a])-MIN(width[a]))/MEAN(width[a])),1),charsize=1.5
continuum=(continuumL+continuumR)/2.d0
PLOT,vtest,continuum,xrange=[-6500,6500],xst=1,tit='!17',charsize=1.5,xtit='Input velocity (m/s)',ytit='Continuum intensity (a.u.)',yrange=[0.88*contin,contin*1.015],yst=1
XYOUTS,-5000,0.885*contin,'Peak-to-peak variation='+STRTRIM(STRING((MAX(continuum[a])-MIN(continuum[a]))/MEAN(continuum[a])),1),charsize=1.5
OPLOT,[-6500,6500],[MEAN(continuum[a]),MEAN(continuum[a])],linestyle=2
XYOUTS,-5000,0.9*contin,'Actual continuum='+STRTRIM(STRING(contin),1),charsize=1.5
PLOT,vtest,depth/continuum,xrange=[-6500,6500],xst=1,tit='!17',charsize=1.5,xtit='Input velocity (m/s)',ytit='Linedepth'
XYOUTS,-5000,0.2,'Actual linedepth='+STRTRIM(STRING((contin-MIN(line2))/contin),1),charsize=1.5
XYOUTS,-5000,0.1,'Peak-to-peak variation='+STRTRIM(STRING((MAX(depth[a])-MIN(depth[a]))/MEAN(depth[a])),1),charsize=1.5
OPLOT,[-6500,6500],[MEAN(depth[a]/continuum[a]),MEAN(depth[a]/continuum[a])],linestyle=2
DEVICE,/close
SET_PLOT,'x'

SET_PLOT,'x'
!P.MULTI=0
READ,pause
endo:



        ;GOTO,endo2
velocity   = FLTARR(40,10)
velocity2  = FLTARR(40,10)
magnetic   = FLTARR(40,10)
magnetic2  = FLTARR(40,10)
widthL     = FLTARR(40,10)
widthR     = FLTARR(40,10)
truelinewidth= FLTARR(40,10)
depthL     = widthL
depthR     = widthR
continuumL = widthL
continuumR = widthR
width0     = FLTARR(40,10)

voffset    = 0

FOR richard2   =0,39 DO BEGIN   ; TO TEST LINEWIDTH
    FOR richard1 = 0,9 DO BEGIN
        
        
                                ; 1ST PASSAGE TO ESTIMATE VELOCITY
        rp         = REFORM(lcp0[*,richard1,richard2]) ; LCP

        line2      = INTERPOL(rp,rlam,lam+voffset*dlamdv)
        a          = WHERE(lam ge MAX(rlam))
        line2[a]   = 1.d0
        a          = WHERE(lam le MIN(rlam))
        line2[a]   = 1.d0
        
        a          = WHERE(line2 LE (1.d0+MIN(line2))/2.d0,na)
        width0[richard2,richard1]=INT_TABULATED(lam,1.d0-line2) ;(1.d0-min(line2))
        PLOT,lam,line2,xst=1,xrange=[-1,1]
        inten2     = filters## line2
        c1         = REFORM(COS(x)##inten2) ; coefficient a1 (for the cosine)
        s1         = REFORM(SIN(x)##inten2) ; coefficient b1 (for the sine)
        c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
        s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
        phi1       = ATAN(-s1,-c1)
        phi2       = ATAN(-s2,-c2)
        vel1x      = phi1*pv1/2.d0/!dpi ; velocity measurement
        vel2x      = phi2*pv2/2.d0/!dpi ; velocity measurement
        res1       = INTERPOL(vtest,vel1a,vel1x)
        vel2x      = (vel2x-vel1x+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1x
        res3       = INTERPOL(vtest,vel2a,vel2x)
        
        
        rp         = REFORM(rcp0[*,richard1,richard2]) ; RCP
        line2      = INTERPOL(rp,rlam,lam+voffset*dlamdv)
        a          = WHERE(lam ge MAX(rlam))
        line2[a]   = 1.d0
        a          = WHERE(lam le MIN(rlam))
        line2[a]   = 1.d0 
        
        OPLOT,lam,line2,linestyle=2
        inten2     = filters## line2
        c1         = REFORM(COS(x)##inten2) ; coefficient a1 (for the cosine)
        s1         = REFORM(SIN(x)##inten2) ; coefficient b1 (for the sine)
        c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
        s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
        phi1       = ATAN(-s1,-c1)
        phi2       = ATAN(-s2,-c2)
        vel1y      = phi1*pv1/2.d0/!dpi ; velocity measurement
        vel2y      = phi2*pv2/2.d0/!dpi ; velocity measurement
        res2       = INTERPOL(vtest,vel1a,vel1y)
        vel2y      = (vel2y-vel1y+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1y
        res4       = INTERPOL(vtest,vel2a,vel2y)
                  
                    ; 2nd PASSAGE TO CORRECT FOR THE VELOCITIES

        rp         = REFORM(lcp0[*,richard1,richard2]) ;LCP
        vLCP       = 0.         ;(res1[0]+res3[0])/2.d0
        line2      = INTERPOL(rp,rlam,lam-vLCP*dlamdv)
        a          = WHERE(lam ge MAX(rlam))
        line2[a]   = 1.d0
        a          = WHERE(lam le MIN(rlam))
        line2[a]   = 1.d0 
        a          = WHERE(line2 LE (1.d0+MIN(line2))/2.d0,na)
        width0[richard2,richard1]=INT_TABULATED(lam,1.d0-line2) ;(1.d0-min(line2))
        OPLOT,lam,line2,col=180
        inten2     = filters## line2
        inten2     = inten2/MEAN(norm)*dlam ; average integral of filters
        c1         = REFORM(COS(x)##inten2) ; coefficient a1 (for the cosine)
        s1         = REFORM(SIN(x)##inten2) ; coefficient b1 (for the sine)
        c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
        s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
        phi1       = ATAN(-s1,-c1)
        phi2       = ATAN(-s2,-c2)
        vel1x      = phi1*pv1/2.d0/!dpi ; velocity measurement
        vel2x      = phi2*pv2/2.d0/!dpi ; velocity measurement
        res1       = INTERPOL(vtest,vel1a,vel1x)
        vel2x      = (vel2x-vel1x+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1x
        res3       = INTERPOL(vtest,vel2a,vel2x)
        
        
        c1         = -c1*dtune/period*2.d0
        s1         = -s1*dtune/period*2.d0
        c2         = -c2*dtune/period*2.d0
        s2         = -s2*dtune/period*2.d0
        
        
        widthL[richard2,richard1] = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(c1^2.d0+s1^2.d0)/SQRT(c2^2.d0+s2^2.d0)))
        depthL[richard2,richard1] = period/2.d0/SQRT(!dpi)*SQRT(c1^2.d0+s1^2.d0)/widthL[richard2,richard1]*exp(!dpi^2.d0*widthL[richard2,richard1]^2.d0/period^2.d0)
        continuumL[richard2,richard1]= TOTAL(inten2)/6.d0+depthL[richard2,richard1]/6.d0*(EXP(-(tune[0,0]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0)+EXP(-(tune[0,1]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0)+EXP(-(tune[0,2]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0)+EXP(-(tune[0,3]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0)+EXP(-(tune[0,4]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0)+EXP(-(tune[0,5]-res1*dlamdv)^2.d0/widthL[richard2,richard1]^2.d0))
        widthL[richard2,richard1] = widthL[richard2,richard1]*2.d0*SQRT(ALOG(2.d0))
        
        rp         = REFORM(rcp0[*,richard1,richard2]) ;RCP
        vRCP       = 0.         ;(res2[0]+res4[0])/2.d0
        line2      = INTERPOL(rp,rlam,lam-vRCP*dlamdv)
        a          = WHERE(lam ge MAX(rlam))
        line2[a]   = 1.d0
        a          = WHERE(lam le MIN(rlam))
        line2[a]   = 1.d0 
        OPLOT,lam,line2,col=180,linestyle=2
        inten2     = filters## line2
        inten2     = inten2/MEAN(norm)*dlam ; average integral of filters
        c1         = REFORM(COS(x)##inten2) ; coefficient a1 (for the cosine)
        s1         = REFORM(SIN(x)##inten2) ; coefficient b1 (for the sine)
        c2         = REFORM(COS(2.d0*x)##inten2) ; coefficient a2 (for the cosine)
        s2         = REFORM(SIN(2.d0*x)##inten2) ; coefficient b2 (for the sine)
        phi1       = ATAN(-s1,-c1)
        phi2       = ATAN(-s2,-c2)
        vel1y      = phi1*pv1/2.d0/!dpi ; velocity measurement
        vel2y      = phi2*pv2/2.d0/!dpi ; velocity measurement
        res2       = INTERPOL(vtest,vel1a,vel1y)
        vel2y      = (vel2y-vel1y+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1y
        res4       = INTERPOL(vtest,vel2a,vel2y)
        
        c1         = -c1*dtune/period*2.d0
        s1         = -s1*dtune/period*2.d0
        c2         = -c2*dtune/period*2.d0
        s2         = -s2*dtune/period*2.d0
        
        widthR[richard2,richard1] = period/!dpi*SQRT(1.d0/3.d0*ALOG(SQRT(c1^2.d0+s1^2.d0)/SQRT(c2^2.d0+s2^2.d0)))
        depthR[richard2,richard1] = period/2.d0/SQRT(!dpi)*SQRT(c1^2.d0+s1^2.d0)/widthR[richard2,richard1]*exp(!dpi^2.d0*widthR[richard2,richard1]^2.d0/period^2.d0)
        continuumR[richard2,richard1]= TOTAL(inten2)/6.d0+depthR[richard2,richard1]/6.d0*(EXP(-(tune[0,0]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0)+EXP(-(tune[0,1]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0)+EXP(-(tune[0,2]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0)+EXP(-(tune[0,3]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0)+EXP(-(tune[0,4]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0)+EXP(-(tune[0,5]-res2*dlamdv)^2.d0/widthR[richard2,richard1]^2.d0))
        
        widthR[richard2,richard1] = widthR[richard2,richard1]*2.d0*SQRT(ALOG(2.d0))
        
        
        velocity [richard2,richard1] = res1 ;(res1+res2)/2.
        magnetic [richard2,richard1] = res2 ;(res2-res1)/dvdb/2.d0
        velocity2[richard2,richard1] = res3 ;(res3+res4)/2.
        magnetic2[richard2,richard1] = res4 ;(res4-res3)/dvdb/2.d0
        
    ENDFOR
ENDFOR

read,pause

!p.multi=0
SET_PLOT,'PS'


DEVICE,file='yo3b.ps',xoffset=0,yoffset=0,xsize=22,ysize=14,/color
loadct,3
plot,(widthL[0:39,7]+widthR[0:39,7])/2.d0,velocity[0:39,7]-voffset,xst=1,yst=1,charsize=1.5,tit='!17',xtit='Estimated Linewidth (A)',ytit='Error on velocity (m/s)',yrange=[-500,200] ;,thick=10,xrange=[0.015,0.165],yrange=[-10,2500]
FOR i=0,6 DO oplot,(widthL[0:39,i]+widthR[0:39,i])/2.d0,velocity[0:39,i]-voffset ;,col=180
xyouts,0.025,-10,'Velocity='+STRTRIM(string(voffset),1),charsize=1.5

plot,(depthL[0:39,7]+depthR[0:39,7])/2.d0,velocity[0:39,7]-voffset,xst=1,yst=1,charsize=1.5,tit='!17',xtit='Estimated Linewidth (A)',ytit='Error on velocity (m/s)' ;,thick=10,xrange=[0.015,0.165],yrange=[-10,2500]
FOR i=0,6 DO oplot,(depthL[0:39,i]+depthR[0:39,i])/2.d0,velocity[0:39,i]-voffset ;,col=180

DEVICE,/CLOSE


DEVICE,file='yo3.ps',xoffset=0,yoffset=0,xsize=22,ysize=28,/color
loadct,4

richard1=0

plot,findgen(40)*100.*COS(richard1*10./180.*!dpi),magnetic[*,richard1],charsize=1.5,tit='!17',xtit='l.o.s. Magnetic Field Strength (Gauss)',ytit='l.o.s. Magnetic Field Strength (Gauss)',xst=1
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),magnetic[i,*,richard1],col=i/(nwidth-1.)*255.
oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),magnetic2[*,richard1],linestyle=3
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),magnetic2[i,*,richard1],col=i/(nwidth-1.)*255.,linestyle=3
oplot,[0,4000],[0,4000],linestyle=2,thick=2
;xyouts,250,3500,'offset velocity
;(m/s)='+STRTRIM(STRING(voffset),1),charsize=1.5
labels=STRARR(nwidth)
FOR i=0,nwidth-1 DO labels[i]='FWHM*'+STRTRIM(STRING((15+nwidth-10-i)/10.0),1)
legend,[labels,'First Fourier coefficient','Second Fourier coefficient'],linestyle=[0,0,0,0,0,0,0,0,0,0,0,3],color=[0,28,57,85,113,142,170,198,226,254,0,0],thick=2,charsize=1.5,/right,/bottom

plot,findgen(40)*100.*COS(richard1*10./180.*!dpi),velocity[*,richard1],charsize=1.5,tit='!17',xtit='l.o.s Magnetic Field Strength (Gauss)',ytit='Doppler velocity (m/s)',xst=1,yrange=[-150+voffset,150+voffset],yst=1
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),velocity[i,*,richard1],col=i/(nwidth-1.)*255.
oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),velocity2[*,richard1],linestyle=3
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),velocity2[i,*,richard1],col=i/(nwidth-1.)*255.,linestyle=3
oplot,[0,4000],[0,0],linestyle=2,thick=2


plot,findgen(40)*100.*COS(richard1*10./180.*!dpi),0.5*(magnetic[*,richard1]+magnetic2[*,richard1]),charsize=1.5,tit='!17',xtit='l.o.s Magnetic Field Strength (Gauss)',ytit='l.o.s Magnetic Field Strength (Gauss)',xst=1
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),0.5*(magnetic[i,*,richard1]+magnetic2[i,*,richard1]),col=i/(nwidth-1.)*255.
oplot,[0,4000],[0,4000],linestyle=2,thick=2

plot,findgen(40)*100.*COS(richard1*10./180.*!dpi),0.5*(velocity[*,richard1]+velocity2[*,richard1]),charsize=1.5,tit='!17',xtit='l.o.s Magnetic Field Strength (Gauss)',ytit='Doppler velocity (m/s)',xst=1,yrange=[-150+voffset,150+voffset],yst=1
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),0.5*(velocity[i,*,richard1]+velocity2[i,*,richard1]),col=i/(nwidth-1.)*255.
oplot,[0,4000],[0,0],linestyle=2,thick=2

plot,findgen(40)*100.*COS(richard1*10./180.*!dpi),(widthL[*,richard1]+widthR[*,richard1])/2.d0,charsize=1.5,tit='!17',xtit='l.o.s. Magnetic Field Strength (Gauss)',ytit='Linewidth (A)',xst=1,yrange=[0.,0.15],yst=1
FOR i=1,nwidth-1 DO oplot,findgen(40)*100.*COS(richard1*10./180.*!dpi),(widthL[i,*,richard1]+widthR[i,*,richard1])/2.d0,col=i/(nwidth-1.)*255.
legend,labels,linestyle=0,color=[0,28,57,85,113,142,170,198,226,254],thick=2,charsize=1.5,/right,/bottom


plot,findgen(40)*100.*COS(0.*10./180.*!dpi),smooth(width0[5,*,0],2),xst=1,yrange=[0.04,0.15],yst=1
oplot,findgen(40)*100.*COS(0.*10./180.*!dpi),smooth((widthL[5,*,0]+widthR[5,*,0])/2.d0,2),linestyle=2
FOR i=1,9 DO oplot,findgen(40)*100.*COS(i*10./180.*!dpi),smooth(width0[5,*,i],2),col=i/9.*255.
FOR i=1,9 DO oplot,findgen(40)*100.*COS(i*10./180.*!dpi),smooth((widthL[5,*,i]+widthR[5,*,i])/2.d0,2),col=i/9.*255.,linestyle=2


DEVICE,/close
SET_PLOT,'x'

SET_PLOT,'x'
!P.MULTI=0
endo2:



        
       ;SAVE,vtest,vel1a,vel2a,FILE='lookuptable'+STRTRIM(STRING(inclination),1)+'_'+STRTRIM(STRING(field),1)+'_2bin' ; 2bin is for RCP, bin is for LCP !!!
       
read,pause

;data=FLTARR(31,ntest,2)
;for j=0,30 do begin & restore,'lookuptable0_'+STRTRIM(STRING(j),1)+'_bin' & vel1a0=vel1a & restore,'lookuptable0_'+STRTRIM(STRING(j),1)+'_2bin' & data[j,*,0]=vel1a[*] & data[j,*,1]=vel1a0[*] & endfor
;set_plot,'ps'
;!p.multi=0
;device,file='yo.ps',/color,xoffset=0,yoffset=0,xsize=20,ysize=24,bits=24
;tvim,data[*,*,0],/scale,aspect=1,tit='!17',xtit='l.o.s magnetic field (Gauss)',xrange=[0,3000],yrange=[vtest[0],vtest[900]],ytit='Input velocity (m/s)',pcharsize=1.5
;tvim,data[*,*,1],/scale,aspect=1,tit='!17',xtit='l.o.s magnetic field (Gauss)',xrange=[0,3000],yrange=[vtest[0],vtest[900]],ytit='Input velocity (m/s)',pcharsize=1.5
;device,/close


;SET_PLOT,'PS'
;!P.MULTI=0
;device,file='yo2.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,/color
;loadct,4
;restore,'lookuptable0_0_bin' 
;vel1a0=vel1a
;vel2a0=vel2a
;restore,'lookuptable0_0_2bin' 
;plot,vtest,(vel2a+vel2a0)/2.d0,xst=1,charsize=1.5,tit='!17',xtit='input velocity (m/s)',ytit='output velocity (m/s)',yrange=[-6000,6000],yst=1,thick=3,xrange=[-7000,7000]
;for i=0,4 do for j=1,30 do begin &  restore,'lookuptable'+STRTRIM(STRING(i*2),1)+'_'+STRTRIM(STRING(j),1)+'_bin' & vel1a0=vel1a & vel2a0=vel2a & restore,'lookuptable'+STRTRIM(STRING(i*2),1)+'_'+STRTRIM(STRING(j),1)+'_2bin' & IF( j MOD 10 EQ 0) THEN oplot,vtest,(vel2a+vel2a0)/2.d0,col=j*5,linestyle=i & endfor
;legend,['0 Gauss','1000 Gauss','2000 Gauss','3000 Gauss','0 degrees','20 degrees','40 degrees','60 degrees','80 degrees'],linestyle=[0,0,0,0,0,1,2,3,4],color=[0,50,100,150,0,0,0,0,0],charsize=1.5,thick=2
;device,/close
;set_plot,'x'

END
