; PROGRAM TO PERFORM A MDI-LIKE ALGORITHM
; USING 2D LOOK-UP TABLES (a1 LCP=f(B,V) AND a1 RCP=f(B,V))
; THIS CODE RETURNS, FOR a1LCP AND a1RCP THE VALUE OF DOPPLER
; VELOCITY AND L.O.S MAGNETIC FIELD

PRO lookup

SET_PLOT,'x'
!P.MULTI=0
WINDOW,0,RETAIN=2

ntest   = 901 ; number of input velocities for the look-up table
nB      = 34  ; number of l.o.s magnetic field strengths (for 0 inclination)
data    = FLTARR(nB,ntest,2) ; look-up tables
lam0    = 6173.3433d0
well    = 200000.d0;150000.d0 ; CCD full well

ntune   = 6      ; Set number of tuning positions (NUMBER OF FILTERS)
ntunep  = 11     ; All positions for the purpose of plotting.
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)
dvdb    = 1./2.13d12*lam0*2.5d0*299792458.d0  ; conversion from magnetic field to Doppler velocity

;-----------------------------------------------------------------------
; BUILD HMI FILTERS
;-----------------------------------------------------------------------

; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
;phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
;contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
phase      = FLTARR(7)
contrast   = FLTARR(7)+1.0

wmich      = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]          ; FSRs of the tunable elements
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements

dtune  = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune


nlamp  = 30061 ; Set wavelength points for line profile
lamp   = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

; SOLAR LINE
RESTORE,'/scr20/richard/hmi/line_profile_library.sav'
rlam          = lambda-6173.3433d0

; GRID WE WANT IN WAVELENGTH
nlam          = 38500.
dlam          = 1.d0/1.75d3
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = lam0/299792458.d0
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

blocker       = 0.57471215*exp(-lam^2.d0/4.5d0^2.d0) ;perfect blocker

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

;-------------------------------------------
; BUILDS 2D LOOK-UP TABLES
; REQUIRES HMIrichard.PRO TO HAVE BEEN RUN BEFORE 
;------------------------------------------

for j=0,nB-1 do begin 
restore,'lookuptable0_'+STRTRIM(STRING(j),1)+'_bin'
vel1a0=vel1a
restore,'lookuptable0_'+STRTRIM(STRING(j),1)+'_2bin'
data[j,*,0]=vel1a[*]
data[j,*,1]=vel1a0[*]
endfor ; vtest is in the .bin

;data[*,*,0] = RCP
;data[*,*,1] = LCP

M = 3000;
L1= FINDGEN(M)/(M-1.)*(MAX(data)-MIN(data))+MIN(data)

dataBL=FLTARR(ntest,M)  ; B=f(v,LCP)
dataRL=FLTARR(ntest,M)  ; RCP=f(v,LCP)
datavL=FLTARR(nB,M)     ; v=f(B,LCP)
dataRL2=FLTARR(nB,M)    ; RCP=f(B,LCP)



FOR i=0,ntest-1 DO BEGIN
    ;mini=MIN(data[*,i,1])
    ;maxi=MAX(data[*,i,1])
    ;a=WHERE(L1 GE mini AND L1 LE maxi)
    dataBL[i,*]=INTERPOL(FINDGEN(nB)*100.d0,data[*,i,1],L1[*],/SPLINE) ;data[*,i,0] is monotonic ;dataBL=f(B,L1)
ENDFOR


FOR i=0,nB-1 DO BEGIN
    ;mini=MIN(data[i,*,1])
    ;maxi=MAX(data[i,*,1])
    ;a=WHERE(L1 GE mini AND L1 LE maxi)
    datavL[i,*]=INTERPOL(vtest,data[i,*,1],L1[*],/SPLINE) ;datavL=f(B,L1)
    ;PLOT,L1,datavL[i,*],xst=1
    ;OPLOT,data[i,*,1],vtest,col=180
    ;PRINT,i
    ;read,pause
ENDFOR


FOR i=0,ntest-1 DO BEGIN
   ;a=WHERE(dataBL[i,*] GE 0. AND dataBL[i,*] LE 3000.)
   temp=dataBL[i,*]
   dataRL[i,*]=INTERPOL(data[*,i,0],FINDGEN(nB)*100.d0,temp[*],/SPLINE) ;dataRL=f(v,L1)
ENDFOR

FOR i=0,nB-1 DO BEGIN
   ;a=WHERE(datavL[i,*] GE vtest[0] AND datavL[i,*] LE vtest[ntest-1])
   temp=REFORM(datavL[i,*])
   dataRL2[i,*]=INTERPOL(data[i,*,0],vtest,temp[*],/SPLINE) ;dataRL=f(v,L1)
ENDFOR


datavLR=FLTARR(M,M) ; v=f(LCP,RCP)
dataBLR=FLTARR(M,M) ; B=f(LCP,RCP)


FOR i=0,M-1 DO BEGIN
    ;mini=MIN(dataRL[*,i])
    ;maxi=MAX(dataRL[*,i])
    ;a=WHERE(L1 GE mini AND L1 LE maxi)
    datavLR[i,*]=INTERPOL(vtest,dataRL[*,i],L1[*],/SPLINE) ;datavLR=f(LCP,RCP)

    ;mini=MIN(dataRL2[*,i])
    ;maxi=MAX(dataRL2[*,i])
    ;a=WHERE(L1 GE mini AND L1 LE maxi)
    dataBLR[i,*]=INTERPOL(FINDGEN(nB)*100.d0,dataRL2[*,i],L1[*],/SPLINE) ;dataBLR=f(LCP,RCP)
ENDFOR

; my two 2D inverse look-up tables are:
; datavLR and dataBLR

voffset=0.d0; v=offset in m/s

errorv=FLTARR(nB,10)
errorB=FLTARR(nB,10)

FOR field=0,nB-1 DO BEGIN ;30
FOR inclination=0,9 DO BEGIN ;9

rp  = REFORM(lcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables
line= INTERPOL(rp,rlam,lam+voffset*dlamdv)
a   = WHERE(lam ge MAX(rlam))
line[a] = 1.d0
a   = WHERE(lam le MIN(rlam))
line[a] = 1.d0
PLOT,lam,line,xrange=[-1,1],xst=1


; INTENSITIES
inten      = filters##line
x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune
c1         = REFORM(COS(x)##inten)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)##inten)  ; coefficient b1 (for the sine)
c2         = REFORM(COS(2.d0*x)##inten)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)##inten)  ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0   ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi; velocity measurement

rp  = REFORM(rcp0[*,inclination,field]) ; solar line used as the reference to build the look-up tables
line= INTERPOL(rp,rlam,lam+voffset*dlamdv)
a   = WHERE(lam ge MAX(rlam))
line[a] = 1.d0
a   = WHERE(lam le MIN(rlam))
line[a] = 1.d0
OPLOT,lam,line,col=180,linestyle=2,thick=3
inten      = filters##line
c1         = REFORM(COS(x)##inten)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)##inten)  ; coefficient b1 (for the sine)
c2         = REFORM(COS(2.d0*x)##inten)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)##inten)  ; coefficient b2 (for the sine)
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel3       = phi1*pv1/2.d0/!dpi; velocity measurement
vel4       = phi2*pv2/2.d0/!dpi; velocity measurement
IX=INTERPOL(FINDGEN(M),L1,vel1,/SPLINE)
IY=INTERPOL(FINDGEN(M),L1,vel3,/SPLINE)
velocity=BILINEAR(datavLR,IX,IY)
magnetic=BILINEAR(dataBLR,IX,IY)
PRINT,inclination*10.,field*100.d0,voffset,velocity,magnetic
errorv[field,inclination]=velocity;(velocity-voffset)
errorB[field,inclination]=magnetic;-field*100.d0*COS(inclination*10.d0/180.d0*!dpi))

ENDFOR
ENDFOR

SET_PLOT,'PS'
!P.MULTI=0
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=22,ysize=28,/color
LOADCT,4
plot,100.d0*findgen(nB),errorB[*,0],xrange=[0,3900],yrange=[0,4000],yst=1,tit='!17',xtit='l.o.s Magnetic Field Strength (Gauss)',ytit='l.o.s Magnetic Field Strength (Gauss)',xst=1,charsize=1.5
FOR i=1,7 DO oplot,findgen(nB)*100.*COS(i*10./180.*!dpi),errorB[*,i],col=i/9.*255.
oplot,[0,4000],[0,4000],linestyle=2,thick=2

plot,100.d0*findgen(nB),errorv[*,0],xrange=[0,3900],yrange=[-150+voffset,150+voffset],yst=1,tit='!17',xtit='l.o.s Magnetic Field Strength (Gauss)',ytit='Doppler velocity (m/s)',xst=1,charsize=1.5
FOR i=1,7 DO oplot,findgen(nB)*100.*COS(i*10./180.*!dpi),errorv[*,i],col=i/9.*255.
oplot,[0,4000],[0,0],linestyle=2,thick=2


DEVICE,/CLOSE
SET_PLOT,'X'

read,pause

END
