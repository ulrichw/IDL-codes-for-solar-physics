PRO figurespapier_HMI


SET_PLOT,'ps'
!P.MULTI=0



; PHASE AND CONTRAST FOR NON-TUNABLE ELEMENTS
;------------------------------------------------------------------------------------------------------------------

nx=256
anglim=965
xcenter     = nx/2
ycenter     = nx/2
distance    = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx

RESTORE,'RESULTS/RESULTS_June09_710660_CAL_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN' ; RECALCULATED IN NOVEMBER 2009
phigCal=phig
BgCal=Bg
;RESTORE,'RESULTS/RESULTS_June09_710600_OBS_'+STRTRIM(STRING(LONG(nx)),1)+'_SIDE.BIN'
;phigObs=phig
;BgObs=Bg
phigObs=phigCal
BgObs=BgCal


a=WHERE(distance GT anglim OR lateconv NE 1,COMPLEMENT=b)

FOR i=3,6 DO BEGIN & temp=REFORM(PhigObs[*,*,i]) & temp[a]=-10000.d0 & phigObs[*,*,i]=temp[*,*] & temp=REFORM(BgObs[*,*,i]) & temp[a]=-10000.d0 & BgObs[*,*,i]=temp[*,*] & phigObs[0:1,*,i]=-10000.d0 & BgObs[0:1,*,i]=-10000.d0 & ENDFOR
FOR i=3,6 DO BEGIN & temp=REFORM(PhigCal[*,*,i]) & temp[a]=-10000.d0 & phigCal[*,*,i]=temp[*,*] & temp=REFORM(BgCal[*,*,i]) & temp[a]=-10000.d0 & BgCal[*,*,i]=temp[*,*] & phigCal[0:1,*,i]=-10000.d0 & BgCal[0:1,*,i]=-10000.d0 & ENDFOR


; ADD 180 DEGREES TO PHASES < -180
FOR i=3,6 DO BEGIN
temp = REFORM(PhigObs[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigObs[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigObs[*,*,i]=temp
    ENDELSE
temp = REFORM(PhigCal[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigCal[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigCal[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moyobs=DBLARR(4)
moybobs=moyobs
moycal=DBLARR(4)
moybcal=moycal
FOR i=3,6 DO BEGIN & temp=REFORM(Phigobs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyobs[i-3]=MEAN(temp[a]) & temp[b]=moyobs[i-3] & Phigobs[*,*,i]=temp & temp=REFORM(Bgobs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moybobs[i-3]=MEAN(temp[a]) & temp[b]=moybobs[i-3] & Bgobs[*,*,i]=temp  & ENDFOR
FOR i=3,6 DO BEGIN & temp=REFORM(Phigcal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moycal[i-3]=MEAN(temp[a]) & temp[b]=moycal[i-3] & Phigcal[*,*,i]=temp & temp=REFORM(Bgcal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moybcal[i-3]=MEAN(temp[a]) & temp[b]=moybcal[i-3] & Bgcal[*,*,i]=temp  & ENDFOR

PRINT,'AVERAGE VALUES FOR NON-TUNABLE ELEMENTS'
PRINT,moycal
PRINT,moybcal


; WE PLOT THE RESULT
!P.MULTI=[0,2,2]
device,file='phasesntObs.ps',bits=24,xoffset=0.5,yoffset=1,xsize=21,ysize=19,/color
LOADCT,4
tvim,phigobs[*,*,3],/scale,tit='!17E2 '+STRTRIM(STRING(moyobs[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigobs[*,*,3]),MAX(phigobs[*,*,3])]
tvim,phigobs[*,*,4],/scale,tit='!17E3 '+STRTRIM(STRING(moyobs[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigobs[*,*,4]),MAX(phigobs[*,*,4])]
tvim,phigobs[*,*,5],/scale,tit='!17E4 '+STRTRIM(STRING(moyobs[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigobs[*,*,5]),MAX(phigobs[*,*,5])]
tvim,phigobs[*,*,6],/scale,tit='!17E5 '+STRTRIM(STRING(moyobs[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigobs[*,*,6]),MAX(phigobs[*,*,6])]
DEVICE,/close

!P.MULTI=[0,2,2]
device,file='contrastsObs.ps',bits=24,xoffset=0.5,yoffset=1,xsize=21,ysize=19,/color
tvim,Bgobs[*,*,3]  ,/scale,tit='!17E2 '+STRTRIM(STRING(moybobs[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgobs[*,*,3]),MAX(Bgobs[*,*,3])]
tvim,Bgobs[*,*,4]  ,/scale,tit='!17E3 '+STRTRIM(STRING(moybobs[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgobs[*,*,4]),MAX(Bgobs[*,*,4])]
tvim,Bgobs[*,*,5]  ,/scale,tit='!17E4 '+STRTRIM(STRING(moybobs[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgobs[*,*,5]),MAX(Bgobs[*,*,5])]
tvim,Bgobs[*,*,6]  ,/scale,tit='!17E5 '+STRTRIM(STRING(moybobs[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgobs[*,*,6]),MAX(Bgobs[*,*,6])]
DEVICE,/close

!P.MULTI=[0,2,2]
device,file='phasesntCal.ps',bits=24,xoffset=0.5,yoffset=1,xsize=21,ysize=19,/color
LOADCT,4
tvim,phigcal[*,*,3],/scale,tit='!17E2 '+STRTRIM(STRING(moycal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigcal[*,*,3]),MAX(phigcal[*,*,3])]
tvim,phigcal[*,*,4],/scale,tit='!17E3 '+STRTRIM(STRING(moycal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigcal[*,*,4]),MAX(phigcal[*,*,4])]
tvim,phigcal[*,*,5],/scale,tit='!17E4 '+STRTRIM(STRING(moycal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigcal[*,*,5]),MAX(phigcal[*,*,5])]
tvim,phigcal[*,*,6],/scale,tit='!17E5 '+STRTRIM(STRING(moycal[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigcal[*,*,6]),MAX(phigcal[*,*,6])]
DEVICE,/close

!P.MULTI=[0,2,2]
device,file='contrastsCal.ps',bits=24,xoffset=0.5,yoffset=1,xsize=21,ysize=19,/color
tvim,Bgcal[*,*,3]  ,/scale,tit='!17E2 '+STRTRIM(STRING(moybcal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgcal[*,*,3]),MAX(Bgcal[*,*,3])]
tvim,Bgcal[*,*,4]  ,/scale,tit='!17E3 '+STRTRIM(STRING(moybcal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgcal[*,*,4]),MAX(Bgcal[*,*,4])]
tvim,Bgcal[*,*,5]  ,/scale,tit='!17E4 '+STRTRIM(STRING(moybcal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgcal[*,*,5]),MAX(Bgcal[*,*,5])]
tvim,Bgcal[*,*,6]  ,/scale,tit='!17E5 '+STRTRIM(STRING(moybcal[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bgcal[*,*,6]),MAX(Bgcal[*,*,6])]
DEVICE,/close


; SPATIALLY AVERAGED NON-TUNABLE TRANSMISSION PROFILE
;------------------------------------------------------------------------------------------------------------------
SET_PLOT,'ps'
nlam        = 36000               ; number of wavelengths
dlam        = 0.25d0/1.d3         ; resolution in wavelength
lam         = dlam*(DINDGEN(nlam)-(nlam-1)/2.)
FSR         = DBLARR(7)    ; FSR in Angstrom
FSR[0]      = 0.172-0.0010576d0; for the narrow-band Michelson
FSR[1]      = 0.344-0.00207683d0; for the broad-band  Michelson
FSR[2]      = 0.693+0.000483467d0  ; for E1
FSR[3]      = 1.407d0;1.405d0       ; for E2
FSR[4]      = 2.779d0               ; for E3
FSR[5]      = 5.682d0               ; for E4
FSR[6]      = 11.354d0              ; for E5
lamref      = 6173.3433d0
transmission = FLTARR(2,401)
OPENR,1,'frontwindow3.txt'
READF,1,transmission
CLOSE,1
wavelength   = REFORM(transmission[0,*])
transmission = REFORM(transmission[1,*])
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)

q            = READFITS('blocker11.fits') ; blocker filter profile from http://www.lmsal.com/~shine/Public/hmi/blocker11.fits
blocker0    = blocker0*INTERPOL(q[*,1]/100.d0,q[*,0]+2.6-lamref,lam) ; MODIF OF JUNE 2009

Bg          = [0.97723184d0,0.99716162d0,0.97434975d0,0.9845d0,0.9643d0,0.9874d0,0.998987d0]
Phig        = [0.d0,0.d0,0.d0,-5.5004062d0,0.34564038d0,-5.8693346d0,0.84729495d0]*!dpi/180.d0 ; from RESULTS/RESULTS_June09_710660_CAL_256_FRONT.BIN
;Phig        = [0.d0,0.d0,0.d0,-5.7121230d0,0.21370590d0,-5.8707862d0,0.69710796d0]*!dpi/180.d0 ; from RESULTS/RESULTS_June09_710660_CAL_256_SIDE.BIN
;Phig        = [0.d0,0.d0,0.d0,-5.5975558d0,0.35280715d0,-5.8704731d0,0.80044277d0]*!dpi/180.d0 ; from RESULTS/RESULTS_June09_710660_CAL_256_FRONT2.BIN
;Phig        = [0.d0,0.d0,0.d0,-5.8216680d0,0.21216013d0,-5.8796427d0,0.64316043d0]*!dpi/180.d0 ; from RESULTS/RESULTS_June09_710660_CAL_256_SIDE2.BIN

; WHEN BLOCKER IS AT +3.59
;Bg          = [0.97723184d0,0.99716162d0,0.97434975d0,0.9874d0,0.9655d0,0.99097d0,1.01344d0]
;Phig        = [0.d0,0.d0,0.d0,-5.997d0,-0.08763d0,-4.46157d0,9.00517d0]*!dpi/180.d0 


profilef    = blocker0
FOR i=3,6 DO profilef = profilef*(1.d0+Bg[i]*COS(2.d0*!dpi*lam/FSR[i]+Phig[i]))/2.d0
profilef=profilef-4.5411594e-05 ; THRESHOLD REBINNED TO (1,1)

a=WHERE(profilef EQ MAX(profilef))
PRINT,'MAX NON-TUNABLE PROFILE AT:',lam[a]

RESTORE,'DATA_NONTUNABLE_CAL.BIN'

!P.MULTI=0
RESTORE,'LyotLampData.bin'
LOADCT,4
DEVICE,file='averageprofile.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color
PLOT,lam,profilef*6.2680245/3.2416076,xrange=[-3.5,3.5],xst=1,tit='!17',xtit='Wavelength (A)',ytit='Normalized transmittance',thick=2,/ylog,yrange=[0.0001,1.1],yst=1
OPLOT,lresp+lam[a[0]],data,col=180,thick=2
OPLOT,lam0,8.0*inten2/27./3.2416076,psym=1,thick=4,col=100
DEVICE,/CLOSE




READ,pause

; ANGULAR DEPENDENCE TESTS; CALMODE DATA; FEB 2007; VACUUM; OVEN AT
; 29.8 C
; FOR E1 ONLY
;------------------------------------------------------------------------------------------------------------------


tr=[-0.69996291,-0.28368565,0.0,0.28368565,0.69996291]; angle of incidence (obtained assuming that small field stop = 175" radius)
tl=[-0.74148435,-0.31073995,0.0,0.31073995,0.74148435]
td=[-0.89334183,-0.40869834,0.0,0.40869834,0.89334183]
tu=[-0.66808730,-0.29886431,0.0,0.29886431,0.66808730]
pr=[-156.5,-141.5,-137.,-141.5,-156.5]+137. ; phases: only rough values (I took the maximum or minimum visible on the plot)
pl=[-140.5,-135.5,-137.,-135.5,-140.5]+137.
pd=[-163.5,-144.75,-137.,-144.75,-163.5]+137.
pu=[-138.75,-134.,-137.,-134.,-138.75]+137.

;res=poly_fit(td,pd,2) ; stronger angular dependence is for "down" positions
res=FLTARR(3)
res[0]=0.d0 ; passes through (0,0)
res[1]=0.d0 ; symmetrical
res[2]=-35.45

res2=-6.366 ; for pl by ignoring the 1st point...

x=findgen(1000)/999.

!P.MULTI=0
DEVICE,file='angdep.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color
LOADCT,3
PLOT,x,-(res[0]+res[1]*x+res[2]*x^2.d0)*693.5d0/360.d0,tit='!17',xtit='Angle of incidence (degrees)',ytit='Wavelength change (mA)',charsize=1.5,thick=2,yst=1,yrange=[0,64]
OPLOT,td,-pd*693.5d0/360.d0,psym=4,thick=3
OPLOT,x,-res2*x^2.d0*693.5d0/360.d0,thick=3,col=180
OPLOT,[tl[4],tl[4]],-[pl[4],pl[4]]*693.5d0/360.d0,psym=4,thick=3,col=180
DEVICE,/CLOSE


; PHASE AND CONTRAST DIFFERENCE BETWEEN OBS AND CAL
;------------------------------------------------------------------------------------------------------------------

nx2=256

RESTORE,'CPT/CPT_laser_front_710600.BIN' ; OBSMODE
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.86541d0+113.287d0,6.0024801d0-5.26063d0,-138.18929d0+138.755d0]*!dpi/180.d0
p1=PHIG0
B1=BG0
RESTORE,'CPT/CPT_laser_side_710600.BIN' ; OBSMODE'
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-113.46681d0+113.287d0,5.9007082d0-5.26063d0,-138.12755d0+138.755d0]
p2=PHIG0
B2=BG0
RESTORE,'CPT/CPT_laser_front_712408.BIN' ; OBSMODE
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.84284d0+113.287d0,5.8180826d0-5.26063d0,-138.89860d0+138.755d0]
p3=PHIG0
B3=BG0
RESTORE,'CPT/CPT_laser_side_712408.BIN' ; OBSMODE
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-113.90070d0+113.287d0,6.0767232d0-5.26063d0,-138.86148d0+138.755d0]
p4=PHIG0
B4=BG0

phigOBS=(p1+p2+p3+p4)/4.d0
bgOBS  =(b1+b2+b3+b4)/4.d0

RESTORE,'CPT/CPT_laser_front_710660.BIN' ; CALMODE (980 arcsec)
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.42624d0+112.853d0,7.0971522d0-6.20089d0,-135.83704d0+135.252d0]
p1=PHIG0
B1=BG0
RESTORE,'CPT/CPT_laser_side_710660.BIN' ; CALMODE (980 arcsec)
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-113.18183d0+112.853d0,7.2059481d0-6.20089d0,-135.30350d0+135.252d0]
p2=PHIG0
B2=BG0
RESTORE,'CPT/CPT_laser_front_712468.BIN' ; CALMODE (980 arcsec)
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.93758d0+112.853d0,7.3798297d0-6.20089d0,-134.51804d0+135.252d0]
p3=PHIG0
B3=BG0
RESTORE,'CPT/CPT_laser_side_712468.BIN' ; CALMODE (980 arcsec)
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-113.53256d0+112.853d0,7.1696120d0-6.20089d0,-134.46911d0+135.252d0]
p4=PHIG0
B4=BG0

phigCAL=(p1+p2+p3+p4)/4.d0
bgCAL  =(b1+b2+b3+b4)/4.d0

distance=dist(256,256)*0.5*16.
distance=shift(distance,128,128)
a=where(distance le 955,complement=b)
distance[a]=1.d0
distance[b]=0.d0

device,file='calobs.ps',xoffset=0,yoffset=0,xsize=20,ysize=48,/color,bits=24
!p.multi=[0,1,3]
LOADCT,4
tvim,(phigCAL[*,*,0]-phigOBS[*,*,0])*distance,/scale,range=[-2.5,10.5],tit='!17',xtit='Pixel number',ytit='Pixel number',xrange=[0,4095],yrange=[0,4095],stit='Phase difference (in degrees)',pcharsize=2.5
tvim,(phigCAL[*,*,1]-phigOBS[*,*,1])*distance,/scale,range=[-2.5,10.5],tit='!17',xtit='Pixel number',ytit='Pixel number',xrange=[0,4095],yrange=[0,4095],stit='Phase difference (in degrees)',pcharsize=2.5
tvim,(phigCAL[*,*,2]-phigOBS[*,*,2])*distance,/scale,range=[-2.5,10.5],tit='!17',xtit='Pixel number',ytit='Pixel number',xrange=[0,4095],yrange=[0,4095],stit='Phase difference (in degrees)',pcharsize=2.5
DEVICE,/CLOSE

temp=REFORM(phigCAL[*,*,0]-phigOBS[*,*,0])
print,mean(temp[a])
temp=REFORM(phigCAL[*,*,1]-phigOBS[*,*,1])
print,mean(temp[a])
temp=REFORM(phigCAL[*,*,2]-phigOBS[*,*,2])
print,mean(temp[a])


; PHASE + CONTRAST MAPS IN LASER LIGHT AND OBSMODE AND CALMODE
;------------------------------------------------------------------------------------------------------------------

anglim=970.
distance=dist(4096,4096)*0.5d0
distance=shift(distance,2048,2048)
distance=REBIN(distance,256,256)
a=WHERE(distance GT anglim,COMPLEMENT=b) ; depends on the size of the target ;963 obsmode; 925 calmode
distance[a]=1.d0
distance[b]=0.d0

FOR i=0,2 DO BEGIN 
    temp=REFORM(PhigObs[*,*,i])
    temp[a]=-10000.d0 
    phigObs[*,*,i]=temp[*,*]
    temp=REFORM(BgObs[*,*,i])
    temp[a]=-10000.d0
    BgObs[*,*,i]=temp[*,*]
    phigObs[0:1,*,i]=-10000.d0
    BgObs[0:1,*,i]=-10000.d0
    temp=REFORM(PhigCal[*,*,i])
    temp[a]=-10000.d0 
    phigCal[*,*,i]=temp[*,*]
    temp=REFORM(BgCal[*,*,i])
    temp[a]=-10000.d0
    BgCal[*,*,i]=temp[*,*]
    phigCal[0:1,*,i]=-10000.d0
    BgCal[0:1,*,i]=-10000.d0
ENDFOR

; ADD 180 DEGREES TO PHASES < -180
FOR i=0,2 DO BEGIN
temp = REFORM(PhigObs[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigObs[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigObs[*,*,i]=temp
    ENDELSE
temp = REFORM(PhigCal[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigCal[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigCal[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(3)
moyb=moy
moycal=DBLARR(3)
moybcal=moycal
FOR i=0,2 DO BEGIN & temp=REFORM(PhigObs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy[i]=MEAN(temp[a]) & temp[b]=moy[i] & PhigObs[*,*,i]=temp & temp=REFORM(BgObs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb[i]=MEAN(temp[a]) & temp[b]=moyb[i] & BgObs[*,*,i]=temp  & ENDFOR
FOR i=0,2 DO BEGIN & temp=REFORM(PhigCal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moycal[i]=MEAN(temp[a]) & temp[b]=moycal[i] & PhigCal[*,*,i]=temp & temp=REFORM(BgCal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moybcal[i]=MEAN(temp[a]) & temp[b]=moybcal[i] & BgCal[*,*,i]=temp  & ENDFOR


; WE PLOT THE RESULT
!P.MULTI=[0,2,3]
device,file='phasesObs.ps',bits=24,xoffset=0.,yoffset=0.,xsize=23,ysize=26,/color
LOADCT,4
tvim,phigObs[*,*,0],/scale,tit='!17NB '+STRTRIM(STRING(moy[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,0]),MAX(phigObs[*,*,0])],pcharsize=2.0
tvim,BgObs[*,*,0]  ,/scale,tit='!17NB '+STRTRIM(STRING(moyb[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,0]),MAX(BgObs[*,*,0])],pcharsize=2.0
tvim,phigObs[*,*,1],/scale,tit='!17WB '+STRTRIM(STRING(moy[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,1]),MAX(phigObs[*,*,1])],pcharsize=2.0
tvim,BgObs[*,*,1]  ,/scale,tit='!17WB '+STRTRIM(STRING(moyb[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,1]),MAX(BgObs[*,*,1])],pcharsize=2.0
tvim,phigObs[*,*,2],/scale,tit='!17E1 '+STRTRIM(STRING(moy[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,2]),MAX(phigObs[*,*,2])],pcharsize=2.0
tvim,BgObs[*,*,2]  ,/scale,tit='!17E1 '+STRTRIM(STRING(moyb[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,2]),MAX(BgObs[*,*,2])],pcharsize=2.0

DEVICE,/close

device,file='phasesCal.ps',bits=24,xoffset=0.,yoffset=0.,xsize=23,ysize=26,/color
LOADCT,4
tvim,phigCal[*,*,0],/scale,tit='!17NB '+STRTRIM(STRING(moycal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,0]),MAX(phigCal[*,*,0])],pcharsize=2.0
tvim,BgCal[*,*,0]  ,/scale,tit='!17NB '+STRTRIM(STRING(moybcal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,0]),MAX(BgCal[*,*,0])],pcharsize=2.0
tvim,phigCal[*,*,1],/scale,tit='!17WB '+STRTRIM(STRING(moycal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,1]),MAX(phigCal[*,*,1])],pcharsize=2.0
tvim,BgCal[*,*,1]  ,/scale,tit='!17WB '+STRTRIM(STRING(moybcal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,1]),MAX(BgCal[*,*,1])],pcharsize=2.0
tvim,phigCal[*,*,2],/scale,tit='!17E1 '+STRTRIM(STRING(moycal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,2]),MAX(phigCal[*,*,2])],pcharsize=2.0
tvim,BgCal[*,*,2]  ,/scale,tit='!17E1 '+STRTRIM(STRING(moybcal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,2]),MAX(BgCal[*,*,2])],pcharsize=2.0


DEVICE,/close


; STABILITY : PHASE + CONTRAST MAPS IN LASER WHEN OVEN TEMPERATURE VARIES
;------------------------------------------------------------------------------------------------------------------
nx2=256
RESTORE,'CPT/CPT_laser_front_june09_90306.BIN' ; OBSMODE TS12=29.86
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.86541d0+113.287d0,6.0024801d0-5.26063d0,-138.18929d0+138.755d0]*!dpi/180.d0
PhigObs=Phig0
BgObs=Bg0
RESTORE,'CPT/CPT_laser_front_june09_90636.BIN' ; OBSMODE TS12=31.305
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.86541d0+113.287d0,6.0024801d0-5.26063d0,-138.18929d0+138.755d0]*!dpi/180.d0
PhigObs1=Phig0
BgObs1=Bg0
RESTORE,'CPT/CPT_laser_front_june09_90846.BIN' ; OBSMODE TS12=32.04
FOR iii=0,nx2-1 DO FOR jjj=0,nx2-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]*180.d0/!dpi;-[-112.86541d0+113.287d0,6.0024801d0-5.26063d0,-138.18929d0+138.755d0]*!dpi/180.d0
PhigCal=Phig0 ; even though both detunes are in Obsmode
BgCal=Bg0

anglim=963.
distance=dist(4096,4096)*0.5d0
distance=shift(distance,2048,2048)
distance=REBIN(distance,256,256)
a=WHERE(distance GT anglim,COMPLEMENT=b) ; depends on the size of the target ;963 obsmode; 925 calmode
distance[a]=1.d0
distance[b]=0.d0

FOR i=0,2 DO BEGIN 
    temp=REFORM(PhigObs[*,*,i])
    temp[a]=-10000.d0 
    phigObs[*,*,i]=temp[*,*]
    temp=REFORM(BgObs[*,*,i])
    temp[a]=-10000.d0
    BgObs[*,*,i]=temp[*,*]
    phigObs[0:1,*,i]=-10000.d0
    BgObs[0:1,*,i]=-10000.d0
    temp=REFORM(PhigCal[*,*,i])
    temp[a]=-10000.d0 
    phigCal[*,*,i]=temp[*,*]
    temp=REFORM(BgCal[*,*,i])
    temp[a]=-10000.d0
    BgCal[*,*,i]=temp[*,*]
    phigCal[0:1,*,i]=-10000.d0
    BgCal[0:1,*,i]=-10000.d0
    temp=REFORM(PhigObs1[*,*,i])
    temp[a]=-10000.d0 
    phigObs1[*,*,i]=temp[*,*]
    temp=REFORM(BgObs1[*,*,i])
    temp[a]=-10000.d0
    BgObs1[*,*,i]=temp[*,*]
    phigObs1[0:1,*,i]=-10000.d0
    BgObs1[0:1,*,i]=-10000.d0
ENDFOR

; ADD 180 DEGREES TO PHASES < -180
FOR i=0,2 DO BEGIN
temp = REFORM(PhigObs[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigObs[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigObs[*,*,i]=temp
    ENDELSE
temp = REFORM(PhigObs1[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigObs1[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigObs1[*,*,i]=temp
    ENDELSE
temp = REFORM(PhigCal[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        PhigCal[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        PhigCal[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(3)
moyb=moy
moycal=DBLARR(3)
moybcal=moycal
moy1=DBLARR(3)
moyb1=moy1
FOR i=0,2 DO BEGIN & temp=REFORM(PhigObs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy[i]=MEAN(temp[a]) & temp[b]=moy[i] & PhigObs[*,*,i]=temp & temp=REFORM(BgObs[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb[i]=MEAN(temp[a]) & temp[b]=moyb[i] & BgObs[*,*,i]=temp  & ENDFOR
FOR i=0,2 DO BEGIN & temp=REFORM(PhigCal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moycal[i]=MEAN(temp[a]) & temp[b]=moycal[i] & PhigCal[*,*,i]=temp & temp=REFORM(BgCal[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moybcal[i]=MEAN(temp[a]) & temp[b]=moybcal[i] & BgCal[*,*,i]=temp  & ENDFOR
FOR i=0,2 DO BEGIN & temp=REFORM(PhigObs1[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy1[i]=MEAN(temp[a]) & temp[b]=moy1[i] & PhigObs1[*,*,i]=temp & temp=REFORM(BgObs1[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb1[i]=MEAN(temp[a]) & temp[b]=moyb1[i] & BgObs1[*,*,i]=temp  & ENDFOR

; WE PLOT THE RESULT
!P.MULTI=[0,3,3]
device,file='stability.ps',bits=24,xoffset=0.,yoffset=0.,xsize=31,ysize=26,/color
LOADCT,4

tvim,phigObs[*,*,0],/scale,tit='!17NB '+STRTRIM(STRING(moy[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,0]),MAX(phigObs[*,*,0])],pcharsize=2.0
tvim,phigObs1[*,*,0],/scale,tit='!17NB '+STRTRIM(STRING(moy1[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs1[*,*,0]),MAX(phigObs1[*,*,0])],pcharsize=2.0
tvim,phigCal[*,*,0],/scale,tit='!17NB '+STRTRIM(STRING(moycal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,0]),MAX(phigCal[*,*,0])],pcharsize=2.0


tvim,phigObs[*,*,1],/scale,tit='!17WB '+STRTRIM(STRING(moy[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,1]),MAX(phigObs[*,*,1])],pcharsize=2.0
tvim,phigObs1[*,*,1],/scale,tit='!17WB '+STRTRIM(STRING(moy1[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs1[*,*,1]),MAX(phigObs1[*,*,1])],pcharsize=2.0

tvim,phigCal[*,*,1],/scale,tit='!17WB '+STRTRIM(STRING(moycal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,1]),MAX(phigCal[*,*,1])],pcharsize=2.0

tvim,phigObs[*,*,2],/scale,tit='!17E1 '+STRTRIM(STRING(moy[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs[*,*,2]),MAX(phigObs[*,*,2])],pcharsize=2.0
tvim,phigObs1[*,*,2],/scale,tit='!17E1 '+STRTRIM(STRING(moy1[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigObs1[*,*,2]),MAX(phigObs1[*,*,2])],pcharsize=2.0
tvim,phigCal[*,*,2],/scale,tit='!17E1 '+STRTRIM(STRING(moycal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phigCal[*,*,2]),MAX(phigCal[*,*,2])],pcharsize=2.0

;tvim,BgObs[*,*,0]  ,/scale,tit='!17NB '+STRTRIM(STRING(moyb[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,0]),MAX(BgObs[*,*,0])],pcharsize=2.0
;tvim,BgCal[*,*,0]  ,/scale,tit='!17NB '+STRTRIM(STRING(moybcal[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,0]),MAX(BgCal[*,*,0])],pcharsize=2.0
;tvim,BgObs[*,*,1]  ,/scale,tit='!17WB '+STRTRIM(STRING(moyb[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,1]),MAX(BgObs[*,*,1])],pcharsize=2.0
;tvim,BgCal[*,*,1]  ,/scale,tit='!17WB '+STRTRIM(STRING(moybcal[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,1]),MAX(BgCal[*,*,1])],pcharsize=2.0
;tvim,BgObs[*,*,2]  ,/scale,tit='!17E1 '+STRTRIM(STRING(moyb[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgObs[*,*,2]),MAX(BgObs[*,*,2])],pcharsize=2.0
;tvim,BgCal[*,*,2]  ,/scale,tit='!17E1 '+STRTRIM(STRING(moybcal[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(BgCal[*,*,2]),MAX(BgCal[*,*,2])],pcharsize=2.0


DEVICE,/close





; PHASE + CONTRAST MAPS IN SUNLIGHT AND CALMODE
;------------------------------------------------------------------------------------------------------------------

nx=256
ny=256
RESTORE,'RESULTS/RESULTS_listSun070905_553222_256.BIN'

a = WHERE(lateconv EQ 1,COMPLEMENT=aa)
IF(a[0] EQ -1) THEN BEGIN
    PRINT,'NO EARLY CONVERGENCES: PROBLEM'
    STOP
ENDIF

PRINT,'AVERAGED PHASES:'
tab = DBLARR(nx,ny)
moy = DBLARR(3)
dis = moy
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR

anglim=970.
distance=dist(256,256)*0.5*16.
distance=shift(distance,128,128)
a=WHERE(distance LE anglim,COMPLEMENT=b)
distance[a]=1.d0
IF(b[0] NE -1) THEN distance[b]=0.d0

; TO SUPRESS THE POTENTIAL STRIPES AT THE LEFT AND UPPER EDGE
;distance[0:7,*]=0.0
;distance[*,120:nx-1]=0.0
;distance[*,0:2]=0.0
;distance[251:nx-1,*]=0.0

FOR i=0,2 DO Phig[*,*,i]=Phig[*,*,i]*distance
TVIM,distance

PRINT,'CORRECTED AVERAGED PHASES'
a = WHERE(Phig[*,*,0] NE 0.d0 AND lateconv EQ 1,COMPLEMENT=aa)
FOR i=0,2 DO BEGIN & tab = Phig[*,*,i] & moy[i]=MEAN(tab[a]) & dis[i]=SIGMA(tab[a]) & PRINT,moy[i],dis[i] & ENDFOR
    
!P.MULTI=[0,2,3]
DEVICE,file='phaseSun.ps',/color,bits=24,xoffset=0,yoffset=0,xsize=23,ysize=26
titl=STRARR(3)
titl[0]='NB '
titl[1]='WB '
titl[2]='E1 '

LOADCT,4
TVIM,Phig[*,*,0],/scale,range=[moy[0]-2.5*dis[0],moy[0]+2.5*dis[0]],tit='!17'+titl[0]+STRTRIM(STRING(moy[0]),1),stit='!17Relative phase (in degrees)',barwidth=0.5,xtit='Pixels',ytit='Pixels',pcharsize=2.0
LOADCT,3
TVIM,linewidth*distance,/scale,stit='!17Linewidth (A)',barwidth=0.5,range=[MIN(linewidth[a]),MAX(linewidth[a])],xtit='Pixels',ytit='Pixels',pcharsize=2.0

LOADCT,4
TVIM,Phig[*,*,1],/scale,range=[moy[1]-2.5*dis[1],moy[1]+2.5*dis[1]],tit='!17'+titl[1]+STRTRIM(STRING(moy[1]),1),stit='!17Relative phase (in degrees)',barwidth=0.5,xtit='Pixels',ytit='Pixels',pcharsize=2.0
LOADCT,3
TVIM,(linedepth/continuum)*distance,/scale,stit='!17Linedepth (a.u.)',barwidth=0.5,range=[MIN(linedepth[a]/continuum[a]),MAX(linedepth[a]/continuum[a])],xtit='Pixels',ytit='Pixels',pcharsize=2.0


LOADCT,4
TVIM,Phig[*,*,2],/scale,range=[moy[2]-2.5*dis[2],moy[2]+2.5*dis[2]],tit='!17'+titl[2]+STRTRIM(STRING(moy[2]),1),stit='!17Relative phase (in degrees)',barwidth=0.5,xtit='Pixels',ytit='Pixels',pcharsize=2.0
LOADCT,3
TVIM,continuum*distance,/scale,stit='!17Continuum intensity (a.u.)',barwidth=0.5,range=[MIN(continuum[a]),MAX(continuum[a])],xtit='Pixels',ytit='Pixels',pcharsize=2.0


DEVICE,/close


; RECONSTRUCTION OF DETUNE SEQUENCE
;------------------------------------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS_listSun070905_553222_reconstruction.BIN'
a=WHERE(Inten EQ -1.0)
Inten[a]=0.0
nseq=27
!P.MULTI=[0,1,2]
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='reconstruction.ps',xsize=20,ysize=26,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,FINDGEN(nseq)+1,REBIN(Inten,1,1,nseq),xst=1,tit='!17',xtit='Position number',ytit='Intensity (a.u.)',charsize=1.5,thick=3,yrange=[0.15,0.3],yst=1;,psym=10
OPLOT,FINDGEN(nseq)+1,REBIN(Inten2,1,1,nseq),linestyle=2,color=180,thick=3;,psym=10
PLOT,FINDGEN(nseq)+1,(REBIN(Inten2,1,1,nseq)-REBIN(Inten,1,1,nseq))/REBIN(Inten,1,1,nseq),xst=1,tit='!17',xtit='Position number',ytit='Relative Difference',charsize=1.5,thick=3,yst=1;,psym=10
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten,1,1,nseq))-REFORM(REBIN(Inten2,1,1,nseq)))^2.0 )/SQRT(TOTAL(REBIN(Inten2,1,1,nseq)^2.0))
device,/close
!P.MULTI=0


; FRINGES IN WHITE LIGHT
;------------------------------------------------------------------------------------------------------------------

nx=256
ny=256
RESTORE,'RESULTS/RESULTS_LAMP_561378_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN' ; CALMODE
proxy1Cal=proxy1
proxy2Cal=proxy2
RESTORE,'RESULTS/RESULTS_LAMP_560718_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN' ; OBSMODE
proxy1Obs=proxy1
proxy2Obs=proxy2

center   = [128,128]  ; !!!WARNING, depends on the images
distance = SHIFT(DIST(nx,ny),center[0],center[1])*0.5d0*4096.d0/nx ; distance in arcseconds from the image center
anglim= 950.
a = WHERE(distance LT anglim, COMPLEMENT=b)
distance[a] = 1.d0
distance[b] = 0.d0

mini=MIN([proxy1Obs[a],proxy1Cal[a],proxy2Obs[a],proxy2Cal[a]])
maxi=MAX([proxy1Obs[a],proxy1Cal[a],proxy2Obs[a],proxy2Cal[a]])


SET_PLOT,'ps'
!p.multi=[0,2,4]
loadct,4
DEVICE,file='fringesObs.ps',xoffset=0,yoffset=0,xsize=20,ysize=27,/color,bits=24
FOR i=0,3 DO BEGIN
    TVIM,proxy1Obs[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(i+1)+'/(4. FSR!DNB!N)',pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
    TVIM,proxy2Obs[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(i+1)+'/(4. FSR!DNB!N)',pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
ENDFOR
DEVICE,/close
!p.multi=[0,2,4]
DEVICE,file='fringesCal.ps',xoffset=0,yoffset=0,xsize=20,ysize=27,/color,bits=24
FOR i=0,3 DO BEGIN
    TVIM,proxy1Cal[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(i+1)+'/(4. FSR!DNB!N)',pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
    TVIM,proxy2Cal[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(i+1)+'/(4. FSR!DNB!N)',pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
ENDFOR
DEVICE,/close

RESTORE,'RESULTS/RESULTS_LAMP2_561378_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
mini=MIN([proxy1[a],proxy2[a]])
maxi=MAX([proxy1[a],proxy2[a]])
lamref = 6173.3433d0
freqs  = ([.003d0,.006d0,.009d0]*1.d10*2.d0*1.516d0+lamref)/lamref^2.d0 ; in A^{-1}

!p.multi=[0,2,3]
loadct,4
DEVICE,file='fringesFW.ps',xoffset=0,yoffset=0,xsize=20,ysize=23,/color,bits=24
FOR i=0,2 DO BEGIN
    TVIM,proxy1[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(freqs[i]),pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
    TVIM,proxy2[*,*,i]*distance,range=[mini,maxi],/scale,tit='!17'+STRING(freqs[i]),pcharsize=2.0,xtit='Pixels',ytit='Pixels',barwidth=0.5,stit='Fourier coefficient'
ENDFOR
DEVICE,/close





; THROUGHTPUT
;------------------------------------------------------------------------------------------------------------------


  

; TEMPERATURE DEPENDENCE OF TUNABLE ELEMENTS
;------------------------------------------------------------------------------------------------------------------


Toven=[27.65,29.62,34.1] ;HMI_TS12_OVN_LYOT & HMI_TS13_OVN_WBM
Tmich=[27.19,29.31,33.83];HMI_TS14_OVN_NBM

; FRONT CAMERA
NBobs1=[-188.08,-79.24,-116.47]
WBobs1=[-16.34,29.20,-3.79]
E1obs1=[12.92,35.68,33.22]

NBcal1=[-191.35,-82.08,-118.65]
WBcal1=[-15.58,30.32,-5.03]
E1cal1=[17.89,41.54,35.55]

NBobs1=-(NBobs1-NBobs1[1])*171.0/360.
WBobs1=-(WBobs1-WBobs1[1])*342.0/360.
E1obs1=-(E1obs1-E1obs1[1])*693.5/360.

NBcal1=-(NBcal1-NBcal1[1])*171.0/360.
WBcal1=-(WBcal1-WBcal1[1])*342.0/360.
E1cal1=-(E1cal1-E1cal1[1])*693.5/360.

;SIDE CAMERA
NBobs2=[-188.85,-80.55,-118.11]
WBobs2=[-16.55,28.89,-4.65]
E1obs2=[12.93,35.24,33.09]

NBcal2=[-191.82,-83.24,-118.39]
WBcal2=[-16.02,29.08,-5.39]
E1cal2=[17.83,41.60,36.50]

NBobs2=-(NBobs2-NBobs2[1])*171.0/360.
WBobs2=-(WBobs2-WBobs2[1])*342.0/360.
E1obs2=-(E1obs2-E1obs2[1])*693.5/360.
                    
NBcal2=-(NBcal2-NBcal2[1])*171.0/360.
WBcal2=-(WBcal2-WBcal2[1])*342.0/360.
E1cal2=-(E1cal2-E1cal2[1])*693.5/360.

NBobs1=(NBobs1+NBobs2)/2.
WBobs1=(WBobs1+WBobs2)/2.
E1obs1=(E1obs1+E1obs2)/2.

NBcal1=(NBcal1+NBcal2)/2.
WBcal1=(WBcal1+WBcal2)/2.
E1cal1=(E1cal1+E1cal2)/2.


DEVICE,file='tempe.ps',xoffset=0,yoffset=0,xsize=20,ysize=31,/color,bits=24
LOADCT,3
!P.MULTI=[0,1,3]

tempe=FINDGEN(1000)/999.*8.+27.

PLOT,Tmich,NBobs1,psym=4,thick=2,yrange=[-20,60],yst=1,xrange=[27,35],xst=1,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength change (mA)',charsize=1.5
res=poly_fit(Tmich,NBobs1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],thick=2
OPLOT,Tmich,NBcal1,psym=4,thick=2,col=180
res=poly_fit(Tmich,NBcal1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],col=180,thick=2

PLOT,Tmich,WBobs1,psym=4,thick=2,yrange=[-20,60],yst=1,xrange=[27,35],xst=1,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength change (mA)',charsize=1.5
res=poly_fit(Tmich,WBobs1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],thick=2
OPLOT,Tmich,WBcal1,psym=4,thick=2,col=180
res=poly_fit(Tmich,WBcal1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],col=180,thick=2

PLOT,Tmich,E1obs1,psym=4,thick=2,yrange=[-20,60],yst=1,xrange=[27,35],xst=1,tit='!17',xtit='Temperature (degrees Celsius)',ytit='Wavelength change (mA)',charsize=1.5
res=poly_fit(Tmich,E1obs1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],thick=2
OPLOT,Tmich,E1cal1,psym=4,thick=2,col=180
res=poly_fit(Tmich,E1cal1,2)
OPLOT,tempe,res[0]+tempe*res[1]+tempe^2.0*res[2],col=180,thick=2
DEVICE,/CLOSE


; I-RIPPLE (SPATIALLY AVERAGED)
;------------------------------------------------------------------------------------------------------------------

;CF=9
;--------------------------------------------------
RESTORE,'CPT/SEQUENCE_listSun070905_560718_256.BIN'

distance=dist(4096,4096)*0.5d0
distance=shift(distance,2048,2048)
distance=REBIN(distance,256,256)
a=where(distance le 980.,na,complement=b) ;980" was the anglim used to produce the series with wavelength_dependence_test_cotunes.pro
front=FLTARR(27)
for i=0,26 do front[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
; I-RIPPLE = 0.0203095
side=FLTARR(27)
for i=0,26 do side[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
; I-RIPPLE = 0.0203895

!P.MULTI=0
DEVICE,file='Iripple.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
LOADCT,3
PLOT,FINDGEN(27)+1,front/MEAN(front),tit='!17',xtit='Position number',ytit='Normalized intensity',xst=1,charsize=1.5,thick=2,psym=10,yrange=[0.988,1.012],yst=1

;CF=17
;-------------------------------------------------
RESTORE,'CPT/SEQUENCE_listSun070905_561378_256.BIN'

front=FLTARR(27)
for i=0,26 do front[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
; I-RIPPLE = 0.0207153
side=FLTARR(27)
for i=0,26 do side[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
; I-RIPPLE =0.0208400


OPLOT,FINDGEN(27)+1,front/MEAN(front),thick=2,col=180,psym=10
DEVICE,/CLOSE

; I-RIPPLE (VARIATION ACROSS APERTURE)
;------------------------------------------------------------------------------------------------------------------

;RESTORE,'CPT/SEQUENCE_listSun070905_560718_256.BIN'
RESTORE,'CPT/SEQUENCE2_listLamp070903_539346_256.BIN'


x0=[50,100,150,200]
y0=[50,100,150,200]

front1=FLTARR(27)
for i=0,26 do front1[i]=REFORM(IMF[x0[0],y0[0],2+i])
front2=FLTARR(27)
for i=0,26 do front2[i]=REFORM(IMF[x0[1],y0[1],2+i])
front3=FLTARR(27)
for i=0,26 do front3[i]=REFORM(IMF[x0[2],y0[2],2+i])
front4=FLTARR(27)
for i=0,26 do front4[i]=REFORM(IMF[x0[3],y0[3],2+i])

!P.MULTI=0
LOADCT,4
DEVICE,file='spatialIripple.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
PLOT,FINDGEN(27)+1 ,front1/MEAN(front1),tit='!17',xtit='Position number',ytit='Normalized intensity',xst=1,charsize=1.5,thick=3,psym=10,yrange=[0.988,1.012],yst=1
OPLOT,FINDGEN(27)+1,front2/MEAN(front2),thick=3,psym=10,col=80
OPLOT,FINDGEN(27)+1,front3/MEAN(front3),thick=3,psym=10,col=160
OPLOT,FINDGEN(27)+1,front4/MEAN(front4),thick=3,psym=10,col=220
DEVICE,/CLOSE


; I-RIPPLE AND STABILITY (CHANGE IN THE FRONT WINDOW TEMPERATURE)
;------------------------------------------------------------------------------------------------------------------

RESTORE,'CPT/SEQUENCE2_listLamp070903_546410_256.BIN'
front1=FLTARR(27)
for i=0,26 do front1[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_546638_256.BIN'
front2=FLTARR(27)
for i=0,26 do front2[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_546986_256.BIN'
front3=FLTARR(27)
for i=0,26 do front3[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)

!P.MULTI=0
LOADCT,4
DEVICE,file='temperatureIripple.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
PLOT,FINDGEN(27)+1 ,front1/MEAN(front1),tit='!17',xtit='Position number',ytit='Normalized intensity',xst=1,charsize=1.5,thick=3,psym=10,yrange=[0.988,1.012],yst=1
OPLOT,FINDGEN(27)+1,front2/MEAN(front2),thick=3,psym=10,col=100
OPLOT,FINDGEN(27)+1,front3/MEAN(front3),thick=3,psym=10,col=200
DEVICE,/CLOSE


; ARTIFACT CHECK AND I-RIPPLE
;------------------------------------------------------------------------------------------------------------------

RESTORE,'CPT/SEQUENCE2_listLamp070903_539006_256.BIN'
distance=dist(4096,4096)*0.5d0
distance=shift(distance,2048,2048)
distance=REBIN(distance,256,256)
a=where(distance le 980.,na,complement=b) ;980" was the anglim used to produce the series with wavelength_dependence_test_cotunes.pro
front1=FLTARR(27)
for i=0,26 do front1[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side1=FLTARR(27)
for i=0,26 do side1[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539074_256.BIN'
front2=FLTARR(27)
for i=0,26 do front2[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side2=FLTARR(27)
for i=0,26 do side2[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539142_256.BIN'
front3=FLTARR(27)
for i=0,26 do front3[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side3=FLTARR(27)
for i=0,26 do side3[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539210_256.BIN'
front4=FLTARR(27)
for i=0,26 do front4[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side4=FLTARR(27)
for i=0,26 do side4[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539278_256.BIN'
front5=FLTARR(27)
for i=0,26 do front5[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side5=FLTARR(27)
for i=0,26 do side5[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539346_256.BIN'
front6=FLTARR(27)
for i=0,26 do front6[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side6=FLTARR(27)
for i=0,26 do side6[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)
RESTORE,'CPT/SEQUENCE2_listLamp070903_539414_256.BIN'
front7=FLTARR(27)
for i=0,26 do front7[i]=TOTAL(IMF[*,*,2+i])/FLOAT(na)
side7=FLTARR(27)
for i=0,26 do side7[i]=TOTAL(IMS[*,*,2+i])/FLOAT(na)

!P.MULTI=0
DEVICE,file='artefact.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
LOADCT,3
PLOT,FINDGEN(27)+1 ,((front2-front1)/front1-MEAN((front2-front1)/front1))*1.d4,tit='!17',xtit='Position number',ytit='Relative difference (x 10!U-4!N)',xst=1,charsize=1.5,thick=2,psym=10,yst=1,yrange=[-1.5,1.5]
;OPLOT,FINDGEN(27)+1,(front2-front1)/front1,thick=2,psym=10
OPLOT,FINDGEN(27)+1,((front3-front1)/front1-MEAN((front3-front1)/front1))*1.d4,thick=2,psym=10
OPLOT,FINDGEN(27)+1,((front4-front1)/front1-MEAN((front4-front1)/front1))*1.d4,thick=2,psym=10
OPLOT,FINDGEN(27)+1,((front5-front1)/front1-MEAN((front5-front1)/front1))*1.d4,thick=2,psym=10
OPLOT,FINDGEN(27)+1,((front6-front1)/front1-MEAN((front6-front1)/front1))*1.d4,thick=2,psym=10
;OPLOT,FINDGEN(27)+1,(front7/MEAN(front7),thick=2,psym=10
DEVICE,/CLOSE

END
