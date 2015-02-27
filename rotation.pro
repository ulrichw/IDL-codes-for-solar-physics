; THIS PROGRAM COMPUTES THE RELATIVE IMAGE ROTATION ANGLE
; FOR THE LOCKHEED-MARTIN HELIOSTAT
; BASED ON THE CODE PROVIDED BY R. SHINE
; e-mail 01/05/2006



FUNCTION julian,iy,im,id

iyy = 1900+LONG(iy)
im  = LONG(im)

IF(im LE 2) THEN BEGIN
    iyy = iyy -1
    im  = im + 12
ENDIF

ia = iyy/100
ib = 2-ia+ia/4
ii = LONG(365.25d0*iyy)
ij = LONG(30.6001d0*(im+1))
jd = ii+ij+DOUBLE(ib)+1720994.5d0+DOUBLE(id)

RETURN,jd

END

FUNCTION paloaltosun,iy,im,id,sday,ha,dec,az,el

jd = JULIAN(iy,im,id) ; Julian date

h0 = (jd-2415020.d0)/36525.d0
h  = (jd-2415020.d0+sday/86400.d0)/36525.d0
hh = h*h

ehel= 0.4093198d0-2.2703d-4*h-2.86d-8*hh
eks = 0.01675104d0-0.0000418d0*h

sml = 279.6967d0+36000.769d0*h
sml = sml MOD 360.d0
sml = sml*!dpi/180.d0

anm = 358.4758d0+35999.0498d0*h-0.00015d0*hh
anm = anm*!dpi/180.d0

cc  = (1.91946d0-0.00479d0*h)*SIN(anm)+0.020d0*SIN(2.d0*anm)
cc  = cc*!dpi/180.d0

sl  = sml+cc

; DECLINATION OF THE SUN
dec = ASIN(SIN(ehel)*SIN(sl))*180.d0/!dpi

; RIGHT ASCENSION OF THE SUN
ra  = ATAN ((COS(ehel)*SIN(sl)),COS(sl))*24.d0/2.d0/!dpi
ra  = ra+24.d0*(ra lt 0.d0)

;PRINT,ATAN ((COS(ehel)*SIN(sl)),COS(sl))

sidg0 = 6.6461d0+2400.0513d0*h0
sidg0 = sidg0 MOD 24.d0
sidg  = sday * 1.002737908d0/3600.d0 +  sidg0

lo    = -122.14913d0 ; longitude of Lockheed Martin

sidl  =  sidg + lo*24.d0/360.d0
; HOUR ANGLE OF THE SUN
ha    = sidl-ra

ha = (ha+12.d0) MOD 24.d0 -12.d0
ha = (ha-12.d0) MOD 24.d0 +12.d0

la = 37.407d0*!dpi/180.d0 ; latitude of Lockheed Martin

dr = dec*!dpi/180.d0
hr = ha*(360.d0/24.d0)*!dpi/180.d0

s1 = SIN(la)*SIN(dr)
c1 = COS(la)*COS(dr)*COS(hr)

xq = s1+c1

;xq = MIN([xq,1.d0])
;xq = MAX([xq,-1.d0])

; ALTITUDE OF THE SUN
el = ASIN(xq)

s1 = SIN(dr)-SIN(la)*xq
c1 = COS(la)*SQRT(1.d0-xq*xq)
xq = s1/c1
;xq = MIN([xq,1.d0])
;xq = MAX([xq,-1.d0])

; AZIMUTH OF THE SUN
az = ACOS(xq)
az = az + (ha gt 0.d0)*(2.d0*!dpi-2.d0*az)
az = az-!dpi

;PRINT,dec,ra,ha,el*180.d0/!dpi,az*180.d0/!dpi
;PRINT,'END !'

RETURN,0

END


;--------------------------------------------------------------
;
; MAIN PROGRAM
;
;--------------------------------------------------------------


FUNCTION rotation,year,month,dom,time

; TIME MUST BE A NUMBER OF SECONDS IN UT
; time = heure * 3600 + minutes * 60 + secondes

year0 = year - 1900
month0= month
dom0  = dom
time0 = time

nh    = 10000 ; number of times for which we compute a rotation angle
;sday = time0 + 8.d0*3600.d0 ; conversion PST to UT
;sday  = time0
;sday = FINDGEN(nh)/nh*23.5*3600.+8.5*3600.
sday = FINDGEN(nh)/nh*13.*3600.+14.*3600.

ha    = DBLARR(nh)
dec   = DBLARR(nh)
az    = DBLARR(nh)
el    = DBLARR(nh)
res   = paloaltosun(year0,month0,dom0,sday,ha,dec,az,el)

lat = 37.407d0*!dpi/180.d0 ; latitude of Lockheed Martin

ra1 = ASIN(COS(lat)*SIN(az)/COS(dec*!dpi/180.d0))
ra  = az+ATAN(COS(el),SIN(el))-ra1
ra  = ra*180.d0/!dpi


altitude = 41.62d0       ; secondary mirror 
azimuth  = -6.7d0        ;(180.-6.7)*!pi/180.

;omega = ATAN( -SIN(!pi/2-altitude)*SIN(azimuth),COS(!pi/2.-altitude)*COS(lat)  - SIN(!pi/2.-altitude)*COS(azimuth)*SIN(lat) )
;rho   = ATAN( COS(!pi/2.-altitude)*SIN(lat)+SIN(!pi/2.-altitude)*COS(azimuth)*COS(lat),-SIN(!pi/2.-altitude)*SIN(azimuth)/SIN(omega) )

altaz2hadec,altitude,azimuth,lat*180.d0/!dpi,omega,rho

PRINT,'omega=',omega;*180./!dpi
PRINT,'rho=',rho;*180./!dpi
PRINT,'dec=',dec[0]

rho   = rho  *!dpi/180.d0
omega = omega*!dpi/180.d0
dec   = dec  *!dpi/180.d0
ha    = ha   *!dpi/12.d0

;a = WHERE(ha LT 0.0)
;ha[a] = 24.d0+ha[a]

K = COS(0.5d0*(rho-dec))/COS(0.5d0*(!dpi-rho-dec))
PRINT,'K=',K[0]
rabis = 2.d0*ATAN(K/TAN(0.5d0*(omega-ha)))*180.d0/!dpi
;PRINT,ra*180.d0/!dpi

a = WHERE(rabis LT 0.0 AND sday/3600. GT 20)
rabis[a]=rabis[a]+360.


; Location of the Zenith image
a = -(!pi/2.-dec)
b =  (!pi/2.-el)
c =  (FLTARR(nh)+!pi/2.-lat)
s = 0.5*(a+b+c)
;k = SQRT(SIN(s-a)*SIN(s-b)*SIN(s-c)/SIN(s))
;zenith = 2.*ATAN(k/SIN(s-c))*180./!pi
sz = SQRT(SIN(s-a)*SIN(s-b)/SIN(a)/SIN(b))
cz = SQRT(SIN(s)  *SIN(s-c)/SIN(a)/SIN(b))

a = WHERE(cz EQ MIN(cz))
zenith = 2.*ATAN(sz,cz) *180./!pi
zenith[a[0]:nh-1] = 360.-zenith[a[0]:nh-1]

;WINDOW,0,retain=2
SET_PLOT,'ps'
!P.MULTI=0
DEVICE,FILE='yo.ps'
;plot,sday/3600.d0,ra,tit='!17',xtit='Hour (UT)',ytit='Rotation Angle',xst=1,yst=1
plot,sday/3600.d0-8,el*180./!dpi,tit='!17Elevation',xst=1;,xrange=[17+21./60.-8,22.+48./60.-8]
plot,sday/3600.d0-8,ha*12.d0/!dpi,tit='!17 Angle horaire',xst=1
oplot,[0,32],[0,0],linestyle=2
;plot,sday/3600.d0,az/(!dpi/180.d0),tit='!17Azimuth'
;plot,sday/3600.d0,ra1/(!dpi/180.d0),tit='!17ra1',xst=1
plot,sday/3600.d0-8,rabis,xst=1,tit='!17Rotation Angle',xtit='Hour (PST)',ytit='Angle (degrees)';,xrange=[17+21./60.-8,22.+48./60.-8]
oplot,[0,40],[0,0],linestyle=2
plot,sday/3600.d0-8,DERIV(sday/3600.d0,rabis),xst=1,tit='!17Derivative of the rotation angle',xtit='Hour (PST)',ytit='d angle/d hour'
oplot,[0,32],[15,15],linestyle=2
plot,sday/3600.d0-8,dec*180./!dpi,xst=1,tit='!17 Declination'
plot,sday/3600.d0-8,zenith,tit='!17 Zenith location'
plot,sday/3600.d0-8,cz,xst=1
oplot,sday/3600.d0-8,sz,linestyle=1
plot,sday/3600.d0-8,rabis-zenith,xst=1,tit='!17 Zenith Location on the image',xtit='Time (PST)',ytit='Angle (NCP D Z) (Degrees)'
DEVICE,/CLOSE

READ,pause

RETURN,ra ; IN DEGREES

END
