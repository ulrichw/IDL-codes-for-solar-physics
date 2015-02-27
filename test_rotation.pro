; PROGRAM TO TEST MY EQUATIONS FOR THE PROJECTION
; OF SOLAR ROTATION VELOCITY ONTO THE HMI CCD

; I RUN RICHARD'S CODE (compile before with: .r solar_rotation)
B0=7.d0 ; B ANGLE IN DEGREES
B0=B0*!dpi/180.d0
P =0.d0 ; P ANGLE IN DEGREES
P = P*!dpi/180.d0
r = 1920.d0 ; solar radius in pixels
time=1.d0
;res=rotation(r,2048,2048,P,B0,time)
centerX=2048.d0
centerY=2048.d0
nx=1024.d0
n=4096./nx
res=solar_rotation(4096,n,r,centerX,centerY,p,b0,time)


x=FLTARR(nx,nx)
y=x
FOR i=0,nx-1 DO x[*,i]=dindgen(nx)*(4096.d0/nx)-centerX ; pixels coordinates (center=(0,0))
FOR j=0,nx-1 DO y[j,*]=dindgen(nx)*(4096.d0/nx)-centerY

; REMOVE THE P-ANGLE FROM THE CARTESIAN FRAME

xprime=x*cos(P)-y*sin(P)
yprime=x*sin(P)+y*cos(P)

; NOW PROJECT ONTO SPHERE

Z = sqrt(r^2.d0-xprime^2.d0-yprime^2.d0)

L = ATAN(xprime,(Z*cos(B0)-Yprime*sin(B0)))
B = ASIN((Yprime*cos(B0)+Z*sin(B0))/r) 

Omega = 452.d0 - 49.d0*sin(B)^2.d0 - 84.d0*sin(B)^4.d0 - 31.7d0 ; calculate rotation rate from richard's code
A=WHERE(FINITE(Omega) Eq 0)
Omega[a]=0.d0
Omega=2.d0*!dpi*Omega*1d-9 ; rotation angular velocity in rad/s

; NOW REMOVE THE B-ANGLE

Vx = cos(L)*cos(B)*Omega*r
Vy = sin(L)*sin(B0)*cos(B)*Omega*r
Vz = -cos(B0)*sin(L)*cos(B)*Omega*6.9894d8
A=WHERE(FINITE(Vx) Eq 0)
Vx[a]=0.d0
A=WHERE(FINITE(Vy) Eq 0)
Vy[a]=0.d0
A=WHERE(FINITE(Vz) Eq 0)
Vz[a]=0.d0

!p.multi=[0,1,2]

TVIM,vz,/scale
TVIM,res[*,*,2],/scale


END
