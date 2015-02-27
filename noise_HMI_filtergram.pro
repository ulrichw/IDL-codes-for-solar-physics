; this program calculates the rms variation in DNs on the HMI
; filtergramas due to photon noise, as a function of the Doppler
; velocity

PRO noise_HMI_filtergram


; READING THE FILTER PROFILES

nx2=128
nlam=12001
dlam          = 0.001d0 ; sampling rate in Angstroms
;filter=fltarr(nlam,6,nx2,nx2)
;openr,1,'HMIfilter.128.out'
;readu,1,filter
;close,1
RESTORE,'filter.bin'
lam=(findgen(nlam)-(nlam-1.)/2.)*dlam

;plot,lam,filter[*,0,nx2/2,nx2/2] ; central pixel profiles

; CREATING THE Fe I line PROFILE
a=0.03d0
dlamdv = 2.059205672212074294d-5
width  = 0.06415d0
depth  = 0.56250d0

nvel=501
dvel=26.d0
velocity=(DINDGEN(nvel)-(nvel-1.d0)/2.d0)*dvel

inten=DBLARR(6,nvel)
line=DBLARR(nlam)

FOR i=0,nvel-1 DO BEGIN
    PRINT,i
    lam0 = dlamdv*velocity[i]
    l    = (lam-lam0)/width
    aaa=WHERE(ABS(l) LE 26.5,COMPLEMENT=bbb)
    gaussian= 0.015d0*exp(-((lam-lam0)+0.225d0)*((lam-lam0)+0.225d0)/0.2d0/0.2d0)-0.004d0*exp(-((lam-lam0)-0.150d0)*((lam-lam0)-0.150d0)/0.22d0/0.22d0) ;GAUSSIAN LINE ADDED TO VOIGT PROFILE TO SIMULATE ASYMMETRY
    line[aaa] = 1.d0-depth*exp(-l[aaa]*l[aaa])*(1.d0-a*2.d0/sqrt(!dpi)*(1.d0/2.d0/l[aaa]/l[aaa])*((4.d0*l[aaa]*l[aaa]+3.d0)*(l[aaa]*l[aaa]+1.d0)*exp(-l[aaa]*l[aaa])-1.d0/l[aaa]/l[aaa]*(2.d0*l[aaa]*l[aaa]+3.d0)*sinh(l[aaa]*l[aaa])))-gaussian[aaa]
    line[bbb]=1.d0-gaussian[bbb]    

    FOR k=0,5 DO inten[k,i]=INT_TABULATED(lam,line*REFORM(filter[*,k]),/DOUBLE)

ENDFOR

inten=inten/MAX(inten)*166000.d0  ; for an exposure time corresponding to 166000.d0 photo-electrons ( a camera gain of 16.0, and an exposure time of 0.12992 sec)

SET_PLOT,'PS'
DEVICE,FILE='noise.ps',xoffset=0,yoffset=0,xsize=33,ysize=15,/color,bits=24
LOADCT,4
!P.MULTI=[0,2,1]
PLOT,velocity ,inten[0,*]/MAX(inten),xst=1,yrange=[0.6,1],yst=1,charsize=1.5,tit='!17',xtit='Velocity (m/s)',ytit='intensity/continuum',thick=3
OPLOT,velocity,inten[1,*]/MAX(inten),color=40,thick=3
OPLOT,velocity,inten[2,*]/MAX(inten),color=80,thick=3
OPLOT,velocity,inten[3,*]/MAX(inten),color=120,thick=3
OPLOT,velocity,inten[4,*]/MAX(inten),color=160,thick=3
OPLOT,velocity,inten[5,*]/MAX(inten),color=200,thick=3

PLOT,velocity,1.d0/SQRT(inten[0,*]),xst=1,yrange=[0.0024,0.0032],yst=1,charsize=1.5,tit='!17',xtit='Velocity (m/s)',ytit='photon noise/intensity',thick=3
OPLOT,velocity,1.d0/SQRT(inten[1,*]),color=40,thick=3
OPLOT,velocity,1.d0/SQRT(inten[2,*]),color=80,thick=3
OPLOT,velocity,1.d0/SQRT(inten[3,*]),color=120,thick=3
OPLOT,velocity,1.d0/SQRT(inten[4,*]),color=160,thick=3
OPLOT,velocity,1.d0/SQRT(inten[5,*]),color=200,thick=3


!P.MULTI=0
DEVICE,/CLOSE

READ,pause


END
