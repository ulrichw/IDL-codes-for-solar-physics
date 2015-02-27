; from the data provided by Juan Borrero from IBIS data
; TO SEE IF THE HMI LINE PROFILES ARE LINEAR IN RV (RADIUS VECTOR)

PRO gaussianfit,X,A,F,pder

F = 1.d0-A[0]*EXP(-(X-A[1])^2.d0/A[2]^2.d0)

pder= [ [-EXP(-(X-A[1])^2.d0/A[2]^2.d0)],[-2.d0*(X-A[1])/A[2]^2.d0*A[0]*EXP(-(X-A[1])^2.d0/A[2]^2.d0)],[-2.d0*(X-A[1])^2.d0/A[2]^3.d0*A[0]*EXP(-(X-A[1])^2.d0/A[2]^2.d0)] ]


END

FUNCTION  myinterpol,V,X,U

nU=N_ELEMENTS(U)
P=FLTARR(nU)
nX=N_ELEMENTS(X)
Xmin=MIN(x)
Xmax=MAX(x)

a=WHERE(U GT Xmax,na)
b=WHERE(U GE Xmin AND U LE Xmax,nb)

IF (b[0] NE -1) THEN BEGIN
    FOR i=0,nb-1 DO BEGIN
        c=WHERE(X LE U[b[i]],nc)
        c2=WHERE(X GE U[b[i]])
        P[b[i]]=(V[c2[0]]-V[c[nc-1]])/(X[c2[0]]-X[c[nc-1]])*(U[b[i]]-X[c[nc-1]])+V[c[nc-1]]
    ENDFOR
ENDIF

IF (a[0] NE -1) THEN BEGIN
    FOR i=0,na-1 DO BEGIN
        P[a[i]]=(V[nX-1]-V[nX-2])/(X[nX-1]-X[nX-2])*(U[a[i]]-X[nX-2])+V[nX-2]
    ENDFOR
ENDIF

RETURN,p

END


PRO HMIlineprofile

RESTORE,'profile_6173_north.sav'
wave_all=wave_all-6173.65d0
dw=FLTARR(14,11)
FOR i=0,10 DO dw[*,i]=wave_all[1:14,i]-wave_all[0:13,i]
dw=MEAN([min(dw),max(dw)])

wmini=MIN(wave_all)
wmaxi=MAX(wave_all)
wave=FINDGEN(15)*dw+wmini
rv=FLOAT(rv)
rv=acos(1.d0-rv);*180./!dpi ; to convert from radius vector to angular distance from disk center rv=1-cos(\mu)
profilei=FLTARR(15,11)


weights=FLTARR(15)+1.d0
linewidth=FLTARR(11)
linedepth=FLTARR(11)
FOR i=0,10 DO BEGIN
    res=POLY_FIT([0,1,2,13,14],profile[[0,1,2,13,14],i],1)
    profile[*,i]=profile[*,i]/(res[0]+res[1]*FINDGEN(15))
    A=[0.4,0.0,.1]
    yfit = CURVEFIT(wave_all[*,i],profile[*,i], weights, A, SIGMA, FUNCTION_NAME='gaussianfit',TOL=1.d-7)  
    profilei[*,i]=INTERPOL(profile[*,i],wave_all[*,i]-A[1],wave) ;+A[1] to center the profile
    linewidth[i]=A[0]
    linedepth[i]=A[2]
    ;PLOT,wave,profilei[*,i],xst=1,yst=1
    ;OPLOT,wave,yfit,col=180
    ;READ,pause
ENDFOR



k=10
rvi=rv[k]
profilerv=FLTARR(15) ; check if it's linear in COS(angle)
FOR i=0,14 DO profilerv[i]=INTERPOL(profilei[i,[0,3,5]],COS(rv[[0,3,5]]),COS(rvi))

; [0,3,5] because those are the angles for which R. Ulrich gave me
; line profiles: 0, 45, and 60 degrees

rvi=1.d0-COS(rv[k])
profilerv2=FLTARR(15) ; check if it's linear in rv
profilerv3=profilerv2
FOR i=0,14 DO profilerv2[i]=INTERPOL(REFORM(profilei[i,[0,3,5]]),1.d0-COS(rv[[0,3,5]]),rvi)
FOR i=0,14 DO profilerv3[i]=MYINTERPOL(REFORM(profilei[i,[0,3,5]]),1.d0-COS(rv[[0,3,5]]),rvi) ; to test the interpolation method


WINDOW,0,RETAIN=2,xsize=900,ysize=900
!P.MULTI=[0,2,2]
plot,wave, profilerv,xst=1,yst=1,yrange=[0,1.05]
OPLOT,wave,profilei[*,k],col=180
plot,wave,(profilerv -profilei[*,k])/profilei[*,k],xst=1,yst=1
plot,wave, profilerv2,xst=1,yst=1,yrange=[0,1.05]
OPLOT,wave,profilerv3,thick=2,linestyle=2
OPLOT,wave,profilei[*,k],col=180
plot,wave,(profilerv2-profilei[*,k])/profilei[*,k],xst=1,yst=1

OPENR,1,'Ulrich_Fe_0.txt'
d0=FLTARR(2,98)
READF,1,d0
CLOSE,1
OPENR,1,'Ulrich_Fe_45.txt'
d45=FLTARR(2,98)
READF,1,d45
CLOSE,1
OPENR,1,'Ulrich_Fe_60.txt'
d60=FLTARR(2,98)
READF,1,d60
CLOSE,1

linewidth2=FLTARR(3)
linedepth2=FLTARR(3)
weights=FLTARR(98)+1.d0
yfit = CURVEFIT(REFORM(d0[0,*]),REFORM(d0[1,*]), weights, A, SIGMA, FUNCTION_NAME='gaussianfit',TOL=1.d-7)  
linewidth2[0]=A[0]
linedepth2[0]=A[2]
yfit = CURVEFIT(REFORM(d45[0,*]),REFORM(d45[1,*]), weights, A, SIGMA, FUNCTION_NAME='gaussianfit',TOL=1.d-7)  
linewidth2[1]=A[0]
linedepth2[1]=A[2]
yfit = CURVEFIT(REFORM(d60[0,*]),REFORM(d60[1,*]), weights, A, SIGMA, FUNCTION_NAME='gaussianfit',TOL=1.d-7)  
linewidth2[2]=A[0]
linedepth2[2]=A[2]




READ,pause


END
