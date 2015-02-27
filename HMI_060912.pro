PRO sinewave,X,A,F,pder

A[0]=ABS(A[0])

F = A[0] + (A[1]*cos(X/2.d0)+A[2]*sin(X/2.d0))^2.d0

pder= [[FLTARR(N_ELEMENTS(X))+1.d0],[2.d0*A[1]*cos(X/2.d0)^2.d0+2.d0*A[2]*cos(X/2.d0)*sin(X/2.d0)],[2.d0*A[2]*sin(X/2.d0)^2.d0+2.d0*A[1]*cos(X/2.d0)*sin(X/2.d0)] ]


END


PRO HMI_060912
; for data from 10/28/2006, 215730-221953: js_wl_tune_fine


nim=53

filename=[                                     $
 ["/SUM3/D25965082/D9517226/S00000/i_061028_215730.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_215756.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_215822.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_215848.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_215914.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_215941.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220007.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220033.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220059.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220124.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220150.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220216.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220242.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220308.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220336.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220402.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220428.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220454.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220519.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220544.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220609.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220634.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220658.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220723.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220749.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220814.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220839.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220907.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220933.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_220958.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221024.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221050.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221116.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221142.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221209.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221235.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221301.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221327.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221353.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221419.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221446.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221512.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221538.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221604.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221630.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221656.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221722.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221748.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221813.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221838.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221903.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221927.fit"] $
,["/SUM3/D25965082/D9517226/S00000/i_061028_221953.fit"] ]
time=[215730,215756,215822,215848,215914,215941,220007,220033,220059,220124,220150,220216,220242,220308,220336,220402,220428,220454,220519,220544,220609,220634,220658,220723,220749,220814,220839,220907,220933,220958,221024,221050,221116,221142,221209,221235,221301,221327,221353,221419,221446,221512,221538,221604,221630,221656,221722,221748,221813,221838,221903,221927,221953]


d=FLTARR(4096,4096,nim)
ix=[2098-2047+indgen(2048),2101+indgen(2048)]
iy=[2067-2047+indgen(2048),2132+indgen(2048)]

times=FLTARR(nim)
HWLPOS=FLTARR(4,nim)

FOR i=0,nim-1 DO BEGIN
    temp=READFITS(filename[i],header,/silent)
    HWLPOS[0,i]=sxpar(header,'HWL1POS')
    HWLPOS[1,i]=sxpar(header,'HWL2POS')
    HWLPOS[2,i]=sxpar(header,'HWL3POS')
    HWLPOS[3,i]=sxpar(header,'HWL4POS')
    PRINT,i,sxpar(header,'HSHIEXP'),HWLPOS[0,i],HWLPOS[1,i],HWLPOS[2,i],HWLPOS[3,i]    
    q=temp[ix,*]
    d[*,*,i]=q[*,iy]+0.0
    hour=FLOOR(time[i]/10000.)  ; hours
    minute=FLOOR( (time[i]-FLOOR(time[i]/10000.)*10000.)/100.)
    second=time[i]-hour*10000.-minute*100.
    times[i]=hour*3600.+minute*60.+second ; time in seconds
ENDFOR

dark=[0,13,26,39,52]

nx=256
d=REBIN(d,nx,nx,nim)
d2=d


darks=REBIN(d[*,*,dark],nx,nx,1)

c1=REBIN(d[0:nx/64-1,0:nx/64-1,*],1,1,nim)
c2=REBIN(d[nx-nx/64:nx-1,0:nx/64-1,*],1,1,nim)
c3=REBIN(d[0:nx/64-1,nx-nx/64:nx-1,*],1,1,nim)
c4=REBIN(d[nx-nx/64:nx-1,nx-nx/64:nx-1,*],1,1,nim)

cav0=REFORM(c1+c2+c3+c4)/4
; Try to separate gain variations and dark current variations
pcav=poly_fit(times,cav0,2)
cavfit=poly(times,pcav)
cavres=cav0-cavfit
cav=cavfit
cav=cav/total(REBIN(cav(dark),1))
    
FOR i=0,nim-1 DO d2[*,*,i]=d[*,*,i]-cav[i]*darks
intin=6*TOTAL(REBIN(d2,1))
    
c1=REBIN(d2[0:nx/64-1,0:nx/64-1,*],1,1,nim)
c2=REBIN(d2[nx-nx/64:nx-1,0:nx/64-1,*],1,1,nim)
c3=REBIN(d2[0:nx/64-1,nx-nx/64:nx-1,*],1,1,nim)
c4=REBIN(d2[nx-nx/64:nx-1,nx-nx/64:nx-1,*],1,1,nim)
q=[[c1,c2],[c3,c4]]
d2=d2-REBIN(q,nx,nx,nim,/sample)


temp=rebin(d2,1,1,53)
PLOT,temp

E1=REFORM(temp[1:12])
WB=REFORM(temp[14:25])
NB=REFORM(temp[40:51])

PRINT,(MAX(E1)-MIN(E1))/MEAN(E1)
PRINT,(MAX(WB)-MIN(WB))/MEAN(WB)
PRINT,(MAX(NB)-MIN(NB))/MEAN(NB)

READ,pause

; FIT OF THE I-RIPPLE
E1   = E1/MEAN(E1)
E1[0]=(E1[0]+E1[6])/2.d0
E1[1]=(E1[1]+E1[7])/2.d0
E1[2]=(E1[2]+E1[8])/2.d0
E1[3]=(E1[3]+E1[9])/2.d0
E1[4]=(E1[4]+E1[10])/2.d0
E1[5]=(E1[5]+E1[11])/2.d0
E1=E1[0:5]
WB   = WB/MEAN(WB)
WB[0]=(WB[0]+WB[6])/2.d0
WB[1]=(WB[1]+WB[7])/2.d0
WB[2]=(WB[2]+WB[8])/2.d0
WB[3]=(WB[3]+WB[9])/2.d0
WB[4]=(WB[4]+WB[10])/2.d0
WB[5]=(WB[5]+WB[11])/2.d0
WB=WB[0:5]
NB   = NB/MEAN(NB)
NB[0]=(NB[0]+NB[6])/2.d0
NB[1]=(NB[1]+NB[7])/2.d0
NB[2]=(NB[2]+NB[8])/2.d0
NB[3]=(NB[3]+NB[9])/2.d0
NB[4]=(NB[4]+NB[10])/2.d0
NB[5]=(NB[5]+NB[11])/2.d0
NB=NB[0:5]

angE1=REFORM(-HWLPOS[0,1:12]*6.d0 MOD 360.d0)*!dpi/180.d0 ; only 6 different phase angles
angE1=angE1[0:5]
angNB=REFORM(HWLPOS[1,14:25]*6.d0 MOD 360.d0)*!dpi/180.d0 ; only 6 different phase angles
angNB=angNB[0:5]


ang2=FINDGEN(10000)/9999.*2.d0*!dpi
ang3=-ang2

SET_PLOT,'PS'
DEVICE,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=25,/color
LOADCT,3

!P.MULTI=[0,1,3]

A=[1.d0,-0.1,0.1]
resp      = CURVEFIT(angE1,E1,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2,yerror=YE)
PRINT,A
PRINT,'STANDARD ERROR=',YE
funct=A[0] + (A[1]*cos(ang3/2.d0)+A[2]*sin(ang3/2.d0))^2.d0
PLOT,angE1,E1,tit='!17',xst=1,yst=1,charsize=1.5,yrange=[MIN(funct)*0.99,MAX(funct)*1.01],psym=2,xrange=[-5.4,0.1],xtit='Tuning angle (rad)',ytit='Normalized intensity'
OPLOT,ang3,funct,col=180
; I-RIPPLE VALUE E1
PRINT,'I-RIPPLE (from fit) =',(max(funct)-MIN(funct))/MEAN(funct)

A=[1.d0,-0.1,0.1]
resp      = CURVEFIT(angNB,WB,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2,yerror=YE)
PRINT,A
PRINT,'STANDARD ERROR=',YE
funct=A[0] + (A[1]*cos(ang2/2.d0)+A[2]*sin(ang2/2.d0))^2.d0
PLOT,angNB,WB,tit='!17',xst=1,yst=1,charsize=1.5,yrange=[MIN(funct)*0.99,MAX(funct)*1.01],psym=2,xrange=[-.1,5.4],xtit='Tuning angle (rad)',ytit='Normalized intensity'
OPLOT,ang2,funct,col=180
; I-RIPPLE VALUE WB
PRINT,'I-RIPPLE (from fit) =',(max(funct)-MIN(funct))/MEAN(funct)

A=[1.d0,-0.1,0.1]
resp      = CURVEFIT(angNB,NB,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2,yerror=YE)
PRINT,A
PRINT,'STANDARD ERROR=',YE
funct=A[0] + (A[1]*cos(ang2/2.d0)+A[2]*sin(ang2/2.d0))^2.d0
PLOT,angNB,NB,tit='!17',xst=1,yst=1,charsize=1.5,yrange=[MIN(funct)*0.99,MAX(funct)*1.01],psym=2,xrange=[-.1,5.4],xtit='Tuning angle (rad)',ytit='Normalized intensity'
OPLOT,ang2,funct,col=180
; I-RIPPLE VALUE NB
PRINT,'I-RIPPLE (from fit) =',(max(funct)-MIN(funct))/MEAN(funct)

DEVICE,/CLOSE
SET_PLOT,'x'
!P.MULTI=0


READ,pause


END
