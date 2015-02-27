PRO profile,x,A,F,pder

;x=wavelength
;A[0]=Free Spectral Range
;A[1]=contrast

A[1]=ABS(A[1])
F=(1.0+A[1]*COS(2.*!dpi/A[0]*(x-A[2])))/2.

IF N_PARAMS() GE 4 THEN BEGIN
pder=[ [!dpi*A[1]*(x-A[2])/A[0]^2.0*SIN(2.*!dpi/A[0]*(x-A[2]))],[1./2.*COS(2.*!dpi/A[0]*(x-A[2]))],[!dpi*A[1]/A[0]*SIN(2.*!dpi/A[0]*(x-A[2]))] ]
ENDIF
                                                                                                                                 
END

PRO gaussian,x,A,F,pder

;A[0]=amplitude
;A[1]=central frequency
;A[2]=FWHM

F=A[0]*exp(-(x-A[1])^2.0*4.*ALOG(2.)/A[2]^2.0)

IF N_PARAMS() GE 4 THEN BEGIN
pder= [ [F/A[0]],[F*(A[1]-x)*8.*ALOG(2.)/A[2]^2.0],[F*8.*ALOG(2.)*(x-A[1])^2.0/A[2]^3.0] ]
ENDIF

END

PRO seb,x,A,F,pder
FSRmich=0.189d0
lam0=6173.d0
wlyot = FSRmich*lam0/6768.0*4.0
cmich=2*!dpi/wlyot
F=A[0]*cos(cmich*x+A[1])
pder=[ [F/A[0]],[-A[0]*sin(cmich*x+A[1])] ]
END


;--------------------------------------------------------------------------------------------
;
; PROGRAM FOR THE LYOT FILTER CALIBRATION
;
;--------------------------------------------------------------------------------------------

PRO Lyotcalib,filename

GOTO,spatial
num=5

; READ THE CONFIGURATION FILE
; this files must contain 6 file names: 5 transmission profiles for
; the Lyot elements, and 1 profile for the complete filter
filenames=STRARR(num)
OPENR,1,filename
READF,1,filenames
CLOSE,1

; I/ ON-AXIS TRANSMISSION PROFILES FOR THE SINGLE ELEMENTS
FOR i=0,num-2 DO BEGIN
RESTORE,filenames[i]
lam=lam+6173.0
nlam=N_ELEMENTS(lam)

A=[1.0,1.0,6172.0]
weights=FLTARR(nlam)+1.0
res=CURVEFIT(lam,transmission,weights,A,FUNCTION_NAME='profile',itmax=600,TOL=1.e-10)
PRINT,'LYOT ELEMENT NUMBER',i
PRINT,'FREE SPECTRAL RANGE=',A[0]
PRINT,'CONTRAST=',A[1]
ENDFOR

; II/ ON-AXIS TRANSMISSION PROFILE FOR THE COMPLETE LYOT FILTER

RESTORE,filenames[num-1]
lam=lam+6173.0

A=[1.0,6173.0,1.0]
weights=FLTARR(nlam)+1.0
res=CURVEFIT(lam,transmission,weights,A,FUNCTION_NAME='gaussian',itmax=600,TOL=1.e-10)
b=WHERE(lam GE A[1]-A[2]/2. AND lam LE A[1]+A[2]/2.)
res=CURVEFIT(lam[b],transmission[b],weights[b],A,FUNCTION_NAME='gaussian',itmax=600,TOL=1.e-15)
PRINT,'AMPLITUDE',A[0]
PRINT,'CENTRAL FREQUENCY',A[1]
PRINT,'FWHM',A[2]

; III/ STABILITY TO TEMPERATURE

; IV/ SPATIAL DEPENDENCE
spatial:
set_plot,'x'
!p.multi=[0,1,3]
window,0,retain=2,ysize=800

largeur   =FLTARR(5); provided by Jesper Schou
largeur[0]=1.62;1.96;2.01    ; radius (in mm) of the beam on the E1 Lyot element aperture (1.62 in CAL, 1.96 in OBS)
largeur[1]=0.85    ; E2
largeur[2]=1.16    ; E3
largeur[3]=1.39    ; E4
largeur[4]=1.59    ; E5

; If we place E2 in E5
;largeur[0]=2.01    ; radius (in mm) of the beam on the E1 Lyot element aperture
;largeur[4]=0.85    ; E5
;largeur[2]=1.16    ; E3
;largeur[3]=1.39    ; E4
;largeur[1]=1.59    ; E2


nx=800.0

FOR ii=0,0 DO BEGIN ;1,4
 e1=MRDFITS('elemOne5oct2006.fits')*!pi/180.0 ; new element E1
;e1=mrdfits('element'+STRTRIM(STRING(ii+1),1)+'phase.fits')*!pi/180.0
;e1=mrdfits('wideA+B_C_30jun2005.fits')*!pi/180.0 ; new element E2
a=WHERE(e1 NE 0.0)
av=MEAN(e1[a])
e1[a]=e1[a]-av

a=WHERE(e1[*,400] NE 0.0)   ; the resolution depends on the picture ?
dx=30.0/FLOAT(N_ELEMENTS(a)) ; the disc has a diameter of 30mm for measurements prior to end of 2005

ny=ROUND(largeur[ii]/dx*2.0)
disc=SHIFT(DIST(ny,ny)*dx,ny/2,ny/2)
a=WHERE(disc le largeur[ii],COMPLEMENT=b)
disc[a]=1.0
na=FLOAT(N_ELEMENTS(a))
disc[b]=0.0
resultat=FLTARR(2,nx,nx)

; the transmission profile is (1+contrast cos(---+2 phase))/2
FOR i=ny/2,nx-1-ny/2 DO BEGIN
;PRINT,i
FOR j=ny/2,nx-1-ny/2 DO BEGIN

IF (ny MOD 2 EQ 0) THEN disc2=disc*e1[i-ny/2:i+ny/2-1,j-ny/2:j+ny/2-1] ELSE disc2=disc*e1[i-ny/2:i+ny/2,j-ny/2:j+ny/2]

youp=TOTAL(exp(COMPLEX(0,1)*disc2[a]))
resultat[0,i,j]=1./na*SQRT(youp*TOTAL(exp(-COMPLEX(0,1)*disc2[a]))) ; contrast
IF(resultat[0,i,j] NE 0.0) THEN resultat[1,i,j]=ALOG(youp/resultat[0,i,j]/na)/COMPLEX(0,1)  ELSE resultat[1,i,j]=0.0 ; phase
ENDFOR
tvim,e1*180./!pi,/scale,tit='!17',xtit='!17x',ytit='!17y',stit='!17Phase (in degrees)'
tvim,resultat[0,*,*],/scale
tvim,resultat[1,*,*]*180./!pi,/scale
ENDFOR

disc=dist(nx,nx)*dx
disc=shift(disc,nx/2,nx/2)
a=where(disc LE 12.7,COMPLEMENT=b)  ; 12.7 mm is the image radius at the exit of the Lyot filter in Calmode (11.9 in Obsmode)
disc[a]=1.0
disc[b]=0.0
resultat[0,*,*]=resultat[0,*,*]*disc
resultat[1,*,*]=resultat[1,*,*]*disc
SAVE,resultat,file='E'+STRTRIM(STRING(ii),1)+'.bin'

PRINT,'minimum phase',MIN(resultat[1,*,*])*180./!pi
temp=REFORM(resultat[0,*,*])
PRINT,'minimum contrast',MIN(temp[a])
temp[WHERE(temp EQ 0.0)]=2.0
temp2=REFORM(resultat[1,*,*])
temp2[WHERE(temp2 EQ 0.0)]=MAX(resultat[1,*,*])+5.0*!pi/180.

set_plot,'ps'
!p.multi=[0,1,3]
LOADCT,4
device,file='E'+STRTRIM(STRING(ii),1)+'.ps',bits=24,xsize=20,xoffset=0,ysize=26,yoffset=0,/color
tvim,e1*180./!pi,/scale,tit='!17Element E'+STRTRIM(STRING(ii+1),1)+': non uniformity map',xtit='!17x',ytit='!17y',stit='!17Phase (in degrees)'
tvim,temp2*180./!pi,/scale,tit='!17Transmission I(!7k!17)=(1+Contrast cos(2!7pk!17/FSR+ Phase))/2',xtit='!17x',ytit='!17y',stit='!17Phase (in degrees)'
tvim,temp,/scale,range=[MIN(temp),MAX(resultat[0,*,*])],tit='!17',xtit='!17x',ytit='!17y',stit='!17Contrast'
device,/close
set_plot,'x'
ENDFOR
END
