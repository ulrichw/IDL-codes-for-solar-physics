; this programs fits for the average phase and contrast of the Lyot
; non tunable elements using  LyotProfiles2006aug3.fits


PRO MyLyot,X,A,F,pder

A[0] = ABS(A[0])
A[1] = ABS(A[1])
A[3] = ABS(A[3])
A[5] = ABS(A[5])
A[7] = ABS(A[7])
A[9] = ABS(A[9])
A[10]= ABS(A[10])
A[13]= ABS(A[13])
A[14]= ABS(A[14])
A[15]= ABS(A[15])

nel = N_ELEMENTS(X)
F   = FLTARR(nel)

RESTORE,'temp.bin' ; contains wavelengths, blocker filter profile, and Iripple

dpi   = 2.d0*!dpi

FSR   = FLTARR(4) ; in A
FSR[0]= A[13]; 1.405d0               ; for E2
FSR[1]= A[14] ;2.779d0               ; for E3
FSR[2]= 5.682d0               ; for E4
FSR[3]= 11.354d0              ; for E5

E1    = (1.d0+A[15] *COS(dpi*l/0.7039+A[16]))/2.d0  ; I try to fit for E1 if light from the source is partly polarized
;E1    = (1.d0+A[15] *COS(dpi*l/0.690))/2.d0
;E1    = 1.0 ; check or uncheck if you fit for E1
E2    = (1.d0+A[1] *COS(dpi*l/FSR[0]+A[2]))/2.d0
E3    = (1.d0+A[3] *COS(dpi*l/FSR[1]+A[4]))/2.d0
E4    = (1.d0+A[5] *COS(dpi*l/FSR[2]+A[6]))/2.d0
E5    = (1.d0+A[7] *COS(dpi*l/FSR[3]+A[8]))/2.d0
line  = (1.d0-A[9] *EXP(-(l-A[11])^2.d0/A[10]^2.d0)); SOLAR LINE

F     = A[0]*blocker0*E1*E2*E3*E4*E5*line+A[12]

dE1d1 = COS(dpi*l/0.7039+A[16])/2.d0
dE1d2 = -A[15]*SIN(dpi*l/0.7039+A[16])/2.d0
;dE1d1 = COS(dpi*l/0.690)/2.d0
dE2d1 = COS(dpi*l/FSR[0]+A[2])/2.d0
dE2d2 = -A[1]*SIN(dpi*l/FSR[0]+A[2])/2.d0
dE3d1 = COS(dpi*l/FSR[1]+A[4])/2.d0
dE3d2 = -A[3]*SIN(dpi*l/FSR[1]+A[4])/2.d0
dE4d1 = COS(dpi*l/FSR[2]+A[6])/2.d0
dE4d2 = -A[5]*SIN(dpi*l/FSR[2]+A[6])/2.d0
dE5d1 = COS(dpi*l/FSR[3]+A[8])/2.d0
dE5d2 = -A[7]*SIN(dpi*l/FSR[3]+A[8])/2.d0
dlined9 = -EXP(-(l-A[11])^2.d0/A[10]^2.d0)
dlined10= -A[9]*2.d0*(l-A[11])^2.d0/A[10]^3.d0*EXP(-(l-A[11])^2.d0/A[10]^2.d0)
dlined11= -A[9]*2.d0*(l-A[11])/A[10]^2.d0*EXP(-(l-A[11])^2.d0/A[10]^2.d0)

pder  = [ [blocker0*E1*E2*E3*E4*E5*line] , [dE2d1*E1*E3*E4*E5*A[0]*blocker0*line] , [dE2d2*E1*E3*E4*E5*A[0]*blocker0*line] , [dE3d1*E1*E2*E4*E5*A[0]*blocker0*line] , [dE3d2*E1*E2*E4*E5*A[0]*blocker0*line] , [dE4d1*E1*E2*E3*E5*A[0]*blocker0*line] , [dE4d2*E1*E2*E3*E5*A[0]*blocker0*line] , [dE5d1*E1*E2*E3*E4*A[0]*blocker0*line] , [dE5d2*E1*E2*E3*E4*A[0]*blocker0*line] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined9] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined10] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined11] , [FLTARR(nel)+1.d0] , [A[0]*blocker0*dpi*l*A[1]/FSR[0]^2.0*SIN(dpi*l/FSR[0]+A[2])/2.d0*E1*E3*E4*E5*line] , [A[0]*blocker0*dpi*l*A[3]/FSR[1]^2.0*SIN(dpi*l/FSR[1]+A[4])/2.d0*E1*E2*E4*E5*line] , [dE1d1*E2*E3*E4*E5*A[0]*blocker0*line] , [dE1d2*E2*E3*E4*E5*A[0]*blocker0*line] ]

xx = SHIFT(DINDGEN(nel)-nel/2.d0,nel/2.d0)
xx[nel/2.d0]=-xx[nel/2.d0]
param=0.0011;0.00175
kernel=1./(1.+param*ABS(xx)+param*xx^2.+param*ABS(xx)^3.)
fkernel = FFT(kernel,/DOUBLE)

F = FLOAT(FFT(FFT(F,/DOUBLE)*fkernel,/INVERSE,/DOUBLE))
FOR i=0,N_ELEMENTS(A)-1 DO pder[*,i] = FLOAT(FFT(FFT(pder[*,i],/DOUBLE)*fkernel,/INVERSE,/DOUBLE))

F = F*nel
pder = pder*nel


END



PRO LyotShine_Sun
; we correct for the potential change in intensity
; we correct for the curvature of the spectrum due to the spectrograph


!P.MULTI=0

;filename = '/scr109a/couvidat/HMI_DATA/LyotSunsansP2006oct16.fits'; the complete Lyot with no blocker but front window. E1 with no polarizer and sunlight. Taken around 1 PDT, with front window
filename = '/scr109a/couvidat/HMI_DATA/lyotsunsetSprofiles2006oct16.fits'; 36 Lyot profiles (tuning a polaroid in front of E1), no blocker, but front window IRIS=29 mm, line shifted by -9.59 mA (13:09 PDT)
;filename = '/scr109a/couvidat/HMI_DATA/lyotsunsetTprofiles2006oct16.fits'; 36 Lyot profiles (tuning a polaroid in front of E1), no blocker but front window, taken at 13:23 PDT, meaning line shifted by -9.14 IRIS=12 mm
;filename = '/scr109a/couvidat/HMI_DATA/lyotsunsetXprofiles2006oct18.fits'; idem but WITH blocker and front window, taken at 8:30 PDT, meaning line shifted by -16.8 mA (-815.95 m/s)
;filename = '/scr109a/couvidat/HMI_DATA/lyotThermaltestProfiles2006oct19.fits' ; WITH blocker and front window, thermal tests spanning a day

d        = READFITS(filename)
blocker  = 0
IF filename EQ '/scr109a/couvidat/HMI_DATA/LyotSunsansP2006oct16.fits' THEN d = FLOAT(REFORM(REBIN(d[*,100:190,1],1320,1)))
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetTprofiles2006oct16.fits' OR filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetSprofiles2006oct16.fits' THEN BEGIN
    inten    = FLTARR(36)       ; TO CHECK FOR AN I-RIPPLE
    FOR i=0,35 DO inten[i]=TOTAL(d[*,i])
    FOR i=0,35 DO d[*,i]=d[*,i]/TOTAL(d[*,i]) ; to compensate for intensity variation along the tuning sequence
   ;d = REFORM(d[*,0]+d[*,6]+d[*,12]+d[*,18]+d[*,24]+d[*,30])/6.0
    d = REFORM(REBIN(d,1320,1))
ENDIF
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetXprofiles2006oct18.fits' THEN BEGIN
    d = REFORM(REBIN(d,1320,1))
    blocker = 1
ENDIF
  IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotThermaltestProfiles2006oct19.fits' THEN BEGIN
    d = REFORM(d[*,0]) ; 12 different temperatures
    blocker = 1
ENDIF  


nelem      = N_ELEMENTS(d)
data       = d
data       = REVERSE(data)

l          = READFITS('ref2006oct16.fits') ; reference solar spectrum taken at that time
spectrum   = REVERSE(REFORM(l[*,1]))

IF filename EQ '/scr109a/couvidat/HMI_DATA/LyotSunsansP2006oct16.fits'        THEN l = REVERSE(REFORM(l[*,0]))-6173.3433-0.04385   ; for LyotSunsansP2006oct16.fits
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetTprofiles2006oct16.fits' THEN l = REVERSE(REFORM(l[*,0]))-6173.3433-0.040625   ; for lyotsunsetTprofiles2006oct16.fits
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetSprofiles2006oct16.fits' THEN l = REVERSE(REFORM(l[*,0]))-6173.3433-0.01175
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotsunsetXprofiles2006oct18.fits' THEN l = REVERSE(REFORM(l[*,0]))-6173.3433+0.156;0.15645   ; for lyotsunsetXprofiles2006oct18.fits
IF filename EQ '/scr109a/couvidat/HMI_DATA/lyotThermaltestProfiles2006oct19.fits' THEN l = REVERSE(REFORM(l[*,0]))-6173.3433+0.02;0.02525

; We correct for the solar spectrum curvature due to the spectrograph
x           = FINDGEN(1320)
data        = data/(0.935+0.000305037*x+4.03807e-08*x^2.0-1.69869e-10*x^3.0)
data        = data/MAX(data) 

IF blocker EQ 1 THEN BEGIN
    q          = READFITS('blocker11.fits')
    blocker0   = SHIFT(q[*,1],2.75/(q[1,0]-q[0,0])) ; careful, the shift depends on the focal distance
    blocker0   = INTERPOL(blocker0/100.d0,q[*,0]-6173.3433,l)
ENDIF ELSE blocker0=1.0
SAVE,l,blocker0,FILE='temp.bin'

nel        = N_ELEMENTS(data)
weight     = FLTARR(nel)+1.d0
d          = 0.0

A          = FLTARR(17)
A[0]       = MAX(data)
A[1]       = 1.0
A[3]       = 1.0
A[5]       = 1.0
A[7]       = 1.0
A[9]       = 0.52 ; linedepth
A[10]      = 0.07 ; linewidth
A[11]      = -0.014 ; central wavelength
A[12]      = 0.08  ; dark level
A[13]      = 1.405 ; FSR E2
A[14]      = 2.779 ; FSR E3
A[15]      = 0.1 ; contrast E1

resp = CURVEFIT(FINDGEN(nel),data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2)

print,chi2
PRINT,'contrasts',A[1],A[3],A[5],A[7],A[15]
PRINT,'amplitude',A[0]
PRINT,'Solar Line',A[9],A[10],A[11]
PRINT,'angles',A[2]*180./!pi,A[4]*180./!pi,A[6]*180./!pi,A[8]*180./!pi,A[16]*180./!pi
PRINT,'Dark Level Offset',A[12]
PRINT,'FSRs',A[13],A[14]

SET_PLOT,'ps'
DEVICE,FILE='yo2.ps',/color,xoffset=0,yoffset=0,xsize=20,ysize=16
LOADCT,4
plot,l,data,xst=1,yst=1,charsize=1.5,xtit='Wavelength (A)',ytit='Normalized Intensity',thick=2,yrange=[0,1.4]
oplot,l,resp,col=180
oplot,[0,0],[0,1.5]
dpi   = 2.d0*!dpi
FSR   = FLTARR(4)
FSR[0]= A[13]; 1.405d0     
FSR[1]= A[14] ;2.779d0           
FSR[2]= 5.682d0          
FSR[3]= 11.354d0    
E1    = (1.d0+A[15]*COS(dpi*l/0.7039+A[16]))/2.d0
E2    = (1.d0+A[1] *COS(dpi*l/FSR[0]+A[2] ))/2.d0
E3    = (1.d0+A[3] *COS(dpi*l/FSR[1]+A[4] ))/2.d0
E4    = (1.d0+A[5] *COS(dpi*l/FSR[2]+A[6] ))/2.d0
E5    = (1.d0+A[7] *COS(dpi*l/FSR[3]+A[8] ))/2.d0
F     = A[0]*blocker0*E1*E2*E3*E4*E5+A[12]
nel   = N_ELEMENTS(F)
xx    = SHIFT(DINDGEN(nel)-nel/2.d0,nel/2.d0)
xx[nel/2.d0]=-xx[nel/2.d0]

; WE CONVOLVE TO SIMULATE THE FINITE RESOLUTION 
param=0.0011;0.00175
kernel=1./(1.+param*ABS(xx)+param*xx^2.+param*ABS(xx)^3.)
fkernel = FFT(kernel,/DOUBLE)
F     = FLOAT(FFT(FFT(F,/DOUBLE)*fkernel,/INVERSE,/DOUBLE))
F     = F*nel
line  = A[0]*(1.d0-A[9] *EXP(-(l-A[11])^2.d0/A[10]^2.d0)); SOLAR LINE
line  = FLOAT(FFT(FFT(line,/DOUBLE)*fkernel,/INVERSE,/DOUBLE))
line = line*nel

oplot,l,F,color=90,thick=2
oplot,l,line;/MAX(line)*MAX(F)
DEVICE,/CLOSE

aaa = WHERE(ABS(F-MAX(F)/2.d0) LT 0.0098)
PRINT,'FWHM NON-TUNABLE LYOT=',l[aaa],l[aaa[N_ELEMENTS(aaa)-1]]-l[aaa[0]]

SET_PLOT,'x'
SAVE,l,F,file='temp.bin'

READ,pause

END
