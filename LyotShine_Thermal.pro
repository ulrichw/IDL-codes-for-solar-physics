; this programs fits for the average phase and contrast of the Lyot
; non tunable elements using  LyotProfiles2006aug3.fits


PRO MyLyot,X,A,F,pder

A[0] = ABS(A[0])
A[1] = ABS(A[1])
A[3] = ABS(A[3])
A[5] = ABS(A[5])
A[7] = ABS(A[7])
A[10]= ABS(A[10])
A[12]= ABS(A[12])
A[13]= ABS(A[13])
A[14]= ABS(A[14])

nel = N_ELEMENTS(X)
F   = FLTARR(nel)

RESTORE,'temp.bin' ; contains wavelengths, blocker filter profile, and Iripple

dpi   = 2.d0*!dpi

FSR   = FLTARR(4) ; in A
FSR[0]= 1.405d0;A[13]; 1.405d0               ; for E2
FSR[1]= 2.779d0;A[14] ;2.779d0               ; for E3
FSR[2]= 5.682d0               ; for E4
FSR[3]= 11.354d0              ; for E5

E1    = (1.d0+A[10]*COS(dpi*l/0.690+A[11]))/2.d0
E2    = (1.d0+A[1] *COS(dpi*l/FSR[0]+A[2]))/2.d0
E3    = (1.d0+A[3] *COS(dpi*l/FSR[1]+A[4]))/2.d0
E4    = (1.d0+A[5] *COS(dpi*l/FSR[2]+A[6]))/2.d0
E5    = (1.d0+A[7] *COS(dpi*l/FSR[3]+A[8]))/2.d0
line  = (1.d0-A[12] *EXP(-(l-A[9])^2.d0/A[13]^2.d0)); SOLAR LINE
;line  = (1.d0-0.65*EXP(-(l-A[9])^2.d0/0.0612^2.d0)); SOLAR LINE


F     = A[0]*blocker0*E1*E2*E3*E4*E5*line+A[14]

dE1d1 = COS(dpi*l/0.690+A[11])/2.d0
dE1d2 = -A[10]*SIN(dpi*l/0.690+A[11])/2.d0
dE2d1 = COS(dpi*l/FSR[0]+A[2])/2.d0
dE2d2 = -A[1]*SIN(dpi*l/FSR[0]+A[2])/2.d0
dE3d1 = COS(dpi*l/FSR[1]+A[4])/2.d0
dE3d2 = -A[3]*SIN(dpi*l/FSR[1]+A[4])/2.d0
dE4d1 = COS(dpi*l/FSR[2]+A[6])/2.d0
dE4d2 = -A[5]*SIN(dpi*l/FSR[2]+A[6])/2.d0
dE5d1 = COS(dpi*l/FSR[3]+A[8])/2.d0
dE5d2 = -A[7]*SIN(dpi*l/FSR[3]+A[8])/2.d0
dlined12= -EXP(-(l-A[9])^2.d0/A[13]^2.d0)
dlined13= -A[12]*2.d0*(l-A[9])^2.d0/A[13]^3.d0*EXP(-(l-A[9])^2.d0/A[13]^2.d0)
dlined9 = -A[12]*2.d0*(l-A[9])/A[13]^2.d0*EXP(-(l-A[9])^2.d0/A[13]^2.d0)

pder  = [ [blocker0*E1*E2*E3*E4*E5*line] , [dE2d1*E1*E3*E4*E5*A[0]*blocker0*line] , [dE2d2*E1*E3*E4*E5*A[0]*blocker0*line] , [dE3d1*E1*E2*E4*E5*A[0]*blocker0*line] , [dE3d2*E1*E2*E4*E5*A[0]*blocker0*line] , [dE4d1*E1*E2*E3*E5*A[0]*blocker0*line] , [dE4d2*E1*E2*E3*E5*A[0]*blocker0*line] , [dE5d1*E1*E2*E3*E4*A[0]*blocker0*line] , [dE5d2*E1*E2*E3*E4*A[0]*blocker0*line] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined9] , [dE1d1*E2*E3*E4*E5*A[0]*blocker0*line] , [dE1d2*E2*E3*E4*E5*A[0]*blocker0*line] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined12] , [E1*E2*E3*E4*E5*A[0]*blocker0*dlined13] , [FLTARR(nel)+1.d0] ]

;xx = SHIFT(DINDGEN(nel)-nel/2.d0,nel/2.d0)
;xx[nel/2.d0]=-xx[nel/2.d0]
;param=0.0024
;kernel=1./(1.+param*ABS(xx)+param*xx^2.+param*ABS(xx)^3.)
;fkernel = FFT(kernel,/DOUBLE)

;F = FLOAT(FFT(FFT(F,/DOUBLE)*fkernel,/INVERSE,/DOUBLE))
;FOR i=0,N_ELEMENTS(A)-1 DO pder[*,i] = FLOAT(FFT(FFT(pder[*,i],/DOUBLE)*fkernel,/INVERSE,/DOUBLE))

;F = F*nel
;pder = pder*nel


END



PRO LyotShine_Thermal,number,E1,E2,E3,E4,E5;,linedepth,linewidth
; we correct for the potential change in intensity
; we correct for the curvature of the spectrum due to the spectrograph


!P.MULTI=0


filename = 'lyotThermaltestProfiles2006oct19.fits' ; WITH blocker and front window, thermal tests spanning a day

d        = READFITS(filename)
d        = REFORM(d[*,number])              ; 12 different temperatures
blocker  = 1


nelem     = N_ELEMENTS(d)
data      = d
data      = data/MAX(data) 
data      = REVERSE(data)
l         = READFITS('ref2006oct16.fits')
spectrum  = REVERSE(REFORM(l[*,1]))

CASE number OF
    0: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.01199 ; 0
    1: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0387 ; 1
    2: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.01399 ; 2
    3: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.00305 ; 3
    4: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.00305 ; 4
    5: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0397 ; 5
    6: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.04615 ; 6
    7: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.03832 ; 7
    8: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0435 ; 8
    9: l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0456 ; 9
    10:l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0521 ; 10
    11:l         = REVERSE(REFORM(l[*,0]))-6173.3433+0.0549 ; 11
ENDCASE



; We correct for the solar spectrum curvature due to the spectrograph
x         = FINDGEN(1320)
data      = data/(0.935+0.000305037*x+4.03807e-08*x^2.0-1.69869e-10*x^3.0)

q         = READFITS('blocker11.fits')
blocker0  = SHIFT(q[*,1],2.75/(q[1,0]-q[0,0])) ; careful, the shift depends on the focal distance
blocker0  = INTERPOL(blocker0/100.d0,q[*,0]-6173.3433,l)
SAVE,l,blocker0,FILE='temp.bin'

nel        = N_ELEMENTS(data)
weight     = FLTARR(nel)+1.d0
d          = 0.0

A          = FLTARR(15)
A[0]       = MAX(data)
A[1]       = 1.0
A[3]       = 1.0
A[5]       = 1.0
A[7]       = 1.0
A[12]       = 0.52 ; linedepth
A[13]      = 0.07 ; linewidth
A[9]      = -0.014 ; central wavelength
A[14]      = 0.08  ; dark level
A[10]      = 1.0 ; contrast E1

resp = CURVEFIT(FINDGEN(nel),data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-9,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2)

print,chi2
PRINT,'contrasts',A[1],A[3],A[5],A[7],A[10]
PRINT,'amplitude',A[0]
PRINT,'Solar Line',A[9],A[12],A[13]
PRINT,'angles',A[2]*180./!pi,A[4]*180./!pi,A[6]*180./!pi,A[8]*180./!pi,A[11]*180./!pi
PRINT,'Dark Level Offset',A[14]

linedepth = A[12]
linewidth = A[13]
E1=A[11]*180./!pi
E2=A[2]*180./!pi
E3=A[4]*180./!pi
E4=A[6]*180./!pi
E5=A[8]*180./!pi

SET_PLOT,'ps'
DEVICE,FILE='yo2.ps',/color,xoffset=0,yoffset=0,xsize=14,ysize=11
loadct,3
plot,l,data,xst=1,yst=1,charsize=1.5,xtit='Wavelength (A)',ytit='Normalized Intensity',thick=2
oplot,l,resp,col=180
oplot,[0,0],[0,1]

DEVICE,/CLOSE

SET_PLOT,'x'

;READ,pause

END

; TO PLOT THE WAVELENGTH SHIFT AS A FUNCTION OF THE TEMPERATURE

;res=poly_fit(T,-E1*0.690/360.*1000.,2,yfit=y1)
;res=poly_fit(T,-E2*1.405/360.*1000.,2,yfit=y2)
;res=poly_fit(T,-E3*2.779/360.*1000.,2,yfit=y3)
;res=poly_fit(T,-E4*5.682/360.*1000.,2,yfit=y4)
;set_plot,'ps'
;device,file='yo.ps',xoffset=0,yoffset=0,xsize=21,ysize=26,/color
;plot,T,-E1*0.690/360.*1000.,xst=1,charsize=1.5,xtit='T (Celsius)',ytit='Wavelength Shift (mA)'
;oplot,T,y1,col=180
;plot,T,-E2*1.405/360.*1000.,xst=1,charsize=1.5,xtit='T (Celsius)',ytit='Wavelength Shift (mA)'
;oplot,T,y2,col=180
;plot,T,-E3*2.779/360.*1000.,xst=1,charsize=1.5,xtit='T (Celsius)',ytit='Wavelength Shift (mA)'
;oplot,T,y3,col=180
;plot,T,-E4*5.682/360.*1000.,xst=1,charsize=1.5,xtit='T (Celsius)',ytit='Wavelength Shift (mA)'
;oplot,T,y4,col=180
;device,/close
