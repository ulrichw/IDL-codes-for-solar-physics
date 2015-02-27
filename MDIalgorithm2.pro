; this code tests how well a MDI-like algorithm can retrieve linewidth,
; linedepth, and continuum intensity for HMI.
; USING T= 5 \Delta\lambda INSTEAD OF 6 \Delta\lambda
; AND USING ACTUAL HMI FILTER PROFILES

PRO MDIalgorithm2

ntune   = 5.d0 ; parameter in front of dvtune to calculate the period (NOT THE NUMBER OF SAMPLING POSITIONS, WHICH IS ALWAYS 6)
FSR     = 0.169 ;FSR NB MICHELSON (USED IN IBIS2.pro and with calibration11)
dtune   = 2.0*FSR/5.0 ; tuning interval for HMI 6 positions
tune    = (DINDGEN(6)-2.5)*dtune
T       = dtune*6.d0 ; actual period, not the assumed one (which is ntune*dvtune)
dlamdv  = 2.059205672212074294d-5

; Fe I PROFILE FROM 1h-AVERAGED IBIS DATA

;RESTORE,'../IBIS/IBIS2.bin'
RESTORE,'../IBIS/IBIS_lineref_94_293.bin' ; line at pixel [94,293] (very low magnetic field) from IBIS2.pro
I0      = linerefavgRCP[0]
Id      = (linerefavgRCP[0]-MIN(linerefavgRCP))/linerefavgRCP[0]
dv      = 0.1315/2./SQRT(ALOG(2.))


; GAUSSIAN PROFILE

;I0      = 1.d0
;Id      = 0.62
;dv      = 0.102/2./SQRT(ALOG(2.)) ; dv is FWHM in mA, values from Yang & Aimee
;linerefavgRCP = I0 - Id*exp(-wavelength^2.d0/dv^2.d0)


RESTORE,'../IBIS/filter.bin' ; HMI filter profiles from the code IBIS2.pro
meanfilt=MEAN([TOTAL(filters[5,*]),TOTAL(filters[4,*]),TOTAL(filters[3,*]),TOTAL(filters[2,*]),TOTAL(filters[1,*]),TOTAL(filters[0,*])])
meanfiltex=MEAN([INT_TABULATED(wavelength,filters[5,*]),INT_TABULATED(wavelength,filters[4,*]),INT_TABULATED(wavelength,filters[3,*]),INT_TABULATED(wavelength,filters[2,*]),INT_TABULATED(wavelength,filters[1,*]),INT_TABULATED(wavelength,filters[0,*])])

selec=WHERE(wavelength GE -T/2. AND wavelength LE T/2.) ; we select the wavelength range of 1 period


F0       = TOTAL(linerefavgRCP*filters[5,*]) ; I5
F1       = TOTAL(linerefavgRCP*filters[4,*]) ; I4
F2       = TOTAL(linerefavgRCP*filters[3,*]) ; I3
F3       = TOTAL(linerefavgRCP*filters[2,*]) ; I2
F4       = TOTAL(linerefavgRCP*filters[1,*]) ; I1
F5       = TOTAL(linerefavgRCP*filters[0,*]) ; I0

F0ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[5,*])) ; I5
F1ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[4,*])) ; I4
F2ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[3,*])) ; I3
F3ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[2,*])) ; I2
F4ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[1,*])) ; I1
F5ex       = INT_TABULATED(wavelength,(linerefavgRCP*filters[0,*])) ; I0



PLOT,wavelength,linerefavgRCP,xst=1
OPLOT,wavelength,filters[0,*]
OPLOT,wavelength,filters[1,*]
OPLOT,wavelength,filters[2,*]
OPLOT,wavelength,filters[3,*]
OPLOT,wavelength,filters[4,*]
OPLOT,wavelength,filters[5,*]

a1exact = -2.0/T*INT_TABULATED(wavelength[selec],linerefavgRCP[selec]*COS(2.d0*!dpi*wavelength[selec]/T))
b1exact = -2.0/T*INT_TABULATED(wavelength[selec],linerefavgRCP[selec]*SIN(2.d0*!dpi*wavelength[selec]/T))
a2exact = -2.0/T*INT_TABULATED(wavelength[selec],linerefavgRCP[selec]*COS(4.d0*!dpi*wavelength[selec]/T))
b2exact = -2.0/T*INT_TABULATED(wavelength[selec],linerefavgRCP[selec]*SIN(4.d0*!dpi*wavelength[selec]/T))
linewidthexact = T/!dpi*SQRT(ALOG(SQRT(a1exact^2.d0+b1exact^2.d0)/SQRT(a2exact^2.d0+b2exact^2.d0))/3.d0)
linedepthexact = SQRT(a1exact^2.d0+b1exact^2.d0)*T/2.d0/SQRT(!dpi)/linewidthexact*exp(!dpi^2.d0*linewidthexact^2.d0/T^2.d0)
continuumexact  = 1.d0/T*INT_TABULATED(wavelength[selec],linerefavgRCP[selec])+linedepthexact/T*INT_TABULATED(wavelength[selec],exp(-wavelength[selec]^2.d0/linewidthexact^2.d0))

a1approx= -2.0/ntune*(F0*COS(2.d0*!dpi*(-2.5)/6.)+F1*COS(2.d0*!dpi*(-1.5)/6.)+F2*COS(2.d0*!dpi*(-0.5)/6.)+F3*COS(2.d0*!dpi*(0.5)/6.)+F4*COS(2.d0*!dpi*(1.5)/6.)+F5*COS(2.d0*!dpi*(2.5)/6.)) 
b1approx= -2.0/ntune*(F0*SIN(2.d0*!dpi*(-2.5)/6.)+F1*SIN(2.d0*!dpi*(-1.5)/6.)+F2*SIN(2.d0*!dpi*(-0.5)/6.)+F3*SIN(2.d0*!dpi*(0.5)/6.)+F4*SIN(2.d0*!dpi*(1.5)/6.)+F5*SIN(2.d0*!dpi*(2.5)/6.))
a2approx= -2.0/ntune*(F0*COS(4.d0*!dpi*(-2.5)/6.)+F1*COS(4.d0*!dpi*(-1.5)/6.)+F2*COS(4.d0*!dpi*(-0.5)/6.)+F3*COS(4.d0*!dpi*(0.5)/6.)+F4*COS(4.d0*!dpi*(1.5)/6.)+F5*COS(4.d0*!dpi*(2.5)/6.))
b2approx= -2.0/ntune*(F0*SIN(4.d0*!dpi*(-2.5)/6.)+F1*SIN(4.d0*!dpi*(-1.5)/6.)+F2*SIN(4.d0*!dpi*(-0.5)/6.)+F3*SIN(4.d0*!dpi*(0.5)/6.)+F4*SIN(4.d0*!dpi*(1.5)/6.)+F5*SIN(4.d0*!dpi*(2.5)/6.))
linewidthapprox = ntune*dtune/!dpi*SQRT(ALOG(SQRT(a1approx^2.d0+b1approx^2.d0)/SQRT(a2approx^2.d0+b2approx^2.d0))/3.d0)
linedepthapprox = SQRT(a1approx^2.d0+b1approx^2.d0)*ntune*dtune/2.d0/SQRT(!dpi)/linewidthapprox*exp(!dpi^2.d0*linewidthapprox^2.d0/(ntune*dtune)^2.d0)
continuumapprox = (1.d0/6.d0*(F0+F1+F2+F3+F4+F5)+linedepthapprox/6.d0*(exp(-tune[0]^2.d0/linewidthapprox^2.d0)+exp(-tune[1]^2.d0/linewidthapprox^2.d0)+exp(-tune[2]^2.d0/linewidthapprox^2.d0)+exp(-tune[3]^2.d0/linewidthapprox^2.d0)+exp(-tune[4]^2.d0/linewidthapprox^2.d0+exp(-tune[5]^2.d0/linewidthapprox^2.d0))))
continuumapprox = continuumapprox/meanfilt
linedepthapprox = linedepthapprox/meanfilt

!P.MULTI=0

PRINT,'EXACT RESULTS',linerefavgRCP[0],dv,Id
PRINT,'RESULTS FROM CONTINUOUS ALGORITHM',continuumexact,linewidthexact,linedepthexact/continuumexact

PRINT,'RESULTS FROM DISCRETE ALGORITHM'
PRINT,'CONTINUUM RELATIVE ERROR:',(continuumapprox-linerefavgRCP[0])/linerefavgRCP[0],continuumapprox
PRINT,'LINEWIDTH RELATIVE ERROR:',(linewidthapprox-dv)/dv,linewidthapprox
linedepthapprox=linedepthapprox/continuumapprox ; !!! WARNING, CHANGES THE RESULT !!!
PRINT,'LINEDEPTH RELATIVE ERROR:',(linedepthapprox-Id)/Id,linedepthapprox

READ,pause
;line2    = I0 - Id*(6.d0/5.d0)*exp(-wavelength^2.d0/(dv*5.d0/6.d0)^2.d0);line    = SHIFT(line,10000)
;print,int_tabulated(wavelength,line)
;print,int_tabulated(wavelength,line2)

END
