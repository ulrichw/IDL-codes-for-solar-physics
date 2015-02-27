; this code tests how well a MDI-like algorithm can retrieve linewidth,
; linedepth, and continuum intensity for HMI.

PRO MDIalgorithm

FSR     = 0.172 ;FSR NB MICHELSON
dtune   = 2.0*FSR/5.0 ; tuning interval for HMI 6 positions
tune    = (DINDGEN(6)-2.5)*dtune
T       = dtune*6.0
nlam    = 38500.
dlam    = T/nlam
lam     = DINDGEN(nlam)*dlam-T/2.

I0      = 1.d0
Id      = 0.62
dv      = 0.102/2./SQRT(ALOG(2.))

line    = I0 - Id*exp(-lam^2.d0/dv^2.d0);line    = SHIFT(line,10000)

a       = WHERE(ABS(lam-tune[0]) EQ MIN(ABS(lam-tune[0])))
F0      = line[a[0]]
a       = WHERE(ABS(lam-tune[1]) EQ MIN(ABS(lam-tune[1])))
F1      = line[a[0]]
a       = WHERE(ABS(lam-tune[2]) EQ MIN(ABS(lam-tune[2])))
F2      = line[a[0]]
a       = WHERE(ABS(lam-tune[3]) EQ MIN(ABS(lam-tune[3])))
F3      = line[a[0]]
a       = WHERE(ABS(lam-tune[4]) EQ MIN(ABS(lam-tune[4])))
F4      = line[a[0]]
a       = WHERE(ABS(lam-tune[5]) EQ MIN(ABS(lam-tune[5])))
F5      = line[a[0]]

PLOT,lam,line,xst=1
OPLOT,[tune[0],tune[0]],[0,F0]
OPLOT,[tune[1],tune[1]],[0,F1]
OPLOT,[tune[2],tune[2]],[0,F2]
OPLOT,[tune[3],tune[3]],[0,F3]
OPLOT,[tune[4],tune[4]],[0,F4]
OPLOT,[tune[5],tune[5]],[0,F5]

a1exact = -2.0/T*INT_TABULATED(lam,line*COS(2.d0*!dpi*lam/T))
a1approx= -2.0/6.*(F0*COS(2.d0*!dpi*(-2.5)/6.)+F1*COS(2.d0*!dpi*(-1.5)/6.)+F2*COS(2.d0*!dpi*(-0.5)/6.)+F3*COS(2.d0*!dpi*(0.5)/6.)+F4*COS(2.d0*!dpi*(1.5)/6.)+F5*COS(2.d0*!dpi*(2.5)/6.))
b1exact = -2.0/T*INT_TABULATED(lam,line*SIN(2.d0*!dpi*lam/T))
b1approx= -2.0/6.*(F0*SIN(2.d0*!dpi*(-2.5)/6.)+F1*SIN(2.d0*!dpi*(-1.5)/6.)+F2*SIN(2.d0*!dpi*(-0.5)/6.)+F3*SIN(2.d0*!dpi*(0.5)/6.)+F4*SIN(2.d0*!dpi*(1.5)/6.)+F5*SIN(2.d0*!dpi*(2.5)/6.))

a2exact = -2.0/T*INT_TABULATED(lam,line*COS(4.d0*!dpi*lam/T))
a2approx= -2.0/6.*(F0*COS(4.d0*!dpi*(-2.5)/6.)+F1*COS(4.d0*!dpi*(-1.5)/6.)+F2*COS(4.d0*!dpi*(-0.5)/6.)+F3*COS(4.d0*!dpi*(0.5)/6.)+F4*COS(4.d0*!dpi*(1.5)/6.)+F5*COS(4.d0*!dpi*(2.5)/6.))
b2exact = -2.0/T*INT_TABULATED(lam,line*SIN(4.d0*!dpi*lam/T))
b2approx= -2.0/6.*(F0*SIN(4.d0*!dpi*(-2.5)/6.)+F1*SIN(4.d0*!dpi*(-1.5)/6.)+F2*SIN(4.d0*!dpi*(-0.5)/6.)+F3*SIN(4.d0*!dpi*(0.5)/6.)+F4*SIN(4.d0*!dpi*(1.5)/6.)+F5*SIN(4.d0*!dpi*(2.5)/6.))


linewidthexact = T/!dpi*SQRT(ALOG(SQRT(a1exact^2.d0+b1exact^2.d0)/SQRT(a2exact^2.d0+b2exact^2.d0))/3.d0)
linedepthexact = SQRT(a1exact^2.d0+b1exact^2.d0)*T/2.d0/SQRT(!dpi)/linewidthexact*exp(!dpi^2.d0*linewidthexact^2.d0/T^2.d0)
linewidthapprox = T/!dpi*SQRT(ALOG(SQRT(a1approx^2.d0+b1approx^2.d0)/SQRT(a2approx^2.d0+b2approx^2.d0))/3.d0)
linedepthapprox = SQRT(a1approx^2.d0+b1approx^2.d0)*T/2.d0/SQRT(!dpi)/linewidthapprox*exp(!dpi^2.d0*linewidthapprox^2.d0/T^2.d0)
continuumexact  = 1.d0/T*INT_TABULATED(lam,line)+linedepthexact/T*INT_TABULATED(lam,exp(-lam^2.d0/linewidthexact^2.d0))
continuumapprox = 1.d0/6.d0*(F0+F1+F2+F3+F4+F5)+1.d0/6.d0*(exp(-tune[0]^2.d0/linewidthapprox^2.d0)+exp(-tune[1]^2.d0/linewidthapprox^2.d0)+exp(-tune[2]^2.d0/linewidthapprox^2.d0)+exp(-tune[3]^2.d0/linewidthapprox^2.d0)+exp(-tune[4]^2.d0/linewidthapprox^2.d0+exp(-tune[5]^2.d0/linewidthapprox^2.d0)))

PRINT,a1exact,a1approx
PRINT,a2exact,a2approx
PRINT,linewidthexact,linewidthapprox,(linewidthapprox-linewidthexact)/linewidthexact,dv
PRINT,linedepthexact,linedepthapprox,(linedepthapprox-linedepthexact)/linedepthexact,Id
PRINT,continuumexact,continuumapprox,(continuumapprox-continuumexact)/continuumexact,I0

END
