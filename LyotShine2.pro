PRO MyLyot,X,A,F,pder

F        = FLTARR(3)

FSR      = 0.7039d0
nl       = 4001
l        = FINDGEN(nl)/FLOAT(nl-1)*16.*FSR-8.*FSR
nangle   = 100
angle    = FINDGEN(nangle)/FLOAT(nangle-1)*2.*!pi
Ep       = COMPLEXARR(nangle,nl)
integral = FLTARR(nangle)
angle0   =   6.*!pi/180.

;A[1]     = ABS(A[1])

FOR k=0,2 DO BEGIN

  ; entrance polarizer
    Amp      = FLTARR(nl)+1.0

  ; first calcite bloc
    Ef       = COS(!pi/4.+angle0) * Amp   * EXP(+COMPLEX(0,1)*!pi*l/FSR/2.) ; fast axis of calcite
    Es       = SIN(!pi/4.+angle0) * Amp   * EXP(-COMPLEX(0,1)*!pi*l/FSR/2.) ; slow axis of calcite

  ; half-wave plate
    Eo       = ( COS(!pi/4.+A[0]-X[k]) * Ef + SIN(!pi/4.+A[0]-X[k]) * Es) * EXP(+COMPLEX(0,1)*!pi/2.)
    Ee       = (-SIN(!pi/4.+A[0]-X[k]) * Ef + COS(!pi/4.+A[0]-X[k]) * Es) * EXP(-COMPLEX(0,1)*!pi/2.)

  ; second calcite bloc
    Ef       = ( COS(!pi/4.-A[0]+X[k]) * Eo + SIN(!pi/4.-A[0]+X[k]) * Ee) * EXP(+COMPLEX(0,1)*!pi*l/FSR/2.)
    Es       = (-SIN(!pi/4.-A[0]+X[k]) * Eo + COS(!pi/4.-A[0]+X[k]) * Ee) * EXP(-COMPLEX(0,1)*!pi*l/FSR/2.)

  ; quarter-wave plate
    Eo14     = ( COS(!pi/4.) * Ef + SIN(!pi/4.) * Es) * EXP(+COMPLEX(0,1)*!pi/4.) ; achromatic 1/4 wave plate
    Ee14     = (-SIN(!pi/4.) * Ef + COS(!pi/4.) * Es) * EXP(-COMPLEX(0,1)*!pi/4.)

    FOR i=0,nangle-1 DO BEGIN
        Ep[i,*]   = SIN(angle[i])*Eo14 - COS(angle[i])*Ee14
        integral[i]=TOTAL(ABS(Ep[i,*])^2.0)*(l[1]-l[0])
    ENDFOR

    F[k] = (max(integral)-min(integral))/mean(integral)
    
ENDFOR

END


PRO LyotShine2

data = [0.315,0.19,0.336]
angle= [0.,7.,15.]*!pi/180.
weight=FLTARR(3)+1.0

A    = [7.*!pi/180.];,.789]

resp = CURVEFIT(angle,data,weight,A,FUNCTION_NAME='MyLyot',TOL=1.d-15,ITMAX=4000,/DOUBLE,STATUS=stat,CHISQ=chi2,/NODERIVATIVE)

print,stat,chi2,A*[180./!pi];,1.0]


END
