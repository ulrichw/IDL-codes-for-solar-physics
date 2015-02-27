; to fit a sine wave for the wobble sequences (js_wl_fine_tune)


PRO sinewave,X,A,F,pder

A[3] = ABS(A[3])
A[0]=ABS(A[0]) ; we want a positive amplitude

F = A[0]*SIN(2.*!pi/A[3]*X+A[1])+A[2]

;pder= [
;[SIN(2.*!pi/T*X+A[1])],[A[0]*COS(2.*!pi/T*X+A[1])],[FLTARR(N_ELEMENTS(X))+1.0],[-A[0]*2.*!pi/T^2.0*COS(2.*!pi/T*X+A[1])]] ; WRONG !!!


END

PRO sinewavefit,data

nd = N_ELEMENTS(data)

weights = FLTARR(nd)+1.0

A = [10,0,5,6.]

resp      = CURVEFIT(FINDGEN(nd),data,weights,A,FUNCTION_NAME='sinewave',TOL=1.d-7,ITMAX=5000,/DOUBLE,CHISQ=chi2,/NODERIVATIVE)

PRINT,A
READ,pause

END
