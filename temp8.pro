PRO myfunction,X,A,F;,pder

angle2 = A[0]*!pi/180. + [0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.,0.,15.,30.,45.,60.,75.,90.]*!pi/180. ; tuning 1/4-wave plate
angle3 = A[1]*!pi/180. - [0.,0.,0.,0.,0.,0.,0.,30.,30.,30.,30.,30.,30.,30.,60.,60.,60.,60.,60.,60.,60.,90.,90.,90.,90.,90.,90.,90.,120.,120.,120.,120.,120.,120.,120.,150.,150.,150.,150.,150.,150.,150.,180.,180.,180.,180.,180.,180.,180.]*!pi/180. ; tuning polarizer

angle1 = A[3]*!pi/180.

; 1/4-wave plate at the entrance of E2
;FOR i=0,48 DO BEGIN
Eo     = A[2] * EXP(+COMPLEX(0,1)*!pi/4.) *  COS(angle1)
Ee     = A[2] * EXP(-COMPLEX(0,1)*!pi/4.) *(-SIN(angle1))
; tuning 1/4-wave plate
Eo2    = (COS(angle2-angle1) * Eo + SIN(angle2-angle1) * Ee) * EXP(+COMPLEX(0,1)*!pi/4.)
Ee2    = (COS(angle2-angle1) * Ee - SIN(angle2-angle1) * Eo) * EXP(-COMPLEX(0,1)*!pi/4.)
; tuning polarizer
Ep     = COS(angle3-angle2)  * Eo2 + SIN(angle3-angle2) * Ee2

F      = ABS(Ep)^2.d0

;ENDFOR

END



PRO temp8

Inten  = [0.036,0.043,0.064,0.1,0.144,0.174,0.180,0.071,0.105,0.138,0.169,0.184,0.174,0.139,0.144,0.170,0.182,0.177,0.149,0.109,0.066,0.181,0.174,0.153,0.118,0.074,0.045,0.036,0.145,0.112,0.079,0.049,0.036,0.046,0.076,0.073,0.046,0.036,0.041,0.07,0.111,0.149,0.036,0.043,0.064,0.101,0.145,0.175,0.180]-0.036 ; dark level removed

ni     = N_ELEMENTS(Inten)

weight = FLTARR(ni)+1.0

;N = 90
;residu = 1.d10
;
;FOR i=0,N-1 DO BEGIN
;    PRINT,i
;    FOR j=0,N-1 DO BEGIN
;        ;FOR k=0,N-1 DO BEGIN
;            A      = [i/FLOAT(N-1)*180.-90.,j/FLOAT(N-1)*180.,0.384708,49.5]
;            myfunction,findgen(49),A,F
;            b      = TOTAL((F-Inten)^2.0)
;            IF(b LT residu) THEN BEGIN
;                residu = b
;                bingo=[i,j];[i,j,k]
;            ENDIF
;        ;ENDFOR
;    ENDFOR
;ENDFOR
;PRINT,residu,[bingo[0]/FLOAT(N-1)*90.,bingo[1]/FLOAT(N-1)*180.,0.384708,49.5]
;myfunction,findgen(49),[bingo[0]/FLOAT(N-1)*90.,bingo[1]/FLOAT(N-1)*180.,0.384708,49.4],F & plot,F


A      = [0.0,45.0,0.3847,45.]
resp   = CURVEFIT(FINDGEN(ni),Inten,weight,A,FUNCTION_NAME='myfunction',TOL=1.d-9,ITMAX=4000,STATUS=stat,CHISQ=chi2,/NODERIVATIVE)

SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=14,ysize=11,/color
loadct,3
plot,Inten,thick=3
oplot,resp,color=180
DEVICE,/CLOSE
SET_PLOT,'x'


PRINT,A,chi2
READ,pause



END
