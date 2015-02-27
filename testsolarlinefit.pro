PRO fit3,X,A,F,pder

A[0]=ABS(A[0])
A[1]=ABS(A[1])

width = 0.061

F   = A[0]*(1.0-A[1]*EXP(-(X-A[2])^2.0/width^2.0))

pder=[ [F/A[0]],[-A[0]*EXP(-(X-A[2])^2.0/width^2.0)], [-A[0]*A[1]*2.0*(X-A[2])/width^2.0*EXP(-(X-A[2])^2.0/width^2.0)] ]

END


PRO fit4,X,A,F,pder

A[0]=ABS(A[0])
A[1]=ABS(A[1])
A[3]=ABS(A[3])

F   = A[0]*(1.0-A[1]*EXP(-(X-A[2])^2.0/A[3]^2.0))

pder=[ [F/A[0]],[-A[0]*EXP(-(X-A[2])^2.0/A[3]^2.0)], [-A[0]*A[1]*2.0*(X-A[2])/A[3]^2.0*EXP(-(X-A[2])^2.0/A[3]^2.0)], [-A[0]*A[1]*2.0*(X-A[2])^2.0/A[3]^3.0*EXP(-(X-A[2])^2.0/A[3]^2.0)] ]

END

PRO testsolarlinefit

WINDOW,0,RETAIN=2

nlamp = 1223 ; Set wavelength points for line profile
lamp  = -0.3D0+DINDGEN(nlamp)*0.6D/(nlamp-1)
OPENR,1,'Ulrich_Fe_0.txt'
roger  = FLTARR(2,98)
READF,1,roger
CLOSE,1
nroger = 98
rlam   = REFORM(roger(0,*))
rp     = REFORM(roger(1,*))
rp1    = rp/INTERPOL([rp(0),rp(nroger-1)],[rlam(0),rlam(nroger-1)],rlam)
line   = INTERPOL(rp,rlam,lamp)
line   = line/INTERPOL([line(0),line(nlamp-1)],[lamp(0),lamp(nlamp-1)],lamp)
nlam   = 4000
dlam   = 1/1d3
lam    = (DINDGEN(nlam)-(nlam-1)/2.)*dlam
dlamdv = 6173.3433d0/3d8
ntest  = 101 
dvtest = 130.
vtest  = dvtest*(DINDGEN(ntest)-(ntest-1)/2)
lines  = DBLARR(nlam,ntest)
dlinesdv  = DBLARR(nlam,ntest)
q1     = [MAX(line),line,MAX(line)]
q2     = [MIN(lamp)-10,lamp,MAX(lamp)+10] 
FOR i  = 0,ntest-1 DO lines(*,i) = INTERPOL(q1,q2,lam+vtest(i)*dlamdv)

erreur = FLTARR(ntest)
FOR i  = 0,ntest-1 DO BEGIN
    weightsp=FLTARR(nlam)+1.0
     A=[1.05,0.55,0.0]
    ;A=[1.05,0.55,0.0,0.06]

    ; from -0.1725 A to +0.1725 A by steps of 0.0685 A
    wavel  = [-0.1725,-0.1035,-0.0345,+0.0345,+0.1035,+0.1725]
    inten  = [lines[1827,i],lines[1896,i],lines[1965,i],lines[2034,i],lines[2103,i],lines[2172,i]]

    resp   = CURVEFIT(wavel,inten,weightsp,A,FUNCTION_NAME='fit3',TOL=1.d-8,ITMAX=800,/DOUBLE,STATUS=stat,CHISQ=chi2)
    PLOT,lam,lines[*,i]
    OPLOT,wavel,resp,thick=2,linestyle=2
    IF stat EQ 0 THEN erreur[i]=(-A[2]/dlamdv-vtest[i]) ELSE erreur[i]=0.0
    PRINT,'STATUS=',stat,i
ENDFOR    

READ,pause

END
