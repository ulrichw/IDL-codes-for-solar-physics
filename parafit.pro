; TO FIT A PARABOLA FOR THE ANGULAR DEPENDENCE TEST

PRO parafit

param = DBLARR(3)
param[0] = 2.d0
param[1] = 2.d0
param[2] = 2.d0
seed = 1l
ntest   = 1000
err = DBLARR(3,ntest)
inten = DBLARR(3)

angle = DINDGEN(1000)/999.-0.5d0
parab = param[0] + param[1]*angle + param[2]*angle^2.d0

plot,angle,parab

FOR i=0,ntest-1 DO BEGIN


inten[0] = parab[250] ; at -0.2507 degrees
inten[1] = parab[500]
inten[2] = parab[749]

inten  = inten * (RANDOMU(seed,3)*0.2+0.9)

res   = POLY_FIT([angle[250],angle[500],angle[749]],inten,2,yfit=y)

oplot,angle,res[0]+res[1]*angle+res[2]*angle^2.d0,linestyle=2

err[0,i] = (res[0] - param[0])/param[0]
err[1,i] = (res[1] - param[1])/param[1]
err[2,i] = (res[2] - param[2])/param[2]

ENDFOR


PRINT,'ERRORS:'
PRINT,MEAN(err[0,*]),SIGMA(err[0,*])
PRINT,MEAN(err[1,*]),SIGMA(err[1,*])
PRINT,MEAN(err[2,*]),SIGMA(err[2,*])

END
