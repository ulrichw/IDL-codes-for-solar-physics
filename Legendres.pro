; Expansion of y=f(x) into a Legendre Polynomial series

PRO legendres,x,y,y2

nL  = 18;Polynomial order
leg = DBLARR(nL)

x2  = x/MAX(x)

FOR i=0,nL-1 DO leg[i]=INT_TABULATED(x2,LEGENDRE(x2,i,/DOUBLE)*y)*(2.d0*i+1.d0)/2.d0

y2  = x2
y2[*]= 0.d0

FOR i=0,nL-1 DO y2 = y2+leg[i]*LEGENDRE(x2,i,/DOUBLE)

END
