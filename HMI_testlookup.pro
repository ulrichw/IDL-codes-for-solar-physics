; program to test the implementation of lookup tables in
; the observable code velocity.c


PRO HMI_testlookup

d=FLOAT(readfits('/SUM4/D8762668/S00000/file.fits')) ; look-up table

map=FLTARR(4096,4096)
map2=map

ratio = 4096l/256l
ntest=501
axist=[1002l,256l,256l]

FOR i=0,4095 DO BEGIN
    PRINT,i
    FOR j=0,4095 DO BEGIN


	      ;bilinear interpolation of the look-up tables at pixel (x,y)
	      ;NB: it depends on how the filtergrams rebinning from 4096*4096 to axist[1]*axist[2] was done in phasemaps.c
	      ;find the 4 neighbors of (iii,jjj) on the grid of the look-up tables
	      ;and deal with boundary problems
	      x0 = (i/ratio);
	      y0 = (j/ratio);
	      x1 = x0+1;
	      y1 = y0+1;
	      x  = double(i MOD ratio)/ double(ratio)+double(x0);
	      y  = double(j MOD ratio)/ double(ratio)+double(y0);
	      if(x1 GE 256) THEN BEGIN
		  x0 = x0-1;
		  x1 = x1-1;
              ENDIF
	      if(y1 GE 256) THEN BEGIN
		
		  y0 = y0-1;
		  y1 = y1-1;
              ENDIF

	      ;perform the bilinear interpolation
	     ;for(v=0,ntest-1) DO BEGIN
	;{
              v=100l

		  RR1=d[v+x0*axist[0]+y0*axist[0]*axist[1]]*(double(x1)-x)+d[v+x1*axist[0]+y0*axist[0]*axist[1]]*(x-double(x0));
 		  RR2=d[v+x0*axist[0]+y1*axist[0]*axist[1]]*(double(x1)-x)+d[v+x1*axist[0]+y1*axist[0]*axist[1]]*(x-double(x0));
		  RR1 =(double(y1)-y)*RR1+(y-double(y0))*RR2;
 		  map[i,j]=RR1;  //for 1st Fourier coefficient
	

		  RR1=d[v+ntest+x0*axist[0]+y0*axist[0]*axist[1]]*(double(x1)-x)+d[v+ntest+x1*axist[0]+y0*axist[0]*axist[1]]*(x-double(x0));
 		  RR2=d[v+ntest+x0*axist[0]+y1*axist[0]*axist[1]]*(double(x1)-x)+d[v+ntest+x1*axist[0]+y1*axist[0]*axist[1]]*(x-double(x0));
		  RR1 =(double(y1)-y)*RR1+(y-double(y0))*RR2;
		  map2[i,j]=RR1; //for 2nd Fourier coefficient
	;}

    ENDFOR
ENDFOR

READ,pause

END
