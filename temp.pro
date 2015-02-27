PRO gaussianfit,X,A,F,pder

; A[0] is central wavelength
; A[1] is width
; A[2] is depth
; A[3] is continuum
F= A[2]*exp(-(X-A[0])^2.d0/A[1]);

pder=[ [A[2]*exp(-(X-A[0])^2.d0/A[1])*2.d0*(X-A[0])/A[1]],[A[2]*(X-A[0])^2.d0/A[1]^2.d0*exp(-(X-A[0])^2.d0/A[1])],[exp(-(X-A[0])^2.d0/A[1])]] 

END

PRO temp

OPENR,1,'/home/couvidat/cvs/JSOC/temp'
d=strarr(16754)
readf,1,d
close,1

yo=fltarr(16754)
yo2=yo
for i=0,16753 do begin
    temp=fitsio_read_image(STRTRIM(d[i]+'/magnetogram.fits',1),h)
    res=histogram(temp[1500:2500,1800:2200],binsize=0.1,location=x,min=-70,max=70)
    A=[0.,10.,1800.]
    weights=fltarr(N_ELEMENTS(x))+1.0
    resp = CURVEFIT(x,res,weights,A,FUNCTION_NAME='gaussianfit',TOL=1.d-9,ITMAX=5000)   
    yo[i]=SQRT(A[1])
    yo2[i]=A[0]
    if(i MOD 10 EQ 0) then print,i,' ',d[i],' ',yo[i]
endfor

plot,yo,xst=1,yst=1
SAVE,yo,yo2,d,file='temp.bin'
READ,pause
END
