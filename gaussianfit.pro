; FIT THE SOLAR LINE BY A GAUSSIAN PROFILE

PRO gaussianfit,X,A,F,pder

; A[0] is central wavelength
; A[1] is width
; A[2] is depth
; A[3] is continuum
F= A[3]-A[2]*exp(-(X-A[0])^2.d0/A[1]);
	
pder=[ [-A[2]*exp(-(X-A[0])^2.d0/A[1])*2.d0*(X-A[0])/A[1]],[-A[2]*(X-A[0])^2.d0/A[1]^2.d0*exp(-(X-A[0])^2.d0/A[1])],[-exp(-(X-A[0])^2.d0/A[1])],[FLTARR(N_ELEMENTS(X))+1.0] ] 


END


PRO doublegaussianfit,X,A,F,pder

; A[0] is central wavelength
; A[1] is width
; A[2] is depth
; A[3] is continuum
F= A[3]-A[2]*exp(-(X-A[0])^2.d0/A[1])-A[6]*exp(-(X-A[4])^2.d0/A[5]);
	
pder=[ [-A[2]*exp(-(X-A[0])^2.d0/A[1])*2.d0*(X-A[0])/A[1]],[-A[2]*(X-A[0])^2.d0/A[1]^2.d0*exp(-(X-A[0])^2.d0/A[1])],[-exp(-(X-A[0])^2.d0/A[1])],[FLTARR(N_ELEMENTS(X))+1.0],[-A[6]*exp(-(X-A[4])^2.d0/A[5])*2.d0*(X-A[4])/A[5]],[-A[6]*(X-A[4])^2.d0/A[5]^2.d0*exp(-(X-A[4])^2.d0/A[5])],[-exp(-(X-A[4])^2.d0/A[5])]] 


END

PRO linefit,X,A,F,pder

; A[0] is central wavelength
; A[1] is width
; A[2] is depth
; A[3] is continuum
F= A[2]*exp(-(X-A[0])^2.d0/A[1]);
	
pder=[ [A[2]*exp(-(X-A[0])^2.d0/A[1])*2.d0*(X-A[0])/A[1]],[A[2]*(X-A[0])^2.d0/A[1]^2.d0*exp(-(X-A[0])^2.d0/A[1])],[exp(-(X-A[0])^2.d0/A[1])]] 


END


;openr,1,"/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceFeLine.bin"
; d=fltarr(7000)
;readu,1,d
; close,1
; openr,1,"/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/ReferenceWavelength.bin"
; l=fltarr(7000)
;readu,1,l
;close,1
;plot,l,d,xst=1
;A=[0.,0.004,0.4,1.0]
;aa=where(l ge -0.172 and l le 0.172,na)
;weights=fltarr(7000)+1.0
; resp      = CURVEFIT(l,d,weights,A,FUNCTION_NAME='gaussianfit',TOL=1.d-9,ITMAX=5000)       
;plot,l,resp,col=180
