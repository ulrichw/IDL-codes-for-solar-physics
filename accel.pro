;restore,'pow.save'
restore,'hrpow.save'

nsample=6 ; Number of samples
dtsample=45D0 ; Time between samples
dtint=0.125D0 ; Time between interpolated samples
ntint=nint(dtsample/dtint)+1 ; Number of points to interpolate to

tint=dtint*dindgen(ntint)
errix=fltarr(ntint)
errjx=fltarr(ntint)
coeffix=dblarr(nsample,ntint)
coeffjx=dblarr(nsample,ntint)
for idt=0,ntint-1 do begin
;dtsample=idt+0.0
tshift=dtsample/2
tshift=idt*dtint
fnyq=1d6/dtsample/2

ntest=2500
dtest=10.d0
ftest=dtest*dindgen(ntest)
omtest=2*!dpi*ftest/1d6
hh=abs(hr12pow-650000)
hh(360)=(hh(359)+hh(361))/2
hh(720)=(hh(719)+hh(721))/2
ptest=interpol(hh(2:*),hr12fx(2:*),ftest,/spline)
; Cut at x times Nyquist
;ptest(where(ftest ge 0.90*fnyq))=0
; White noise to x times Nyquist
;ptest=0*ptest
;ptest(where(ftest lt 1.00*fnyq))=1d10
fc=1.1*fnyq
qc=ftest/fc
qn=ftest/fnyq+1d-300
;ptest=0*ptest
w=where(qc lt 1)
; White before PSF
;ptest(w)=1d10
;restore,'lapalma.save'
; La Palma radial
;ptest(w)=1e3*interpol(pb,kb,qn(w))
; La Palma at ky=0
;wx=where(kx0 ge 0)
;ptest(w)=1e3*interpol(px(wx),kx0(wx),qn(w))
; 1d PSF
;ptest(w)=ptest(w)*((1-qc(w))>0)^2
; 2d PSF @ ky=0
;ptest(w)=ptest(w)*(2/!pi*(acos(qc(w))-qc(w)*sqrt(1-qc(w)^2)))^2
; sinc^2 from pixel integration. Add 1d-300 to get sinc(0)=1
;ptest=ptest*(sin(qn*!dpi/2)/(qn*!dpi/2))^2
; Two point average
;ptest=ptest*cos(omtest*dtsample/2)^2

; Calculate autocorrelation
dta=1D0
nta=10000
fa=1d6/(dta*nta)*((indgen(nta)+nta/2) mod nta-nta/2)
pa=dblarr(nta)
w=where(abs(fa) le max(ftest))
pa(w)=interpol(ptest,ftest,abs(fa(w)))
q=double(fft(pa,1))
ta=dta*(lindgen(nta)-nta/2)
acor=q((nta+lindgen(nta)-nta/2) mod nta)/q(0)

c=dtest*2/n12x*12/1d6

errp=dblarr(101)
errs=dblarr(101)
erri=dblarr(101)
errj=dblarr(101)
errc=dblarr(101)
errp(0)=total(sqrt(c*total(ptest)))
errs(0)=total(sqrt(c*total(ptest)))
erri(0)=total(sqrt(c*total(ptest)))
errj(0)=total(sqrt(c*total(ptest)))
errc(0)=total(sqrt(c*total(ptest)))
;plot_io,ftest,ptest,yrange=[1d2,1d10]
;plot_io,ftest,ptest,/nodata,xrange=[0,1.5d4],yrange=[1e-10,10],ystyle=1
;for nsample=6,6,2 do begin
  if ((nsample mod 2) eq 0) then begin
    tsample=dtsample*(dindgen(nsample)-nsample/2+1)
  endif else begin
    tsample=dtsample*(dindgen(nsample)-nsample/2)
  end
  cx=nsample/2D0
  x=tsample/dtsample/cx
  x0=tshift/dtsample/cx
  
  a=dblarr(nsample,nsample)
  rh=dblarr(nsample)
  for i=0,nsample-1 do begin
    a(i,*)=x^i
    rh(i)=x0^i
  end
  a(0,*)=1/sqrt(nsample)
  rh(0,*)=1/sqrt(nsample)
  for i=1,nsample-1 do begin
    a1=x*a(i-1,*)
    rh1=x0*rh(i-1)
;   for j=(i-2)>0,i-1 do begin
    for j=0,i-1 do begin
      c1=total(a1*a(j,*))
      a1=a1-a(j,*)*c1
      rh1=rh1-rh(j)*c1
    end
    cn=1/sqrt(total(a1^2))
    a(i,*)=a1*cn
    rh(i)=rh1*cn
  end
  coeffp=invert(a)#rh
  
  cerr=cos(tshift*omtest)
  for i=0,nsample-1 do cerr=cerr-coeffp(i)*cos(tsample(i)*omtest)
  serr=sin(tshift*omtest)
  for i=0,nsample-1 do serr=serr-coeffp(i)*sin(tsample(i)*omtest)
  ptransp=cerr^2+serr^2
  errp(nsample)=total(sqrt(c*total(ptransp*ptest)))
; print,nsample,errp(nsample),total(coeffp^2)
  
  om0=dblarr(99,99)
  om0(0,0)=[3250,0]
  om0(1,0:1)=[2190,3670]
  om0(2,0:2)=[1200,3300,4000]
  om0(3,0:3)=[2000,3500,5000,6500]
  om0(3,0:3)=[1000,2000,3000,4000]
  om0(4,0:4)=[1000,2000,3000,4000,5000]
; om0(4,0:4)=[1000,2800,3500,3900,4100]
  om0(5,0:5)=[500,1000,2000,3000,4000,5000]
  for i=0,98 do om0(i,0:i)=(1+dindgen(i+1))/(i+2)*1d6/dtsample/2
  om0=om0*2*!dpi/1d6
  a=dblarr(nsample,nsample)
  rh=dblarr(nsample)
  a(nsample-1,*)=1
  rh(nsample-1)=1
  if ((nsample mod 2) eq 0) then begin
    a(nsample-2,*)=(tsample/dtsample)^1
    rh(nsample-2)=(tshift/dtsample)^1
  end
  for i=0,(nsample-1)/2-1 do begin
    a(2*i,*)=cos(tsample*om0((nsample-1)/2-1,i))
    a(2*i+1,*)=sin(tsample*om0((nsample-1)/2-1,i))
    rh(2*i)=cos(tshift*om0((nsample-1)/2-1,i))
    rh(2*i+1)=sin(tshift*om0((nsample-1)/2-1,i))
  end
  coeffs=invert(a)#rh
  
  cerr=cos(tshift*omtest)
  for i=0,nsample-1 do cerr=cerr-coeffs(i)*cos(tsample(i)*omtest)
  serr=sin(tshift*omtest)
  for i=0,nsample-1 do serr=serr-coeffs(i)*sin(tsample(i)*omtest)
  ptranss=cerr^2+serr^2
  errs(nsample)=total(sqrt(c*total(ptranss*ptest)))
; print,nsample,errs(nsample),total(coeffs^2)
  
; Do ideal case based on power spectrum
  a=dblarr(nsample,nsample)
  rh=dblarr(nsample)
  for i=0,nsample-1 do begin
    a(*,i)=interpol(acor,ta,tsample-tsample(i))
    rh(i)=interpol(acor,ta,tshift-tsample(i))
  end
  coeffi=invert(a)#rh

  cerr=cos(tshift*omtest)
  for i=0,nsample-1 do cerr=cerr-coeffi(i)*cos(tsample(i)*omtest)
  serr=sin(tshift*omtest)
  for i=0,nsample-1 do serr=serr-coeffi(i)*sin(tsample(i)*omtest)
  ptransi=cerr^2+serr^2
  erri(nsample)=total(sqrt(c*total(ptransi*ptest)))
  errix(idt)=total(sqrt(c*total(ptransi*ptest)))
; print,nsample,erri(nsample),total(coeffi^2)
; print,dtsample,tshift,nsample,erri(nsample),total(coeffi^2)

; Ideal case but forcing constants to be preserved
  ai=a
  rhi=rh
  a=dblarr(nsample+1,nsample+1)
  rh=dblarr(nsample+1)
  a(0:nsample-1,0:nsample-1)=ai
  a(0:nsample-1,nsample)=1
  a(nsample,0:nsample-1)=1
  rh(0:nsample-1)=rhi
  rh(nsample)=1
  coeffj=invert(a)#rh
  coeffj=coeffj(0:nsample-1)
; coeffj=[0.11010416,-0.22947822,0.61937406,0.61937406,-0.22947822,0.11010416]

  cerr=cos(tshift*omtest)
  for i=0,nsample-1 do cerr=cerr-coeffj(i)*cos(tsample(i)*omtest)
  serr=sin(tshift*omtest)
  for i=0,nsample-1 do serr=serr-coeffj(i)*sin(tsample(i)*omtest)
  ptransj=cerr^2+serr^2
  errj(nsample)=total(sqrt(c*total(ptransj*ptest)))
  errjx(idt)=total(sqrt(c*total(ptransj*ptest)))
; print,dtsample,tshift,nsample,errj(nsample),total(coeffj^2)

; Truncated sinc. Pretty crappy for 2d. Good for white.
  q=(tsample-tshift)/dtsample+1d-300
  coeffc=sin(q*!dpi)/(q*!dpi)
; Truncate by sinc. Pretty good for 2d! At least when power to Nyq.
; coeffc=coeffc*sin(q*!dpi*2/nsample)/(q*!dpi*2/nsample)
; Force sum=1. No diff for nsample even. Worse for nsample odd (tshift=1/2).
  coeffc=coeffc/total(coeffc)
  cerr=cos(tshift*omtest)
  for i=0,nsample-1 do cerr=cerr-coeffc(i)*cos(tsample(i)*omtest)
  serr=sin(tshift*omtest)
  for i=0,nsample-1 do serr=serr-coeffc(i)*sin(tsample(i)*omtest)
  ptransc=cerr^2+serr^2
  errc(nsample)=total(sqrt(c*total(ptransc*ptest)))
; print,nsample,errc(nsample),total(coeffc^2)

; oplot,ftest,ptransp*ptest,linestyle=(nsample-1)/2
; oplot,ftest,ptransc*ptest,linestyle=(nsample-1)/2
; oplot,ftest,ptransj,linestyle=(nsample-1)
; oplot,ftest,ptransp*ptest,linestyle=0
; oplot,ftest,ptranss*ptest,linestyle=1
; oplot,ftest,ptransi*ptest,linestyle=2
; oplot,ftest,ptransj*ptest,linestyle=2
;end ; nsample

coeffix(*,idt)=coeffi
coeffjx(*,idt)=coeffj
end ; idt

SAVE,coeffjx,tint,FILE='interpolation_coefficients.bin'

g=dblarr(ntest)
g(where(ftest ge 1d6/dtsample/2))=1
;print,sqrt(c*total(ptest)),sqrt(c*total(ptest*g))

;print,coeffj(0:5)

end
