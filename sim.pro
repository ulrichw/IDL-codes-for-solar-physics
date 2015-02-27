; Program to set filter profiles and calculate expected noise performance.

; Set line to look at (0=6768, 1=6173)
iline=1

if (iline eq 0) then begin

; Wavelength of line in A
  lam0=6768.0

; Effective g
  geff=1.426

; Line depth and width relative to 6768
  fdepth=1.0
  fwidth=1.0

  title='Ni6768'

end

if (iline eq 1) then begin

; Wavelength of line in A
  lam0=6173.0

; Effective g
  geff=2.500

; Line depth and width relative to 6768 from Yang's web page
  fdepth=0.623/0.53
; fdepth=0.623/0.53*0.55/0.7
  fwidth=0.10/0.12

  title='Fe6173'

end

; Conversion factor from field to velocity
dvdb=1.35*(lam0/6768.)*(geff/1.426)

; Set full well depth (or actually the maximum exposure level) in electrons
well=1.25d5

; Set number of tuning positions
ntune=6
; All positions for the purpose of plotting.
ntunep=11
; Set number of tuning positions over wnarrow. (1/spacing)
; Don't mess with this unless you know what you are doing.
inttune=5./2

; Set period of narrowest Michelson
wnarrow=0.189D0*lam0/6768.0
; Set widths of tunable elements in units of wnarrow
; [1,2] is present MDI. [1,2,4] is with a third tunable element
mlist=[1,2,4]
nmich=n_elements(mlist)
wmich=mlist*wnarrow
;wmich(*)=1.05*wmich(*)
; Set mich contrasts
michcont=dblarr(nmich)+1.0D0
;michcont(*)=0.95D0
; Set tuning positions
dtune=wnarrow/inttune
tune=dblarr(nmich,ntune)
for i=0,nmich-1 do tune(i,*)=(-(ntune-1)/2.+dindgen(ntune))*dtune
;i=2 & tune(i,*)=tune(i,*)+0.167*wmich(i)
;tune(2,*)=0
tunep=(-(ntunep-1)/2.+dindgen(ntunep))*dtune

; Set period of narrowest Lyot element
lyotw0=8*wnarrow
;nlyot=4
;lyotw=dblarr(nlyot)
;lyotw(0)=lyotw0
;for i=1,nlyot-1 do lyotw(i)=lyotw(i-1)*2
lyotw=lyotw0*[1,2,4,8]
;lyotw=lyotw0*[1,1,2,2,4,8] ; MDI
nlyot=n_elements(lyotw)
;lyotw(*)=1.05*lyotw(*)
; Set lyot contrasts
lyotcont=dblarr(nlyot)+1.0D0
;lyotcont=[0.9396,0.9609,0.9764,1]
;lyotcont(*)=0.95D0
;lyotcont(3)=0.95D0
lyottune=dblarr(nlyot)
;lyottune=[0.52,4.58,-0.07,11.22]/360*lyotw
;lyottune(*)=0.05*lyotw(*)

; Identical exposures?
fixexp=1

; Number of line parameters
npar=3
 
; Set wavelength points for line profile
nlamp=1223
lamp=-0.3D0+dindgen(nlamp)*0.6D/(nlamp-1)

if 0 then begin
; The following lines read the line profile from Yang's file
  restore,'/scr305/schou/iquv.sav'
  profile=iquv0000(*,0:nlamp-1,*)
  line=reform(profile(0,*,0))
end

if 0 then begin
; Rogers files
  roger=readx('mu.60.Fe.txt')
  nroger=n_elements(roger)/2
  rlam=reform(roger(0,*))
  rp=reform(roger(1,*))
; Make continuum 1
  rp1=rp/interpol([rp(0),rp(nroger-1)],[rlam(0),rlam(nroger-1)],rlam)
  line=interpol(rp,rlam,lamp)
  line=line/interpol([line(0),line(nlamp-1)],[lamp(0),lamp(nlamp-1)],lamp)
  dlinedd=1-line
  dlinedi=line
end

if 1 then begin
; The following lines use Ted's approximate equation
  t=-2.7726*(lamp/0.125/fwidth)^2
  t = -25. + (t+25.)*(t gt -25.)
  line=1-0.6*exp(t)*fdepth
  dlinedd=exp(t)*fdepth
  dlinedi=line
end

; Set wavelength mesh for filter calculations
nlam=2880*5*2
dlam=1/1d3
lam=(dindgen(nlam)-(nlam-1)/2.)*dlam

; Convert wavelength to velocity
dlamdv=lam0/3d8
dvdlam=1/dlamdv
vel=lam*dvdlam
dvel=dlam*dvdlam
dvtune=dvdlam*dtune
vtune=dvdlam*tune

; Number of test filter profiles to test sensitivities
ntestx=1
iseed=1000l
tunex=dblarr(nmich,ntune,ntestx)
for i=0,ntune-1 do for j=0,nmich-1 do tunex(j,i,*)=tune(j,i)
; Make 1mA errors
; Do all elements
;for i=0,ntune-1 do for j=0,nmich-1 do tunex(j,i,*)=tune(j,i)+randomn(iseed,ntestx)*0.001
; Do one element
;for i=0,ntune-1 do tunex(2,i,*)=tune(j,i)+randomn(iseed,ntestx)*0.001

; Make 1 arcmin errors. 1 arcmin is 1/60/90 of 90 degrees
;for i=0,ntune-1 do for j=0,nmich-1 do tunex(j,i,*)=tune(j,i)+randomn(iseed,ntestx)*wmich(j)/60/90
; One element at a time
;j=2
;for i=0,ntune-1 do tunex(j,i,*)=tune(j,i)+randomn(iseed,ntestx)*wmich(j)/60/90

blocker=dblarr(nlam)
; Simple boxcar blocker to arbitrary null in Lyot
blocker(where(abs(lam) le lyotw0))=1

bfac=1.6
q=readx('smooth.txt')
;q=readx('nomatch2.txt')
blocker=interpol(q(1,*)/100,bfac*(q(0,*)-6170),lam,/spline)
;bfac=1.
;b11=readfits('/home2/schou/calib/sdo/inverse/blocker11.fits')
;blocker11=interpol(b11(*,1)/max(b11(*,1)),bfac*(b11(*,0)-6169.7),lam)
;blocker=blocker11
blockerx=dblarr(nlam,ntestx)
for ix=0,ntestx-1 do blockerx(*,ix)=blocker
; Put in blocker fringes
for ix=0,ntestx/2-1 do begin
  blockerx(*,2*ix)=blocker*(1+0.01*cos(2*!dpi*lam*ix/40.0))
  blockerx(*,2*ix+1)=blocker*(1+0.01*sin(2*!dpi*lam*ix/40.0))
end

clyot=2*!dpi/lyotw
; Set Lyot filter profile. No Doubled elements.
lyot=blocker
for i=0,nlyot-1 do lyot=lyot*(1+lyotcont(i)*cos(2*!dpi*(lam+lyottune(i))/lyotw(i)))/2
lyotx=blockerx
;for ix=0,ntestx-1 do for i=0,nlyot-1 do lyotx(*,ix)=lyotx(*,ix)*(1+lyotcont(i)*cos(2*!dpi*(lam+lyottune(i))/lyotw(i)))/2
for ix=0,ntestx-1 do for i=0,nlyot-1 do lyotx(*,ix)=lyotx(*,ix)*(1+lyotcont(i)*cos(clyot(i)*(lam+lyottune(i))))/2

; MDI Lyot. Assume 1,1,2,2,4,8 periods.
;t1=!dpi*lam/lyotw0
;lyotm=(cos(t1)^2*cos(t1/2)^2*cos(t1/4)*cos(t1/8))^2

cmich=2*!dpi/wmich
filters=dblarr(nlam,ntune)
filtersx=dblarr(nlam,ntune,ntestx)
for itune=0,ntune-1 do begin
  filters(*,itune)=lyot
  for ix=0,ntestx-1 do filtersx(*,itune,ix)=lyotx(*,ix)
; for i=0,nmich-1 do filters(*,itune)=filters(*,itune)*(1+michcont(i)*cos(2*!dpi*(lam+tune(itune))/(wmich(i))))/2
  for i=0,nmich-1 do filters(*,itune)=filters(*,itune)*(1+michcont(i)*cos(cmich(i)*(lam+tune(i,itune))))/2
  for ix=0,ntestx-1 do for i=0,nmich-1 do filtersx(*,itune,ix)=filtersx(*,itune,ix)*(1+michcont(i)*cos(cmich(i)*(lam+tunex(i,itune,ix))))
end
filtersx=filtersx/(2^nmich)

; Adjust filters if variable exposures are allowed
if (fixexp eq 0) then begin
; Find total filter transmissions
  trans=dblarr(ntune)
  for itune=0,ntune-1 do trans(itune)=total(filters(*,itune))
; Keep max transmission at roughly 1
  trans=trans/max(trans)
  for itune=0,ntune-1 do filters(*,itune)=filters(*,itune)/trans(itune)
  for itune=0,ntune-1 do filtersx(*,itune,*)=filtersx(*,itune,*)/trans(itune)
end

; Set profiles for plotting.
filtersp=dblarr(nlam,ntunep)
for itune=0,ntunep-1 do begin
  filtersp(*,itune)=lyot
; for i=0,nmich-1 do filtersp(*,itune)=filtersp(*,itune)*cos(!dpi*(lam+tunep(itune))/(wmich(i)))^2
  for i=0,nmich-1 do filtersp(*,itune)=filtersp(*,itune)*(1+michcont(i)*cos(2*!dpi*(lam+tunep(itune))/(wmich(i))))/2
end

; Adjust filters if variable exposures are allowed
if (fixexp eq 0) then begin
; Find total filter transmissions
  transp=dblarr(ntunep)
  for itune=0,ntunep-1 do transp(itune)=total(filtersp(*,itune))
; Keep max transmission at roughly 1
  transp=transp/max(transp)
  for itune=0,ntunep-1 do filtersp(*,itune)=filtersp(*,itune)/transp(itune)
end

; Set test velocities
ntest=101
;ntest=33
dvtest=500.
vtest=dvtest*(dindgen(ntest)-(ntest-1)/2)
lines=dblarr(nlam,ntest)
dlinesdv=dblarr(nlam,ntest)
; Calculate Dopler shifted line profiles
; Following two lines take care of continuum outside of where profile is given
q1=[max(line),line,max(line)]
q2=[min(lamp)-10,lamp,max(lamp)+10]
for i=0,ntest-1 do lines(*,i)=interpol(q1,q2,lam+vtest(i)*dlamdv)
; Get derivative with respect to v
dv=100D0
for i=0,ntest-1 do dlinesdv(*,i)=(interpol(q1,q2,lam+(vtest(i)+dv)*dlamdv)-interpol(q1,q2,lam+(vtest(i)-dv)*dlamdv))/2/dv
; Get derivative with respect to depth
dlinesdd=dblarr(nlam,ntest)
q1=[max(line),dlinedd,max(line)]
for i=0,ntest-1 do dlinesdd(*,i)=interpol(q1,q2,lam+vtest(i)*dlamdv)
; Get derivative with respect to overall intensity
dlinesdi=dblarr(nlam,ntest)
q1=[max(line),dlinedi,max(line)]
for i=0,ntest-1 do dlinesdi(*,i)=interpol(q1,q2,lam+vtest(i)*dlamdv)

; Calculate sigma for a single velocity determination from linearized
; inverse without fitting for depth or intensity.
; Calculate intensities
inten=transpose(filters)#lines
; And derivatives
dintendv=transpose(filters)#dlinesdv
; Calculate "exposure" time
texp=well/max(inten)
sigma=dblarr(ntest)
for itest=0,ntest-1 do begin
; Calculate intensities
  i0=texp*inten(*,itest)
; Calculate derivative with respect to velocity
  sens=texp*dintendv(*,itest)
; Calculate least squares estimate of noise. (This is just a weighted mean.)
  sigma(itest)=1/sqrt(total(sens^2/i0))
end

intenx=dblarr(ntune,ntest,ntestx)
for ix=0,ntestx-1 do intenx(*,*,ix)=transpose(filtersx(*,*,ix))#lines

; Calculate sigma(1) for a single velocity determination from linearized
; inverse also fitting for depth and intensity.
dintendd=transpose(filters)#dlinesdd
dintendi=transpose(filters)#dlinesdi
sigma1=dblarr(ntest)
a=dblarr(ntune,npar)
for itest=0,ntest-1 do begin
; Calculate intensities
  i0=texp*inten(*,itest)
; Calculate inverse sigmas
  si0=1/sqrt(i0)
; Calculate derivative with respect to velocity
  a(*,0)=texp*dintendv(*,itest)*si0
  a(*,1)=texp*dintendi(*,itest)*si0
  a(*,2)=texp*dintendd(*,itest)*si0
  cov=invert(transpose(a)#a)
; Calculate least squares estimate of noise. (This is just a weighted mean.)
  sigma1(itest)=sqrt(cov(0,0))
end

; Calculate sigma(v) for simultaneous fit of LCP and RCP assuming same
; intensity and depth.
vt1=dblarr(ntest,ntest)
vt2=dblarr(ntest,ntest)
sigmas=dblarr(npar+1,ntest,ntest)
a=dblarr(2*ntune,npar+1)
for it1=0,ntest-1 do begin
  for it2=0,ntest-1 do begin
    vt1(it1,it2)=vtest(it1)
    vt2(it1,it2)=vtest(it2)
;   Calculate intensities
    i0=texp*[inten(*,it1),inten(*,it2)]
;   Calculate inverse sigmas
    si0=1/sqrt(i0)
;   Calculate derivative with respect to velocity
    a(*,0)=texp*[dintendv(*,it1),0*dintendv(*,it2)]*si0
    a(*,1)=texp*[0*dintendv(*,it1),dintendv(*,it2)]*si0
    a(*,2)=texp*[dintendi(*,it1),dintendi(*,it2)]*si0
    a(*,3)=texp*[dintendd(*,it1),dintendd(*,it2)]*si0
    cov=invert(transpose(a)#a)
;   Calculate least squares estimate of noise. (This is just a weighted mean.)
    for i=0,npar do sigmas(i,it1,it2)=sqrt(cov(i,i))
  end
end
vta=(vt1+vt2)/2
vts=(vt2-vt1)/2
vt1r=reform(vt1,ntest^2)
vt2r=reform(vt2,ntest^2)
vtar=reform(vta,ntest^2)
vtsr=reform(vts,ntest^2)
sigmasr=reform(sigmas,npar+1,ntest^2)
sigmavr=sqrt(sigmasr(0,*)^2+sigmasr(1,*)^2)/2
sigmav=reform(sigmavr,ntest,ntest)

x=2*!dpi*(-(ntune-1)/2.+dindgen(ntune))/ntune
c0=reform(cos(0*x)#inten)
c1=reform(cos(1*x)#inten)
s1=reform(sin(1*x)#inten)
c2=reform(cos(2*x)#inten)
s2=reform(sin(2*x)#inten)
pv1=dvtune*inttune*2
pv2=dvtune*inttune
phi1=atan(-s1,-c1)
phi2=atan(-s2,-c2)
vel1=phi1*pv1/2/!dpi
vel2=phi2*pv2/2/!dpi
vel1a=(vel1-vtest+10.5*pv1) mod pv1-pv1/2+vtest
vel2a=(vel2-vtest+10.5*pv2) mod pv2-pv2/2+vtest

nt=10000
sigmat1=fltarr(ntest)
sigmat2=fltarr(ntest)
sigmat1a=fltarr(ntest) ; Sigma from first phase velocity
sigmat2a=fltarr(ntest) ; Sigma form second phase velocity
sigmat3a=fltarr(ntest) ; Sigma from crudely weighted average velocity
sigmat4a=fltarr(ntest) ; Sigma from properly weighted average velocity
itest=50
for itest=0,ntest-1 do begin
  vt=vtest(itest)
  i0=inten(*,itest)*texp
  iseed=1000l
  intent=dblarr(ntune,nt)
  for i=0,ntune-1 do intent(i,*)=i0(i)+sqrt(i0(i))*randomn(iseed,nt)
  c0=reform(cos(0*x)#intent)
  c1=reform(cos(1*x)#intent)
  s1=reform(sin(1*x)#intent)
  c2=reform(cos(2*x)#intent)
  s2=reform(sin(2*x)#intent)
  
  phi1t=atan(-s1,-c1)
  phi2t=atan(-s2,-c2)
  vel1t=phi1t*pv1/2/!dpi
  vel2t=phi2t*pv2/2/!dpi
  vel1at=(vel1t-vt+10.5*pv1) mod pv1-pv1/2+vt ; Fix phase ambiguity
  vel2at=(vel2t-vt+10.5*pv2) mod pv2-pv2/2+vt
  sigmat1(itest)=sqrt(rebin((vel1at-vel1a(itest))^2,1))
  sigmat2(itest)=sqrt(rebin((vel2at-vel2a(itest))^2,1))
  v1t=interpol(vtest,vel1a,vel1at,/spline) ; Put through lookup table
  v2t=interpol(vtest,vel2a,vel2at,/spline)
  rest1=v1t-vtest(itest) ; Calculate residual
  rest2=v2t-vtest(itest)
  vart1=mean(rest1^2) ; Calculate variance
  vart2=mean(rest2^2)
  covt12=mean(rest1*rest2) ; Calculate covariance
  sigmat1a(itest)=sqrt(vart1)
  sigmat2a(itest)=sqrt(vart2)
; Weigths for weighted average assuming uncorrelated errors
  ct1=vart2/(vart1+vart2)
  ct2=1-ct1
  rest3=ct1*rest1+ct2*rest2 ; Make linear combination
  vart3=mean(rest3^2)
  sigmat3a(itest)=sqrt(vart3)
; Weigths for weighted average assuming ncorrelated errors
  ct1x=(vart2-covt12)/(vart1+vart2-2*covt12)
  ct2x=1-ct1x
  rest4=ct1x*rest1+ct2x*rest2 ; Make linear combination
  vart4=mean(rest4^2)
  sigmat4a(itest)=sqrt(vart4)
end

;Need to finish use of filters1 to calculate errors

vel1x=dblarr(ntest,ntestx)
vel2x=dblarr(ntest,ntestx)
for ix=0,ntestx-1 do begin
  intent=intenx(*,*,ix)
  c0=reform(cos(0*x)#intent)
  c1=reform(cos(1*x)#intent)
  s1=reform(sin(1*x)#intent)
  c2=reform(cos(2*x)#intent)
  s2=reform(sin(2*x)#intent)
  
  vt=vtest
  phi1t=atan(-s1,-c1)
  phi2t=atan(-s2,-c2)
  vel1t=phi1t*pv1/2/!dpi
  vel2t=phi2t*pv2/2/!dpi
  vel1at=(vel1t-vt+10.5*pv1) mod pv1-pv1/2+vt
  vel2at=(vel2t-vt+10.5*pv2) mod pv2-pv2/2+vt
  vel1x(*,ix)=interpol(vtest,vel1a,vel1at)
  vel2x(*,ix)=interpol(vtest,vel2a,vel2at)
end
err1x=vel1x-rebin(reform(vtest,ntest,1),ntest,ntestx)
err2x=vel2x-rebin(reform(vtest,ntest,1),ntest,ntestx)
sig1x=sqrt(rebin((err1x)^2,ntest,1))
sig2x=sqrt(rebin((err2x)^2,ntest,1))
w5=where(abs(vtest) le 5000)
rms1x=sqrt(rebin(err1x(w5,*)^2,1,ntestx))

vi0=findgen(ntest)
; This is just the diagonal elements. ii=lindgen(ntest) & sig0=sigmav(ii,ii)
sig0=interpolate(sigmav,vi0,vi0)
vix=vi0-dvdb*1000./dvtest
viy=vi0+dvdb*1000./dvtest
sig1000=interpolate(sigmav,vix,viy)
vix=vi0-dvdb*2000./dvtest
viy=vi0+dvdb*2000./dvtest
sig2000=interpolate(sigmav,vix,viy)
vix=vi0-dvdb*3000./dvtest
viy=vi0+dvdb*3000./dvtest
sig3000=interpolate(sigmav,vix,viy)
vix=vi0-dvdb*4000./dvtest
viy=vi0+dvdb*4000./dvtest
sig4000=interpolate(sigmav,vix,viy)
stop

ps,file='sigmav.ps'
!p.charsize=1.5
xyouts,0,0,'!17',/device

plot,vtest,sig0,xrange=[0,8000],yrange=[0,30],xtitle='V (m/s)',ytitle='!7r!17!S!DV!R (m/s)',title=title,xstyle=1
;oplot,vtest,sig2000,linestyle=2
;oplot,vtar(w3),sigmavr(w3),linestyle=2
;oplot,vtest,sigma1/sqrt(2),linestyle=2 ; Divide by sqrt(2) for back and forth.
;oplot,vtest,1/sqrt((1/sigmat1a^2+1/sigmat2a^2))/sqrt(2),linestyle=3
oplot,vtest,sigmat1a/sqrt(2),linestyle=1
oplot,vtest,sigmat2a/sqrt(2),linestyle=2
oplot,vtest,sigmat4a/sqrt(2),linestyle=3


psend1
xwin
stop

psc,file='filtersc.ps'
!p.charsize=1.5
xyouts,0,0,'!17',/device
jscolors

xr=[-0.5,0.5]
colx=[1,2,0,3,5]
plot,lam,filters(*,0),xrange=xr,xstyle=9,/nodata,xtitle='!7dk!17(A)',ytitle='Transmission',yrange=[0,1.1],ystyle=1,ymargin=[4,4]
axis,xaxis=1,xrange=xr*dvdlam/1000,xstyle=1,xtitle='!7d!17V (km/s)'
for i=0,4 do oplot,lam,filters(*,i),color=colx(i)
oplot,lam,lines(*,ntest/2),linestyle=1
oplot,lam,filtersp(*,0),color=colx(2),linestyle=2

psend1
xwin

end


