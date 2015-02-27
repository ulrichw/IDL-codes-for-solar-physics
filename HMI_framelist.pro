; THIS PROGRAM IS DESIGNED TO CHECK WHICH FILTERGRAM ORDER (FRAMELIST)
; PROVIDES THE LEAST NOISE ON THE OBSERVABLES
; THE NOISE IS ASSUMED TO COME FROM 2 MAIN SOURCES:
; 1)  FROM ONE FILTERGRAM TO THE OTHER THE SOLAR SURFACE HAS ROTATED A
; LITTLE BIT
; 2) THE FILTERGRAMS NEED TO BE INTERPOLATED IN TIME
; 
; FOR THE TEMPORAL INTERPOLATION, I USE accel.pro PROVIDED BY JESPER SCHOU
; AND WHICH DOES A WIENER INTERPOLATION
; FOR THE SOLAR ROTATION I USE ~richard/idlwork/hmi/level1/rotation.pro
; PROVIDED BY RICHARD WACHTER
;-------------------------------------------------------------------------


; SHIFT IN THE FOURIER DOMAIN
FUNCTION decalage,f,lag,dx

; FUNCTION f MUST BE PERIODICAL (f[0]=f[N-1])
; f MUST HAVE AN EVEN NUMBER OF POINTS

N=N_ELEMENTS(f)
k=FINDGEN(N)/FLOAT(N)/dx
k[N/2+1:N-1]=-REVERSE(k[1:N/2-1])
f2=FFT((FFT(f)*exp(-2.0*!dpi*COMPLEX(0,1)*k*lag)),/INVERSE)

RETURN,FLOAT(f2)

END

; FOURIER INTERPOLATION
FUNCTION i_f, C, fac
; C is REAL (no imaginary part)
; C must be periodic (C[0]=C[N-1])
; C must have an even number of points
; the sampling rate must be uniform
; and initial and interpolated domain must cover THE SAME RANGE (NT=N0T0)                
                                                                
nt0 = DOUBLE(N_ELEMENTS(C))
nt  = nt0*fac
C_i = DBLARR(nt)
G   = DCOMPLEXARR(nt)
                                                                                
C2              = FFT(C,/DOUBLE)
G[0:nt0/2]      = C2[0:nt0/2]
G[nt0/2+1:nt/2] = DCOMPLEX(0,0)
G[nt/2+1:nt-1]  = CONJ(REVERSE(G[1:nt/2-1])) ; function is real
C_i=FFT(G,/INVERSE,/DOUBLE)
RETURN,DOUBLE(C_i)
                                                                                
END

; CODE FROM RICHARD WACHTER TO SIMULATE SOLAR ROTATION
FUNCTION rotation, radius, cent_x, cent_y, p0, b0, t, shift2=shift2, dist=dist

; radius is the radius of the sun in HMI pixels units (assuming 4096*4096 pixels)
; usually it is around 1970 pixels
; cent_x and cent_y are the coordinates of the center of the Sun
; in HMI pixel units (usually 2048,2048)
; p0 and B0 are the solar angles (B0 varies by +/- 7.2 degrees, p0
; is always close to 0
; t is the time at which the shift is computed
; the result is shift: an array giving the vector field
; the shift units are in HMI pixel sizes 
; to speed up calculation, n < size, but the result is always
; given for the HMI pixel size (4096*4096 array)
; if n < size that only means that the array is rebinned to
; perform calculations faster, but then the result is re-sized
; again to the full resolution

n=1024
size=4096
fac=size/n

dist=1.0 ; in A.U. distance from SDO to the Sun
 
shift=dblarr(n, n, 2)
shift2=shift

au=214.94272d0*dist             ;(1 astonomical unit in solar radii)
maxphi=asin(1.d0/au)

; read crop table
openr,1,'/home/richard/ccode/hmicomp/crop11'
readf,1,a
readf,1,a
readf,1,a
crp=fltarr(2,size)
readf,1,crp
close,1

pang_matrix=double([[1.0, 0.0, 0.0], [0.0, cos(p0), sin(p0)], [0.0, -sin(p0), cos(p0)]]) ; p-angle rotation matrix (rotate around x)
bang_matrix=double([[cos(b0), 0.0, -sin(b0)], [0.0, 1.0, 0.0], [sin(b0), 0.0, cos(b0)]]) ; b-angle rotation matrix (rotate around y)
  
pang_2d=[[cos(p0), sin(p0)], [-sin(p0), cos(p0)]] ;p-angle rotation
for i=0, n-1 do begin  
    for j=0, n-1 do begin  

        if ((j*fac+fac/2) ge crp[0,i*fac+fac/2] and (j*fac+fac/2) lt (crp[0,i*fac+fac/2]+crp[1,i*fac+fac/2])) then begin
            
            xy=transpose([((double(i)+0.5)*fac)-cent_x, ((double(j)+0.5)*fac-cent_y)]/radius) ; center of pixel coordinate
            inr=sqrt(xy[0]^2+xy[1]^2) 
            
            
            if keyword_set(dist) then begin
                
                if inr ge 1.0 then begin
                    xy=xy/inr
                    inr=1.0d0
                    phie=maxphi
                    ind=sin(-phie+!pi/2.) 
                endif else begin

                    phie=inr*maxphi
                    ind=sin(-phie+asin(au*sin(phie)))
                endelse
 
                xyd=xy/inr*ind
                xyz=transpose([sqrt(1.d0-ind^2), xyd[0], xyd[1]]) 
              
            endif else begin
                if inr ge 1.0 then begin
                    xy=xy/inr
                    inr=1.0 
                endif
                xyz=transpose([sqrt(1.-inr^2), xy[0], xy[1]])
            endelse
 

                                ; place point on surface of a sphere
  
 
                                ; p-angle rotation    
            xyp=pang_matrix ## xyz
                                ; b-angle rotation
            xyzb=bang_matrix ## xyp
      

      ; now rotate sphere according to rotation rate:
            sinlat=xyzb[2]
            Omega = 452.0 - 49.0*sinlat^2. - 84.0*sinlat^4. - 31.7 ; calculate rotation rate
            rotang=2.*!pi*Omega*1e-9*t ; rotation angle
            
            rot_matrix=double([[cos(rotang), -sin(rotang), 0.0], [sin(rotang), cos(rotang), 0.0], [0.0, 0.0, 1.0]]) ; rotate around z
            
            xyzr=rot_matrix ## xyzb

; now rotate back b-angle, p-angle, and project into plane

            xyzrbi=invert(bang_matrix, /double) ## xyzr ; b-angle
            xyzrpi=invert(pang_matrix, /double) ## xyzrbi ; p-angle

;;;;;;;; projection in plane from finite distance
            if keyword_set(dist) then begin
    
                singam=sqrt(xyzrpi[1]^2+xyzrpi[2]^2)
                cosgam=xyzrpi[0]
                
                phi=atan(singam/(au-cosgam))
                ing=phi/maxphi
                xyr=xyzrpi[[1,2]]/singam*ing
 
            endif else begin
                xyr=transpose(xyzrpi[[1,2]]) ;projection
            endelse


            shift[i,j,*]=(xyr-xy)*radius ; calculate difference

; alternative way: calculate derivative directly


            xys=pang_2d ## xyz[1:2]
            xym=transpose([xyz[0], xys[0], xys[1]])
            
            slat=(sin(b0)*xym[0]+cos(b0)*xym[2])
            omeg=2*!pi*1e-9*(452. - 49.0*slat^2 -84.0*slat^4 - 31.7)
            svec=transpose([-omeg*xym[1]*cos(b0), omeg*(xym[0]*cos(b0)-xym[2]*sin(b0)), omeg*xym[1]*sin(b0)]*t) ; !! change sign of b0
            
            if keyword_set(dist) then begin
   ;finite distance proj.
                singam=sqrt(svec[1]^2+svec[2]^2)
                cosgam=svec[0]
                
                phi=atan(singam/(au-cosgam))
                ing=phi/maxphi
                
                svecp=svec[[1,2]]/singam*ing*radius
            endif else begin
                svecp=svec[[1,2]]*radius
            endelse

            shift2[i,j,*]=invert(pang_2d) ## svecp

        endif
        

    endfor
endfor

return, shift

END


;-------------------------------------------------------------------------
;
;
; MAIN PROGRAM
;
;-------------------------------------------------------------------------


PRO HMI_framelist

lam0    = 6173.3433d0
ntune   = 6           ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0   ; Set number of tuning positions over wnarrow. (1/spacing)
ntime   = 600.*3.        ; duration of observation in minutes
nframe  = ntune*2     ; total number of filtergrams taken per observable (each filter is used twice
                      ; because of LCP and RCP) 
cadence = 45.d0
dt      = cadence/nframe ; time interval between 2 filtergrams, or integration time (?)
LANDE   = 1.0/(2.0*4.67e-5*0.000061733433*2.5*299792458.0) ; conversion from Doppler velocity to magnetic field
; see page 25 IPD (B=LANDE*(VLCP-VRCP))
dtime   = 60.d0/dt
nt      = ntime*dtime        ; total number of filtergrams
nt2     = nt/nframe   ; total number of observables
time    = FINDGEN(nt)*dt

; Parameters for the solar line
;------------------------------------------------------------------------

nlam       = 18000.;21500.
dlam       = 1.d0/1.75d3
lam        = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam


; Fe I LINE PROFILE FROM ROGER ULRICH'S WEBSITE
;-------------------------------------------------------------------------

OPENR,1,'Ulrich_Fe_0.txt'
roger  = FLTARR(2,98)
READF,1,roger
CLOSE,1
rlam   = REFORM(roger[0,*])
rp     = REFORM(roger[1,*])
rlam   = [-10.d0,rlam,-10.d0]
rp     = [1.d0,rp,1.d0]
line   = INTERPOL(rp,rlam,lam)
line   = line/INTERPOL([line[0],line[nlam-1]],[lam[0],lam[nlam-1]],lam)

; MODEL OF SOLAR LINE
;-------------------------------------------------------------------------
width  = 0.102d0/2.d0/SQRT(ALOG(2.d0)) ; w in A^2 (part added by seb)
depth  = 0.62d0
line  = 1.d0-depth*EXP(-lam^2.d0/width^2.d0)
dlinedv= depth*2.d0*lam/width^2.d0*EXP(-lam^2.d0/width^2.d0)

; to draw the derivative of the solar line
;SET_PLOT,'ps'
;device,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
;plot,lam,dlinedv,xrange=[-0.2,0.2],xst=1,charsize=1.5,tit='!17',xtit='Wavelength (A)',thick=2
;OPLOT,[-0.171,-0.171],[-10,10]
;OPLOT,[-0.103,-0.103],[-10,10]
;OPLOT,[-0.034,-0.034],[-10,10]
;OPLOT,[0.171,0.171],[-10,10]
;OPLOT,[0.103,0.103],[-10,10]
;OPLOT,[0.034,0.034],[-10,10]
;DEVICE,/CLOSE


; HMI FILTER TRANSMISSION PROFILES
;--------------------------------------------------------------------------
 
;phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
phase      = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
contrast   = [1.,1.,1.,1.,1.,1.,1.]
wmich      = [0.172d0-0.0010576d0,0.344d0-0.00207683d0,0.693d0+0.000483467d0]  
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0]

dtune      = wmich[0]/inttune  ; Set tuning position interval
tune       = DBLARR(3,ntune)
FOR  i     = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

dlamdv     = lam0/299792458.d0
dvdlam     = 1.d0/dlamdv
vel        = lam*dvdlam   ; DOPPLER shift in cm/s
dvel       = dlam*dvdlam
dvtune     = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune      = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s


; FRONT WINDOW AND BLOCKER FILTER TRANSMISSION PROFILES
RESTORE,'frontwindow.bin' ; front window
blocker    = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
;q          = READFITS('blocker11.fits')
;blocker    = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) 

; theoretical blocker filter
blocker    = blocker*EXP(-lam^2.d0/2.d0/(8.d0/2.354820d0)^2.d0)

; NON-TUNABLE PROFILE
lyot       = blocker
FOR i = 0,3 DO lyot = lyot*(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0

; TUNABLE PROFILE
cmich      = 2.d0*!dpi/wmich
filters    = DBLARR(nlam,ntune)
FOR itune  = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i  = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
ENDFOR

read,pause

; TO PLOT THE FILTERS
;SET_PLOT,'ps'
;device,file='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=14,/color,bits=24
;plot,lam,filters[*,0],xrange=[-0.75,0.75],xst=1,charsize=1.5,tit='!17',yrange=[0,1.1],yst=1,xtit='Wavelength (A)',thick=2
;xyouts,0.165,0.85,'F0',charsize=1.5
;oplot,lam,filters[*,1],col=40,thick=2
;xyouts,0.097,0.95,'F1',charsize=1.5,col=40
;oplot,lam,filters[*,2],col=80,thick=2
;xyouts,0.026,1.025,'F2',charsize=1.5,col=40
;oplot,lam,filters[*,3],col=120,thick=2
;xyouts,-0.060,1.025,'F3',charsize=1.5,col=120
;oplot,lam,filters[*,4],col=160,thick=2
; xyouts,-0.130,0.95,'F4',charsize=1.5,col=160
;oplot,lam,filters[*,5],col=200,thick=2
;xyouts,-0.200,0.85,'F5',charsize=1.5,col=200
;device,/CLOSE


; WE DEFINE THE FRAMELIST
;------------------------------------------------------------------------------------------

; the filters are in this order:
;--------------------------------------------------------
;filters[*,0] is F0 centered on +171 mA
;filters[*,1] is F1 centered on +103 mA
;filters[*,2] is F2 centered on +034 mA
;filters[*,3] is F3 centered on -034 mA
;filters[*,4] is F4 centered on -103 mA
;filters[*,5] is F5 centered on -171 mA
;--------------------------------------------------------


;tframelist='F0L,F0R,F1R,F1L,F2L,F2R,F3R,F3L,F4L,F4R,F5R,F5L'
;framelist=[0,6,7,1,2,8,9,3,4,10,11,5]

;tframelist='F0L,F0R,F1L,F1R,F2L,F2R,F3L,F3R,F4L,F4R,F5L,F5R'
;framelist=[0,6,1,7,2,8,3,9,4,10,5,11]

;tframelist='F1L,F1R,F2L,F2R,F3L,F3R,F4L,F4R,F5L,F5R,F0L,F0R'
;framelist=[1,7,2,8,3,9,4,10,5,11,0,6]

;tframelist='F2L,F2R,F3L,F3R,F4L,F4R,F5L,F5R,F0L,F0R,F1L,F1R'
;framelist=[2,8,3,9,4,10,5,11,0,6,1,7]

;tframelist='F3L,F3R,F4L,F4R,F5L,F5R,F0L,F0R,F1L,F1R,F2L,F2R'
;framelist=[3,9,4,10,5,11,0,6,1,7,2,8]

;tframelist='F3L,F3R,F4R,F4L,F5L,F5R,F0R,F0L,F1L,F1R,F2R,F2L'
;framelist=[3,9,10,4,5,11,6,0,1,7,8,2]

tframelist='F3L,F3R,F4R,F4L,F0L,F0R,F5R,F5L,F1L,F1R,F2R,F2L'
framelist=[3,9,10,4,0,6,11,5,1,7,8,2]

;tframelist='F3L,F3R,F4L,F4R,F0L,F0R,F5R,F5L,F1R,F1L,F2R,F2L'
;framelist=[3,9,4,10,0,6,11,5,7,1,8,2]

;tframelist='F4L,F4R,F5L,F5R, F0L,F0R,F1L,F1R,F2L,F2R,F3L,F3R'
;framelist=[4,10,5,11,0,6,1,7,2,8,3,9]

;tframelist='F5L,F5R,F0L,F0R,F1L,F1R,F2L,F2R,F3L,F3R,F4L,F4R'
;framelist=[5,11,0,6,1,7,2,8,3,9,4,10]

;tframelist='F2L,F2R,F1L,F1R,F0L,F0R,F5L,F5R,F4L,F4R,F3L,F3R'
;framelist=[2,8,1,7,0,6,5,11,4,10,3,9]

;tframelist='F2L,F2R,F3R,F3L,F4L,F4R,F1R,F1L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,4,10,7,1,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F1L,F1R,F4R,F4L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,1,7,10,4,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F5L,F5R,F0R,F0L,F1L,F1R,F4R,F4L'
;framelist=[2,8,9,3,5,11,6,0,1,7,10,4]

;tframelist='F0L,F0R,F5R,F5L,F1L,F1R,F4R,F4L,F2L,F2R,F3R,F3L'
;framelist=[0,6,11,5,1,7,10,4,2,8,9,3]

;tframelist='F0L,F0R,F2R,F2L,F4L,F4R,F1R,F1L,F3L,F3R,F5R,F5L'
;framelist=[0,6,8,2,4,10,7,1,3,9,11,5]

;tframelist='F2L,F2R,F3R,F3L,F1L,F1R,F4R,F4L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,1,7,10,4,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F4L,F4R,F1R,F1L,F0L,F0R,F5R,F5L'
;framelist=[2,8,9,3,4,10,7,1,0,6,11,5]

;tframelist='F2L,F2R,F3R,F3L,F0L,F0R,F5R,F5L,F1L,F1R,F4R,F4L'
;framelist=[2,8,9,3,0,6,11,5,1,7,10,4]

;tframelist='F3L,F3R,F2R,F2L,F0L,F0R,F5R,F5L,F1L,F1R,F4R,F4L'
;framelist=[3,9,8,2,0,6,11,5,1,7,10,4]

;tframelist='F3L,F3R,F2R,F2L,F0L,F0R,F5R,F5L,F4L,F4R,F1R,F1L'
;framelist=[3,9,8,2,0,6,11,5,4,10,7,1]

;tframelist='F3L,F3R,F2R,F2L,F5L,F5R,F0R,F0L,F4L,F4R,F1R,F1L'
;framelist=[3,9,8,2,5,11,6,0,4,10,7,1]

;tframelist='F3L,F3R,F2L,F2R,F5L,F5R,F0L,F0R,F4L,F4R,F1L,F1R'
;framelist=[3,9,2,8,5,11,0,6,4,10,1,7]

;tframelist='F3L,F3R,F2R,F2L,F5L,F5R,F0R,F0L,F4L,F4R,F1R,F1L'
;framelist=[3,9,8,2,5,11,6,0,4,10,7,1]

;tframelist='F2L,F2R,F4R,F4L,F1L,F1R,F3R,F3L,F0L,F0R,F5R,F5L'
;framelist=[2,8,10,4,1,7,9,3,0,6,11,5]

;tframelist='F3L,F3R,F0R,F0L,F2L,F2R,F5R,F5L,F1L,F1R,F4R,F4L'
;framelist=[3,9,6,0,2,8,11,5,1,7,10,4]

;tframelist='F3L,F3R,F2R,F2L,F4L,F4R,F1R,F1L,F5L,F5R,F0R,F0L'
;framelist=[3,9,8,2,4,10,7,1,5,11,6,0]

;tframelist='F3L,F3R,F5L,F5R,F4L,F4R,F1L,F1R,F0L,F0R,F2L,F2R'
;framelist=[3,9,5,11,4,10,1,7,0,6,2,8]
;framelist=[3,9,11,5,4,10,7,1,0,6,8,2]

;tframelist='F3L,F3R,F1L,F1R,F4L,F4R,F0L,F0R,F5L,F5R,F2L,F2R'
;framelist=[3,9,1,7,4,10,0,6,5,11,2,8]

;tframelist='F3L,F3R,F1R,F1L,F4L,F4R,F0R,F0L,F5L,F5R,F2R,F2L'
;framelist=[3,9,7,1,4,10,6,0,5,11,8,2]

;tframelist='F3L,F3R,F1L,F1R,F5L,F5R,F0L,F0R,F4L,F4R,F2L,F2R'
;framelist=[3,9,1,7,5,11,0,6,4,10,2,8]
;framelist=[3,9,7,1,5,11,6,0,4,10,8,2]

;tframelist='F4L,F4R,F3R,F3L,F0L,F0R,F5R,F5L,F2L,F2R,F1R,F1L'
;framelist=[4,10,9,3,0,6,11,5,2,8,7,1]

;tframelist='F3L,F3R,F0R,F0L,F1L,F1R,F4R,F4L,F5L,F5R,F2R,F2L'
;framelist=[3,9,6,0,1,7,10,4,5,11,8,2]

;tframelist='F3L,F3R,F0R,F0L,F1L,F1R,F5R,F5L,F4L,F4R,F2R,F2L'
;framelist=[3,9,6,0,1,7,11,5,4,10,8,2]

;tframelist='F3L,F3R,F0R,F0L,F4L,F4R,F1R,F1L,F5L,F5R,F2R,F2L'
;framelist=[3,9,6,0,4,10,7,1,5,11,8,2]

;tframelist='F3L,F3R,F0R,F0L,F4L,F4R,F5R,F5L,F1L,F1R,F2R,F2L'
;framelist=[3,9,6,0,4,10,11,5,1,7,8,2]

;tframelist='LR,RL,LR,RL,LR,RL'
;framelist=[2,8,9,3,0,6,11,5,1,7,10,4]


; LOOK-UP TABLES WITH MDI-LIKE ALGORITHM
;-----------------------------------------------------------------------------

ntest      = 501 ; Set test velocities
dvtest     = dlam/dlamdv;100.d0 ; in m/s
vtest      = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)
lines      = DBLARR(nlam,ntest)
FOR i      = 0,ntest-1 DO lines[*,i] = INTERPOL(line,lam,lam+(vtest[i])*dlamdv)
inten      = TRANSPOSE(filters)#lines

x          = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune
c1         = REFORM(COS(x)#inten)  ; coefficient a1 (for the cosine)
s1         = REFORM(SIN(x)#inten)  ; coefficient b1 (for the sine)
linedepths = sqrt(c1^2.d0+s1^2.d0)
c2         = REFORM(COS(2.d0*x)#inten)  ; coefficient a2
s2         = REFORM(SIN(2.d0*x)#inten)  ; coefficient b2 (for the sine)
pv1        = dvtune*inttune*2.d0   ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2        = dvtune*inttune
phi1       = ATAN(-s1,-c1)
phi2       = ATAN(-s2,-c2)
vel1       = phi1*pv1/2.d0/!dpi    ; velocity measurement
vel2       = phi2*pv2/2.d0/!dpi    ; velocity measurement
vel1a      = (vel1-vtest+10.5d0*pv1) MOD pv1-pv1/2.d0+vtest ; the lookup tables
vel2a      = (vel2-vel1+10.5d0*pv2)  MOD pv2-pv2/2.d0+vel1


; WE LOOP OVER THE TIME AND SPACE COORDINATES
;------------------------------------------------------------------------------

;RESTORE,'velocity_im3dV2.bin' ; Doppler velocity measured by MDI at location [100,50] of datacube im3dV2.fits, dt=60 sec, nt=776
RESTORE,'velocity_90015.bin' ; Doppler velocity + l.o.s magnetic field at locations [*,172,*] and [*,22,*] of datacube *90015*, nt=600, dt=60 secs
; Spline interpolation in space to have HMI-like full-disk resolution:
; dx=695.99 [Mm]/1970 [pixels]=0.35329442 Mm instead of dx=0.826018 Mm
; from file /scratch/109a/couvidat/SUNSPOTS/hr_Vm_700x700_01h.90015_lat06.2Nlon082.2.fits

; INTERPOLATION COEFFICIENTS
; interpolation coefficients were obtained from accel.pro
RESTORE,'interpolation_coefficients.bin'
dtint= tint[1]-tint[0]
nint = N_ELEMENTS(tint)

spatialinterpolation=0 ; take into account the spatial drift due to solar rotation (yes=1, no=0) ?

; HORIZONTAL SHIFT DUE TO SOLAR ROTATION (ASSUMING P0=0 AND B0=0)
depart   = 0 ; x position at which we start on the Doppler velocity and magnetic field maps
fin      = nt+depart-1
PRINT,'BILINEAR INTERPOLATION'
shiftx=0.01960173d0 ; pixels per 45./12.=3.75 seconds for a solar radius of 1970 pixels
IX=findgen(30405)*shiftx
IY=findgen(600*60./dt)*dt/60.d0

; IN THE REST, WE ASSUME: B=(V_LCP-V_RCP) AND V=(V_LCP+V_RCP)/2
; HERE ARE THE SERIES OF REFERENCE: velocity AND magnetic
IF(spatialinterpolation EQ 1) THEN BEGIN
    velocity = BILINEAR(velocity1int,IX,IY)
    magnetic = BILINEAR(magnetic1int,IX,IY)
ENDIF ELSE BEGIN
;    velocity = INTERPOL(REFORM(velocity1int[depart,*]),FINDGEN(600),IY,/QUADRATIC)
;    magnetic = INTERPOL(REFORM(magnetic1int[depart,*]),FINDGEN(600),IY,/QUADRATIC)
    velocity=[I_F(REFORM(velocity1int[depart,*]),16),I_F(REFORM(velocity1int[depart+100,*]),16),I_F(REFORM(velocity1int[depart+200,*]),16)]
    magnetic=[I_F(REFORM(magnetic1int[depart,*]),16),I_F(REFORM(magnetic1int[depart+100,*]),16),I_F(REFORM(magnetic1int[depart+200,*]),16)]
ENDELSE

PRINT,'DONE'

inten    = DBLARR(nt2,nframe)
inten2   = inten

PRINT,'INTERPOLATE SOLAR LINE'


; WE ADD THE SDO ORBITAL VELOCITY TO THE DOPPLER VELOCITY
; SDO orbital velocity ; between -3.1 to +3.1 km/s in 12 hours (I
; assume a circular orbit)
phase=0.d0
vSDO=3079.89d0*COS(2.d0*!dpi*time/24.d0/3600.d0+phase)
;vSDO=vSDO-(FINDGEN(nt)/(nt-1.d0)*6840.22-3420.11d0)
;vSDO=-(FINDGEN(nt)/(nt-1.d0)*13000.20-6500.d0)

; we add the solar rotation velocity
offset=0.d0;at disk center
IF(spatialinterpolation EQ 1) THEN BEGIN
    velocity[depart:fin,*]=(velocity[depart:fin,*]-MEAN(velocity[depart:fin,*]))+vSDO+offset
ENDIF ELSE BEGIN
    velocity = (velocity-MEAN(velocity))+vSDO+offset
ENDELSE


IF(spatialinterpolation EQ 1) THEN BEGIN
    FOR t=0,nt2-1 DO BEGIN      ; observable index
        
        FOR j=0,nframe-1 DO BEGIN
                                ; INTERPOL THE SOLAR LINE SPECTRUM AT THE DOPPLER VELOCITY
            IF(framelist[j] LE nframe/2-1) THEN BEGIN
                LCP   = velocity[t*nframe+j+depart,t*nframe+j]+magnetic[t*nframe+j+depart,t*nframe+j]/2.d0/LANDE
                velo=-LCP*dlamdv
                line2 = decalage(line,velo,dlam)
                inten[t,j]  = TOTAL(filters[*,framelist[j]]*line2)
            ENDIF ELSE BEGIN
                RCP   = velocity[t*nframe+j+depart,t*nframe+j]-magnetic[t*nframe+j+depart,t*nframe+j]/2.d0/LANDE
                velo=-RCP*dlamdv
                line2 = decalage(line,velo,dlam)
                inten[t,j]  = TOTAL(filters[*,framelist[j]-nframe/2]*line2)
            ENDELSE
        ENDFOR
        
    ENDFOR
ENDIF ELSE BEGIN
    FOR t=0,nt2-1 DO BEGIN      ; observable index
        
        FOR j=0,nframe-1 DO BEGIN
                                ; INTERPOL THE SOLAR LINE SPECTRUM AT THE DOPPLER VELOCITY
            IF(framelist[j] LE nframe/2-1) THEN BEGIN
                LCP   = velocity[t*nframe+j]+magnetic[t*nframe+j]/2.d0/LANDE 
                velo=-LCP*dlamdv
                line2 = decalage(line,velo,dlam)
                inten[t,j]  = TOTAL(filters[*,framelist[j]]*line2)
            ENDIF ELSE BEGIN
                RCP   = velocity[t*nframe+j]-magnetic[t*nframe+j]/2.d0/LANDE
                velo=-RCP*dlamdv
                line2 = decalage(line,velo,dlam)
                inten[t,j]  = TOTAL(filters[*,framelist[j]-nframe/2]*line2)
            ENDELSE
        ENDFOR
        
    ENDFOR
ENDELSE


PRINT,'DONE'

v1LCP= DBLARR(nt2)
v2LCP= v1LCP
v1RCP= v1LCP
v2RCP= v1LCP

PRINT,'COMPUTE OBSERVABLES'

FOR t=3,nt2-3 DO BEGIN ; we start at 2 because with the interpolation algorithm we need 2 before and 3 after

    PRINT,t

; WE CORRECT THE FILTERGRAMS FOR THE FACT THAT THEY ARE NOT TAKEN AT
; THE EXACT SAME TIME AND THAT BECAUSE OF SOLAR ROTATION THEY ARE NOT
; TAKEN AT THE EXACT SAME LOCATION
;------------------------------------------------------------------------------------

; WE INTERPOLATE IN TIME
; the first frame (filtergram) is taken at relative time t=0
; the second frame is taken at t=dt
; the third at t=2*dt, and so on....


    FOR j=0,nframe-1 DO BEGIN
         dtime=j*dt
         dtime=ROUND(dtime/dtint)
         inten2[t,j] = coeffjx[0,nint-1-dtime]*inten[t-3,j]+coeffjx[1,nint-1-dtime]*inten[t-2,j]+coeffjx[2,nint-1-dtime]*inten[t-1,j]+coeffjx[3,nint-1-dtime]*inten[t,j]+coeffjx[4,nint-1-dtime]*inten[t+1,j]+coeffjx[5,nint-1-dtime]*inten[t+2,j]
;        inten2[t,j]=INTERPOL(inten[*,j],FINDGEN(nt2)*cadence+dt*j,t*cadence)
    ENDFOR

; WE INTERPOLATE IN SPACE
;    FOR j=0,nframe-1 DO BEGIN
;        dshiftx=FINDGEN(nt2)+shiftx*j
;        inten2[t,j] = INTERPOL(inten2[*,j],dshiftx,t,/SPLINE)
;    ENDFOR



; MDI-LIKE ALGORITHM TO OBTAIN THE DOPPLER VELOCITY
;-----------------------------------------------------------------------------------

    temp=SORT(framelist)

    ; LCP
    intent=REFORM(inten2[t,temp[0:nframe/2-1]]) ; depends on framelist definition
    c1   = TOTAL(COS(x)*intent)
    s1   = TOTAL(SIN(x)*intent)
    c2   = TOTAL(COS(2.d0*x)*intent)
    s2   = TOTAL(SIN(2.d0*x)*intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel1at=vel1t
    vel2t = phi2t*pv2/2.d0/!dpi
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1LCP[t] = INTERPOL(vtest,vel1a,vel1at) ; velocity obtained from look-up table
    v2LCP[t] = INTERPOL(vtest,vel2a,vel2at) ; velocity obtained from look-up table

    ; RCP
    intent=REFORM(inten2[t,temp[nframe/2:nframe-1]]) ; depends on framelist definition
    c1   = TOTAL(COS(x)*intent)
    s1   = TOTAL(SIN(x)*intent)
    c2   = TOTAL(COS(2.d0*x)*intent)
    s2   = TOTAL(SIN(2.d0*x)*intent)
    phi1t = ATAN(-s1,-c1)
    phi2t = ATAN(-s2,-c2)
    vel1t = phi1t*pv1/2.d0/!dpi
    vel1at=vel1t
    vel2t = phi2t*pv2/2.d0/!dpi
    vel2at= (vel2t-vel1t+10.5d0*pv2) MOD pv2-pv2/2.d0+vel1t
    v1RCP[t] = INTERPOL(vtest,vel1a,vel1at) ; velocity obtained from look-up table
    v2RCP[t] = INTERPOL(vtest,vel2a,vel2at) ; velocity obtained from look-up table


    
ENDFOR
PRINT,'DONE'

vel=FLTARR(nt2)
mag=vel

IF(spatialinterpolation EQ 1) THEN BEGIN
    velocity=REFORM(REBIN(velocity[depart:fin,*],nt2,nt2))
    magnetic=REFORM(REBIN(magnetic[depart:fin,*],nt2,nt2))
    FOR i=0,nt2-1 DO BEGIN
        vel[i]=velocity[i,i]
        mag[i]=magnetic[i,i]
    ENDFOR
ENDIF ELSE BEGIN
   ;velocity=REFORM(REBIN(velocity[depart,*],1,nt2))
   ;magnetic=REFORM(REBIN(magnetic[depart,*],1,nt2))
    velocity=velocity[FINDGEN(nt2)*nframe]
    magnetic=magnetic[FINDGEN(nt2)*nframe]
    FOR i=0,nt2-1 DO BEGIN
        vel[i]=velocity[i]
        mag[i]=magnetic[i]
    ENDFOR
ENDELSE

velocity=vel[3:nt2-3]
magnetic=mag[3:nt2-3]


estvel=(v1LCP+v1RCP+v2LCP+v2RCP)/4.d0 ; estimated velocity
estmag=((v1LCP-v1RCP)+(v2LCP-v2RCP))*LANDE/2.d0 ; estimated magnetic field
estvel=estvel[3:nt2-3]
estmag=estmag[3:nt2-3]

SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xoffset=1,yoffset=1,xsize=20,ysize=26,/color,bits=24
LOADCT,3
!P.MULTI=[0,2,3]

;a=WHERE(ABS(estvel-velocity) GT 20.d0,na,COMPLEMENT=b)
;PRINT,'ERRORS=',na,a
;estvel[a]=velocity[a]+MEAN(estvel[b]-velocity[b])
;estmag[a]=magnetic[a]+MEAN(estmag[b]-magnetic[b])

aa=[795,796,797,798,799,1596,1597,1598,1599]
a=WHERE(ABS(estvel-velocity) GT 10.d0,na,COMPLEMENT=b)
estvel[aa]=velocity[aa]+MEAN(estvel[b]-velocity[b])
estmag[aa]=magnetic[aa]+MEAN(estmag[b]-magnetic[b])

PLOT,velocity,estvel-velocity,psym=1,xst=1,tit='!17'+tframelist,xtit='Velocity (m/s)',ytit='Error (m/s)',charsize=1.,yrange=[-15,15],yst=1

res=HISTOGRAM(estvel-velocity,BINSIZE=.25,LOCATIONS=xh,MAX=25,MIN=-25)
;histo,estvel-velocity,-50,50,1.
PLOT,xh,res,charsize=1.,tit='!17'+STRING(MIN(estvel-velocity))+STRING(MAX(estvel-velocity))+STRING(SIGMA(estvel-velocity))+STRING(MEAN(estvel-velocity)),xtit='Error on velocity (m/s)',ytit='Histogram',xst=1


PLOT,FINDGEN(nt2-5)*cadence/60.d0,velocity,tit='!17',xtit='Time (min)',ytit='Velocity (m/s)',charsize=1.,yst=1,xst=1

PLOT,FINDGEN(nt2-5)*cadence/60.d0,magnetic,tit='!17',xtit='Time (min)',ytit='L.o.s. magnetic field strength (B)',charsize=1.,yst=1,xst=1

;t=sort(velocity)
;PLOT,velocity[t],estvel[t]-velocity[t],xst=1,tit='!17'+tframelist,xtit='Velocity (m/s)',ytit='Error (m/s)',charsize=1.5,yrange=[-15,15],yst=1


;histo,(estvel-velocity)/velocity,-0.04,0.04,.001
;res=HISTOGRAM((estvel-velocity)/velocity,BINSIZE=.00025,LOCATIONS=xh,MAX=0.04,MIN=-0.04)
;PLOT,xh,res,charsize=1.,tit='!17'+STRING(SIGMA(estvel-velocity))+STRING(MEAN(estvel-velocity)),xtit='Relative error (m/s)',ytit='Histogram',xst=1

PLOT,magnetic,estmag-magnetic,psym=1,xst=1,tit='!17'+tframelist,xtit='L.o.s. magnetic field strength (B)',ytit='Error (B)',charsize=1.,yrange=[-5,5],yst=1

res=HISTOGRAM(estmag-magnetic,BINSIZE=.1,LOCATIONS=xh,MAX=4,MIN=-4)
PLOT,xh,res,charsize=1.,tit='!17'+STRING(MIN(estmag-magnetic))+STRING(MAX(estmag-magnetic))+STRING(SIGMA(estmag-magnetic))+STRING(MEAN(estmag-magnetic)),xtit='Error on field (B)',ytit='Histogram',xst=1


DEVICE,/CLOSE

READ,pause


END
