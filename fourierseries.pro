
PRO gaussianfit,X,A,F,pder

F=1.d0-A[0]*exp(-(X-A[2])^2.d0/A[1]^2.d0)

pder=[ [-exp(-(X-A[2])^2.d0/A[1]^2.d0)],[-2.d0*A[0]*(X-A[2])^2.d0/A[1]^3.d0*exp(-(X-A[2])^2.d0/A[1]^2.d0)], [-A[0]/A[1]^2.d0*2.d0*(X-A[2])*exp(-(X-A[2])^2.d0/A[1]^2.d0)]]



END

PRO fourierseries



;lam=findgen(1000)/999.-0.5
;width=0.134
;depth=0.8
;l0=-0.12
;y=1.d0-depth*exp(-(lam-l0)^2.d0/width^2.d0)
;period=lam[999]-lam[0]
;
;
;c1=2.d0/period*int_tabulated(lam,y*cos(2.d0*!dpi*lam/period))
;s1=2.d0/period*int_tabulated(lam,y*sin(2.d0*!dpi*lam/period))
;c2=2.d0/period*int_tabulated(lam,y*cos(4.d0*!dpi*lam/period))
;s2=2.d0/period*int_tabulated(lam,y*sin(4.d0*!dpi*lam/period))
;
;plot,lam,y,xst=1
;print,c1,s1,c2,s2
;print,period/!dpi*sqrt(1.d0/3.d0*alog((sqrt(c1^2.d0+s1^2.d0))/(sqrt(c2^2.d0+s2^2.d0)))),width
;print,period/2.d0/!dpi*atan(s1/c1),l0
;print,sqrt(c1^2.d0+s1^2.d0)/2.d0*period/sqrt(!dpi)/width*exp(!dpi^2.d0*width^2.d0/period^2.d0),depth


;-----------------------------------------------------------------------------------------------------------------------
; NOW HMI CASE
;-----------------------------------------------------------------------------------------------------------------------

lam0    = 6173.3433d0

ntune   = 6      ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)


; PARAMETERS
;-----------------------------------------------------------------------
; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
;phase      = FLTARR(7)
;contrast   = FLTARR(7)+1.0

wmich      = [0.172d0,0.344d0,0.693d0]          ; FSRs of the tunable elements
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements

dtune  = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
period = (ntune-1.d0)*dtune
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

nlamp  = 30061 ; Set wavelength points for line profile
lamp   = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

; GRID WE WANT IN WAVELENGTH
nlam          = 60000.
dlam          = 1.d0/1.75d3/5.d0
lam           = (DINDGEN(nlam)-(nlam-1.d0)/2.d0)*dlam

dlamdv        = 2.059205672212074294d-5
dvdlam        = 1.d0/dlamdv
vel           = lam*dvdlam   ; DOPPLER shift in cm/s
dvel          = dlam*dvdlam
dvtune        = dvdlam*dtune ; DOPPLER shift between tuning positions in m/s
vtune         = dvdlam*tune  ; DOPPLER shift of the tuning positions in m/s

; BLOCKER FILTER
RESTORE,'frontwindow.bin' ; front window
blocker       = INTERPOL(transmission/100.d0,wavelength*10.d0-6173.3433d0,lam)
q             = READFITS('blocker11.fits')
blocker       = blocker*INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-6173.3433d0,lam) ; blocker used for look-up tables

;wavelength = [-0.475470,-0.465460,-0.455450,-0.445440,-0.435430,-0.425420,-0.415410,-0.405400,-0.395390,-0.385380,-0.375370,-0.365360,-0.355350,-0.345340,-0.335330,-0.325320,-0.315310,-0.305300,-0.295290,-0.285280,-0.275270,-0.265260,-0.255250,-0.245240,-0.235230,-0.225220,-0.215210,-0.205200,-0.195190,-0.185180,-0.175170,-0.165160,-0.155150,-0.145140,-0.135130,-0.125120,-0.115110,-0.105100,-0.0950900,-0.0850800,-0.0750700,-0.0650600,-0.0550500,-0.0450400,-0.0350300,-0.0250200,-0.0150100,-0.00500000,0.00501000,0.0150200,0.0250300,0.0350400,0.0450500,0.0550600,0.0650700,0.0750800,0.0850900,0.0951000,0.105110,0.115120,0.125130,0.135140,0.145150,0.155160,0.165170,0.175180,0.185190,0.195200,0.205210,0.215220,0.225230,0.235240,0.245250,0.255260,0.265270,0.275280,0.285290,0.295300,0.305310,0.315320,0.325330,0.335340,0.345350,0.355360,0.365370,0.375380,0.385390,0.395400,0.405410,0.415420,0.425430,0.435440,0.445450,0.455460,0.465470,0.475480,0.485490,0.495500]

;lineprofile0 =[1.00447,1.00444,1.00501,1.00508,1.00511,1.00511,1.00509,1.00508,1.00473,1.00401,1.00340,1.00269,1.00174,1.00009,0.998082,0.996317,0.994873,0.993968,0.993783,0.993825,0.994025,0.994263,0.994643,0.994890,0.995349,0.995453,0.995142,0.994554,0.994067,0.993847,0.993929,0.993886,0.994019,0.993814,0.993427,0.992296,0.989101,0.980927,0.964652,0.935828,0.889673,0.823512,0.739270,0.643219,0.544847,0.456429,0.388372,0.347983,0.342219,0.375538,0.447979,0.550264,0.662702,0.764941,0.843935,0.896826,0.929316,0.948714,0.960595,0.968336,0.973091,0.975876,0.977372,0.978072,0.978475,0.979077,0.979950,0.980730,0.982280,0.983586,0.985572,0.987222,0.988758,0.989460,0.989678,0.989728,0.989511,0.989442,0.989508,0.989419,0.989547,0.991192,0.991897,0.991990,0.991844,0.991826,0.992328,0.992060,0.991896,0.992043,0.992026,0.991417,0.991130,0.990842,0.991048,0.990855,0.989203,0.989201]

wavelength = [-0.295290,-0.285280,-0.275270,-0.265260,-0.255250,-0.245240,-0.235230,-0.225220,-0.215210,-0.205200,-0.195190,-0.185180,-0.175170,-0.165160,-0.155150,-0.145140,-0.135130,-0.125120,-0.115110,-0.105100,-0.0950900,-0.0850800,-0.0750700,-0.0650600,-0.0550500,-0.0450400,-0.0350300,-0.0250200,-0.0150100,-0.00500000,0.00501000,0.0150200,0.0250300,0.0350400,0.0450500,0.0550600,0.0650700,0.0750800,0.0850900,0.0951000,0.105110,0.115120,0.125130,0.135140,0.145150,0.155160,0.165170,0.175180,0.185190,0.195200,0.205210,0.215220,0.225230,0.235240,0.245250,0.255260,0.265270,0.275280,0.285290,0.295300,0.305310,0.315320,0.325330,0.335340,0.345350,0.355360,0.365370,0.375380,0.385390]

lineprofile0 =[0.993783,0.993825,0.994025,0.994263,0.994643,0.994890,0.995349,0.995453,0.995142,0.994554,0.994067,0.993847,0.993929,0.993886,0.994019,0.993814,0.993427,0.992296,0.989101,0.980927,0.964652,0.935828,0.889673,0.823512,0.739270,0.643219,0.544847,0.456429,0.388372,0.347983,0.342219,0.375538,0.447979,0.550264,0.662702,0.764941,0.843935,0.896826,0.929316,0.948714,0.960595,0.968336,0.973091,0.975876,0.977372,0.978072,0.978475,0.979077,0.979950,0.980730,0.982280,0.983586,0.985572,0.987222,0.988758,0.989460,0.989678,0.989728,0.989511,0.989442,0.989508,0.989419,0.989547,0.991192,0.991897,0.991990,0.991844,0.991826,0.992328]+0.0062170029

;lineprofile0=1.d0-0.632*exp(-wavelength^2.d0/0.04342^2.d0) ; perfect Gaussian profile
;print,"actual parameters",0.04342d0*2.d0*SQRT(ALOG(2.d0)),0.632d0

; INPUT VELOCITIES
ntest         = 501 ; Set test velocities
dvtest        = 30. ; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)

; NON-TUNABLE PROFILE
lyot          = blocker
FOR i = 0,3 DO BEGIN
    lyot      = lyot   *(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0
ENDFOR

cmich     = 2.d0*!dpi/wmich

pv1      = dvtune*inttune*2.d0 ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2      = dvtune*inttune

line     = INTERPOL(lineprofile0,wavelength,lamp)
a        = WHERE(lamp ge MAX(wavelength))
line[a]  = 1.d0
a        = WHERE(lamp le MIN(wavelength))
line[a]  = 1.d0
q2       = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1       = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles

x        = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune

vel1a=fltarr(ntest)
vel2a=vel1a
widtha=vel1a
deptha=vel1a
contin=vel1a
contin2=vel1a
c1=vel1a
c2=vel1a
s1=vel1a
s2=vel1a
a1=vel1a
a2=vel1a
b1=vel1a
b2=vel2a
width=vel1a
depth=vel1a
vel=vel1a
contin3=vel1a
aaaa=where(lam ge -period/2.d0 and lam le period/2.d0)

filters   = DBLARR(nlam,ntune)
FOR itune = 0,ntune-1 DO BEGIN
    filters[*,itune] = lyot
    FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+phase[i]))/2.d0
ENDFOR
norm=TOTAL(filters)/6.d0

lines   = INTERPOL(q1,q2,lam)   
PLOT,lam,lines,xst=1,xrange=[-0.7,0.7]
W=FLTARR(nlam)+1.d0
A=[1.d0,0.07,0.0]
res=curvefit(lam[aaaa],lines[aaaa],W,A,FUNCTION_NAME='gaussianfit')
oplot,lam[aaaa],res,col=180
print,A[0],A[1]*2.d0*SQRT(ALOG(2.d0)),A[2]
lines   = INTERPOL(q1,q2,lam-6500.*dlamdv) 
OPLOT,lam,lines,linestyle=2
lines   = INTERPOL(q1,q2,lam+6500.*dlamdv) 
OPLOT,lam,lines,linestyle=3
for i=0,5 do OPLOT,lam,filters[*,i],thick=2


FOR kkk=0,ntest-1 DO BEGIN

    lines   = INTERPOL(q1,q2,lam+vtest[kkk]*dlamdv)    
          ; INTENSITIES
            inten      = TRANSPOSE(filters)#lines
            c1[kkk]         = COS(x[0])*inten[0]+COS(x[1])*inten[1]+COS(x[2])*inten[2]+COS(x[3])*inten[3]+COS(x[4])*inten[4]+COS(x[5])*inten[5]
            s1[kkk]         = SIN(x[0])*inten[0]+SIN(x[1])*inten[1]+SIN(x[2])*inten[2]+SIN(x[3])*inten[3]+SIN(x[4])*inten[4]+SIN(x[5])*inten[5]
            c2[kkk]         = COS(2.d0*x[0])*inten[0]+COS(2.d0*x[1])*inten[1]+COS(2.d0*x[2])*inten[2]+COS(2.d0*x[3])*inten[3]+COS(2.d0*x[4])*inten[4]+COS(2.d0*x[5])*inten[5]
            s2[kkk]         = SIN(2.d0*x[0])*inten[0]+SIN(2.d0*x[1])*inten[1]+SIN(2.d0*x[2])*inten[2]+SIN(2.d0*x[3])*inten[3]+SIN(2.d0*x[4])*inten[4]+SIN(2.d0*x[5])*inten[5]

            a1[kkk]=2.d0/period*int_tabulated(lam[aaaa],lines[aaaa]*cos(2.d0*!dpi*lam[aaaa]/period))
            b1[kkk]=2.d0/period*int_tabulated(lam[aaaa],lines[aaaa]*sin(2.d0*!dpi*lam[aaaa]/period))
            a2[kkk]=2.d0/period*int_tabulated(lam[aaaa],lines[aaaa]*cos(4.d0*!dpi*lam[aaaa]/period))
            b2[kkk]=2.d0/period*int_tabulated(lam[aaaa],lines[aaaa]*sin(4.d0*!dpi*lam[aaaa]/period))

            phi1       = ATAN(-s1[kkk],-c1[kkk])
            phi2       = ATAN(-s2[kkk],-c2[kkk])
            vel1a[kkk] = phi1*pv1/2.d0/!dpi ; velocity measurement
            vel2       = phi2*pv2/2.d0/!dpi ; velocity measurement
            vel2a[kkk] = ((vel2-vel1a[kkk]+10.5d0*pv2) MOD pv2)-pv2/2.d0+vel1a[kkk]
            temp       = period/!dpi*sqrt(alog(sqrt(c1[kkk]^2.d0+s1[kkk]^2.d0)/sqrt(c2[kkk]^2.d0+s2[kkk]^2.d0))/3.d0)
            widtha[kkk]= temp*2.d0*sqrt(alog(2.d0))
            deptha[kkk]= (dtune*2.d0/period)*sqrt(c1[kkk]^2.d0+s1[kkk]^2.d0)*period/2.d0/sqrt(!dpi)/temp*exp(!dpi^2.d0*temp^2.d0/period^2.d0)
            contin[kkk]= ((inten[0]+deptha[kkk]*exp(-(tune[0,0]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[1]+deptha[kkk]*exp(-(tune[0,1]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[2]+deptha[kkk]*exp(-(tune[0,2]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[3]+deptha[kkk]*exp(-(tune[0,3]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[4]+deptha[kkk]*exp(-(tune[0,4]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[5]+deptha[kkk]*exp(-(tune[0,5]-vel1a[kkk]*dlamdv)^2.d0/temp^2.d0)))/6.d0
            contin3[kkk]= ((inten[0]+deptha[kkk]*exp(-(tune[0,0]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[1]+deptha[kkk]*exp(-(tune[0,1]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[2]+deptha[kkk]*exp(-(tune[0,2]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[3]+deptha[kkk]*exp(-(tune[0,3]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[4]+deptha[kkk]*exp(-(tune[0,4]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0))+(inten[5]+deptha[kkk]*exp(-(tune[0,5]-vtest[kkk]*dlamdv)^2.d0/temp^2.d0)))/6.d0
            ;if(vel1a[kkk] GT 0) then contin2[kkk]=inten[0]/norm else contin2[kkk]=inten[5]/norm
            contin2[kkk]= (((vel1a[kkk]+pv1/2.)/pv1*inten[0]+(-vel1a[kkk]+pv1/2.)/pv1*inten[5]))


            vel[kkk]=-period/2.d0/!dpi*ATAN(-b1[kkk],-a1[kkk])/dlamdv
            temp=period/!dpi*sqrt(alog((a1[kkk]^2.d0+b1[kkk]^2.d0)/(a2[kkk]^2.d0+b2[kkk]^2.d0))/6.d0)
            width[kkk]= temp*2.d0*sqrt(alog(2.d0))
            depth[kkk]= sqrt(a1[kkk]^2.d0+b1[kkk]^2.d0)*period/2.d0/sqrt(!dpi)/temp*exp(!dpi^2.d0*temp^2.d0/period^2.d0)

ENDFOR

;plot,vtest,vel1a,xst=1
;plot,vtest,widtha,xst=1
print,width[250],depth[250]
print,widtha[250],deptha[250]
widthe=fltarr(501)+0.093888341
plot,vtest,(widtha-widthe)/widthe,xst=1,xrange=[-6500,6500]


read,pause


END
