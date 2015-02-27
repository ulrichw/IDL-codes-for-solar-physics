; idl equivalent to the C code lookup.c used to produce
; lookup tables for HMI

PRO lookup

lam0    = 6173.3433d0

ntune   = 6      ; Set number of tuning positions (NUMBER OF FILTERS)
inttune = 5.d0/2.d0  ; Set number of tuning positions over wnarrow. (1/spacing)

nx      = 256
ny      = 256

; PARAMETERS
;-----------------------------------------------------------------------
; ACTUAL PHASES AND CONTRASTS OF THE FILTER ELEMENTS 
;phase      = [0.0d0,0.d0,0.d0,-5.0700685d0,-0.4029859d0,-4.2042572d0,-8.1685274d0]*!dpi/180.d0 
;contrast   = [0.98,0.99,0.98,0.952,0.964,0.987,1.0]  
phase      = FLTARR(7)
contrast   = FLTARR(7)+1.0

wmich      = [0.172d0,0.344d0,0.693d0]          ; FSRs of the tunable elements
lyotw      = [1.407d0,2.779d0,5.682d0,11.354d0] ; FSRs non-tunable elements

dtune  = wmich[0]/inttune  ; Set tuning positions (dtune = interval between two tuning positions, in Angstrom)
period = (ntune-1.d0)*dtune
tune   = DBLARR(3,ntune)
FOR  i = 0,2 DO tune[i,*]=(-(ntune-1)/2.d0+dindgen(ntune))*dtune

nlamp  = 30061 ; Set wavelength points for line profile
lamp   = -8.D0+DINDGEN(nlamp)*16.D0/(nlamp-1)

; GRID WE WANT IN WAVELENGTH
nlam          = 20000.
dlam          = 1.d0/1.75d3
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

wavelength = [-0.475470,-0.465460,-0.455450,-0.445440,-0.435430,-0.425420,-0.415410,-0.405400,-0.395390,-0.385380,-0.375370,-0.365360,-0.355350,-0.345340,-0.335330,-0.325320,-0.315310,-0.305300,-0.295290,-0.285280,-0.275270,-0.265260,-0.255250,-0.245240,-0.235230,-0.225220,-0.215210,-0.205200,-0.195190,-0.185180,-0.175170,-0.165160,-0.155150,-0.145140,-0.135130,-0.125120,-0.115110,-0.105100,-0.0950900,-0.0850800,-0.0750700,-0.0650600,-0.0550500,-0.0450400,-0.0350300,-0.0250200,-0.0150100,-0.00500000,0.00501000,0.0150200,0.0250300,0.0350400,0.0450500,0.0550600,0.0650700,0.0750800,0.0850900,0.0951000,0.105110,0.115120,0.125130,0.135140,0.145150,0.155160,0.165170,0.175180,0.185190,0.195200,0.205210,0.215220,0.225230,0.235240,0.245250,0.255260,0.265270,0.275280,0.285290,0.295300,0.305310,0.315320,0.325330,0.335340,0.345350,0.355360,0.365370,0.375380,0.385390,0.395400,0.405410,0.415420,0.425430,0.435440,0.445450,0.455460,0.465470,0.475480,0.485490,0.495500]

lineprofile0 =[1.00447,1.00444,1.00501,1.00508,1.00511,1.00511,1.00509,1.00508,1.00473,1.00401,1.00340,1.00269,1.00174,1.00009,0.998082,0.996317,0.994873,0.993968,0.993783,0.993825,0.994025,0.994263,0.994643,0.994890,0.995349,0.995453,0.995142,0.994554,0.994067,0.993847,0.993929,0.993886,0.994019,0.993814,0.993427,0.992296,0.989101,0.980927,0.964652,0.935828,0.889673,0.823512,0.739270,0.643219,0.544847,0.456429,0.388372,0.347983,0.342219,0.375538,0.447979,0.550264,0.662702,0.764941,0.843935,0.896826,0.929316,0.948714,0.960595,0.968336,0.973091,0.975876,0.977372,0.978072,0.978475,0.979077,0.979950,0.980730,0.982280,0.983586,0.985572,0.987222,0.988758,0.989460,0.989678,0.989728,0.989511,0.989442,0.989508,0.989419,0.989547,0.991192,0.991897,0.991990,0.991844,0.991826,0.992328,0.992060,0.991896,0.992043,0.992026,0.991417,0.991130,0.990842,0.991048,0.990855,0.989203,0.989201]

lineprofile45 = [1.00391,1.00448,1.00537,1.00534,1.00538,1.00553,1.00567,1.00515,1.00484,1.00458,1.00384,1.00291,1.00148,1.00058,0.999926,0.998906,0.997526,0.996379,0.995958,0.995723,0.996076,0.996137,0.996310,0.995974,0.995705,0.995267,0.995004,0.994715,0.994393,0.993895,0.993235,0.992605,0.992392,0.990756,0.987331,0.980631,0.969781,0.952117,0.926195,0.889543,0.839953,0.777853,0.706090,0.629078,0.554224,0.490055,0.442618,0.417714,0.419243,0.447437,0.500968,0.572437,0.652744,0.732437,0.802353,0.856120,0.896472,0.926850,0.947073,0.960377,0.968769,0.974024,0.977197,0.978859,0.979464,0.979665,0.980069,0.980933,0.982326,0.983639,0.984811,0.986422,0.987440,0.988704,0.989813,0.990314,0.990591,0.990749,0.991216,0.991710,0.992247,0.992572,0.992794,0.993023,0.994427,0.994544,0.994526,0.994284,0.994316,0.993821,0.993642,0.993490,0.992832,0.992035,0.991620,0.990701,0.990598,0.991006]

lineprofile60 = [1.00385,1.00464,1.00515,1.00554,1.00511,1.00547,1.00529,1.00485,1.00449,1.00401,1.00308,1.00209,1.00092,0.999794,0.999266,0.998074,0.997084,0.995469,0.995278,0.995446,0.995702,0.995720,0.995649,0.995616,0.994924,0.994406,0.994285,0.994000,0.993230,0.992485,0.991387,0.989238,0.986354,0.981636,0.973618,0.961482,0.943974,0.920111,0.888135,0.847520,0.797651,0.739986,0.677951,0.615385,0.558416,0.511432,0.477672,0.460486,0.460973,0.480296,0.518138,0.571508,0.635503,0.703607,0.768989,0.827190,0.874685,0.910921,0.936821,0.954723,0.966197,0.972945,0.977044,0.979295,0.979871,0.980277,0.980732,0.981306,0.982846,0.984197,0.985302,0.986847,0.987916,0.989406,0.990105,0.990938,0.991718,0.991600,0.992072,0.992811,0.993348,0.993541,0.993844,0.993738,0.995480,0.995176,0.995551,0.995160,0.994802,0.994830,0.994589,0.993867,0.992855,0.992846,0.991849,0.990751,0.990451,0.991012]

; INPUT VELOCITIES
ntest         = 901 ; Set test velocities
dvtest        = dlam/dlamdv ; in m/s
vtest         = dvtest*(DINDGEN(ntest)-(ntest-1.d0)/2.d0)

; NON-TUNABLE PROFILE
lyot          = blocker
FOR i = 0,3 DO BEGIN
    lyot      = lyot   *(1.d0+contrast[i+3]*COS(2.d0*!dpi/lyotw[i]*lam+phase[i+3]))/2.d0
ENDFOR

cmich     = 2.d0*!dpi/wmich

RESTORE,'CPT/CPT_laser_front.BIN' ; phase and contrast maps

distance = SHIFT(DIST(nx,nx),nx/2,nx/2)*0.5d0*4096.d0/DOUBLE(nx)
aaa      = WHERE(distance le 960.,complement=bbb)

pv1      = dvtune*inttune*2.d0 ; dynamic range of the tuning positions in m/s (about 16.7 km/s)
pv2      = dvtune*inttune

line     = INTERPOL(lineprofile0,wavelength,lamp)
a        = WHERE(lamp ge MAX(wavelength))
line[a]  = 1.d0
a        = WHERE(lamp le MIN(wavelength))
line[a]  = 1.d0
q2       = [MIN(lamp)-10.d0,lamp,MAX(lamp)+10.d0]
q1       = [MAX(line),line,MAX(line)] ; Calculate Doppler shifted line profiles

vel1a    = DBLARR(nx,ny)
vel2a    = vel1a

x        = 2.d0*!dpi*(-(ntune-1)/2.d0+DINDGEN(ntune))/ntune ; phase of the measurements the formula is 2\pi tune/periode of the ray. periode=5+1 intervals * dtune

SET_PLOT,'x'
!P.multi = [0,2,2]
WINDOW,0,RETAIN=2,XSIZE=1100,YSIZE=900

degree=4
fit1=DBLARR(ntest,degree+1,degree+1)
fit2=DBLARR(ntest,degree+1,degree+1)


vela=FLTARR(ntest)
velb=vela
velc=vela

FOR kkk=0,ntest-1 DO BEGIN

    lines   = INTERPOL(q1,q2,lam+vtest[kkk]*dlamdv) ; CONVENTION: VELOCITIES > 0 FOR BLUESHIFTS (SHIFT TOWARD LOWER WAVELENGTHS), INVERSE OF MDI CONVENTION
    
    filters   = DBLARR(nlam,ntune)
    FOR itune = 0,ntune-1 DO BEGIN
        filters[*,itune] = lyot
        FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+Phig0[nx/2,ny/2,i]))/2.d0
    ENDFOR
    
    
    
                                ; INTENSITIES
    inten      = TRANSPOSE(filters)#lines
    
                                ; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd COEFFICIENT)
    c1         = COS(x[0])*inten[0]+COS(x[1])*inten[1]+COS(x[2])*inten[2]+COS(x[3])*inten[3]+COS(x[4])*inten[4]+COS(x[5])*inten[5] ; coefficient a1 (for the cosine)
    s1         = SIN(x[0])*inten[0]+SIN(x[1])*inten[1]+SIN(x[2])*inten[2]+SIN(x[3])*inten[3]+SIN(x[4])*inten[4]+SIN(x[5])*inten[5] ; coefficient b1 (for the sine)
    c2         = COS(2.d0*x[0])*inten[0]+COS(2.d0*x[1])*inten[1]+COS(2.d0*x[2])*inten[2]+COS(2.d0*x[3])*inten[3]+COS(2.d0*x[4])*inten[4]+COS(2.d0*x[5])*inten[5] ; coefficient a2
    s2         = SIN(2.d0*x[0])*inten[0]+SIN(2.d0*x[1])*inten[1]+SIN(2.d0*x[2])*inten[2]+SIN(2.d0*x[3])*inten[3]+SIN(2.d0*x[4])*inten[4]+SIN(2.d0*x[5])*inten[5] ; coefficient b2 (for the sine)
    
    phi1       = ATAN(-s1,-c1)
    phi2       = ATAN(-s2,-c2)
    vel1       = phi1*pv1/2.d0/!dpi ; velocity measurement
    vel2       = phi2*pv2/2.d0/!dpi ; velocity measurement
    velc[kkk]  = vel1
    vela[kkk]  = ((vel1-vtest[kkk]+10.5d0*pv1) MOD pv1)-pv1/2.d0+vtest[kkk] ; look-up table for 1st Fourier coefficient
    velb[kkk]  = ((vel2-vela[kkk]+10.5d0*pv2) MOD pv2)-pv2/2.d0+vela[kkk] ; look-up table for 2nd Fourier coefficient
    
    
ENDFOR


READ,pause

FOR kkk=100,ntest-1 DO BEGIN

    PRINT,kkk
    lines   = INTERPOL(q1,q2,lam+vtest[kkk]*dlamdv) ; CONVENTION: VELOCITIES > 0 FOR BLUESHIFTS (SHIFT TOWARD LOWER WAVELENGTHS), INVERSE OF MDI CONVENTION

    FOR iii=0,nx-1 DO BEGIN
        PRINT,iii
        FOR jjj=0,ny-1 DO BEGIN
           
            IF(distance[iii,jjj] LE 960.d0) THEN BEGIN
          ; TUNABLE PROFILE
            filters   = DBLARR(nlam,ntune)
            FOR itune = 0,ntune-1 DO BEGIN
                filters[*,itune] = lyot
                FOR i = 0,2 DO filters[*,itune] = filters[*,itune]*(1.d0+contrast[i]*COS(cmich[i]*(lam+tune[i,itune])+Phig0[iii,jjj,i]))/2.d0
            ENDFOR
                    
 
    
          ; INTENSITIES
            inten      = TRANSPOSE(filters)#lines
        
          ; BUILD LOOK-UP TABLES (1 FOR 1st FOURIER COEFFICIENT, 1 FOR 2nd COEFFICIENT)
            c1         = COS(x[0])*inten[0]+COS(x[1])*inten[1]+COS(x[2])*inten[2]+COS(x[3])*inten[3]+COS(x[4])*inten[4]+COS(x[5])*inten[5]; coefficient a1 (for the cosine)
            s1         = SIN(x[0])*inten[0]+SIN(x[1])*inten[1]+SIN(x[2])*inten[2]+SIN(x[3])*inten[3]+SIN(x[4])*inten[4]+SIN(x[5])*inten[5] ; coefficient b1 (for the sine)
            c2         = COS(2.d0*x[0])*inten[0]+COS(2.d0*x[1])*inten[1]+COS(2.d0*x[2])*inten[2]+COS(2.d0*x[3])*inten[3]+COS(2.d0*x[4])*inten[4]+COS(2.d0*x[5])*inten[5] ; coefficient a2
            s2         = SIN(2.d0*x[0])*inten[0]+SIN(2.d0*x[1])*inten[1]+SIN(2.d0*x[2])*inten[2]+SIN(2.d0*x[3])*inten[3]+SIN(2.d0*x[4])*inten[4]+SIN(2.d0*x[5])*inten[5] ; coefficient b2 (for the sine)

            phi1       = ATAN(-s1,-c1)
            phi2       = ATAN(-s2,-c2)
            vel1       = phi1*pv1/2.d0/!dpi ; velocity measurement
            vel2       = phi2*pv2/2.d0/!dpi ; velocity measurement
            vel1a[iii,jjj] = ((vel1-vtest[kkk]+10.5d0*pv1) MOD pv1)-pv1/2.d0+vtest[kkk] ; look-up table for 1st Fourier coefficient
            vel2a[iii,jjj] = ((vel2-vel1a[iii,jjj]+10.5d0*pv2) MOD pv2)-pv2/2.d0+vel1a[iii,jjj] ; look-up table for 2nd Fourier coefficient
            ENDIF
                
        ENDFOR
    ENDFOR

mean1=MEAN(vel1a[aaa])
mean2=MEAN(vel2a[aaa])
vel1a[bbb]= mean1
vel2a[bbb]= mean2
Result1 = SFIT(vel1a-mean1,degree,kx=res1)
TVIM,vel1a,/scale
TVIM,Result1,/scale
Result2 = SFIT(vel2a-mean2,degree,kx=res2)
TVIM,vel2a,/scale
TVIM,Result2,/scale

fit1[kkk,*,*]=res1
fit2[kkk,*,*]=res2

READ,pause

ENDFOR
    
END
    
