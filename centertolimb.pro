; center-to-limb variation of Fe I lineprofile
; from Balthasar, H. (1988)

PRO gaussianfit,X,A,F,pder

F=A[0]*exp(-(X-A[2])^2.d0/A[1]^2.d0)

pder=[ [F/A[0]],[2.d0*A[0]*(X-A[2])^2.d0/A[1]^3.d0*exp(-(X-A[2])^2.d0/A[1]^2.d0)], [A[0]/A[1]^2.d0*2.d0*(X-A[2])*exp(-(X-A[2])^2.d0/A[1]^2.d0)]]



END

; fit the linedepth using eq 32 of Unno, W. (1956)

PRO linedepth,X,A,F,pder

F=A[0]*A[1]*X/(1.d0+A[1]*X)

pder= [[F/A[0]],[A[0]*X/(1.d0+A[1]*X)-A[0]*A[1]*X^2.d0/(1.d0+A[1]*X)^2.d0] ]

END


PRO centertolimb,angle

;angle (in degrees) is the center to disk angular distance for which you want a
;solar line profile to be interpolated

wavelength=[-0.475470,-0.465460,-0.455450,-0.445440,-0.435430,-0.425420,-0.415410,-0.405400,-0.395390,-0.385380,-0.375370,-0.365360,-0.355350,-0.345340,-0.335330,-0.325320,-0.315310,-0.305300,-0.295290,-0.285280,-0.275270,-0.265260,-0.255250,-0.245240,-0.235230,-0.225220,-0.215210,-0.205200,-0.195190,-0.185180,-0.175170,-0.165160,-0.155150,-0.145140,-0.135130,-0.125120,-0.115110,-0.105100,-0.0950900,-0.0850800,-0.0750700,-0.0650600,-0.0550500,-0.0450400,-0.0350300,-0.0250200,-0.0150100,-0.00500000,0.00501000,0.0150200,0.0250300,0.0350400,0.0450500,0.0550600,0.0650700,0.0750800,0.0850900,0.0951000,0.105110,0.115120,0.125130,0.135140,0.145150,0.155160,0.165170,0.175180,0.185190,0.195200,0.205210,0.215220,0.225230,0.235240,0.245250,0.255260,0.265270,0.275280,0.285290,0.295300,0.305310,0.315320,0.325330,0.335340,0.345350,0.355360,0.365370,0.375380,0.385390,0.395400,0.405410,0.415420,0.425430,0.435440,0.445450,0.455460,0.465470,0.475480,0.485490,0.495500]

lineprofile0=[1.00447,1.00444,1.00501,1.00508,1.00511,1.00511,1.00509,1.00508,1.00473,1.00401,1.00340,1.00269,1.00174,1.00009,0.998082,0.996317,0.994873,0.993968,0.993783,0.993825,0.994025,0.994263,0.994643,0.994890,0.995349,0.995453,0.995142,0.994554,0.994067,0.993847,0.993929,0.993886,0.994019,0.993814,0.993427,0.992296,0.989101,0.980927,0.964652,0.935828,0.889673,0.823512,0.739270,0.643219,0.544847,0.456429,0.388372,0.347983,0.342219,0.375538,0.447979,0.550264,0.662702,0.764941,0.843935,0.896826,0.929316,0.948714,0.960595,0.968336,0.973091,0.975876,0.977372,0.978072,0.978475,0.979077,0.979950,0.980730,0.982280,0.983586,0.985572,0.987222,0.988758,0.989460,0.989678,0.989728,0.989511,0.989442,0.989508,0.989419,0.989547,0.991192,0.991897,0.991990,0.991844,0.991826,0.992328,0.992060,0.991896,0.992043,0.992026,0.991417,0.991130,0.990842,0.991048,0.990855,0.989203,0.989201]
lineprofile45=[1.00391,1.00448,1.00537,1.00534,1.00538,1.00553,1.00567,1.00515,1.00484,1.00458,1.00384,1.00291,1.00148,1.00058,0.999926,0.998906,0.997526,0.996379,0.995958,0.995723,0.996076,0.996137,0.996310,0.995974,0.995705,0.995267,0.995004,0.994715,0.994393,0.993895,0.993235,0.992605,0.992392,0.990756,0.987331,0.980631,0.969781,0.952117,0.926195,0.889543,0.839953,0.777853,0.706090,0.629078,0.554224,0.490055,0.442618,0.417714,0.419243,0.447437,0.500968,0.572437,0.652744,0.732437,0.802353,0.856120,0.896472,0.926850,0.947073,0.960377,0.968769,0.974024,0.977197,0.978859,0.979464,0.979665,0.980069,0.980933,0.982326,0.983639,0.984811,0.986422,0.987440,0.988704,0.989813,0.990314,0.990591,0.990749,0.991216,0.991710,0.992247,0.992572,0.992794,0.993023,0.994427,0.994544,0.994526,0.994284,0.994316,0.993821,0.993642,0.993490,0.992832,0.992035,0.991620,0.990701,0.990598,0.991006]
lineprofile60=[1.00385,1.00464,1.00515,1.00554,1.00511,1.00547,1.00529,1.00485,1.00449,1.00401,1.00308,1.00209,1.00092,0.999794,0.999266,0.998074,0.997084,0.995469,0.995278,0.995446,0.995702,0.995720,0.995649,0.995616,0.994924,0.994406,0.994285,0.994000,0.993230,0.992485,0.991387,0.989238,0.986354,0.981636,0.973618,0.961482,0.943974,0.920111,0.888135,0.847520,0.797651,0.739986,0.677951,0.615385,0.558416,0.511432,0.477672,0.460486,0.460973,0.480296,0.518138,0.571508,0.635503,0.703607,0.768989,0.827190,0.874685,0.910921,0.936821,0.954723,0.966197,0.972945,0.977044,0.979295,0.979871,0.980277,0.980732,0.981306,0.982846,0.984197,0.985302,0.986847,0.987916,0.989406,0.990105,0.990938,0.991718,0.991600,0.992072,0.992811,0.993348,0.993541,0.993844,0.993738,0.995480,0.995176,0.995551,0.995160,0.994802,0.994830,0.994589,0.993867,0.992855,0.992846,0.991849,0.990751,0.990451,0.991012]

cosa=COS(angle/180.d0*!dpi)
cos1=45.;COS(45.d0/180.d0*!dpi)
cos2=60.;COS(60.d0/180.d0*!dpi)


;if(angle GE 60.) then lineprofile=(cosa-1.0)/(cos2-1.0)*(lineprofile60-lineprofile0)+lineprofile0
;if(angle GE 0 AND angle LT 60.) THEN BEGIN
;    if(angle LT 45.d0 and angle GE 0.) then lineprofile=(cosa-1.0)/(cos1-1.0)*(lineprofile45-lineprofile0)+lineprofile0
;    if(angle LT 60.d0 AND angle GE 45.) then lineprofile=(cosa-cos1)/(cos2-cos1)*(lineprofile60-lineprofile45)+lineprofile45
;ENDIF

;C=((lineprofile60-lineprofile0)-cos2/cos1*(lineprofile45-lineprofile0))/(cos2^2.d0-cos1*cos2)
;B=(lineprofile45-lineprofile0)/cos1-C*cos1
;lineprofile=lineprofile0+B*cosa+C*cosa^2.d0
;print,cosa

x=wavelength
cost2=cos([0,45.,60.]*!dpi/180.d0)
;FWHM2= [92.76,107.03,122.13]; lims=0.075
FWHM2=[93.878368,107.41855,123.57606] ;lims=0.055
;R2   = [0.6651,0.58618,0.543871] ;lims=0.075
R2   = [0.662964,0.585863,0.542768] ;lims=0.055
;res  = interpol(FWHM2,cost2,cosa)
;res  = 150.03272-58.075195*cosa; lims=0.075
res  = 151.34559-58.521771*cosa; lims=0.055
;res2 = interpol(R2,cost2,cosa)
;res2  = 0.41868914+0.24424895*cosa; lims=0.075
res2  = 0.41922611+0.24190794*cosa; lims=0.055
x    = x*res/FWHM2[0]
lineprofile=1.d0-interpol((1.d0-lineprofile0)*res2/R2[0],x,wavelength)

plot,wavelength,lineprofile,xst=1
oplot,wavelength,lineprofile0,col=180
oplot,wavelength,lineprofile60,col=180,thick=2,linestyle=3
oplot,wavelength,lineprofile45,col=180,thick=2,linestyle=2

lims=0.075;0.055
yo=where(wavelength gt -lims and wavelength lt lims,na)

W=FLTARR(na)+1.0
A=[0.6,0.09,0.0]
res=curvefit(wavelength[yo],1.-lineprofile0[yo],W,A,FUNCTION_NAME='gaussianfit')
print,A[0],A[1]*2.d0*SQRT(ALOG(2.d0)),A[2]
plot,wavelength,1.-lineprofile0,xst=1
ziggy=findgen(1000)/999.-0.5
yo=where(ziggy gt -lims and ziggy lt lims,na)
oplot,ziggy[yo],A[0]*exp(-(ziggy[yo]-A[2])^2.d0/A[1]^2.d0),col=180,thick=2
oplot,ziggy,A[0]*exp(-(ziggy-A[2])^2.d0/A[1]^2.d0),col=180,linestyle=2
oplot,wavelength,1.d0-lineprofile0,linestyle=2
read,pause

;values for Balthasar (1988)
cost=[1.0,0.922,0.806,0.742,0.671,0.592,0.5,0.387,0.316,0.224,0.194,0.158,0.112];cos(theta)
R=[0.363,0.392,0.413,0.425,0.438,0.440,0.456,0.462,0.475,0.485,0.493,0.505,0.516];intensity of minimum
dl=[0.,-0.7,-0.5,-2.1,-2.4,-1.6,0.1,3.7,1.4,6.5,6.3,7.3,9.4];limb shift
Wl=[76.4,76.9,78.1,79.2,79.5,82.0,82.9,84.7,85.5,85.7,85.4,84.1,82.1];equivalent width
FWHM=Wl*2.d0*sqrt(alog(2.d0)/!dpi)/(1.d0-R)

;for Roger Ulrich's line profiles (values obtained by gaussianfit
;between -0.075 and 0.075 A
cost2=cos([0,45.,60.]*!dpi/180.d0)
FWHM2=[92.76,107.03,122.13]
R2=1.d0-[0.6651,0.58618,0.543871]


;result of linear fit:
; FWHM2= 150.03272-58.075195*cost (it's linear in cost)
; R2   = 0.58131086-0.24424895*cost

; result of linear interpolation of Ulrich's lines in cost
cost3=cos([0,10.,20.,30.,40.,50.,60.,70.,80.,90.]*!dpi/180.d0)
angle3=[0,10.,20.,30.,40.,50.,60.,70.,80.,90.]
FWHM3=[0.092761627,0.093364045,0.095227387,0.098532534,0.10363087,0.11118493,0.12213213,0.13828101,0.16383893,0.21108853]
R3=[0.665100,0.660925,0.648573,0.628568,0.601770,0.572800,0.543871,0.513161,0.481989,0.451561]

FWHM3=[0.092761627,0.093364045,0.095227387,0.098532534,0.10363087,0.11118493,0.12213213,0.13853351,0.16531924,0.21804844]
R3=[0.665100,0.660925,0.648573,0.628568,0.601770,0.572800,0.543871,0.508196,0.471774,0.436037]

FWHM3=[0.092761627,0.093579153,0.095982400,0.099829882,0.10490900,0.11269270,0.12391474,0.13599022,0.14924792,0.16279879]
R3=[0.665100,0.660661,0.647621,0.626732,0.598998,0.569375,0.539534,0.507195,0.472629,0.437281]

cost3=cos([0,10.,20.,30.,40.,45.,50.,60.,70.,80.,90.]*!dpi/180.d0)
FWHM3=[0.092825288,0.093713361,0.096478083,0.10103028,0.10727729,0.11100868,0.11513222,0.12440321,0.13434564,0.14451631,0.15542210]
R3=[0.660933,0.657427,0.645774,0.626866,0.601615,0.586905,0.570968,0.535983,0.497925,0.457393,0.415691]

END
