; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE CO-TUNE SEQUENCES
; OBJECTIVE: TO DERIVE THE NON-TUNABLE TRANSMISSION PROFILE

; LIKE HMI_LASER2.PRO BUT WE ONLY FIT FOR THE PHASES
; BECAUSE THE CONTRASTS ARE SO DIFFICULT TO OBTAIN

;----------------------------------------------------------------------

FUNCTION reconstruct,moy,moyb

; moy = average phases of non-tunable elements in degrees
; moyb= average contrasts

lamref      = 6173.3433d0
dpi         = 2.d0*!dpi
FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.172d0-0.0010576d0           ; for the narrow-band Michelson
FSR[1]      = 0.344d0-0.00207683d0             ; for the broad-band  Michelson
FSR[2]      = 0.693d0+0.000483467d0            ; for E1
FSR[3]      = 1.407d0            ; for E2
FSR[4]      = 2.779d0            ; for E3
FSR[5]      = 5.682d0            ; for E4
FSR[6]      = 11.354d0           ; for E5

; RECONSTRUCT THE SPATIALLY-AVERAGED TRANSMISSION PROFILE OF THE NON TUNABLE ELEMENTS:
lam         = FINDGEN(20000)/19999.*30.-15.
RESTORE,'frontwindow.bin'
blocker0    = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q           = READFITS('blocker11.fits')
blocker0    = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.20854-lamref,lam) ; I center the profile ;+3.61328

profileg    = blocker0
FOR i=3,6 DO profileg = profileg * (1.d0+moyb[i-3]*COS(dpi/FSR[i]*lam+moy[i-3]*!pi/180.))/2.d0
;SET_PLOT,'x'
;WINDOW,0,RETAIN=2
!P.MULTI=0
SET_PLOT,'PS'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=20,ysize=16,/color
LOADCT,3
a=where(lam ge -2 and lam le 2)
PLOT,lam,profileg/max(profileg[a]),xst=1,xrange=[-5.25,4.],yrange=[0,1],yst=1,charsize=1.5,tit='!17',xtit='Wavelength (A)',ytit='Transmittance'
oplot,lam,profileg/max(profileg[a]),col=180
oplot,[-7.5,7.5],[0.05,0.05],linestyle=2
;oplot,[-0.5,-0.5],[0,1],linestyle=2
;oplot,[0.5,0.5],[0,1],linestyle=2
oplot,[0,0],[0,1],linestyle=2
oplot,[-7.5,7.5],[0.5,0.5],linestyle=2
DEVICE,/CLOSE
SET_PLOT,'x'
res  = FLTARR(2,20000)
res[0,*] = lam
res[1,*] = profileg;/max(profileg[a])

RETURN,res

END


;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser2_nocontrast,draw


fitE5       = 1 ;fit for contrasts (0=no, 1=yes) ?

;SET_PLOT,'X'
;WINDOW,0,RETAIN=2,xsize=900,ysize=900
;!P.MULTI    = [0,2,2]
time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0        ; Fe I 6173 line central wavelength in air

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 64               ; number of rows
ny          = 64               ; number of columns
factor      = 10000.
nseq        = 27               ; number of positions in the co-tune sequence
nl0         = 24.;16.              ; number of wavelength positions for the dye laser
nparam      = 10                ; number of filter parameters we fit for with the cotune: 4 phases, 4 contrasts, laser intensity
IF (fitE5 EQ 0) THEN nparam = 6
anglim      = 950
xcenter     = nx/2
ycenter     = nx/2
distance    = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;------------------------------------------------------

FSR         = DBLARR(7)          ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]      = 0.1710098d0;0.172d0-0.0010576d0            ; for the narrow-band Michelson
FSR[1]      = 0.3421506d0;0.344d0-0.00207683d0            ; for the broad-band  Michelson
FSR[2]      = 0.6943613d0;0.693d0+0.000483467d0           ; for E1
FSR[3]      = 1.407d0;*0.993d0            ; for E2
FSR[4]      = 2.779d0;*0.993d0           ; for E3
FSR[5]      = 5.682d0;*0.993d0            ; for E4
FSR[6]      = 11.354d0;*0.993d0           ; for E5


; BLOCKER FILTER + FRONT WINDOW AVERAGED PROFILE
; from R. Shine's website (http://www.lmsal.com/~shine/Public/hmi/)
;------------------------------------------------------

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
RESTORE,'frontwindow.bin'
blocker0     = transmission/100.d0;INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam);,/LSQUADRATIC)
q            = READFITS('blocker11.fits')
lam          = wavelength*10.d0-lamref
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.5975d0-lamref,lam) ; I center the profile ;+3.20854

; MEASURED INTENSITIES
;-------------------------------------------------------

; list containing the names of filtergrams
; USE readfile3_128 TO PRODUCE THE INTENSITIES


list    = STRARR(nl0)
expo    = FLTARR(nl0)
lam0    = FLTARR(nl0)
relI    = FLTARR(nl0)

;list[0  ] = 'listLaser071015_710420'
list[1 -1 ] = 'listLaser071015_710600'
list[2 -1 ] = 'listLaser071015_710720'
;list[3 -1 ] = 'listLaser071015_710900'
list[4 -2 ] = 'listLaser071015_711020'
list[5 -2 ] = 'listLaser071015_711140'
list[6 -2 ] = 'listLaser071015_711260'
list[7 -2 ] = 'listLaser071015_711380'
list[8 -2 ] = 'listLaser071015_711500'
list[9 -2 ] = 'listLaser071015_711620'
list[10-2 ] = 'listLaser071015_711740'
list[11-2 ] = 'listLaser071015_711860'
;list[12-1 ] = 'listLaser071015_711980'
;list[13-1 ] = 'listLaser071015_712100'
list[14-4 ] = 'listLaser071015_712220'
list[15-4 ] = 'listLaser071015_712408'
list[16-4 ] = 'listLaser071015_712528'
list[17-4 ] = 'listLaser071015_712708'
list[18-4 ] = 'listLaser071015_712828'
list[19-4 ] = 'listLaser071015_712948'
;list[20-4 ] = 'listLaser071015_713068'
list[21-5 ] = 'listLaser071015_713188'
list[22-5 ] = 'listLaser071015_713308'
list[23-5 ] = 'listLaser071015_713428'
list[24-5 ] = 'listLaser071015_713608'
list[25-5 ] = 'listLaser071015_713728'
;list[26-2 ] = 'listLaser071015_713848'
list[27-6 ] = 'listLaser071015_713968'
list[28-6 ] = 'listLaser071015_714088'
list[29-6 ] = 'listLaser071015_714208'
;lam0[0]  = 6173.3222;       6173.3090    0.0149381      9.89776
lam0[1-1]  = 6173.4612;       6173.4590    0.0163509      6.32749
lam0[2-1]  = 6173.7014;       6173.6991    0.0299019      2.60589
;lam0[3-1]  = 6173.8458;       6173.8387    0.0702081      6.90940
lam0[4-2]  = 6174.1000;       6174.0989    0.0429304     0.496411
lam0[5-2]  = 6174.2971;       6174.2985    0.0391197      1.74875
lam0[6-2]  = 6174.6419;       6174.6388    0.0261406     0.688459
lam0[7-2]  = 6174.9020;       6174.8985    0.0216780     0.718730
lam0[8-2]  = 6175.8128;       6175.8084    0.0369006     0.564670
lam0[9-2]  = 6176.0073;       6176.0085    0.0200389     0.331574
lam0[10-2] = 6176.3070;       6176.3082    0.0455064    0.0417319
lam0[11-2] = 6176.5073;       6176.5084    0.0683433    0.0449396
;lam0[12-1] = 6178.6021;       6178.6074    0.0317385   0.00963807
;lam0[13-1] = 6178.9987;       6179.0077    0.0564096   0.00461222
lam0[14-4] = 6175.4288;       6175.4286    0.0525492    0.0505462
lam0[15-4] = 6173.2976;       6173.2988    0.0166755      5.17158
lam0[16-4] = 6173.1124;       6173.1089    0.0219495      2.81007
lam0[17-4] = 6172.9099;       6172.9093    0.0217191      4.92253
lam0[18-4] = 6172.7360;       6172.7394    0.0585246      1.06826
lam0[19-4] = 6172.5099;       6172.5095    0.0663397      1.30987
;lam0[20-4] = 6172.3069;       6172.2991    0.0263406      2.51027
lam0[21-5] = 6172.0090;       6172.0096    0.0169810     0.535824
lam0[22-5] = 6171.7069;       6171.7094    0.0140644     0.609684
lam0[23-5] = 6170.8131;       6170.8097    0.0292688     0.569126
lam0[24-5] = 6170.6130;       6170.6096    0.0186601     0.470184
lam0[25-5] = 6170.3073;       6170.3099    0.0382485     0.389437
;lam0[26-2] = 6168.0092;       6168.0107    0.0532376     0.196522
lam0[27-6] = 6167.6059;       6167.6104    0.0376108    0.0713017
lam0[28-6] = 6171.3060;       6171.3096    0.0505170    0.0685530
lam0[29-6] = 6173.6151;       6173.6088    0.0261074      3.86625

;lam0[0] = 6173.3210;6173.3090d0
;lam0[1] = 6173.4606;6173.4590d0
;lam0[2] = 6173.7022;6173.6991d0
;lam0[3] = 6173.8474;6173.8387d0
;lam0[4] = 6174.1031;6174.0989d0
;lam0[5] = 6174.3010;6174.2985d0
;lam0[6] = 6174.6477;6174.6388d0
;lam0[7] = 6174.9091;6174.8985d0
;lam0[8] = 6175.8247;6175.8084d0
;lam0[9] = 6176.0205;6176.0085d0
;lam0[10]= 6176.3217;6176.3082d0
;lam0[11]= 6176.5232;6176.5084d0
;lam0[12]= 6178.6292;6178.6074d0
;lam0[13]= 6179.0277;6179.0077d0
;lam0[14]= 6175.4389;6175.4286d0
;lam0[15]= 6173.2965;6173.2988d0
;lam0[16]= 6173.1101;6173.1089d0
;lam0[17]= 6172.9065;6172.9093d0
;lam0[18]= 6172.7318;6172.7394d0
;lam0[19]= 6172.5041;6172.5095d0
;lam0[20]= 6172.3004;6172.2991d0
;lam0[21]= 6172.0009;6172.0096d0
;lam0[22]= 6171.6968;6171.7094d0
;lam0[23]= 6170.7986;6170.8097d0
;lam0[24]= 6170.5974;6170.6096d0
;lam0[25]= 6170.2900;6170.3099d0
;lam0[26]= 6167.9792;6168.0107d0
;lam0[27]= 6167.5744;6167.6104d0
;lam0[28]= 6171.2943;6171.3096d0
;lam0[29]= 6173.6152;6173.6088d0
;lam0=[6173.3090d0,6173.4590d0,6173.6991d0,6173.8387d0,6174.0989d0,6174.2985d0,6174.6388d0,6174.8985d0,6175.8084d0,6176.0085d0,6176.3082d0,6176.5084d0,6178.6074d0,6179.0077d0,6175.4286d0,6173.2988d0,6173.1089d0,6172.9093d0,6172.7394d0,6172.5095d0,6172.2991d0,6172.0096d0,6171.7094d0,6170.8097d0,6170.6096d0,6170.3099d0,6168.0107d0,6167.6104d0,6171.3096d0,6173.6088d0]
;lam0=[6173.4590d0,6173.6991d0,6173.8387d0,6174.0989d0,6174.2985d0,6174.6388d0,6174.8985d0,6175.8084d0,6176.0085d0,6176.3082d0,6176.5084d0,6178.6074d0,6179.0077d0,6175.4286d0,6173.2988d0,6173.1089d0,6172.9093d0,6172.7394d0,6172.5095d0,6172.2991d0,6172.0096d0,6171.7094d0,6170.8097d0,6170.6096d0,6170.3099d0,6168.0107d0,6167.6104d0,6171.3096d0,6173.6088d0]
;lam0=[6173.4590d0,6173.6991d0,6173.8387d0,6174.0989d0,6174.2985d0,6174.6388d0,6174.8985d0,6175.8084d0,6176.0085d0,6176.3082d0,6176.5084d0,6178.6074d0,6179.0077d0,6175.4286d0,6173.1089d0,6172.9093d0,6172.7394d0,6172.5095d0,6172.2991d0,6172.0096d0,6171.7094d0,6170.8097d0,6170.6096d0,6170.3099d0,6168.0107d0,6167.6104d0,6171.3096d0,6173.6088d0]

lam0    = lam0-lamref
;relI[0] = 4.0
;relI[1] = 5.6
;relI[2] = 4.8
;relI[3] = 4.1
;relI[4] = 4.3
;relI[5] = 5.15
;relI[6] = 4.85
;relI[7] = 4.4
;relI[8] = 4.7
;relI[9] = 6.0
;relI[10]= 6.2
;relI[11]= 5.2
;relI[12]= 7.05
;relI[13]= 6.6
;relI[14]= 5.1
;relI[15]= 7.3
;relI[16]= 8.1
;relI[17]= 6.2
;relI[18]= 6.5
;relI[19]= 7.3
;relI[20]= 5.4
;relI[21]= 8.0
;relI[22]= 8.05
;relI[23]= 9.3
;relI[24]= 6.5
;relI[25]= 6.0
;relI[26]= 5.4
;relI[27]= 3.6
;relI[28]= 6.0
;relI[29]= 5.35
;expo[0] = 705.
expo[1 -1] = 505.
expo[2 -1] = 505.
;expo[3 -1] = 3755.
expo[4 -2] = 3755.
expo[5 -2] = 3755.
expo[6 -2] = 3755.
expo[7 -2] = 3755.
expo[8 -2] = 3755.
expo[9 -2] = 3755.
expo[10-2]= 3755.
expo[11-2]= 3755.
;expo[12-1]= 3755.
;expo[13-1]= 3755.
expo[14-4]= 3755.
expo[15-4]= 405.
expo[16-4]= 405.
expo[17-4]= 2005.
expo[18-4]= 3755.
expo[19-4]= 3755.
;expo[20-4]= 3755.
expo[21-5]= 3755.
expo[22-5]= 3755.
expo[23-5]= 3755.
expo[24-5]= 3755.
expo[25-5]= 3755.
;expo[26-2]= 3755.
expo[27-6]= 3755.
expo[28-6]= 3755.
expo[29-6]= 505.


;;list[0] = 'listLaser071015_710480'
;list[1-1] = 'listLaser071015_710660'
;list[2-1] = 'listLaser071015_710780'
;list[3-1] = 'listLaser071015_710960'
;list[4-1] = 'listLaser071015_711080'
;list[5-1] = 'listLaser071015_711200'
;list[6-1] = 'listLaser071015_711320'
;list[7-1] = 'listLaser071015_711440'
;list[8-1] = 'listLaser071015_711560'
;list[9-1] = 'listLaser071015_711680'
;list[10-1]= 'listLaser071015_711800'
;list[11-1]= 'listLaser071015_711920'
;list[12-1]= 'listLaser071015_712040'
;list[13-1]= 'listLaser071015_712160'
;list[14-1]= 'listLaser071015_712280'
;list[15-1]= 'listLaser071015_712468'
;list[16-1]= 'listLaser071015_712588'
;list[17-1]= 'listLaser071015_712768'
;list[18-1]= 'listLaser071015_712888'
;list[19-1]= 'listLaser071015_713008'
;list[20-1]= 'listLaser071015_713128'
;list[21-1]= 'listLaser071015_713248'
;list[22-1]= 'listLaser071015_713368'
;list[23-1]= 'listLaser071015_713488'
;list[24-1]= 'listLaser071015_713668'
;list[25-1]= 'listLaser071015_713788'
;list[26-1]= 'listLaser071015_713908'
;list[27-1]= 'listLaser071015_714028'
;list[28-1]= 'listLaser071015_714148'
;list[29-1]= 'listLaser071015_714268'
;lam0=[6173.4590d0,6173.6991d0,6173.8387d0,6174.0989d0,6174.2985d0,6174.6388d0,6174.8985d0,6175.8084d0,6176.0085d0,6176.3082d0,6176.5084d0,6178.6074d0,6179.0077d0,6175.4286d0,6173.2988d0,6173.1089d0,6172.9093d0,6172.7394d0,6172.5095d0,6172.2991d0,6172.0096d0,6171.7094d0,6170.8097d0,6170.6096d0,6170.3099d0,6168.0107d0,6167.6104d0,6171.3096d0,6173.6088d0]
;lam0    = lam0-lamref
;;expo[0] = 705.
;expo[1-1] = 505.
;expo[2-1] = 1005.
;expo[3-1] = 3755.
;expo[4-1] = 3755.
;expo[5-1] = 3755.
;expo[6-1] = 3755.
;expo[7-1] = 3755.
;expo[8-1] = 3755.
;expo[9-1] = 3755.
;expo[10-1]= 3755.
;expo[11-1]= 3755.
;expo[12-1]= 3755.
;expo[13-1]= 3755.
;expo[14-1]= 3755.
;expo[15-1]= 405.
;expo[16-1]= 405.
;expo[17-1]= 2005.
;expo[18-1]= 3755.
;expo[19-1]= 3755.
;expo[20-1]= 3755.
;expo[21-1]= 3755.
;expo[22-1]= 3755.
;expo[23-1]= 3755.
;expo[24-1]= 3755.
;expo[25-1]= 3755.
;expo[26-1]= 3755.
;expo[27-1]= 3755.
;expo[28-1]= 3755.
;expo[29-1]= 605.


;list[0 ]= 'listLaser070905_549188'
;list[1 ]= 'listLaser070905_549324'
;list[2 ]= 'listLaser070905_549596'
;list[3 ]= 'listLaser070905_549732'
;list[4 ]= 'listLaser070905_549868'
;list[5 ]= 'listLaser070905_550004'
;list[6 ]= 'listLaser070905_550412'
;list[7 ]= 'listLaser070905_550548'
;list[8 ]= 'listLaser070905_550820'
;list[9 ]= 'listLaser070905_550956'
;list[10]= 'listLaser070905_551092'
;list[11]= 'listLaser070905_551364'
;list[12]= 'listLaser070905_551500'
;list[13]= 'listLaser070905_551636'
;list[14]= 'listLaser070905_551772'
;list[15]= 'listLaser070905_551908'
;lam0[0 ] = 6173.0491;6172.3566
;lam0[1 ] = 6172.9141;6172.2250
;lam0[2 ] = 6172.4215;6171.7294
;lam0[3 ] = 6172.1727;6171.4854
;lam0[4 ] = 6171.6754;6170.9934
;lam0[5 ] = 6171.1905;6170.5000
;lam0[6 ] = 6169.6988;6169.0295
;lam0[7 ] = 6173.4235;6172.7279
;lam0[8 ] = 6173.6538;6172.9519
;lam0[9 ] = 6173.7778;6173.0754
;lam0[10] = 6173.9083;6173.2033
;lam0[11] = 6174.4087;6173.6985
;lam0[12] = 6174.6634;6173.9552
;lam0[13] = 6174.8915;6174.1862
;lam0[14] = 6175.1511;6174.4400
;lam0[15] = 6175.6438;6174.9383
;lam0     = lam0-lamref;+0.693d0  !!!WARNING!!!
;relI[0 ] = 4.56688 ;4.60826
;relI[1 ] = 2.95802; 2.96647
;relI[2 ] = 3.43708; 3.43317
;relI[3 ] = 4.29308; 4.31731
;relI[4 ] = 3.13255; 3.12851
;relI[5 ] = 2.49317; 2.49997
;relI[6 ] = 4.04419; 4.05959
;relI[7 ] = 3.97266; 3.98424
;relI[8 ] = 3.1608; 3.15832
;relI[9 ] = 3.26378; 3.27738
;relI[10] = 2.8725; 2.86727
;relI[11] = 2.6067; 2.57554
;relI[12] = 2.96883; 2.95981
;relI[13] = 3.15652; 3.13853
;relI[14] = 2.93554; 2.93016
;relI[15] = 2.81568; 2.80612
;expo[0]  = 1005.
;expo[1]  = 1005.
;expo[2]  = 3755.
;expo[3]  = 3755.
;expo[4]  = 3755.
;expo[5]  = 3755.
;expo[6]  = 3755.
;expo[7]  = 1005.
;expo[8]  = 1005.
;expo[9]  = 1005.
;expo[10] = 3755.
;expo[11] = 3755.
;expo[12] = 3755.
;expo[13] = 3755.
;expo[14] = 3755.
;expo[15] = 3755.


;list[0] = 'listLaser070905_549256'
;list[1] = 'listLaser070905_549392'
;list[2] = 'listLaser070905_549664'
;list[3] = 'listLaser070905_549800'
;list[4] = 'listLaser070905_549936'
;list[5] = 'listLaser070905_550072'
;list[6] = 'listLaser070905_550480'
;list[7] = 'listLaser070905_550616'
;list[8 ]= 'listLaser070905_550888'
;list[9 ]= 'listLaser070905_551024'
;list[10]= 'listLaser070905_551160'
;list[11]= 'listLaser070905_551432'
;list[12]= 'listLaser070905_551568'
;list[13]= 'listLaser070905_551704'
;list[14]= 'listLaser070905_551840'
;list[15]= 'listLaser070905_551976'
;lam0[0] = 6173.0489            ; 0.0245983
;lam0[1] = 6172.9138            ; 0.0301045
;lam0[2] = 6172.4212            ; 0.0246679
;lam0[3] = 6172.1727            ; 0.0227419
;lam0[4] = 6171.6749            ; 0.0671509
;lam0[5] = 6171.1903            ; 0.0570870
;lam0[6] = 6169.6986            ;  0.102417
;lam0[7] = 6173.4235            ; 0.0336450
;lam0[8] = 6173.6535            ; 0.0538418
;lam0[9] = 6173.7776            ; 0.0525744
;lam0[10]= 6173.9083            ; 0.0870561
;lam0[11]= 6174.4092            ; 0.0642272
;lam0[12]= 6174.6632            ; 0.0770916
;lam0[13]= 6174.8915            ; 0.0963618
;lam0[14]= 6175.1514            ;  0.115435
;lam0[15]= 6175.6437            ;  0.115777
;lam0     = lam0-lamref
;relI[*] = 0.d0
;expo[0]  = 1005.
;expo[1]  = 1005.
;expo[2]  = 3755.
;expo[3]  = 3755.
;expo[4]  = 3755.
;expo[5]  = 3755.
;expo[6]  = 3755.
;expo[7]  = 1005.
;expo[8]  = 1005.
;expo[9]  = 1005.
;expo[10] = 3755.
;expo[11] = 3755.
;expo[12] = 3755.
;expo[13] = 3755.
;expo[14] = 3755.
;expo[15] = 3755.


nx2     = 256
ny2     = 256
Inten   = DBLARR(nx2,ny2,nseq*nl0) ; measured output intensities (measured on a HMI CCD) with the cotune
FOR i=0,nl0-1 DO BEGIN
    RESTORE,'CPT/SEQUENCE_'+STRTRIM(list[i],1)+'_'+STRTRIM(STRING(LONG(nx2)),1)+'.BIN'
    Inten[*,*,nseq*i:nseq*(i+1)-1] = imx[*,*,2:28]/factor
ENDFOR
Intena   = REBIN(Inten,nx,ny,nseq*nl0)
;RESTORE,'CORRECTION_INTENSITY_549188.BIN'; OBTAINED BY intensity_correction.pro RECOMPUTE EACH TIME !!!!
;RESTORE,'CORRECTION_INTENSITY_710420.BIN' ;OBSMODE
 RESTORE,'CORRECTION_INTENSITY_710600.BIN' ;OBSMODE WITH FIRST DETUNE REMOVED
;RESTORE,'CORRECTION_INTENSITY_710720.BIN' ;OBSMODE WITH FIRST 2 DETUNES REMOVED
;RESTORE,'CORRECTION_INTENSITY_710480.BIN' ;CALMODE
;RESTORE,'CORRECTION_INTENSITY_710660.BIN' ;CALMODE WITH FIRST DETUNE REMOVED
INTENSITYdetunes=INTENSITYdetunes-0.031d-9 ; background intensity


; WE NORMALIZE THE INTENSITIES
FOR i=0,nl0-1 DO BEGIN
    FOR j=0,nseq-1 DO BEGIN
       ;Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/relI[i]*MEAN(relI)/expo[i]*MEAN(expo)
        Inten[*,*,i*nseq+j]=Inten[*,*,i*nseq+j]/INTENSITYdetunes[i*nseq+j]*MEAN(INTENSITYdetunes)/expo[i]*MEAN(expo)
    ENDFOR
ENDFOR

; WE REQUIRE INTENSITIES TO BE LARGER THAN 0 !!! BIAS !!!
;temp = REFORM(REBIN(inten,1,1,nseq*nl0))
;IF(MIN(temp) LT 0.0) THEN inten= inten-MIN(temp)

prof=fltarr(64,64,nl0)
profa=prof
for i=0,nl0-1 do begin 
    prof[*,*,i]=rebin(Inten[*,*,nseq*i:nseq*(i+1)-1],64,64,1)
    profa[*,*,i]=rebin(Intena[*,*,nseq*i:nseq*(i+1)-1],64,64,1)
endfor    
prof=prof[*,*,sort(lam0)]
profa=profa[*,*,sort(lam0)]
moncul=lam0[sort(lam0)]
set_plot,'x'
!p.multi=0
window,0,retain=2
plot,moncul+lamref,REBIN(prof,1,1,nl0),xst=1,psym=4
read,pause

Inten   = REBIN(Inten,nx,ny,nseq*nl0)

; VARIABLE DEFINITION
;--------------------------------------------------------

Bg          = DBLARR(nx,ny,7)  ; contrasts of the Lyot and Michelson elements
Phig        = DBLARR(nx,ny,7)  ; relative phases
I0g         = DBLARR(nx,ny)    ; laser "intensity"
threshold   = DBLARR(nx,ny)
Residual    = DBLARR(nseq*nl0) ; residual for the least-squares fit
dIntendB3   = DBLARR(nseq*nl0) ; derivatives to compute the Jacobian matrix
dIntendB4   = DBLARR(nseq*nl0)
dIntendB5   = DBLARR(nseq*nl0)
dIntendB6   = DBLARR(nseq*nl0)
dIntendPhi3 = DBLARR(nseq*nl0)
dIntendPhi4 = DBLARR(nseq*nl0)
dIntendPhi5 = DBLARR(nseq*nl0)
dIntendPhi6 = DBLARR(nseq*nl0)
dIntendInt  = DBLARR(nseq*nl0)
dIntendThr  = DBLARR(nseq*nl0)
errorbar    = DBLARR(nseq*nl0)
lateconv    = INTARR(nx,ny)    ; to discard spurious results
maxsteps    = 5555
Jac         = DBLARR(nparam,nseq*nl0) ; Jacobian matrix of the least-squares fit
contrasts   = 0.95d0


; DEFINITION OF THE CO-TUNE SEQUENCE
; co-tuning in the range [-2 FSR[E1],+2 FSR[E1]]
; TABLE PROVIDED BY JESPER
;----------------------------------------------------------

tuningJS       = DBLARR(nseq,3)
tuningJS[0 ,*] = [         0.d0,          0.d0,          0.d0]
tuningJS[1 ,*] = [        80.d0,          0.d0,          0.d0]
tuningJS[2 ,*] = [       160.d0,          0.d0,          0.d0]
tuningJS[3 ,*] = [         0.d0,         80.d0,          0.d0]
tuningJS[4 ,*] = [        80.d0,         80.d0,          0.d0]
tuningJS[5 ,*] = [       160.d0,         80.d0,          0.d0]
tuningJS[6 ,*] = [         0.d0,        160.d0,          0.d0]
tuningJS[7 ,*] = [        80.d0,        160.d0,          0.d0]
tuningJS[8 ,*] = [       160.d0,        160.d0,          0.d0]
tuningJS[9 ,*] = [         0.d0,          0.d0,         80.d0]
tuningJS[10,*] = [        80.d0,          0.d0,         80.d0]
tuningJS[11,*] = [       160.d0,          0.d0,         80.d0]
tuningJS[12,*] = [         0.d0,         80.d0,         80.d0] 
tuningJS[13,*] = [        80.d0,         80.d0,         80.d0]
tuningJS[14,*] = [       160.d0,         80.d0,         80.d0]
tuningJS[15,*] = [         0.d0,        160.d0,         80.d0]
tuningJS[16,*] = [        80.d0,        160.d0,         80.d0]
tuningJS[17,*] = [       160.d0,        160.d0,         80.d0]
tuningJS[18,*] = [         0.d0,          0.d0,        160.d0]
tuningJS[19,*] = [        80.d0,          0.d0,        160.d0]
tuningJS[20,*] = [       160.d0,          0.d0,        160.d0]
tuningJS[21,*] = [         0.d0,         80.d0,        160.d0]
tuningJS[22,*] = [        80.d0,         80.d0,        160.d0]
tuningJS[23,*] = [       160.d0,         80.d0,        160.d0]
tuningJS[24,*] = [         0.d0,        160.d0,        160.d0]
tuningJS[25,*] = [        80.d0,        160.d0,        160.d0]
tuningJS[26,*] = [       160.d0,        160.d0,        160.d0]

FOR i=0,nseq-1 DO tuningJS[i,*]  = tuningJS[i,*] * [6.d0,6.d0,-6.d0]*!dpi/180.d0


IF draw EQ 1 THEN GOTO,draw


; PHASE AND CONTRAST OF TUNING ELEMENTS
;-----------------------------------------------------

;RESTORE,'CPT/CPT_laser_front_710600.BIN' ; OBSMODE
;RESTORE,'CPT/CPT_laser_front_712408.BIN' ; OBSMODE
;RESTORE,'CPT/CPT_laser_front_712408_NEWFSR.BIN' ; OBSMODE NEW FSRs
;RESTORE,'CPT/CPT_laser_front_712468_NEWFSR.BIN' ; CALMODE NEW FSRs
;RESTORE,'CPT/CPT_laser_front_549732.BIN' ; OBSMODE
;RESTORE,'CPT/CPT_laser_front_549800.BIN' ; CALMODE
RESTORE,'CPT/CPT_laser_front_AVERAGE_NEWFSR.BIN' ; created by average_phase_contrast.pro

a=WHERE(FINITE(Bg0) EQ 0)
IF (a[0] NE -1) THEN Bg0[a]  =1.d0
a=WHERE(FINITE(Phig0) EQ 0)
IF (a[0] NE -1) THEN Phig0[a]=0.d0

;FOR iii=0,nx-1 DO FOR jjj=0,nx-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]+[-106.26-169.22735,42.35+4.2600189,-140.32+164.03892]*!pi/180.d0 ; 79336
;FOR iii=0,nx-1 DO FOR jjj=0,nx-1 DO Phig0[iii,jjj,*]=Phig0[iii,jjj,*]+[-106.26-170.95412,42.35+2.8555888,-140.32+158.35755]*!pi/180.d0 ; 79916


Bg0 = REBIN(Bg0,nx,nx,3)
Phig0=REBIN(Phig0,nx,nx,3)

Bg[*,*,0:2]   = Bg0[*,*,*]         ; contrasts and phases of the tuning elements
Phig[*,*,0:2] = Phig0[*,*,*]

;----------------------------------------------------

;RESTORE,'I0g_549188.BIN' ; WARNING WARNING WARNING WARNING WARNING !!!!!!!!!!!!! for I0g in Obsmode
;I0g=I0g/100.d0 ; for factor=100
;RESTORE,'temp2.bin' ; WARNING WARNING WARNING WARNING WARNING !!!!!!!!!!!!! for I0g in Calmode

;In=FLTARR(nparam,nparam)
;FOR i=0,nparam-1 DO In[i,i]=1.d0

FOR jjj=0,ny-1 DO BEGIN

   ;TVIM,Phig[*,*,3],/scale 
   ;TVIM,Phig[*,*,4],/scale
   ;TVIM,Phig[*,*,5],/scale
   ;TVIM,Phig[*,*,6],/scale
    TVIM,lateconv

    FOR iii=0,nx-1 DO BEGIN

   ; BEGINNING OF THE LEAST-SQUARES FIT
   ;-----------------------------------


        IF(distance[iii,jjj] LE anglim) THEN BEGIN

           ;Bg[iii,jjj,3:6]  = [0.9,0.99,1.3,0.6] ; contrasts
           ;Phig[iii,jjj,3:6]= [-6.7,0.08,-1.45,177.]/180.d0*!dpi 
            Phig[iii,jjj,3:6]= [0.,0.,0.,0.]/180.d0*!dpi
            Bg[iii,jjj,3:6]  = [1.0,1.0,1.0,1.0]
            IF (fitE5 EQ 0) THEN BEGIN
                Phig[iii,jjj,3:6]= [0.0,0.0,0.0,0.0]/180.d0*!dpi 
                Bg[iii,jjj,3:6]  = [1.0,1.0,1.0,1.0]
            ENDIF

            ; guess of the initial intensity
           ;sum = 0.d0
           ;FOR j=0,nseq-1 DO BEGIN
           ;    FOR ii=0,nl0-1 DO BEGIN
           ;        profileg  = INTERPOL(blocker0,lam,lam0[ii])
           ;        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
           ;        FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0
           ;        sum = sum+profileg
           ;    ENDFOR
           ;ENDFOR
           ;thresh       = TOTAL(Inten[iii,jjj,*])/sum
           ;I0g[iii,jjj] = thresh


            ;sum = FLTARR(nl0)
            ;FOR ii=0,nl0-1 DO BEGIN
            ;ii=10
            ;    profileg  = INTERPOL(blocker0,lam,lam0[ii])
            ;    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]))/2.d0
            ;    sum = 8.d0*TOTAL(Inten[iii,jjj,*])/27.d0/profileg
           ;ENDFOR
            thresh       = 7.;1.4;sum/nl0
           ;thresh=I0g[iii,jjj];1050.;4000.
            I0g[iii,jjj] = thresh
            ;PRINT,'yo',SIGMA(sum)


            converg = 0
            jj      = 0

            WHILE converg EQ 0 DO BEGIN

                FOR j=0,nseq-1 DO BEGIN
                    FOR ii=0,nl0-1 DO BEGIN

                        profileg  = INTERPOL(blocker0,lam,lam0[ii],/QUADRATIC)
                        FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]+tuningJS[j,i]))/2.d0
                        profilegf = profileg
                        FOR i=3,6 DO profileg = profileg * (1.d0-Bg[iii,jjj,i]*SIN(!dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]/2.d0)^2.d0)
                        errorbar[ii*nseq+j]  = 1.d0;rmsINTENSITY[ii]*Inten[iii,jjj,ii*nseq+j]/MEAN(INTENSITYdetunes[ii*nseq:(ii+1)*nseq-1])+.001d0

                        Residual[ii*nseq+j]  = (Inten[iii,jjj,ii*nseq+j] - (profileg*I0g[iii,jjj]+threshold[iii,jjj]))/errorbar[ii*nseq+j]
            
                        IF(fitE5 EQ 1) THEN BEGIN
                            dIntendB3[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB4[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]$
                       *SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB5[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                            dIntendB6[ii*nseq+j] = profilegf*(-SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,4]$
                       *SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)$
                       *(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]
                        ENDIF

                        dIntendPhi3[ii*nseq+j]= profilegf*(-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)*COS(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi4[ii*nseq+j]= profilegf*(-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)*COS(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0))$
                       *(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi5[ii*nseq+j]= profilegf*(-Bg[iii,jjj,5]*SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)*COS(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]$
                       *SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendPhi6[ii*nseq+j]= profilegf*(-Bg[iii,jjj,6]*SIN(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0)*COS(!dpi/FSR[6]*lam0[ii]+Phig[iii,jjj,6]/2.d0))$
                       *(1.d0-Bg[iii,jjj,4]*SIN(!dpi/FSR[4]*lam0[ii]+Phig[iii,jjj,4]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,5]$
                       *SIN(!dpi/FSR[5]*lam0[ii]+Phig[iii,jjj,5]/2.d0)^2.d0)*(1.d0-Bg[iii,jjj,3]*SIN(!dpi/FSR[3]*lam0[ii]+Phig[iii,jjj,3]/2.d0)^2.d0)*I0g[iii,jjj]/errorbar[ii*nseq+j]

                        dIntendInt[ii*nseq+j] = profileg/errorbar[ii*nseq+j]

                        dIntendThr[ii*nseq+j] = 1.d0

                    ENDFOR
                ENDFOR

              ; Computation of the Jacobian matrix
                FOR i= 0,nseq*nl0-1 DO BEGIN
                    IF(fitE5 NE 0) THEN BEGIN
                        Jac[0,i]  = dIntendB3[i]
                        Jac[1,i]  = dIntendB4[i]
                        Jac[2,i]  = dIntendB5[i]
                        Jac[3,i]  = dIntendB6[i]
                        Jac[4,i]  = dIntendPhi3[i]
                        Jac[5,i]  = dIntendPhi4[i]
                        Jac[6,i]  = dIntendPhi5[i]
                        Jac[7,i]  = dIntendPhi6[i]
                        Jac[8,i]  = dIntendInt[i]
                        Jac[9,i]  = dIntendThr[i]
                    ENDIF ELSE BEGIN
                        Jac[0,i]  = dIntendPhi3[i]
                        Jac[1,i]  = dIntendPhi4[i]
                        Jac[2,i]  = dIntendPhi5[i]
                        Jac[3,i]  = dIntendPhi6[i]
                        Jac[4,i]  = dIntendInt[i]
                        Jac[5,i]  = dIntendThr[i]
                    ENDELSE
                ENDFOR               

                TJac    = TRANSPOSE(Jac)
                Hessian = TJac##Jac;+lambda*In

               ;LA_SVD,Jac,W,U,V,/DOUBLE ; Singular Value Decomposition using the LAPACK algorithm
               ;LA_SVD,Hessian,W,U,V,/DOUBLE 
                Hessiani=LA_INVERT(Hessian,/DOUBLE)
                
               ;Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual
               ;Dx     = (V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U))##TJac##Residual
                Dx     =  Hessiani##TJac##Residual

                IF (fitE5 NE 0) THEN err = TRANSPOSE(Dx)/[REFORM(Bg[iii,jjj,3]),REFORM(Bg[iii,jjj,4]),REFORM(Bg[iii,jjj,5]),REFORM(Bg[iii,jjj,6]),REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj],threshold[iii,jjj]] ELSE err = TRANSPOSE(Dx)/[REFORM(Phig[iii,jjj,3]),REFORM(Phig[iii,jjj,4]),REFORM(Phig[iii,jjj,5]),REFORM(Phig[iii,jjj,6]),I0g[iii,jjj],threshold[iii,jjj]]

                err = MAX(ABS(err))

                IF(jj EQ maxsteps-1 AND err GT 1.d-3 )      THEN          converg         = 2 ;d-7
                IF(err LE 1.d-3)                            THEN          converg         = 1

                IF( fitE5 NE 0) THEN BEGIN
                    Bg[iii,jjj,3]   = Bg[iii,jjj,3]  +Dx[0]
                    IF(Bg[iii,jjj,3] LT 0.d0 OR Bg[iii,jjj,3] GT 10.4d0) THEN Bg[iii,jjj,3]   = contrasts ;10.4
                    Bg[iii,jjj,4]   = Bg[iii,jjj,4]  +Dx[1]
                    IF(Bg[iii,jjj,4] LT 0.d0 OR Bg[iii,jjj,4] GT 10.4d0) THEN Bg[iii,jjj,4]   = contrasts
                    Bg[iii,jjj,5]   = Bg[iii,jjj,5]  +Dx[2]
                    IF(Bg[iii,jjj,5] LT 0.d0 OR Bg[iii,jjj,5] GT 10.4d0) THEN Bg[iii,jjj,5]   = contrasts
                    Bg[iii,jjj,6]   = Bg[iii,jjj,6]  +Dx[3]
                    IF(Bg[iii,jjj,6] LT 0.d0 OR Bg[iii,jjj,6] GT 10.4d0) THEN Bg[iii,jjj,6]   = contrasts
                    Phig[iii,jjj,3] = Phig[iii,jjj,3]+Dx[4]
                    IF(ABS(Phig[iii,jjj,3]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,3] = 0.d0
                    Phig[iii,jjj,4] = Phig[iii,jjj,4]+Dx[5]
                    IF(ABS(Phig[iii,jjj,4]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,4] = 0.d0
                    Phig[iii,jjj,5] = Phig[iii,jjj,5]+Dx[6]
                    IF(ABS(Phig[iii,jjj,5]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,5] = 0.d0
                    Phig[iii,jjj,6] = Phig[iii,jjj,6]+Dx[7]
                    IF(ABS(Phig[iii,jjj,6]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,6] = 0.d0
                    I0g[iii,jjj]    = I0g[iii,jjj]   +Dx[8]
                    IF(I0g[iii,jjj] LT 0.2*thresh OR I0g[iii,jjj] GT 4.0*thresh) THEN I0g[iii,jjj] = thresh
                    threshold[iii,jjj]=threshold[iii,jjj]+Dx[9]
                    IF(threshold[iii,jjj] LT -0.5) THEN threshold[iii,jjj] = 0.0d0
                ENDIF ELSE BEGIN
                    Phig[iii,jjj,3] = Phig[iii,jjj,3]+Dx[0]
                    IF(ABS(Phig[iii,jjj,3]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,3] = 0.d0
                    Phig[iii,jjj,4] = Phig[iii,jjj,4]+Dx[1]
                    IF(ABS(Phig[iii,jjj,4]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,4] = 0.d0
                    Phig[iii,jjj,5] = Phig[iii,jjj,5]+Dx[2]
                    IF(ABS(Phig[iii,jjj,5]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,5] = 0.d0
                    Phig[iii,jjj,6] = Phig[iii,jjj,6]+Dx[3]
                    IF(ABS(Phig[iii,jjj,6]) GT 1.20d0*!dpi)              THEN Phig[iii,jjj,6] = 0.d0
                    I0g[iii,jjj]    = I0g[iii,jjj]   +Dx[4]
                    IF(I0g[iii,jjj] LT 0.2*thresh OR I0g[iii,jjj] GT 4.0*thresh) THEN I0g[iii,jjj] = thresh
                    threshold[iii,jjj]=threshold[iii,jjj]+Dx[5]
                    IF(threshold[iii,jjj] LT -0.5) THEN threshold[iii,jjj] = 0.0d0
                ENDELSE


                lateconv[iii,jjj] = converg
                jj = jj+1

            ENDWHILE

            PRINT,' '
            PRINT,'PIXEL',iii,jjj
            PRINT,'PARAMETERS='
            PRINT,'CONTRASTS'
            Bg[iii,jjj,3:6]=Bg[iii,jjj,3:6]/(2.d0-Bg[iii,jjj,3:6])
            FOR i=3,6 DO PRINT,'Bgi =',Bg[iii,jjj,i]
            PRINT,'PHASES'
            Phig[iii,jjj,*] = Phig[iii,jjj,*]*180.d0/!dpi  ; in degrees
            FOR i=3,6 DO PRINT,'Phigi=',Phig[iii,jjj,i]
            PRINT,'LASER INTENSITY'
            I0g[iii,jjj]=I0g[iii,jjj]*2.d0/(1.d0+Bg[iii,jjj,3])*2.d0/(1.d0+Bg[iii,jjj,4])*2.d0/(1.d0+Bg[iii,jjj,5])*2.d0/(1.d0+Bg[iii,jjj,6])
            PRINT,'I0g=',I0g[iii,jjj]
            PRINT,'threshold=',threshold[iii,jjj]
            PRINT,'lateconv=',converg            
            
        ENDIF

    ENDFOR
ENDFOR

SAVE,Bg,Phig,I0g,lateconv,threshold,FILE='RESULTS/RESULTS2_May08_710420_712100_removed'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
PRINT,'ELAPSED TIME=',SYSTIME(1)-time0

draw:

;-----------------------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS2_May08_710420_712100_removed'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

a=WHERE(distance GT anglim OR lateconv NE 1,COMPLEMENT=b) ; depends on the size of the target ;963 obsmode; 925 calmode

FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & temp[a]=-10000.d0 & phig[*,*,i]=temp[*,*] & temp=REFORM(Bg[*,*,i]) & temp[a]=-10000.d0 & Bg[*,*,i]=temp[*,*] & phig[0:1,*,i]=-10000.d0 & Bg[0:1,*,i]=-10000.d0 & ENDFOR

a = WHERE(FINITE(Phig) EQ 0)
IF(a[0] NE -1) THEN Phig[a]= -10000.d0
a = WHERE(FINITE(Bg) EQ 0)
IF(a[0] NE -1) THEN Bg[a]  = -10000.d0


; ADD 180 DEGREES TO PHASES < -180
FOR i=3,6 DO BEGIN
temp = REFORM(Phig[*,*,i])
    a= WHERE(temp NE -10000.d0)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-360.d0
        Phig[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -150.)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+360.d0
        Phig[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(4)
moyb=moy
FOR i=3,6 DO BEGIN & temp=REFORM(Phig[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moy[i-3]=MEAN(temp[a]) & temp[b]=moy[i-3] & Phig[*,*,i]=temp & temp=REFORM(Bg[*,*,i]) & a=WHERE(temp NE -10000.d0,COMPLEMENT=b) & moyb[i-3]=MEAN(temp[a]) & temp[b]=moyb[i-3] & Bg[*,*,i]=temp  & ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
!P.MULTI=[0,2,2]
device,file='yo.ps',bits=24,xoffset=-0.5,yoffset=1,xsize=22.5,ysize=19,/color
LOADCT,4
tvim,phig[*,*,3],/scale,tit='!17E2 '+STRTRIM(STRING(moy[0]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,3]),MAX(phig[*,*,3])]
tvim,phig[*,*,4],/scale,tit='!17E3 '+STRTRIM(STRING(moy[1]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,4]),MAX(phig[*,*,4])]
tvim,phig[*,*,5],/scale,tit='!17E4 '+STRTRIM(STRING(moy[2]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,5]),MAX(phig[*,*,5])]
tvim,phig[*,*,6],/scale,tit='!17E5 '+STRTRIM(STRING(moy[3]),1),xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Relative Phase (in degrees)',range=[MIN(phig[*,*,6]),MAX(phig[*,*,6])]
tvim,Bg[*,*,3]  ,/scale,tit='!17E2',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,3]),MAX(Bg[*,*,3])]
tvim,Bg[*,*,4]  ,/scale,tit='!17E3',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,4]),MAX(Bg[*,*,4])]
tvim,Bg[*,*,5]  ,/scale,tit='!17E4',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,5]),MAX(Bg[*,*,5])]
tvim,Bg[*,*,6]  ,/scale,tit='!17E5',xtit='!17pixels',ytit='!17pixels',barwidth=0.5,stit='!17Contrast',range=[MIN(Bg[*,*,6]),MAX(Bg[*,*,6])]
tvim,I0g,/scale,tit='!17Intensity',barwidth=0.5,range=[MIN(I0g),MAX(I0g)]
DEVICE,/close
PRINT,'PHASE AVERAGES'
FOR i=0,3 DO PRINT,moy[i]
PRINT,'CONTRAST AVERAGES'
FOR i=0,3 DO PRINT,moyb[i]

; RECONSTRUCTION OF THE SEQUENCE

Intenrecons = FLTARR(nx,ny,nseq*nl0)
a=WHERE(distance GT anglim,COMPLEMENT=b)
distance[a] = 0.0
distance[b] = 1.0
FOR i=0,nseq*nl0-1 DO Inten[*,*,i] = Inten[*,*,i]*distance

FOR jjj=0,ny-1 DO BEGIN
    PRINT,jjj
    FOR iii=0,nx-1 DO BEGIN

        IF(distance[iii,jjj] NE 0.0) THEN BEGIN
            
            FOR j=0,nseq-1 DO BEGIN
                FOR ii=0,nl0-1 DO BEGIN
                    
                    profileg  = INTERPOL(blocker0,lam,lam0[ii])
                    FOR i=0,2 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.+tuningJS[j,i]))/2.d0
                    profilegf = profileg
                    FOR i=3,6 DO profileg = profileg * (1.d0+Bg[iii,jjj,i]*COS(dpi/FSR[i]*lam0[ii]+Phig[iii,jjj,i]*!pi/180.))/2.d0
                    
                    Intenrecons[iii,jjj,ii*nseq+j]  = profileg*I0g[iii,jjj]
                    
                ENDFOR
            ENDFOR
            
        ENDIF
        
    ENDFOR
ENDFOR

!P.MULTI=0
; WE PLOT THE SPATIALY AVERAGED RECONSTRUCTION ERROR
device,file='yo3.ps',xsize=35,ysize=25,xoffset=0,yoffset=0,/color
LOADCT,3
PLOT,REBIN(Inten,1,1,nseq*nl0),xst=1,tit='!17',xtit='position number',ytit='intensity',charsize=1.5,yst=1,thick=3;,/ylog
;OPLOT,REBIN(Inten,1,1,nseq*nl0),psym=2
OPLOT,REBIN(Intenrecons,1,1,nseq*nl0),color=150,linestyle=2,thick=2
device,/close
PRINT,'RESIDUAL=',TOTAL( (REFORM(REBIN(Inten,1,1,nseq*nl0))-REFORM(REBIN(Intenrecons,1,1,nseq*nl0)))^2.0 )


a           = WHERE(Inten EQ 0.0)
IF (a[0] NE -1) THEN Inten[a] = -1.d0
DEVICE,file='yo2.ps',xoffset=0,yoffset=0,xsize=20,ysize=26,bits=24,/color
LOADCT,3
!P.multi=[0,2,5]
FOR i=0,nseq*nl0-1 DO BEGIN & temp=REFORM((inten[*,*,i]-intenrecons[*,*,i])/inten[*,*,i])*distance & a=WHERE(FINITE(temp) EQ 0) & IF(a[0] NE -1) THEN temp[a]=0.0 & TVIM,temp,/scale,barwidth=0.5,range=[-0.1,0.1] & ENDFOR
DEVICE,/CLOSE




READ,pause

END

