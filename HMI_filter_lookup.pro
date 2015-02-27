
; idl code that calculate the filter transmission profiles following what is done in vfisv
; MODIFIED TO ACCOUNT FOR CALIBRATION 11
; VERSION OF APRIL 2011: ALSO CALCULATES LOOK-UP TABLE
; AND ESTIMATE PHOTON NOISE MAP AT 3 VELOCITIES
;---------------------------------------------------------------------------------------


PRO HMI_filter_lookup,row,column

; PARAMETERS
 
calibration11 = 1 ; use the calibration 11 or not (yes=1, no=0)?
row=LONG(row)
column=LONG(column)
n=column+row*4096l
Num_lambda=60000l
Lambda_Min=-15000.d0
Delta_Lambda=.5d0
Num_lambda_filter=6

FSR=DBLARR(7)
IF(calibration11 EQ 1) THEN BEGIN
    FSR[0]=0.1689d0  ; //FSR in Angstroms, NB Michelson
    FSR[1]=0.33685d0 ; //WB Michelson
    FSR[2]=0.695d0   ; //Lyot element E1
    FSR[3]=1.417d0   ; //E2
    FSR[4]=2.779d0   ; //E3
    FSR[5]=5.682d0   ; //E4
    FSR[6]=11.354d0  ; //E5
ENDIF ELSE BEGIN
    FSR[0]=0.169d0   ; //FSR in Angstroms, NB Michelson
    FSR[1]=0.337d0   ; //WB Michelson
    FSR[2]=0.695d0   ; //Lyot element E1
    FSR[3]=1.407d0   ; //E2
    FSR[4]=2.779d0   ; //E3
    FSR[5]=5.682d0   ; //E4
    FSR[6]=11.354d0  ; //E5
ENDELSE


nx2=256l
ny2=256l
nelemPHASENT  =4l*nx2*ny2
nelemCONTRASTT=3l*nx2*nx2

filephasesnontunable="/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/non_tunable_phases_710660_June09_cal_256_2.bin";  //binary file containing the phases of the 4 non-tunable elements, format A[element][row][column],AVERAGE FRONT + SIDE CAMERA, 1020", CALMODE, 256x256, with the phases of the tunable elements obtained from an average of the laser detunes
filecontrastsnontunable="/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/non_tunable_contrasts_710660_June09_cal_256_2.bin"; //name of the binary file containing the contrasts of the 4 non-tunable elements A[element][row][column], AVERAGE FRONT+ SIDE CAMERA, 1020", CALMODE, 256x256
filecontraststunable="/home/couvidat/cvs/JSOC/proj/lev1.5_hmi/apps/tunable_contrasts_710660_June09_cal_256.bin";  //binary file containing the contrasts of the 3 tunable elements A[element][row][column], AVERAGE FRONT + SIDE CAMERA, 1020", CALMODE, 256x256
;for these 3 files, there must be values up to a radius = solarradiusmax, and values=zero at larger radii
  
phaseNT=FLTARR(nelemPHASENT)
contrastNT=phaseNT
contrastT=FLTARR(nelemCONTRASTT)

OPENR,1,filephasesnontunable
READU,1,phaseNT
CLOSE,1
OPENR,1,filecontrastsnontunable
READU,1,contrastNT
CLOSE,1
OPENR,1,filecontraststunable
READU,1,contrastT
CLOSE,1

phaseNT=phaseNT*!dpi/180.d0

nblocker=201
IF(calibration11 EQ 1) THEN centerblocker=2.7d0 ELSE centerblocker=2.6d0
nfront=401

 wavelengthbd=[6150.00,6150.20,6150.40,6150.60,6150.80,6151.00,6151.20,6151.40,6151.60,6151.80,6152.00,6152.20,6152.40,6152.60,6152.80,6153.00,6153.20,6153.40,6153.60,6153.80,6154.00,6154.20,6154.40,6154.60,6154.80,6155.00,6155.20,6155.40,6155.60,6155.80,6156.00,6156.20,6156.40,6156.60,6156.80,6157.00,6157.20,6157.40,6157.60,6157.80,6158.00,6158.20,6158.40,6158.60,6158.80,6159.00,6159.20,6159.40,6159.60,6159.80,6160.00,6160.20,6160.40,6160.60,6160.80,6161.00,6161.20,6161.40,6161.60,6161.80,6162.00,6162.20,6162.40,6162.60,6162.80,6163.00,6163.20,6163.40,6163.60,6163.80,6164.00,6164.20,6164.40,6164.60,6164.80,6165.00,6165.20,6165.40,6165.60,6165.80,6166.00,6166.20,6166.40,6166.60,6166.80,6167.00,6167.20,6167.40,6167.60,6167.80,6168.00,6168.20,6168.40,6168.60,6168.80,6169.00,6169.20,6169.40,6169.60,6169.80,6170.00,6170.20,6170.40,6170.60,6170.80,6171.00,6171.20,6171.40,6171.60,6171.80,6172.00,6172.20,6172.40,6172.60,6172.80,6173.00,6173.20,6173.40,6173.60,6173.80,6174.00,6174.20,6174.40,6174.60,6174.80,6175.00,6175.20,6175.40,6175.60,6175.80,6176.00,6176.20,6176.40,6176.60,6176.80,6177.00,6177.20,6177.40,6177.60,6177.80,6178.00,6178.20,6178.40,6178.60,6178.80,6179.00,6179.20,6179.40,6179.60,6179.80,6180.00,6180.20,6180.40,6180.60,6180.80,6181.00,6181.20,6181.40,6181.60,6181.80,6182.00,6182.20,6182.40,6182.60,6182.80,6183.00,6183.20,6183.40,6183.60,6183.80,6184.00,6184.20,6184.40,6184.60,6184.80,6185.00,6185.20,6185.40,6185.60,6185.80,6186.00,6186.20,6186.40,6186.60,6186.80,6187.00,6187.20,6187.40,6187.60,6187.80,6188.00,6188.20,6188.40,6188.60,6188.80,6189.00,6189.20,6189.40,6189.60,6189.80,6190.00]
 blockerbd=[0.0701790,0.0723149,0.0747684,0.0713996,0.0782131,0.0758304,0.0789970,0.0762436,0.0806648,0.0828427,0.0801553,0.0830996,0.0882834,0.0885202,0.0869452,0.0877748,0.0974292,0.0942963,0.0968998,0.0961026,0.100459,0.104028,0.102757,0.107549,0.111349,0.120277,0.117723,0.127142,0.125108,0.135901,0.146540,0.148481,0.151049,0.161267,0.173912,0.191953,0.204322,0.227430,0.239466,0.255259,0.272536,0.311694,0.341673,0.356651,0.409127,0.452214,0.535866,0.614547,0.667113,0.740491,0.847670,0.958023,1.05927,1.19029,1.33457,1.51771,1.90178,2.28149,2.49949,2.80167,3.20520,3.85124,4.35895,4.98798,5.73421,6.53362,8.32412,9.85849,10.6749,12.2367,13.5532,16.0578,17.5336,19.9408,21.4035,26.3633,28.4878,31.9405,33.4455,36.0767,38.4715,42.1947,44.3560,46.8881,49.1468,51.3640,54.9618,57.0772,58.2497,59.3955,60.7570,62.1459,62.7333,63.6812,64.1301,64.7157,65.1849,65.6286,65.4660,65.5828,65.4650,65.4184,65.0766,64.7526,64.0896,63.4954,62.1029,60.4464,59.8266,58.1582,57.0750,54.6480,53.1968,51.1148,49.0429,46.5159,42.3256,38.8035,37.2384,34.5975,32.0762,28.5152,26.2661,23.6695,21.7173,17.3033,15.3031,13.0296,11.8265,10.3604,9.23128,7.53426,6.70699,5.79359,4.97535,4.35803,3.27063,2.70147,2.42119,2.11748,1.83017,1.53329,1.34299,1.19612,1.02633,0.915612,0.711036,0.622389,0.575185,0.507517,0.450262,0.401576,0.365934,0.322949,0.286864,0.272241,0.232461,0.207537,0.189114,0.173546,0.161978,0.152099,0.134795,0.123677,0.110288,0.108344,0.0948288,0.0818621,0.0804488,0.0753219,0.0693417,0.0643225,0.0620898,0.0559437,0.0540745,0.0485797,0.0486797,0.0432530,0.0439143,0.0401164,0.0367754,0.0359879,0.0343058,0.0336281,0.0330711,0.0339798,0.0271329,0.0281424,0.0299408,0.0264017,0.0278133,0.0250958,0.0248676,0.0223389,0.0238825,0.0259792,0.0226330,0.0204282,0.0209307,0.0207487,0.0209464]
wavelengthd=[607.300,607.350,607.400,607.450,607.500,607.550,607.600,607.650,607.700,607.750,607.800,607.850,607.900,607.950,608.000,608.050,608.100,608.150,608.200,608.250,608.300,608.350,608.400,608.450,608.500,608.550,608.600,608.650,608.700,608.750,608.800,608.850,608.900,608.950,609.000,609.050,609.100,609.150,609.200,609.250,609.300,609.350,609.400,609.450,609.500,609.550,609.600,609.650,609.700,609.750,609.800,609.850,609.900,609.950,610.000,610.050,610.100,610.150,610.200,610.250,610.300,610.350,610.400,610.450,610.500,610.550,610.600,610.650,610.700,610.750,610.800,610.850,610.900,610.950,611.000,611.050,611.100,611.150,611.200,611.250,611.300,611.350,611.400,611.450,611.500,611.550,611.600,611.650,611.700,611.750,611.800,611.850,611.900,611.950,612.000,612.050,612.100,612.150,612.200,612.250,612.300,612.350,612.400,612.450,612.500,612.550,612.600,612.650,612.700,612.750,612.800,612.850,612.900,612.950,613.000,613.050,613.100,613.150,613.200,613.250,613.300,613.350,613.400,613.450,613.500,613.550,613.600,613.650,613.700,613.750,613.800,613.850,613.900,613.950,614.000,614.050,614.100,614.150,614.200,614.250,614.300,614.350,614.400,614.450,614.500,614.550,614.600,614.650,614.700,614.750,614.800,614.850,614.900,614.950,615.000,615.050,615.100,615.150,615.200,615.250,615.300,615.350,615.400,615.450,615.500,615.550,615.600,615.650,615.700,615.750,615.800,615.850,615.900,615.950,616.000,616.050,616.100,616.150,616.200,616.250,616.300,616.350,616.400,616.450,616.500,616.550,616.600,616.650,616.700,616.750,616.800,616.850,616.900,616.950,617.000,617.050,617.100,617.150,617.200,617.250,617.300,617.350,617.400,617.450,617.500,617.550,617.600,617.650,617.700,617.750,617.800,617.850,617.900,617.950,618.000,618.050,618.100,618.150,618.200,618.250,618.300,618.350,618.400,618.450,618.500,618.550,618.600,618.650,618.700,618.750,618.800,618.850,618.900,618.950,619.000,619.050,619.100,619.150,619.200,619.250,619.300,619.350,619.400,619.450,619.500,619.550,619.600,619.650,619.700,619.750,619.800,619.850,619.900,619.950,620.000,620.050,620.100,620.150,620.200,620.250,620.300,620.350,620.400,620.450,620.500,620.550,620.600,620.650,620.700,620.750,620.800,620.850,620.900,620.950,621.000,621.050,621.100,621.150,621.200,621.250,621.300,621.350,621.400,621.450,621.500,621.550,621.600,621.650,621.700,621.750,621.800,621.850,621.900,621.950,622.000,622.050,622.100,622.150,622.200,622.250,622.300,622.350,622.400,622.450,622.500,622.550,622.600,622.650,622.700,622.750,622.800,622.850,622.900,622.950,623.000,623.050,623.100,623.150,623.200,623.250,623.300,623.350,623.400,623.450,623.500,623.550,623.600,623.650,623.700,623.750,623.800,623.850,623.900,623.950,624.000,624.050,624.100,624.150,624.200,624.250,624.300,624.350,624.400,624.450,624.500,624.550,624.600,624.650,624.700,624.750,624.800,624.850,624.900,624.950,625.000,625.050,625.100,625.150,625.200,625.250,625.300,625.350,625.400,625.450,625.500,625.550,625.600,625.650,625.700,625.750,625.800,625.850,625.900,625.950,626.000,626.050,626.100,626.150,626.200,626.250,626.300,626.350,626.400,626.450,626.500,626.550,626.600,626.650,626.700,626.750,626.800,626.850,626.900,626.950,627.000,627.050,627.100,627.150,627.200,627.250,627.300]
 frontwindowd=[0.0824,0.0939,0.0787,0.1035,0.0712,0.0794,0.1009,0.0841,0.0889,0.0439,0.0884,0.1057,0.1161,0.0858,0.0717,0.1085,0.0939,0.1030,0.1027,0.0805,0.0831,0.1088,0.1005,0.0767,0.0920,0.0755,0.0862,0.1045,0.1125,0.1068,0.0964,0.1106,0.1134,0.0919,0.1147,0.1105,0.1171,0.1166,0.0975,0.1113,0.1094,0.1317,0.1379,0.1335,0.1444,0.1323,0.1103,0.1487,0.1564,0.1551,0.1550,0.1705,0.1672,0.1258,0.1330,0.1428,0.1620,0.1656,0.1843,0.1940,0.1811,0.1908,0.1628,0.1760,0.2049,0.1959,0.2186,0.2341,0.2504,0.2322,0.2292,0.2409,0.2355,0.2490,0.2984,0.2840,0.2854,0.2927,0.2823,0.3015,0.3269,0.3337,0.3787,0.3949,0.3853,0.3833,0.4298,0.4670,0.4692,0.4981,0.5129,0.5428,0.5865,0.6077,0.6290,0.6324,0.7459,0.7480,0.8150,0.8657,0.8920,0.9285,1.0100,1.0616,1.1343,1.2479,1.2852,1.3584,1.4890,1.5699,1.6990,1.8128,1.9565,2.1121,2.2822,2.4738,2.6485,2.9025,3.1073,3.3880,3.6499,3.9427,4.3108,4.6652,4.9989,5.5283,5.9813,6.5033,7.1314,7.7739,8.3984,9.1871,10.0736,10.8755,11.9132,13.0402,13.9753,15.3119,16.5884,17.9287,19.5970,21.4054,22.8304,24.9109,27.0054,28.6829,30.9226,33.2495,35.3898,37.8594,40.3831,42.7200,44.9671,47.5569,49.9864,52.2030,54.7069,57.0780,58.8828,60.8594,62.6620,64.4119,66.2320,68.2963,69.4221,70.8045,72.2034,73.4397,74.4796,75.5385,76.0362,76.8120,77.7408,78.4357,78.7422,79.3253,79.8358,80.3452,80.8023,80.7372,81.0416,81.5797,81.9136,82.0703,82.4872,82.4652,82.5696,82.6470,82.7640,82.8444,82.7465,82.4404,82.6908,82.8926,83.1894,83.3542,83.3569,83.2381,83.3289,83.2983,83.3072,83.4495,83.3260,83.3259,83.2807,83.1566,83.0003,82.8881,82.8700,83.2338,83.6685,83.5310,83.3660,83.4200,83.4753,83.5535,83.4297,83.5024,83.3255,83.3739,83.2229,83.0138,82.7990,83.0881,82.8699,82.5380,82.4950,82.2752,81.6813,81.0162,80.4859,79.8935,79.0734,78.1775,77.1055,75.6790,73.7593,72.1025,70.1430,67.8169,65.3109,62.7376,59.8977,56.6526,53.3823,50.7530,47.4383,43.8974,40.5806,37.9460,34.8128,31.9024,29.1972,26.9859,24.4400,22.1973,20.2288,18.4441,16.6123,14.9844,13.5245,12.4264,11.2253,10.1547,9.2421,8.5064,7.7160,6.9712,6.3549,5.8530,5.3515,4.8613,4.4704,4.1431,3.7857,3.4459,3.1822,2.9684,2.7397,2.5416,2.3371,2.1395,2.0216,1.8965,1.7565,1.6370,1.5414,1.4241,1.3505,1.2548,1.1747,1.1431,1.0586,0.9672,0.9421,0.8876,0.8538,0.8146,0.7764,0.7325,0.6940,0.6494,0.6119,0.6018,0.5367,0.5330,0.5568,0.5291,0.4903,0.4759,0.4601,0.4422,0.3874,0.3670,0.3589,0.3565,0.3655,0.3273,0.3272,0.3035,0.3008,0.3135,0.2618,0.2653,0.2760,0.2466,0.2407,0.2260,0.2107,0.2255,0.2064,0.2066,0.1937,0.1648,0.1614,0.1906,0.1738,0.1403,0.1535,0.1480,0.1581,0.1243,0.1326,0.1140,0.1280,0.1621,0.1404,0.1348,0.1110,0.1075,0.1022,0.1140,0.1186,0.1072,0.1250,0.1132,0.1193,0.0794,0.0808,0.0930,0.0886,0.0693,0.0934,0.0827,0.0644,0.0723,0.0732,0.0574,0.0606,0.0555,0.0536,0.0540,0.0625,0.0444,0.0616,0.0663,0.0506,0.0506,0.0492,0.0294,0.0661,0.0511,0.0556,0.0188,0.0257,0.0414,0.0484,0.0167,0.0341,0.0216,0.0269,0.0308,0.0451,0.0700,0.0326,0.0110,0.0288,0.0414,0.0225,0.0119,0.0509]


; ORIGINAL PHASE MAP (END OF MARCH 2010)
;phaseT=fitsio_read_image("/SUM0/D109538835/D46717947/S00000/phases.fits") ; original phasemap, in 256x256: phaseT[3,256,256]
;HCME1=36
;HCMNB=82
;HCMWB=58
;nel=3.

; PHASE MAPS OF END OF AUGUST 4, 2010, FOR FRONT CAMERA
phaseT=fitsio_read_image("/SUM2/D134135089/S00000/phases.fits") ; phaseT[5,128,128]
phaseT=REBIN(phaseT,5,256,256)
HCME1=36
HCMNB=82
HCMWB=58
nel=5.

; PHASE MAPS OF END OF NOVEMBER 2010
;phaseT=fitsio_read_image("/SUM18/D112714488/S00000/phases.fits") ;phaseT[5,128,128], 2010/11/24
;phaseT=REBIN(phaseT,5,256,256)
;HCME1=36
;HCMNB=82
;HCMWB=58
;nel=5.


phaseT= phaseT*!dpi/180.d0
lam0  = 6173.3433d0

HCME1phase=DBLARR(Num_lambda_filter)
HCMWBphase=DBLARR(Num_lambda_filter)
HCMNBphase=DBLARR(Num_lambda_filter)

IF(Num_lambda_filter EQ 6) THEN BEGIN
    HCME1phase[5]= double( ((HCME1+15)*6 MOD 360))*!dpi/180.d0 ; //I0
    HCME1phase[4]= double( ((HCME1+9 )*6 MOD 360))*!dpi/180.d0 ; //I1
    HCME1phase[3]= double( ((HCME1+3 )*6 MOD 360))*!dpi/180.d0 ; //I2
    HCME1phase[2]= double( ((HCME1-3 )*6 MOD 360))*!dpi/180.d0 ; //I3
    HCME1phase[1]= double( ((HCME1-9 )*6 MOD 360))*!dpi/180.d0 ; //I4
    HCME1phase[0]= double( ((HCME1-15)*6 MOD 360))*!dpi/180.d0 ; //I5
    HCMWBphase[5]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[4]= double( ((HCMWB-18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[3]= double( ((HCMWB-6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[2]= double( ((HCMWB+6 )*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[1]= double( ((HCMWB+18)*6 MOD 360))*!dpi/180.d0 ;
    HCMWBphase[0]= double( ((HCMWB-30)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[5]= double( ((HCMNB-0 )*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[4]= double( ((HCMNB+24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[3]= double( ((HCMNB-12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[2]= double( ((HCMNB+12)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[1]= double( ((HCMNB-24)*6 MOD 360))*!dpi/180.d0 ;
    HCMNBphase[0]= double( ((HCMNB+0 )*6 MOD 360))*!dpi/180.d0 ;
ENDIF

wavelength=DBLARR(Num_lambda)

FOR i=0l,Num_lambda-1 DO  wavelength[i]= Lambda_Min/1000.d0 + double(i)*Delta_Lambda/1000.d0; //wavelength is the wavelength grid in Angstroms

wavelengthd2=DBLARR(nfront)
frontwindowd2=wavelengthd2

FOR i=0,nfront-1 DO BEGIN
    wavelengthd2[i]=wavelengthd[i]*10.d0-lam0 ;
    frontwindowd2[i]=frontwindowd[i]/100.d0 ;
ENDFOR

wavelengthbd2=DBLARR(nblocker)
blockerbd2=wavelengthbd2

FOR i=0,nblocker-1 DO BEGIN
    wavelengthbd2[i]=wavelengthbd[i]+centerblocker-lam0 ;
    blockerbd2[i]=blockerbd[i]/100.d0 ;
ENDFOR

frontwindowint=INTERPOL(frontwindowd2,wavelengthd2,wavelength)
blockerint=INTERPOL(blockerbd2,wavelengthbd2,wavelength)

blockerint=blockerint*frontwindowint

FWHMCoeffs=[151.34559d0,-58.521771d0]
minimumCoeffs=[0.41922611d0,0.24190794d0]
X0=2048
Y0=2048
row=2048
column=2048
rsun_ref=696000000.d0
dsun_obs=147206569539.31d0
cdeltx=0.504273d0
Rsun=double(asin(rsun_ref/dsun_obs)/3.14159265358979d0*180.d0*3600.d0/cdeltx)
distance      = sqrt((double(row)-Y0)*(double(row)-Y0)+(double(column)-X0)*(double(column)-X0)) ; //distance in pixels
distance      = cos(asin(distance/Rsun)) ;                                                        //cosine of angular distance from disk center

lyot=DBLARR(Num_lambda)
filters=DBLARR(Num_lambda_filter,Num_lambda)


; bilinear interpolation
cols=4096l
column        =n MOD cols         ; //column from 0 to 4095
row           =n / cols           ; //row from 0 to 4095
ratio         =cols/nx2           ; //should be 16
PRINT,'column =',column
PRINT,'row =',row

y0 = (row/ratio)
x0 = (column/ratio)

loc1=x0+y0*nx2
contrastNTi=DBLARR(4)
phaseNTi=DBLARR(4)

contrastNTi[0]=contrastNT[loc1+0.*nx2*ny2]
contrastNTi[1]=contrastNT[loc1+1.*nx2*ny2]
contrastNTi[2]=contrastNT[loc1+2.*nx2*ny2]
contrastNTi[3]=contrastNT[loc1+3.*nx2*ny2]
phaseNTi[0]=phaseNT[loc1+0.*nx2*ny2]
phaseNTi[1]=phaseNT[loc1+1.*nx2*ny2]
phaseNTi[2]=phaseNT[loc1+2.*nx2*ny2]
phaseNTi[3]=phaseNT[loc1+3.*nx2*ny2]
contrastNB=contrastT[loc1+0.*nx2*ny2]
contrastWB=contrastT[loc1+1.*nx2*ny2]
contrastE1=contrastT[loc1+2.*nx2*ny2]
phaseNB=phaseT[loc1*nel]
phaseWB=phaseT[loc1*nel+1]
phaseE1=phaseT[loc1*nel+2]

IF(calibration11 EQ 0) THEN BEGIN

FOR j=0l,Num_lambda-1 DO lyot[j]=blockerint[j]*(1.+contrastNTi[0]*cos(2.0*!dpi/FSR[3]*wavelength[j]+phaseNTi[0]))/2.*(1.+contrastNTi[1]*cos(2.0*!dpi/FSR[4]*wavelength[j]+phaseNTi[1]))/2.*(1.+contrastNTi[2]*cos(2.0*!dpi/FSR[5]*wavelength[j]+phaseNTi[2]))/2.d0*(1.+contrastNTi[3]*cos(2.0*!dpi/FSR[6]*wavelength[j]+phaseNTi[3]))/2. ;
FOR i=0l,Num_lambda_filter-1 DO FOR j=0l,Num_lambda-1 DO filters[i,j] = lyot[j]*(1.d0+contrastNB*cos(2.d0*!dpi/FSR[0]*wavelength[j]+HCMNBphase[i]+phaseNB))/2.*(1.+contrastWB*cos(2.d0*!dpi/FSR[1]*wavelength[j]+HCMWBphase[i]+phaseWB))/2.*(1.d0+contrastE1*cos(2.d0*!dpi/FSR[2]*wavelength[j]-HCME1phase[i]+phaseE1))/2.d0 ;

ENDIF ELSE BEGIN

; WITH CALIBRATION 11
;--------------------------------
    
FOR j=0l,Num_lambda-1 DO lyot[j]=blockerint[j]*(1.+contrastNTi[0]*cos(2.0*!dpi/FSR[3]*wavelength[j]+phaseNTi[0]))/2.*(1.+contrastNTi[1]*cos(2.0*!dpi/FSR[4]*wavelength[j]+phaseNTi[1]))/2.*(1.+contrastNTi[2]*cos(2.0*!dpi/FSR[5]*wavelength[j]+phaseNTi[2]+0.4))/2.d0*(1.+contrastNTi[3]*cos(2.0*!dpi/FSR[6]*wavelength[j]+phaseNTi[3]-1.1))/2. ;
FOR i=0l,Num_lambda_filter-1 DO FOR j=0l,Num_lambda-1 DO filters[i,j] = lyot[j]*(1.d0+contrastNB*cos(2.d0*!dpi/FSR[0]*wavelength[j]+HCMNBphase[i]+phaseNB+1.573*!dpi/180.d0))/2.*(1.+contrastWB*cos(2.d0*!dpi/FSR[1]*wavelength[j]+HCMWBphase[i]+phaseWB+0.59*!dpi/180.d0))/2.*(1.d0+contrastE1*cos(2.d0*!dpi/FSR[2]*wavelength[j]-HCME1phase[i]+phaseE1-3.27*!dpi/180.d0))/2.d0 ;

ENDELSE

; NORMALIZATION

filters = filters/MAX(filters)
;SAVE,wavelength,filters,file='filters.bin'


;-----------------------------------------------------------------------------
;
; CALCULATE THE LINE PROFILE
;
;-----------------------------------------------------------------------------



;nlam=10000
w=wavelength;DINDGEN(nlam)/(nlam-1.d0)*10.d0-5.d0
nlam=Num_lambda

; values with calibration 11
Icg=1.0d0;
fdepthg=0.5625d0;//0.63;//0.56;//0.642710;//0.651742;
fwidthg=0.06415d0;//0.069334;//0.071124;//0.0687;//0.069334;//0.069598;
a=0.03d0; //damping
l=w/fwidthg
line=DINDGEN(nlam)

aa=WHERE(ABS(l) LE 26.5,COMPLEMENT=bb)

line[aa] = Icg-fdepthg*exp(-l[aa]*l[aa])*(1.d0-a*2.d0/sqrt(!dpi)*(1.d0/2.d0/l[aa]/l[aa])*((4.d0*l[aa]*l[aa]+3.d0)*(l[aa]*l[aa]+1.d0)*exp(-l[aa]*l[aa])-1.d0/l[aa]/l[aa]*(2.d0*l[aa]*l[aa]+3.d0)*sinh(l[aa]*l[aa])))-0.015d0*exp(-(w[aa]+0.225d0)*(w[aa]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[aa]-0.15d0)*(w[aa]-0.15d0)/0.22d0/0.22d0)
line[bb] = Icg-0.015d0*exp(-(w[bb]+0.225d0)*(w[bb]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[bb]-0.15d0)*(w[bb]-0.15d0)/0.22d0/0.22d0)
; there is a NaN at wavelength=0
aa=WHERE(FINITE(line) EQ 0,COMPLEMENT=bb)
line[aa[0]]=INTERPOL(line[bb],wavelength[bb],wavelength[aa[0]])

READ,pause

;-----------------------------------------------------------------------------
;
; MULTIPLY THE LINE SHIFTED BY VELOCITY WITH THE PROFILES
;
;-----------------------------------------------------------------------------

dlamdv  = 2.059205672212074294d-5
ntest   = 301
dlam    = (wavelength[1]-wavelength[0])
dvtest  = dlam/dlamdv ; m/s
vtest   = (DINDGEN(ntest)-DOUBLE(ntest-1)/2.d0)*dvtest
L0      = FLTARR(ntest)
L1      = L0
L2      = L0
L3      = L0
L4      = L0
L5      = L0

maxshiftlam = ROUND(vtest[ntest-1]/dvtest)

FOR  i=0,ntest-1 DO BEGIN
    
    PRINT,i
    shiftlam=ROUND(vtest[i]/dvtest) ;//SIGN CONVENTION: POSITIVE VELOCITIES CORRESPOND TO REDSHIFT
    
    L0[i] = 0.d0
    L1[i] = 0.d0
    L2[i] = 0.d0
    L3[i] = 0.d0
    L4[i] = 0.d0
    L5[i] = 0.d0
    
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L0[i] = L0[i]+filters[5,k]*line[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L1[i] = L1[i]+filters[4,k]*line[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L2[i] = L2[i]+filters[3,k]*line[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L3[i] = L3[i]+filters[2,k]*line[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L4[i] = L4[i]+filters[1,k]*line[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L5[i] = L5[i]+filters[0,k]*line[k-shiftlam]
    
ENDFOR
maxi=MAX(L0)
L0=L0/maxi
L1=L1/maxi
L2=L2/maxi
L3=L3/maxi
L4=L4/maxi
L5=L5/maxi

RESTORE,'lineatdiskcenter.bin' ; from HMI_line.pro

lam   =-vtest*dlamdv
plot ,lam+0.172,L0,xrange=[-0.3,0.3],xst=1
oplot,lam+0.1032,L1
oplot,lam+0.0344,L2
oplot,lam-0.0344,L3
oplot,lam-0.1032,L4
oplot,lam-0.172,L5
;oplot,-vel*dlamdv+0.1720,L0b/MAX(L0b)+0.01,col=180
;oplot,-vel*dlamdv+0.1032,L1b/MAX(L0b)+0.01,col=180
;oplot,-vel*dlamdv+0.0344,L2b/MAX(L0b)+0.01,col=180
;oplot,-vel*dlamdv-0.0344,L3b/MAX(L0b)+0.01,col=180
;oplot,-vel*dlamdv-0.1032,L4b/MAX(L0b)+0.01,col=180
;oplot,-vel*dlamdv-0.1720,L5b/MAX(L0b)+0.01,col=180

;GOTO,noimprov
;-------------------------------------------------------------
; ATTEMPT AT IMPROVING THE LINE IN OBSMODE:
;-------------------------------------------------------------

IF(calibration11 EQ 0) THEN BEGIN
; NB: values WITHOUT calibration 11
    Icg=1.0d0                   ;
    fdepthg=0.65                ;0.785d0;0.5625
    fwidthg=0.058d0             ;0.06415
    a=0.03d0                    ; //damping
ENDIF ELSE BEGIN
; NB: values WITH calibration 11
    Icg=1.0d0                   ;
    fdepthg=0.7                 ;0.785d0;0.5625
    fwidthg=0.058d0             ;0.06415
    a=0.03d0                    ; //damping
    l=w/fwidthg
    line2=DINDGEN(nlam)
ENDELSE

aa=WHERE(ABS(l) LE 26.5,COMPLEMENT=bb)

line2[aa] = Icg-fdepthg*exp(-l[aa]*l[aa])*(1.d0-a*2.d0/sqrt(!dpi)*(1.d0/2.d0/l[aa]/l[aa])*((4.d0*l[aa]*l[aa]+3.d0)*(l[aa]*l[aa]+1.d0)*exp(-l[aa]*l[aa])-1.d0/l[aa]/l[aa]*(2.d0*l[aa]*l[aa]+3.d0)*sinh(l[aa]*l[aa])))-0.015d0*exp(-(w[aa]+0.225d0)*(w[aa]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[aa]-0.15d0)*(w[aa]-0.15d0)/0.22d0/0.22d0)
line2[bb] = Icg-0.015d0*exp(-(w[bb]+0.225d0)*(w[bb]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[bb]-0.15d0)*(w[bb]-0.15d0)/0.22d0/0.22d0)
; there is a NaN at wavelength=0
aa=WHERE(FINITE(line2) EQ 0,COMPLEMENT=bb)
line2[aa[0]]=INTERPOL(line2[bb],wavelength[bb],wavelength[aa[0]])
L0c      = FLTARR(ntest)
L1c      = L0c
L2c      = L0c
L3c      = L0c
L4c      = L0c
L5c      = L0c

FOR  i=0,ntest-1 DO BEGIN
    
    PRINT,i
    shiftlam=ROUND(vtest[i]/dvtest) ;//SIGN CONVENTION: POSITIVE VELOCITIES CORRESPOND TO REDSHIFT
    
    L0c[i] = 0.d0
    L1c[i] = 0.d0
    L2c[i] = 0.d0
    L3c[i] = 0.d0
    L4c[i] = 0.d0
    L5c[i] = 0.d0
    
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L0c[i] = L0c[i]+filters[5,k]*line2[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L1c[i] = L1c[i]+filters[4,k]*line2[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L2c[i] = L2c[i]+filters[3,k]*line2[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L3c[i] = L3c[i]+filters[2,k]*line2[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L4c[i] = L4c[i]+filters[1,k]*line2[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L5c[i] = L5c[i]+filters[0,k]*line2[k-shiftlam]
    
ENDFOR
maxi=MAX(L0c)
L0c=L0c/maxi
L1c=L1c/maxi
L2c=L2c/maxi
L3c=L3c/maxi
L4c=L4c/maxi
L5c=L5c/maxi
oplot,lam+0.172,L0c,col=180,thick=2
oplot,lam+0.1032,L1c,col=180,thick=2
oplot,lam+0.0344,L2c,col=180,thick=2
oplot,lam-0.0344,L3c,col=180,thick=2
oplot,lam-0.1032,L4c,col=180,thick=2
oplot,lam-0.172,L5c,col=180,thick=2

line=[L0c,L1c,L2c,L3c,L4c,L5c]
wave=lam
SAVE,line,wave,file='improvedline.bin'
READ,pause

noimprov:
;--------------------------------------------------------------------------------------------
; LINE AND FILTERS AS DERIVED FROM PHASEMAPS_TEST_VOIGT:
;--------------------------------------------------------------------------------------------

; WITH CALIBRATION 11
;--------------------------------
    

;TUNABLE TRANSMISSION PROFILE
    FOR i=0l,Num_lambda_filter-1 DO FOR j=0l,Num_lambda-1 DO filters[i,j] = lyot[j]*(1.d0+contrastNB*cos(2.d0*!dpi/FSR[0]*wavelength[j]+HCMNBphase[i]+phaseNB))/2.*(1.+contrastWB*cos(2.d0*!dpi/FSR[1]*wavelength[j]+HCMWBphase[i]+phaseWB))/2.*(1.d0+contrastE1*cos(2.d0*!dpi/FSR[2]*wavelength[j]-HCME1phase[i]+phaseE1))/2.d0 ;
filters = filters/MAX(filters)

; values from detunes of May 27, 2010
Icg=1.0d0;
fdepthg=0.679164;
fwidthg=0.0676577
a=0.03d0; //damping
l=w/fwidthg
line3=DINDGEN(nlam)

aa=WHERE(ABS(l) LE 26.5,COMPLEMENT=bb)

line3[aa] = Icg-fdepthg*exp(-l[aa]*l[aa])*(1.d0-a*2.d0/sqrt(!dpi)*(1.d0/2.d0/l[aa]/l[aa])*((4.d0*l[aa]*l[aa]+3.d0)*(l[aa]*l[aa]+1.d0)*exp(-l[aa]*l[aa])-1.d0/l[aa]/l[aa]*(2.d0*l[aa]*l[aa]+3.d0)*sinh(l[aa]*l[aa])))-0.015d0*exp(-(w[aa]+0.225d0)*(w[aa]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[aa]-0.15d0)*(w[aa]-0.15d0)/0.22d0/0.22d0)
line3[bb] = Icg-0.015d0*exp(-(w[bb]+0.225d0)*(w[bb]+0.225d0)/0.2d0/0.2d0)+0.004d0*exp(-(w[bb]-0.15d0)*(w[bb]-0.15d0)/0.22d0/0.22d0)
; there is a NaN at wavelength=0
aa=WHERE(FINITE(line3) EQ 0,COMPLEMENT=bb)
line3[aa[0]]=INTERPOL(line3[bb],wavelength[bb],wavelength[aa[0]])
L0d      = FLTARR(ntest)
L1d      = L0d
L2d      = L0d
L3d      = L0d
L4d      = L0d
L5d      = L0d

FOR  i=0,ntest-1 DO BEGIN
    
    PRINT,i
    shiftlam=ROUND(vtest[i]/dvtest) ;//SIGN CONVENTION: POSITIVE VELOCITIES CORRESPOND TO REDSHIFT
    
    L0d[i] = 0.d0
    L1d[i] = 0.d0
    L2d[i] = 0.d0
    L3d[i] = 0.d0
    L4d[i] = 0.d0
    L5d[i] = 0.d0
    
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L0d[i] = L0d[i]+filters[5,k]*line3[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L1d[i] = L1d[i]+filters[4,k]*line3[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L2d[i] = L2d[i]+filters[3,k]*line3[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L3d[i] = L3d[i]+filters[2,k]*line3[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L4d[i] = L4d[i]+filters[1,k]*line3[k-shiftlam]
    FOR k=maxshiftlam,nlam-maxshiftlam-1  DO L5d[i] = L5d[i]+filters[0,k]*line3[k-shiftlam]
    
ENDFOR
maxi=MAX(L0d)
L0d=L0d/maxi
L1d=L1d/maxi
L2d=L2d/maxi
L3d=L3d/maxi
L4d=L4d/maxi
L5d=L5d/maxi
oplot,lam+0.172 ,L0d,col=18000,thick=2
oplot,lam+0.1032,L1d,col=18000,thick=2
oplot,lam+0.0344,L2d,col=18000,thick=2
oplot,lam-0.0344,L3d,col=18000,thick=2
oplot,lam-0.1032,L4d,col=18000,thick=2
oplot,lam-0.172,L5d,col=18000,thick=2


READ,pause

END
