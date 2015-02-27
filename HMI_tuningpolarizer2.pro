; this program checks if the polarizer is working fine (for the CPT)
; by using a detune sequence + 4 polarizer positions and fitting for
; the phases and contrasts of E1, NB, and WB

PRO PROFILES,TRANSMISSION,WAVELENGTH,Q

TRANSMISSION = [      0.079815385,     0.072207692,     0.069323077,     0.058315385,     0.055630769,     0.052015385, $
     0.048530769,     0.047576923,     0.047238462,     0.046415385,     0.047284615,     0.046800000, $
     0.046030769,     0.046446154,     0.045676923,     0.045692308,     0.048015385,     0.050761538, $
     0.051869231,     0.051561538,     0.052215385,     0.051884615,     0.051138462,     0.049823077, $
     0.049876923,     0.049892308,     0.049261538,     0.049953846,     0.051700000,     0.052007692, $
     0.054230769,     0.055784615,     0.056138462,     0.056161538,     0.057823077,     0.059738462, $
     0.060784615,     0.058953846,     0.058761538,     0.060084615,     0.060330769,     0.060715385, $
     0.061800000,     0.063361538,     0.064707692,     0.064823077,     0.066576923,     0.069923077, $
     0.071753846,     0.072553846,     0.073284615,     0.073515385,     0.074753846,     0.075430769, $
     0.076546154,     0.078461538,     0.079207692,     0.080115385,     0.080238462,     0.080169231, $
     0.082061538,     0.082523077,     0.083884615,     0.086030769,     0.088830769,     0.091607692, $
     0.094200000,     0.095576923,     0.098538462,      0.10114615,      0.10226923,      0.10220000, $
      0.10633846,      0.10782308,      0.10943846,      0.11089231,      0.11333077,      0.11494615, $
      0.11745385,      0.11990769,      0.12781538,      0.13343846,      0.13826923,      0.14320769, $
      0.14522308,      0.14779231,      0.15229231,      0.15508462,      0.15857692,      0.16247692, $
      0.16628462,      0.17056923,      0.17987692,      0.19141538,      0.20266154,      0.21387692, $
      0.22080769,      0.22725385,      0.23677692,      0.24432308,      0.25499231,      0.26380769, $
      0.27519231,      0.29003077,      0.30076923,      0.31460769,      0.33273077,      0.34810000, $
      0.36460769,      0.38185385,      0.40283846,      0.42906154,      0.45216154,      0.48198462, $
      0.51460000,      0.54716154,      0.58129231,      0.61933077,      0.66056154,      0.70758462, $
      0.75252308,      0.80444615,      0.86027692,      0.92153846,      0.98621538,       1.0681692, $
       1.1572538,       1.2520231,       1.3602385,       1.4803846,       1.6050615,       1.7420231, $
       1.8848846,       2.0489846,       2.2282846,       2.4297000,       2.6674231,       2.9247769, $
       3.2160615,       3.5376923,       3.8981923,       4.3001231,       4.7465154,       5.2386692, $
       5.7897769,       6.4001462,       7.0971308,       7.8794462,       8.7954000,       9.8169846, $
       10.965323,       12.223577,       13.611100,       15.156415,       16.851254,       18.713315, $
       20.850185,       23.173192,       25.735600,       28.489569,       31.444992,       34.619385, $
       37.940538,       41.397585,       44.961131,       48.587208,       52.290969,       55.923031, $
       59.496015,       62.987469,       66.334431,       69.272046,       71.916238,       74.326546, $
       76.373531,       78.146285,       79.729638,       81.039100,       82.171877,       83.102231, $
       83.936746,       84.586523,       85.097815,       85.637400,       85.994777,       86.183654, $
       86.378377,       86.542108,       86.674346,       86.885777,       87.121323,       87.259969, $
       87.377431,       87.450431,       87.363938,       87.331554,       87.336692,       87.446862, $
       87.515754,       87.553754,       87.631623,       87.563023,       87.491292,       87.471685, $
       87.347977,       87.227700,       87.150177,       86.975031,       86.789754,       86.696085, $
       86.620046,       86.460154,       86.486446,       86.416885,       86.273031,       86.163892, $
       86.029785,       85.936846,       85.790477,       85.673715,       85.692008,       85.638908, $
       85.724615,       85.797985,       85.718015,       85.823392,       85.765908,       85.706177, $
       85.672400,       85.717954,       85.763508,       85.547077,       85.388023,       85.227246, $
       84.836192,       84.541238,       84.288969,       84.074738,       83.662138,       83.064223, $
       82.405046,       81.488623,       80.506915,       79.370508,       78.058131,       76.822562, $
       75.335377,       73.649415,       71.917085,       70.045354,       68.107669,       65.930877, $
       63.701931,       61.504354,       59.133338,       56.696838,       54.297985,       51.768377, $
       49.225531,       46.659469,       44.169608,       41.672885,       39.254346,       36.920238, $
       34.641931,       32.347315,       30.191492,       28.105292,       26.082338,       24.180915, $
       22.410685,       20.751923,       19.179546,       17.722208,       16.391938,       15.121777, $
       13.950746,       12.835738,       11.799569,       10.844308,       9.9634769,       9.1653231, $
       8.4242846,       7.7351385,       7.1281308,       6.5400231,       5.9925923,       5.5200692, $
       5.0875000,       4.6810923,       4.3013231,       3.9545769,       3.6453154,       3.3531846, $
       3.0907692,       2.8614000,       2.6480154,       2.4472615,       2.2679615,       2.0975462, $
       1.9435615,       1.7984000,       1.6625077,       1.5419846,       1.4286538,       1.3257000, $
       1.2343615,       1.1485923,       1.0701077,      0.99720769,      0.92714615,      0.86455385, $
      0.80332308,      0.75009231,      0.70207692,      0.65858462,      0.61944615,      0.58358462, $
      0.54936923,      0.51558462,      0.48363846,      0.45746923,      0.43423077,      0.40927692, $
      0.38575385,      0.36711538,      0.34720000,      0.32970000,      0.31573077,      0.30209231, $
      0.28845385,      0.27267692,      0.26069231,      0.25087692,      0.23942308,      0.23076923, $
      0.22176923,      0.21216154,      0.20382308,      0.19621538,      0.18970000,      0.18182308, $
      0.17176154,      0.16465385,      0.15631538,      0.14856154,      0.14359231,      0.13928462, $
      0.13404615,      0.13016154,      0.12525385,      0.12276923,      0.12067692,      0.11860769, $
      0.11747692,      0.11587692,      0.11295385,      0.10832308,      0.10488462,      0.10130769, $
     0.099369231,     0.098115385,     0.097300000,     0.097053846,     0.095630769,     0.091553846, $
     0.088576923,     0.085561538,     0.083146154,     0.081700000,     0.079984615,     0.078653846, $
     0.076576923,     0.073969231,     0.070961538,     0.070415385,     0.069584615,     0.069323077, $
     0.068461538,     0.069861538,     0.069161538,     0.068446154,     0.067376923,     0.067069231, $
     0.065907692,     0.064115385,     0.062453846,     0.063861538,     0.061915385,     0.061361538, $
     0.060107692,     0.058769231,     0.057830769,     0.056038462,     0.056161538,     0.057076923, $
     0.056953846,     0.059800000,     0.055792308,     0.055792308,     0.055792308 ]


WAVELENGTH = [          607.30000,       607.35000,       607.40000,       607.45000,       607.50000,       607.55000,  $
       607.60000,       607.65000,       607.70000,       607.75000,       607.80000,       607.85000, $
       607.90000,       607.95000,       608.00000,       608.05000,       608.10000,       608.15000, $
       608.20000,       608.25000,       608.30000,       608.35000,       608.40000,       608.45000, $
       608.50000,       608.55000,       608.60000,       608.65000,       608.70000,       608.75000, $
       608.80000,       608.85000,       608.90000,       608.95000,       609.00000,       609.05000, $
       609.10000,       609.15000,       609.20000,       609.25000,       609.30000,       609.35000, $
       609.40000,       609.45000,       609.50000,       609.55000,       609.60000,       609.65000, $
       609.70000,       609.75000,       609.80000,       609.85000,       609.90000,       609.95000, $
       610.00000,       610.05000,       610.10000,       610.15000,       610.20000,       610.25000, $
       610.30000,       610.35000,       610.40000,       610.45000,       610.50000,       610.55000, $
       610.60000,       610.65000,       610.70000,       610.75000,       610.80000,       610.85000, $
       610.90000,       610.95000,       611.00000,       611.05000,       611.10000,       611.15000, $
       611.20000,       611.25000,       611.30000,       611.35000,       611.40000,       611.45000, $
       611.50000,       611.55000,       611.60000,       611.65000,       611.70000,       611.75000, $
       611.80000,       611.85000,       611.90000,       611.95000,       612.00000,       612.05000, $
       612.10000,       612.15000,       612.20000,       612.25000,       612.30000,       612.35000, $
       612.40000,       612.45000,       612.50000,       612.55000,       612.60000,       612.65000, $
       612.70000,       612.75000,       612.80000,       612.85000,       612.90000,       612.95000, $
       613.00000,       613.05000,       613.10000,       613.15000,       613.20000,       613.25000, $
       613.30000,       613.35000,       613.40000,       613.45000,       613.50000,       613.55000, $
       613.60000,       613.65000,       613.70000,       613.75000,       613.80000,       613.85000, $
       613.90000,       613.95000,       614.00000,       614.05000,       614.10000,       614.15000, $
       614.20000,       614.25000,       614.30000,       614.35000,       614.40000,       614.45000, $
       614.50000,       614.55000,       614.60000,       614.65000,       614.70000,       614.75000, $
       614.80000,       614.85000,       614.90000,       614.95000,       615.00000,       615.05000, $
       615.10000,       615.15000,       615.20000,       615.25000,       615.30000,       615.35000, $
       615.40000,       615.45000,       615.50000,       615.55000,       615.60000,       615.65000, $
       615.70000,       615.75000,       615.80000,       615.85000,       615.90000,       615.95000, $
       616.00000,       616.05000,       616.10000,       616.15000,       616.20000,       616.25000, $
       616.30000,       616.35000,       616.40000,       616.45000,       616.50000,       616.55000, $
       616.60000,       616.65000,       616.70000,       616.75000,       616.80000,       616.85000, $
       616.90000,       616.95000,       617.00000,       617.05000,       617.10000,       617.15000, $
       617.20000,       617.25000,       617.30000,       617.35000,       617.40000,       617.45000, $
       617.50000,       617.55000,       617.60000,       617.65000,       617.70000,       617.75000, $
       617.80000,       617.85000,       617.90000,       617.95000,       618.00000,       618.05000, $
       618.10000,       618.15000,       618.20000,       618.25000,       618.30000,       618.35000, $
       618.40000,       618.45000,       618.50000,       618.55000,       618.60000,       618.65000, $
       618.70000,       618.75000,       618.80000,       618.85000,       618.90000,       618.95000, $
       619.00000,       619.05000,       619.10000,       619.15000,       619.20000,       619.25000, $
       619.30000,       619.35000,       619.40000,       619.45000,       619.50000,       619.55000, $
       619.60000,       619.65000,       619.70000,       619.75000,       619.80000,       619.85000, $
       619.90000,       619.95000,       620.00000,       620.05000,       620.10000,       620.15000, $
       620.20000,       620.25000,       620.30000,       620.35000,       620.40000,       620.45000, $
       620.50000,       620.55000,       620.60000,       620.65000,       620.70000,       620.75000, $
       620.80000,       620.85000,       620.90000,       620.95000,       621.00000,       621.05000, $
       621.10000,       621.15000,       621.20000,       621.25000,       621.30000,       621.35000, $
       621.40000,       621.45000,       621.50000,       621.55000,       621.60000,       621.65000, $
       621.70000,       621.75000,       621.80000,       621.85000,       621.90000,       621.95000, $
       622.00000,       622.05000,       622.10000,       622.15000,       622.20000,       622.25000, $
       622.30000,       622.35000,       622.40000,       622.45000,       622.50000,       622.55000, $
       622.60000,       622.65000,       622.70000,       622.75000,       622.80000,       622.85000, $
       622.90000,       622.95000,       623.00000,       623.05000,       623.10000,       623.15000, $
       623.20000,       623.25000,       623.30000,       623.35000,       623.40000,       623.45000, $
       623.50000,       623.55000,       623.60000,       623.65000,       623.70000,       623.75000, $
       623.80000,       623.85000,       623.90000,       623.95000,       624.00000,       624.05000, $
       624.10000,       624.15000,       624.20000,       624.25000,       624.30000,       624.35000, $
       624.40000,       624.45000,       624.50000,       624.55000,       624.60000,       624.65000, $
       624.70000,       624.75000,       624.80000,       624.85000,       624.90000,       624.95000, $
       625.00000,       625.05000,       625.10000,       625.15000,       625.20000,       625.25000, $
       625.30000,       625.35000,       625.40000,       625.45000,       625.50000,       625.55000, $
       625.60000,       625.65000,       625.70000,       625.75000,       625.80000,       625.85000, $
       625.90000,       625.95000,       626.00000,       626.05000,       626.10000,       626.15000, $
       626.20000,       626.25000,       626.30000,       626.35000,       626.40000,       626.45000, $
       626.50000,       626.55000,       626.60000,       626.65000,       626.70000,       626.75000, $
       626.80000,       626.85000,       626.90000,       626.95000,       627.00000,       627.05000, $
       627.10000,       627.15000,       627.20000,       627.25000,       627.30000 ]

Q = [[6190.00,      6189.80,      6189.60,      6189.40,      6189.20,      6189.00,      6188.80,      6188.60, $
      6188.40,      6188.20,      6188.00,      6187.80,      6187.60,      6187.40,      6187.20,      6187.00, $
      6186.80,      6186.60,      6186.40,      6186.20,      6186.00,      6185.80,      6185.60,      6185.40, $
      6185.20,      6185.00,      6184.80,      6184.60,      6184.40,      6184.20,      6184.00,      6183.80, $
      6183.60,      6183.40,      6183.20,      6183.00,      6182.80,      6182.60,      6182.40,      6182.20, $
      6182.00,      6181.80,      6181.60,      6181.40,      6181.20,      6181.00,      6180.80,      6180.60, $
      6180.40,      6180.20,      6180.00,      6179.80,      6179.60,      6179.40,      6179.20,      6179.00, $
      6178.80,      6178.60,      6178.40,      6178.20,      6178.00,      6177.80,      6177.60,      6177.40, $
      6177.20,      6177.00,      6176.80,      6176.60,      6176.40,      6176.20,      6176.00,      6175.80, $
      6175.60,      6175.40,      6175.20,      6175.00,      6174.80,      6174.60,      6174.40,      6174.20, $
      6174.00,      6173.80,      6173.60,      6173.40,      6173.20,      6173.00,      6172.80,      6172.60, $
      6172.40,      6172.20,      6172.00,      6171.80,      6171.60,      6171.40,      6171.20,      6171.00, $
      6170.80,      6170.60,      6170.40,      6170.20,      6170.00,      6169.80,      6169.60,      6169.40, $
      6169.20,      6169.00,      6168.80,      6168.60,      6168.40,      6168.20,      6168.00,      6167.80, $
      6167.60,      6167.40,      6167.20,      6167.00,      6166.80,      6166.60,      6166.40,      6166.20, $
      6166.00,      6165.80,      6165.60,      6165.40,      6165.20,      6165.00,      6164.80,      6164.60, $
      6164.40,      6164.20,      6164.00,      6163.80,      6163.60,      6163.40,      6163.20,      6163.00, $
      6162.80,      6162.60,      6162.40,      6162.20,      6162.00,      6161.80,      6161.60,      6161.40, $
      6161.20,      6161.00,      6160.80,      6160.60,      6160.40,      6160.20,      6160.00,      6159.80, $
      6159.60,      6159.40,      6159.20,      6159.00,      6158.80,      6158.60,      6158.40,      6158.20, $
      6158.00,      6157.80,      6157.60,      6157.40,      6157.20,      6157.00,      6156.80,      6156.60, $
      6156.40,      6156.20,      6156.00,      6155.80,      6155.60,      6155.40,      6155.20,      6155.00, $
      6154.80,      6154.60,      6154.40,      6154.20,      6154.00,      6153.80,      6153.60,      6153.40, $
      6153.20,      6153.00,      6152.80,      6152.60,      6152.40,      6152.20,      6152.00,      6151.80, $
      6151.60,      6151.40,      6151.20,      6151.00,      6150.80,      6150.60,      6150.40,      6150.20, $
      6150.00],$
   [0.0209464,    0.0207487,    0.0209307,    0.0204282,    0.0226330,    0.0259792,    0.0238825,    0.0223389, $
    0.0248676,    0.0250958,    0.0278133,    0.0264017,    0.0299408,    0.0281424,    0.0271329,    0.0339798, $
    0.0330711,    0.0336281,    0.0343058,    0.0359879,    0.0367754,    0.0401164,    0.0439143,    0.0432530, $
    0.0486797,    0.0485797,    0.0540745,    0.0559437,    0.0620898,    0.0643225,    0.0693417,    0.0753219, $
    0.0804488,    0.0818621,    0.0948288,     0.108344,     0.110288,     0.123677,     0.134795,     0.152099, $
     0.161978,     0.173546,     0.189114,     0.207537,     0.232461,     0.272241,     0.286864,     0.322949, $
     0.365934,     0.401576,     0.450262,     0.507517,     0.575185,     0.622389,     0.711036,     0.915612, $
      1.02633,      1.19612,      1.34299,      1.53329,      1.83017,      2.11748,      2.42119,      2.70147, $
      3.27063,      4.35803,      4.97535,      5.79359,      6.70699,      7.53426,      9.23128,      10.3604, $
      11.8265,      13.0296,      15.3031,      17.3033,      21.7173,      23.6695,      26.2661,      28.5152, $
      32.0762,      34.5975,      37.2384,      38.8035,      42.3256,      46.5159,      49.0429,      51.1148, $
      53.1968,      54.6480,      57.0750,      58.1582,      59.8266,      60.4464,      62.1029,      63.4954, $
      64.0896,      64.7526,      65.0766,      65.4184,      65.4650,      65.5828,      65.4660,      65.6286, $
      65.1849,      64.7157,      64.1301,      63.6812,      62.7333,      62.1459,      60.7570,      59.3955, $
      58.2497,      57.0772,      54.9618,      51.3640,      49.1468,      46.8881,      44.3560,      42.1947, $
      38.4715,      36.0767,      33.4455,      31.9405,      28.4878,      26.3633,      21.4035,      19.9408, $
      17.5336,      16.0578,      13.5532,      12.2367,      10.6749,      9.85849,      8.32412,      6.53362, $
      5.73421,      4.98798,      4.35895,      3.85124,      3.20520,      2.80167,      2.49949,      2.28149, $
      1.90178,      1.51771,      1.33457,      1.19029,      1.05927,     0.958023,     0.847670,     0.740491, $
     0.667113,     0.614547,     0.535866,     0.452214,     0.409127,     0.356651,     0.341673,     0.311694, $
     0.272536,     0.255259,     0.239466,     0.227430,     0.204322,     0.191953,     0.173912,     0.161267, $
     0.151049,     0.148481,     0.146540,     0.135901,     0.125108,     0.127142,     0.117723,     0.120277, $
     0.111349,     0.107549,     0.102757,     0.104028,     0.100459,    0.0961026,    0.0968998,    0.0942963, $
    0.0974292,    0.0877748,    0.0869452,    0.0885202,    0.0882834,    0.0830996,    0.0801553,    0.0828427, $
    0.0806648,    0.0762436,    0.0789970,    0.0758304,    0.0782131,    0.0713996,    0.0747684,    0.0723149, $
    0.0701790]]

END

PRO HMI_tuningpolarizer2

nseq    = 31
nparam  = 7
lamref  = 6173.3433d0
lam     = 6173.405d0-lamref

FSR     = DBLARR(7)    ; Free Spectral Ranges in Angstrom of the Lyot and Michelson elements
FSR[0]  = 0.172d0      ; for the narrow-band Michelson in Angstrom
FSR[1]  = 0.344d0      ; for the broad-band  Michelson
FSR[2]  = 0.693d0      ; for E1
FSR[3]  = 1.405d0      ; for E2
FSR[4]  = 2.779d0      ; for E3
FSR[5]  = 5.682d0      ; for E4
FSR[6]  = 11.354d0     ; for E5

sign    = [1,-1]       ; ACCORDING TO JESPER: POLARIZER ROTATES IN THE SAME DIRECTION AS TUNING WAVEPLATE OF NB, AND OPPOSITE DIRECTION AS WAVEPLATE OF WB

tuning       = DBLARR(4,nseq) ; last column is the polarizer
tuning[*,0]  = [         0.d0,          0.d0 ,         0.d0,    0.d0]
tuning[*,1]  = [        80.d0,          0.d0 ,         0.d0,    0.d0]
tuning[*,2]  = [       160.d0,          0.d0 ,         0.d0,    0.d0]
tuning[*,3]  = [         0.d0,         80.d0,          0.d0,    0.d0]
tuning[*,4]  = [        80.d0,         80.d0,          0.d0,    0.d0]
tuning[*,5]  = [       160.d0,         80.d0,          0.d0,    0.d0]
tuning[*,6]  = [         0.d0,        160.d0,          0.d0,    0.d0]
tuning[*,7]  = [        80.d0,        160.d0,          0.d0,    0.d0]
tuning[*,8]  = [       160.d0,        160.d0,          0.d0,    0.d0]
tuning[*,9 ] = [         0.d0,          0.d0,         80.d0,    0.d0]
tuning[*,10] = [        80.d0,          0.d0,         80.d0,    0.d0]
tuning[*,11] = [       160.d0,          0.d0,         80.d0,    0.d0]
tuning[*,12] = [         0.d0,         80.d0,         80.d0,    0.d0] 
tuning[*,13] = [        80.d0,         80.d0,         80.d0,    0.d0]
tuning[*,14] = [       160.d0,         80.d0,         80.d0,    0.d0]
tuning[*,15] = [         0.d0,        160.d0,         80.d0,    0.d0]
tuning[*,16] = [        80.d0,        160.d0,         80.d0,    0.d0]
tuning[*,17] = [       160.d0,        160.d0,         80.d0,    0.d0]
tuning[*,18] = [         0.d0,          0.d0,        160.d0,    0.d0]
tuning[*,19] = [        80.d0,          0.d0,        160.d0,    0.d0]
tuning[*,20] = [       160.d0,          0.d0,        160.d0,    0.d0]
tuning[*,21] = [         0.d0,         80.d0,        160.d0,    0.d0]
tuning[*,22] = [        80.d0,         80.d0,        160.d0,    0.d0]
tuning[*,23] = [       160.d0,         80.d0,        160.d0,    0.d0]
tuning[*,24] = [         0.d0,        160.d0,        160.d0,    0.d0]
tuning[*,25] = [        80.d0,        160.d0,        160.d0,    0.d0]
tuning[*,26] = [       160.d0,        160.d0,        160.d0,    0.d0]
tuning[*,27] = [       160.d0,        160.d0,          0.d0,   20.d0]
tuning[*,28] = [       160.d0,        160.d0,         80.d0,   40.d0]
tuning[*,29] = [       160.d0,        160.d0,        160.d0,   60.d0]
tuning[*,30] = [       160.d0,        160.d0,         40.d0,   80.d0]

FOR i=0,nseq-1 DO tuning[*,i]  = tuning[*,i] * [1.5d0,1.5d0,-1.5d0,3.d0]*!dpi/180.d0

B   = [0.985,0.975,0.95,0.97,0.96,0.96,0.96]
Phi = [-41.,83.,76.,-12.,+13.,-10.,+17]*!dpi/180.d0
I0  = 13432.

; ACTUAL SPATIALLY-AVERAGED FRONT WINDOW FOR HMI
PROFILES,transmission,wavelength,q
blocker0     = INTERPOL(transmission/100.d0,wavelength*10.d0-lamref,lam)
blocker0     = blocker0 * INTERPOL(q[*,1]/100.d0,q[*,0]+3.61328-lamref,lam); I center the profile

Inten   = DBLARR(nseq)
profile = blocker0*I0
FOR i=3,6 DO profile = profile*(1.d0+B[i]*COS(2.d0*!dpi*lam/FSR[i]+Phi[i]))/2.d0
FOR j=0,nseq-1 DO BEGIN
    Inten[j] = profile
    Inten[j] = Inten[j]*(1.d0+B[2]*COS(2.d0*!dpi*lam/FSR[2]+Phi[2]+tuning[2,j]))/2.d0
    FOR i=0,1 DO Inten[j] = Inten[j]*(1.d0+B[i]*COS(2.d0*!dpi*lam/FSR[i]+Phi[i]+tuning[i,j]+sign[i]*tuning[3,j]))/2.d0
ENDFOR

SAVE,inten,file='temp.bin'

; FIT
;-------------------------------------------------------------------------------------

;sign        = [-1,1]

maxsteps    = 105
Residual    = DBLARR(nseq)
history     = DBLARR(maxsteps,nparam)
dIntendI0   = DBLARR(nseq)
dIntendPhi0 = DBLARR(nseq)
dIntendPhi1 = DBLARR(nseq)
dIntendPhi2 = DBLARR(nseq)
dIntendB0   = DBLARR(nseq)
dIntendB1   = DBLARR(nseq)
dIntendB2   = DBLARR(nseq)
err         = DBLARR(maxsteps)

converg     = 0
jj          = 0

Bg          = [0.97,0.97,0.97,0.97,0.97,0.97,0.97]
Phig        = [0.,0.,0.,0.,0.,0.,0.]*!dpi/180.d0

profilef=blocker0
FOR i=3,6 DO profilef = profilef*(1.d0+Bg[i]*COS(2.d0*!dpi*lam/FSR[i]+Phig[i]))/2.d0


I0g = 8.d0/DOUBLE(27)*TOTAL(Inten[0:26])/profilef
PRINT,'GUESS I0=',I0g

WHILE (converg EQ 0) DO BEGIN

    FOR j=0,nseq-1 DO BEGIN

        profileg = 1.d0
        profileg = profileg*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0
        FOR i=0,1 DO profileg = profileg*(1.d0+Bg[i]*COS(2.d0*!dpi*lam/FSR[i]+Phig[i]+tuning[i,j]+sign[i]*tuning[3,j]))/2.d0


        Residual[j]    = Inten[j] - I0g*profilef*profileg

        dIntendI0[j]   =  profilef*profileg
        dIntendPhi0[j] = -Bg[0]*SIN(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j])/2.d0*profilef*(1.d0+Bg[1]*COS(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[1]*tuning[3,j]))/2.d0*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0*I0g
        dIntendPhi1[j] = -Bg[1]*SIN(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[0]*tuning[3,j])/2.d0*profilef*(1.d0+Bg[0]*COS(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j]))/2.d0*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0*I0g
        dIntendPhi2[j] = -Bg[2]*SIN(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j])/2.d0*profilef*(1.d0+Bg[1]*COS(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[1]*tuning[3,j]))/2.d0*(1.d0+Bg[0]*COS(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j]))/2.d0*I0g
        dIntendB0[j]   =  COS(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j])/2.d0*profilef*(1.d0+Bg[1]*COS(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[1]*tuning[3,j]))/2.d0*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0*I0g
        dIntendB1[j]   =  COS(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[1]*tuning[3,j])/2.d0*profilef*(1.d0+Bg[0]*COS(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j]))/2.d0*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0*I0g
        dIntendB2[j]   =  COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j])/2.d0*profilef*(1.d0+Bg[1]*COS(2.d0*!dpi*lam/FSR[1]+Phig[1]+tuning[1,j]+sign[1]*tuning[3,j]))/2.d0*(1.d0+Bg[0]*COS(2.d0*!dpi*lam/FSR[0]+Phig[0]+tuning[0,j]+sign[0]*tuning[3,j]))/2.d0*I0g
        
    ENDFOR

    Jac  = DBLARR(nparam,nseq)
    FOR i= 0,nseq-1 DO BEGIN
        
        Jac[0,i]  = dIntendI0[i]
        Jac[1,i]  = dIntendPhi0[i]
        Jac[2,i]  = dIntendPhi1[i]
        Jac[3,i]  = dIntendPhi2[i]
        Jac[4,i]  = dIntendB0[i]
        Jac[5,i]  = dIntendB1[i]
        Jac[6,i]  = dIntendB2[i]
        
    ENDFOR

    LA_SVD,Jac,W,U,V,/DOUBLE
    Dx     = V##DIAG_MATRIX(1.d0/W)##TRANSPOSE(U)##Residual

    temp   = TRANSPOSE(Dx)/[I0g,Phig[0],Phig[1],Phig[2],Bg[0],Bg[1],Bg[2]]
                err[jj]= MAX(ABS(temp))
                
              ; IF(FINITE(err[jj]) EQ 0 AND jj GT 0)          THEN mu = 2.5d0*mu ELSE mu = mu0
                IF(jj EQ maxsteps-1 AND err[jj] GT 1.d-7 )    THEN converg = 2 ; no convergence 
                IF(err[jj] LE 1.d-7)                          THEN converg = 1 ; convergence
                
                IF(converg EQ 2) THEN BEGIN
                    j       = WHERE(err EQ MIN(err))
                    I0g     = history[j[0],0]
                    Phig[0] = history[j[0],1]
                    Phig[1] = history[j[0],2]
                    Phig[2] = history[j[0],3]
                    Bg[0]   = history[j[0],4]
                    Bg[1]   = history[j[0],5]
                    Bg[2]   = history[j[0],6]
                ENDIF
                
                IF(converg NE 2) THEN BEGIN

                    I0g = I0g+Dx[0]
                    Phig[0] = Phig[0]+Dx[1]
                    IF(ABS(Phig[0]) GT 1.2d0*!dpi) THEN Phig[0] = 0.
                    Phig[1] = Phig[1]+Dx[2]
                    IF(ABS(Phig[1]) GT 1.2d0*!dpi) THEN Phig[1] = 0. 
                    Phig[2] = Phig[2]+Dx[3]
                    IF(ABS(Phig[2]) GT 1.2d0*!dpi) THEN Phig[2] = 0.
                    Bg[0]   = Bg[0]  +Dx[4]
                    IF(Bg[0] LT 0.0) THEN Bg[0] = 0.97
                    Bg[1]   = Bg[1]  +Dx[5]
                    IF(Bg[1] LT 0.0) THEN Bg[1] = 0.97
                    Bg[2]   = Bg[2]  +Dx[6]
                    IF(Bg[2] LT 0.0) THEN Bg[2] = 0.97

                ENDIF
              

               history[jj,*]     = [I0g,Phig[0],Phig[1],Phig[2],Bg[0],Bg[1],Bg[2]]

               jj = jj + 1


ENDWHILE

Intenrecons = FLTARR(nseq)
profilef=blocker0
FOR i=3,6 DO profilef = profilef*(1.d0+Bg[i]*COS(2.d0*!dpi*lam/FSR[i]+Phig[i]))/2.d0
FOR j=0,nseq-1 DO BEGIN
    profileg = 1.d0
    profileg = profileg*(1.d0+Bg[2]*COS(2.d0*!dpi*lam/FSR[2]+Phig[2]+tuning[2,j]))/2.d0
    FOR i=0,1 DO profileg = profileg*(1.d0+Bg[i]*COS(2.d0*!dpi*lam/FSR[i]+Phig[i]+tuning[i,j]+sign[i]*tuning[3,j]))/2.d0
    Intenrecons[j]=profileg*profilef*I0g
ENDFOR


PRINT,'ACTUAL VALUES=',I0,Phi,B
PRINT,'FITTED VALUES=',I0g,Phig,Bg

SET_PLOT,'x'
WINDOW,0,RETAIN=2
PLOT,Inten,xst=1,psym=10
oplot,Intenrecons,col=180,linestyle=2,thick=2,psym=10
READ,pause




END
