PRO clean,imin,imout,headers
; Takes images in imin and headers and subtracts interpolated darks.
; Result is in imout.

nim=n_elements(headers)
print,nim,n_elements(imin),sqrt(n_elements(imin)/nim)
nbin=sqrt(n_elements(imin)/nim)

; Find darks, modA and modC
exposure=getpar(headers,'HSHIEXP')

wdark=where(exposure eq 0,ndark)

;time=gettime(headers) ; time from the keywords
;names=STRARR(53)
names=STRARR(19)

;names[0]='i_061028_213224.fit'
;names[1]='i_061028_213252.fit'
;names[2]='i_061028_213318.fit'
;names[3]='i_061028_213344.fit'
;names[4]='i_061028_213410.fit'
;names[5]='i_061028_213437.fit'
;names[6]='i_061028_213502.fit'
;names[7]='i_061028_213528.fit'
;names[8]='i_061028_213554.fit'
;names[9]='i_061028_213620.fit'
;names[10]='i_061028_213646.fit'
;names[11]='i_061028_213712.fit'
;names[12]='i_061028_213737.fit'
;names[13]='i_061028_213802.fit'
;names[14]='i_061028_213827.fit'
;names[15]='i_061028_213852.fit'
;names[16]='i_061028_213917.fit'
;names[17]='i_061028_213942.fit'
;names[18]='i_061028_214008.fit'
;names[19]='i_061028_214033.fit'
;names[20]='i_061028_214059.fit'
;names[21]='i_061028_214125.fit'
;names[22]='i_061028_214151.fit'
;names[23]='i_061028_214217.fit'
;names[24]='i_061028_214244.fit'
;names[25]='i_061028_214310.fit'
;names[26]='i_061028_214336.fit'
;names[27]='i_061028_214402.fit'
;names[28]='i_061028_214428.fit'
;names[0]='i_061028_214454.fit'
;names[1]='i_061028_214522.fit'
;names[2]='i_061028_214548.fit'
;names[3]='i_061028_214614.fit'
;names[4]='i_061028_214640.fit'
;names[5]='i_061028_214707.fit'
;names[6]='i_061028_214732.fit'
;names[7]='i_061028_214758.fit'
;names[8]='i_061028_214824.fit'
;names[9]='i_061028_214850.fit'
;names[10]='i_061028_214916.fit'
;names[11]='i_061028_214942.fit'
;names[12]='i_061028_215008.fit'
;names[13]='i_061028_215034.fit'
;names[14]='i_061028_215101.fit'
;names[15]='i_061028_215126.fit'
;names[16]='i_061028_215152.fit'
;names[17]='i_061028_215218.fit'
;names[18]='i_061028_215244.fit'
;names[19]='i_061028_215310.fit'
;names[20]='i_061028_215336.fit'
;names[21]='i_061028_215402.fit'
;names[22]='i_061028_215427.fit'
;names[23]='i_061028_215451.fit'
;names[24]='i_061028_215517.fit'
;names[25]='i_061028_215542.fit'
;names[26]='i_061028_215607.fit'
;names[27]='i_061028_215633.fit'
;names[28]='i_061028_215658.fit'
;names[0]='i_061028_215730.fit'
;names[1]='i_061028_215756.fit'
;names[2]='i_061028_215822.fit'
;names[3]='i_061028_215848.fit'
;names[4]='i_061028_215914.fit'
;names[5]='i_061028_215941.fit'
;names[6]='i_061028_220007.fit'
;names[7]='i_061028_220033.fit'
;names[8]='i_061028_220059.fit'
;names[9]='i_061028_220124.fit'
;names[10]='i_061028_220150.fit'
;names[11]='i_061028_220216.fit'
;names[12]='i_061028_220242.fit'
;names[13]='i_061028_220308.fit'
;names[14]='i_061028_220336.fit'
;names[15]='i_061028_220402.fit'
;names[16]='i_061028_220428.fit'
;names[17]='i_061028_220454.fit'
;names[18]='i_061028_220519.fit'
;names[19]='i_061028_220544.fit'
;names[20]='i_061028_220609.fit'
;names[21]='i_061028_220634.fit'
;names[22]='i_061028_220658.fit'
;names[23]='i_061028_220723.fit'
;names[24]='i_061028_220749.fit'
;names[25]='i_061028_220814.fit'
;names[26]='i_061028_220839.fit'
;names[27]='i_061028_220907.fit'
;names[28]='i_061028_220933.fit'
;names[29]='i_061028_220958.fit'
;names[30]='i_061028_221024.fit'
;names[31]='i_061028_221050.fit'
;names[32]='i_061028_221116.fit'
;names[33]='i_061028_221142.fit'
;names[34]='i_061028_221209.fit'
;names[35]='i_061028_221235.fit'
;names[36]='i_061028_221301.fit'
;names[37]='i_061028_221327.fit'
;names[38]='i_061028_221353.fit'
;names[39]='i_061028_221419.fit'
;names[40]='i_061028_221446.fit'
;names[41]='i_061028_221512.fit'
;names[42]='i_061028_221538.fit'
;names[43]='i_061028_221604.fit'
;names[44]='i_061028_221630.fit'
;names[45]='i_061028_221656.fit'
;names[46]='i_061028_221722.fit'
;names[47]='i_061028_221748.fit'
;names[48]='i_061028_221813.fit'
;names[49]='i_061028_221838.fit'
;names[50]='i_061028_221903.fit'
;names[51]='i_061028_221927.fit'
;names[52]='i_061028_221953.fit'
;names[0]='i_061029_005328.fit'
;names[1]='i_061029_005354.fit'
;names[2]='i_061029_005418.fit'
;names[3]='i_061029_005442.fit'
;names[4]='i_061029_005506.fit'
;names[5]='i_061029_005530.fit'
;names[6]='i_061029_005554.fit'
;names[7]='i_061029_005618.fit'
;names[8]='i_061029_005641.fit'
;names[9]='i_061029_005704.fit'
;names[10]='i_061029_005727.fit'
;names[11]='i_061029_005749.fit'
;names[12]='i_061029_005813.fit'
;names[13]='i_061029_005837.fit'
;names[14]='i_061029_005902.fit'
;names[15]='i_061029_005927.fit'
;names[16]='i_061029_005951.fit'
;names[17]='i_061029_010015.fit'
;names[18]='i_061029_010038.fit'
;names[19]='i_061029_010102.fit'
;names[20]='i_061029_010126.fit'
;names[21]='i_061029_010150.fit'
;names[22]='i_061029_010214.fit'
;names[23]='i_061029_010238.fit'
;names[24]='i_061029_010302.fit'
;names[25]='i_061029_010327.fit'
;names[26]='i_061029_010351.fit'
;names[27]='i_061029_010416.fit'
;names[28]='i_061029_010439.fit'
;names[29]='i_061029_010502.fit'
;names[30]='i_061029_010525.fit'
;names[31]='i_061029_010548.fit'
;names[32]='i_061029_010611.fit'
;names[33]='i_061029_010634.fit'
;names[34]='i_061029_010656.fit'
;names[35]='i_061029_010719.fit'
;names[36]='i_061029_010744.fit'
;names[37]='i_061029_010807.fit'
;names[38]='i_061029_010831.fit'
;names[39]='i_061029_010856.fit'
;names[40]='i_061029_010921.fit'
;names[41]='i_061029_010945.fit'
;names[42]='i_061029_011008.fit'
;names[43]='i_061029_011031.fit'
;names[44]='i_061029_011053.fit'
;names[45]='i_061029_011116.fit'
;names[46]='i_061029_011139.fit'
;names[47]='i_061029_011203.fit'
;names[48]='i_061029_011226.fit'
;names[49]='i_061029_011250.fit'
;names[50]='i_061029_011314.fit'
;names[51]='i_061029_011338.fit'
;names[52]='i_061029_011403.fit'
;names[0]='i_061031_205938.fit'
;names[1]='i_061031_210002.fit'
;names[2]='i_061031_210025.fit'
;names[3]='i_061031_210048.fit'
;names[4]='i_061031_210111.fit'
;names[5]='i_061031_210134.fit'
;names[6]='i_061031_210157.fit'
;names[7]='i_061031_210220.fit'
;names[8]='i_061031_210244.fit'
;names[9]='i_061031_210307.fit'
;names[10]='i_061031_210330.fit'
;names[11]='i_061031_210353.fit'
;names[12]='i_061031_210416.fit'
;names[13]='i_061031_210439.fit'
;names[14]='i_061031_210502.fit'
;names[15]='i_061031_210524.fit'
;names[16]='i_061031_210546.fit'
;names[17]='i_061031_210608.fit'
;names[18]='i_061031_210629.fit'
;names[19]='i_061031_210652.fit'
;names[20]='i_061031_210714.fit'
;names[21]='i_061031_210736.fit'
;names[22]='i_061031_210758.fit'
;names[23]='i_061031_210820.fit'
;names[24]='i_061031_210842.fit'
;names[25]='i_061031_210904.fit'
;names[26]='i_061031_210927.fit'
;names[27]='i_061031_210950.fit'
;names[28]='i_061031_211013.fit'
;names[0]='i_061031_200735.fit'
;names[1]='i_061031_200759.fit'
;names[2]='i_061031_200824.fit'
;names[3]='i_061031_200846.fit'
;names[4]='i_061031_200909.fit'
;names[5]='i_061031_200932.fit'
;names[6]='i_061031_200955.fit'
;names[7]='i_061031_201018.fit'
;names[8]='i_061031_201041.fit'
;names[9]='i_061031_201104.fit'
;names[10]='i_061031_201127.fit'
;names[11]='i_061031_201150.fit'
;names[12]='i_061031_201213.fit'
;names[13]='i_061031_201237.fit'
;names[14]='i_061031_201300.fit'
;names[15]='i_061031_201322.fit'
;names[16]='i_061031_201344.fit'
;names[17]='i_061031_201406.fit'
;names[18]='i_061031_201428.fit'
;names[19]='i_061031_201450.fit'
;names[20]='i_061031_201513.fit'
;names[21]='i_061031_201535.fit'
;names[22]='i_061031_201558.fit'
;names[23]='i_061031_201621.fit'
;names[24]='i_061031_201644.fit'
;names[25]='i_061031_201707.fit'
;names[26]='i_061031_201730.fit'
;names[27]='i_061031_201753.fit'
;names[28]='i_061031_201817.fit'
names[0] ='i_061215_233352.fit'
names[1] ='i_061215_233418.fit'
names[2] ='i_061215_233444.fit'
names[3] ='i_061215_233510.fit'
names[4] ='i_061215_233535.fit'
names[5] ='i_061215_233559.fit'
names[6] ='i_061215_233624.fit'
names[7] ='i_061215_233650.fit'
names[8] ='i_061215_233716.fit'
names[9] ='i_061215_233743.fit'
names[10]='i_061215_233809.fit'
names[11]='i_061215_233835.fit'
names[12]='i_061215_233900.fit'
names[13]='i_061215_233926.fit'
names[14]='i_061215_233952.fit'
names[15]='i_061215_234018.fit'
names[16]='i_061215_234044.fit'
names[17]='i_061215_234110.fit'
names[18]='i_061215_234137.fit'
names[0]='i_061215_234720.fit'
names[1]='i_061215_234747.fit'
names[2]='i_061215_234813.fit'
names[3]='i_061215_234839.fit'
names[4]='i_061215_234905.fit'
names[5]='i_061215_234931.fit'
names[6]='i_061215_234956.fit'
names[7]='i_061215_235023.fit'
names[8]='i_061215_235048.fit'
names[9]='i_061215_235113.fit'
names[10]='i_061215_235138.fit'
names[11]='i_061215_235203.fit'
names[12]='i_061215_235228.fit'
names[13]='i_061215_235254.fit'
names[14]='i_061215_235320.fit'
names[15]='i_061215_235346.fit'
names[16]='i_061215_235413.fit'
names[17]='i_061215_235439.fit'
names[18]='i_061215_235505.fit'
names[0]='i_061216_000222.fit'
names[1]='i_061216_000250.fit'
names[2]='i_061216_000316.fit'
names[3]='i_061216_000342.fit'
names[4]='i_061216_000407.fit'
names[5]='i_061216_000433.fit'
names[6]='i_061216_000459.fit'
names[7]='i_061216_000525.fit'
names[8]='i_061216_000551.fit'
names[9]='i_061216_000617.fit'
names[10]='i_061216_000641.fit'
names[11]='i_061216_000707.fit'
names[12]='i_061216_000733.fit'
names[13]='i_061216_000759.fit'
names[14]='i_061216_000825.fit'
names[15]='i_061216_000850.fit'
names[16]='i_061216_000916.fit'
names[17]='i_061216_000942.fit'
names[18]='i_061216_001008.fit'
names[0]='i_061216_002942.fit'
names[1]='i_061216_003009.fit'
names[2]='i_061216_003034.fit'
names[3]='i_061216_003059.fit'
names[4]='i_061216_003125.fit'
names[5]='i_061216_003152.fit'
names[6]='i_061216_003218.fit'
names[7]='i_061216_003244.fit'
names[8]='i_061216_003310.fit'
names[9]='i_061216_003335.fit'
names[10]='i_061216_003401.fit'
names[11]='i_061216_003427.fit'
names[12]='i_061216_003453.fit'
names[13]='i_061216_003519.fit'
names[14]='i_061216_003545.fit'
names[15]='i_061216_003611.fit'
names[16]='i_061216_003636.fit'
names[17]='i_061216_003701.fit'
names[18]='i_061216_003727.fit'
names[0]='i_061216_004722.fit'
names[1]='i_061216_004749.fit'
names[2]='i_061216_004815.fit'
names[3]='i_061216_004841.fit'
names[4]='i_061216_004907.fit'
names[5]='i_061216_004933.fit'
names[6]='i_061216_004959.fit'
names[7]='i_061216_005024.fit'
names[8]='i_061216_005049.fit'
names[9]='i_061216_005113.fit'
names[10]='i_061216_005139.fit'
names[11]='i_061216_005205.fit'
names[12]='i_061216_005231.fit'
names[13]='i_061216_005258.fit'
names[14]='i_061216_005324.fit'
names[15]='i_061216_005350.fit'
names[16]='i_061216_005415.fit'
names[17]='i_061216_005442.fit'
names[18]='i_061216_005507.fit'
names[0]='i_061216_005902.fit'
names[1]='i_061216_005929.fit'
names[2]='i_061216_005955.fit'
names[3]='i_061216_010021.fit'
names[4]='i_061216_010047.fit'
names[5]='i_061216_010114.fit'
names[6]='i_061216_010140.fit'
names[7]='i_061216_010205.fit'
names[8]='i_061216_010231.fit'
names[9]='i_061216_010257.fit'
names[10]='i_061216_010323.fit'
names[11]='i_061216_010349.fit'
names[12]='i_061216_010415.fit'
names[13]='i_061216_010441.fit'
names[14]='i_061216_010507.fit'
names[15]='i_061216_010533.fit'
names[16]='i_061216_010558.fit'
names[17]='i_061216_010623.fit'
names[18]='i_061216_010648.fit'
names[0]='i_061216_011004.fit'
names[1]='i_061216_011031.fit'
names[2]='i_061216_011057.fit'
names[3]='i_061216_011123.fit'
names[4]='i_061216_011148.fit'
names[5]='i_061216_011214.fit'
names[6]='i_061216_011240.fit'
names[7]='i_061216_011306.fit'
names[8]='i_061216_011332.fit'
names[9]='i_061216_011358.fit'
names[10]='i_061216_011424.fit'
names[11]='i_061216_011450.fit'
names[12]='i_061216_011515.fit'
names[13]='i_061216_011540.fit'
names[14]='i_061216_011605.fit'
names[15]='i_061216_011629.fit'
names[16]='i_061216_011654.fit'
names[17]='i_061216_011720.fit'
names[18]='i_061216_011746.fit'
names[0]='i_061216_012326.fit'
names[1]='i_061216_012352.fit'
names[2]='i_061216_012418.fit'
names[3]='i_061216_012443.fit'
names[4]='i_061216_012507.fit'
names[5]='i_061216_012532.fit'
names[6]='i_061216_012557.fit'
names[7]='i_061216_012624.fit'
names[8]='i_061216_012650.fit'
names[9]='i_061216_012716.fit'
names[10]='i_061216_012742.fit'
names[11]='i_061216_012808.fit'
names[12]='i_061216_012833.fit'
names[13]='i_061216_012859.fit'
names[14]='i_061216_012925.fit'
names[15]='i_061216_012951.fit'
names[16]='i_061216_013018.fit'
names[17]='i_061216_013043.fit'
names[18]='i_061216_013110.fit'


PRINT,names

time=getftime(names) ; time from the filename
time=time-total(rebin(time,1))

dark=rebin(imin(*,*,wdark),nbin,nbin,1)

c1=rebin(imin(0:nbin/64-1,0:nbin/64-1,*),1,1,nim)
c2=rebin(imin(nbin-nbin/64:nbin-1,0:nbin/64-1,*),1,1,nim)
c3=rebin(imin(0:nbin/64-1,nbin-nbin/64:nbin-1,*),1,1,nim)
c4=rebin(imin(nbin-nbin/64:nbin-1,nbin-nbin/64:nbin-1,*),1,1,nim)

cav0=reform(c1+c2+c3+c4)/4
; Try to separate gain variations and dark current variations
pcav=poly_fit(time,cav0,2)
cavfit=poly(time,pcav)
cavres=cav0-cavfit
cav=cavfit
cav=cav/total(rebin(cav(wdark),1))

imout=imin
for i=0,nim-1 do imout(*,*,i)=imin(*,*,i)-cav(i)*dark
intin=6*total(rebin(imout,1))

c1=rebin(imout(0:nbin/64-1,0:nbin/64-1,*),1,1,nim)
c2=rebin(imout(nbin-nbin/64:nbin-1,0:nbin/64-1,*),1,1,nim)
c3=rebin(imout(0:nbin/64-1,nbin-nbin/64:nbin-1,*),1,1,nim)
c4=rebin(imout(nbin-nbin/64:nbin-1,nbin-nbin/64:nbin-1,*),1,1,nim)
q=[[c1,c2],[c3,c4]]
imout=imout-rebin(q,nbin,nbin,nim,/sample)

END

;-------------------------------------------------------------------------




; Read images
;nbin=128l
;readimages,names0,imb0,headers0,nbin=nbin,/getnames,iover=iover0
;nim0=n_elements(headers0)
;
;badstr='qqq'
;; Read manually found bad images. Enter or -1 means all good.
;read,badstr,prompt='Bad points: '
;badlist=fix(strsplit(badstr+' -1',/extract))
;w=where(badlist ge 0,n)
;badmask=intarr(nim0)
;if (n ge 1) then badmask(badlist(w))=1
;
;corner0=reform(rebin(imb0(0:nbin/64-1,0:nbin/64-1,*),1,1,nim0))
;
;wgood=getgood(corner0,iover0)
;ngood=n_elements(wgood)
;wgood=wgood(where(badmask(wgood) eq 0,ngood))
;
;; Throw out bad images
;names=names0(wgood)
;imb=imb0(*,*,wgood)+0.0D0
;headers=headers0(wgood)
;iover=iover0(wgood)
;
;nim=ngood
;
;clean,imb,imd,headers
;
;; imd now contains the images with the darks subtracted.
;
;; Get header values.
;exposure=getpar(headers,'HSHIEXP')
