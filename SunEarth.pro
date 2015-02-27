year = 2006
month= INDGEN(12)+1
day  = INDGEN(30)+1
heure= 20
minute=10
seconde=0


res1 = DBLARR(12*30)
spin = DBLARR(24*60)

time = DBLARR(12*30)
time2= FLTARR(24*60)

FOR i=0,11 DO BEGIN
    PRINT,i
    FOR j=0,29 DO BEGIN
        command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month[i]),1)+'.'+STRTRIM(STRING(day[j]),1)+'.'+STRTRIM(STRING(heure),1)+'.'+STRTRIM(STRING(minute),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
        SPAWN,command,res
        res1[i*30+j] = DOUBLE(res) 
        time[i*30+j] = 86400.d0*30.d0*i+j*86400.d0

    ENDFOR
ENDFOR

month = 10
day   = 16
heure = INDGEN(23)+1
minute= INDGEN(59)+1

FOR i=0,22 DO BEGIN
    PRINT,i
    FOR j=0,58 DO BEGIN
        command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day),1)+'.'+STRTRIM(STRING(heure[i]),1)+'.'+STRTRIM(STRING(minute[j]),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -1 | cut -b14-28'
        SPAWN,command,res
        spin[i*60.+j]   = DOUBLE(res)
        time2[i*60.+j]  = 3600.*i+j*60. 
    ENDFOR
ENDFOR


res1= res1*149597870691.d0

window,0,retain=2
!p.multi=[0,1,2]
plot,time,res1,xst=1
plot,time/(86400.d0*30.d0)+1,DERIV(time,res1),xst=1

window,1,retain=2
!p.multi=0
plot,time2/3600.+1,spin,xst=1,psym=3
oplot,[0,23],[0,0]


read,pause

END
