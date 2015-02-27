PRO solarradius

year=2009

radius=FLTARR(12,28)

FOR i=1,12 DO BEGIN
    FOR j=1,28 DO BEGIN
        command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(i),1)+'.'+STRTRIM(STRING(j),1)+'.'+STRTRIM(STRING(12),1)+'.'+STRTRIM(STRING(0),1)+'.'+STRTRIM(STRING(0),1)+'| tail -14 | head -1 | cut -b11-22'
        SPAWN,command,res1
        PRINT,res1
        radius[i-1,j-1] = DOUBLE(res1)/2.d0
    ENDFOR
ENDFOR

READ,pause

END
