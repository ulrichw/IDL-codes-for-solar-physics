PRO temp2

year = 2006
month = 2
seconde = 0

day1 = 15
heure1 = 15
minute1 = 37 

    command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
    SPAWN,command,res1
    res1 = DOUBLE(res1)

day1 = 15
heure1 = 15
minute1 = 49

    command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
    SPAWN,command,res2
    res2 = DOUBLE(res2)
    
  ; radial velocity in m/s (v>0 FOR MOVEMENTS AWAY FROM SUN)
    radial = ( res2 - res1 )*149597870691.d0/12./60.d0 ; Sun-Earth distance in AU taken 6-minute apart

print,radial


END
