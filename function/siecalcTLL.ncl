undef("siecalcTLL")
function siecalcTLL(wrk[*][*][*]:float)
local re,min_lats,max_lats,criterion,rad,time,lat,lon,nt,ny,nx,dy,dx,wcos,ext
begin
   re = 6.37122e3
   rad = get_d2r("float")
   min_lats = 40.5
   max_lats = 89.5
   criterion = 0.15

   wrk!0 = "time"
   wrk!1 = "lat"
   wrk!2 = "lon"
   time = wrk&time
   lat = wrk&lat
   lon = wrk&lon
   nt = dimsizes(time)
   ny = dimsizes(lat)
   nx = dimsizes(lon)

   dy = new(ny,"float")
   dx = new(nx,"float")
   dy = abs(lat(1)-lat(0))*rad
   dx = abs(lon(1)-lon(0))*rad
   wcos = abs(cos(rad*lat))

   ext = new((/nt,ny,nx/),"float")
   ext = re^2*conform(ext,wcos,1)*conform(ext,dy,1)*conform(ext,dx,2)
      copy_VarCoords(wrk,ext)

   ext = where(wrk.gt.criterion,ext,getFillValue(ext))
   ext := dim_sum_n_Wrap(ext(:,{min_lats:max_lats},:),(/1,2/))/10^6
      ext!0 = "time"
      ext&time = time
      ext@units = "millions of km^2"

   return(ext)
end
