tyr = YYEAR
tmm = MMON 
tym = tyr*100+tmm
optMSLP = OOPTION

;----------------------------

fpath = "./data/HadISST_ice.nc"
in = addfile(fpath,"r")
sic = in->sic(:,{40.5:89.5},:)
   sic!1 = "lat"
   sic!2 = "lon"

sic = where(ismissing(sic),default_fillvalue("float"),sic)
   sic@_FillValue = default_fillvalue("float")
   sic@missing_value = default_fillvalue("float")

ym = cd_calendar(sic&time,-1)
iym = ind(ym/100.ge.1981.and.ym/100.le.2010)
sic_clim = clmMonTLL(sic(iym,:,:))
   sic_clim!0 = "mon"
   sic_clim&mon = ispan(1,12,1)
   delete(iym)

iym = ind(ym/100.ge.1979.and.ym/100.le.2022)
sic := sic(iym,:,:)
   delete([/ym,iym/])

ym = cd_calendar(sic&time,-1)
iym = ind(ym.eq.tym)
sic_clim_apd = array_append_record(sic_clim,sic_clim,0)
sic_for_pre = sic
sic_for_pre(iym+1:iym+12,:,:) = (/sic_clim_apd(tmm:tmm+11,:,:)/)
sic_anom = calcMonAnomTLL(sic_for_pre,sic_clim)
   delete(iym)

clat = sqrt(cos(get_d2r(sic_anom&lat)))
sic_anom_reshape = reshape(sic_anom,(/dimsizes(sic_anom&time)/12,12,dimsizes(sic_anom&lat),dimsizes(sic_anom&lon)/))
sic_anom_reshape_wgt = sic_anom_reshape*conform(sic_anom_reshape,clat,2)
   sic_anom_reshape_wgt!0 = "year"
   sic_anom_reshape_wgt!1 = "mon"
   sic_anom_reshape_wgt!2 = "lat"
   sic_anom_reshape_wgt!3 = "lon"
   sic_anom_reshape_wgt&year = ispan(1979,2022,1)
   sic_anom_reshape_wgt&mon = sic_clim&mon
   sic_anom_reshape_wgt&lat = sic_anom&lat
   sic_anom_reshape_wgt&lon = sic_anom&lon
   copy_VarAtts(sic_anom,sic_anom_reshape_wgt)

neof = dimsizes(sic_anom_reshape_wgt&year)-1	; dof
sic_eof = eofunc_n_Wrap(sic_anom_reshape_wgt,neof,True,0)
sic_ets = eofunc_ts_n_Wrap(sic_anom_reshape_wgt,sic_eof,False,0)

sic_eof_tmp = sic_eof/conform(sic_eof,clat,2)
sic_ets_tmp = sic_ets+conform(sic_ets,sic_ets@ts_mean,0)
   copy_VarMeta(sic_eof,sic_eof_tmp)
   copy_VarMeta(sic_ets,sic_ets_tmp)
   delete(clat)

sic_eof := sic_eof_tmp
sic_ets := sic_ets_tmp

load "./function/siecalcTLL.ncl"
sie = siecalcTLL(sic_for_pre)

tlg = ispan(0,12,1)
tld = ispan(0,11,1)
sie_cor = new((/dimsizes(tlg),dimsizes(tld)/),"float")

iym = ind(ym%100.eq.tmm)

do g = 0, dimsizes(tlg)-1
do d = 0, dimsizes(tld)-1
   iym_tlg = iym+g
   iym_tld = iym-d
   iym_tlg := iym_tlg(ind(iym_tlg.le.(dimsizes(ym)-1)))
   iym_tld := iym_tld(ind(iym_tld.ge.0))
   if (dimsizes(iym_tlg).gt.dimsizes(iym_tld)) then
      iym_tlg := iym_tlg(dimsizes(iym_tlg)-dimsizes(iym_tld):)
   else if (dimsizes(iym_tlg).lt.dimsizes(iym_tld)) then
      iym_tld := iym_tld(:dimsizes(iym_tlg)-1)
   end if
   end if
   sie_tlg = sie(iym_tlg)
   sie_tld = sie(iym_tld)
   sie_cor(g,d) = escorc(sie_tlg,sie_tld)
      delete([/iym_tlg,iym_tld,sie_tlg,sie_tld/])
end do
   sie_cor(g,:) = sie_cor(g,:)/avg(abs(sie_cor(g,:)))
end do
   sie_cor!0 = "lag"
   sie_cor!1 = "lead"
   sie_cor&lag = tlg
   sie_cor&lead = tld
    delete([/ym,iym/])

sie_cor := sie_cor(lag|:,lead|::-1)

if (optMSLP) then
   fpath = "./data/ERA5m_msl.nc"
   in = addfile(fpath,"r")
   msl = in->msl(:,{50.:90.},:)
      msl!1 = "lat"
      msl!2 = "lon"
   
   msl = msl/100.
   msl = where(ismissing(msl),default_fillvalue("float"),msl)
      msl@units = "hPa"
      msl@_FillValue = default_fillvalue("float")
      msl@missing_value = default_fillvalue("float")
   
   ym = cd_calendar(msl&time,-1)
   iym = ind(ym/100.ge.1981.and.ym/100.le.2010)
   msl_clim = clmMonTLL(msl(iym,:,:))
      msl_clim!0 = "mon"
      msl_clim&mon = ispan(1,12,1)
      delete(iym)
   
   iym = ind(ym/100.ge.1979.and.ym/100.le.2022)
   msl := msl(iym,:,:)
      delete([/ym,iym/])
   
   ym = cd_calendar(msl&time,-1)
   iym = ind(ym.eq.tym)
   msl_clim_apd = array_append_record(msl_clim,msl_clim,0)
   msl_for_pre = msl
   msl_for_pre(iym+1:iym+12,:,:) = (/msl_clim_apd(tmm:tmm+11,:,:)/)
      delete(iym)
   
   iym = ind(ym%100.eq.tmm)
   msl_anom = msl_for_pre(iym,:,:)-conform(msl_for_pre(iym,:,:),msl_clim({tmm},:,:),(/1,2/))
   
   clat = sqrt(cos(get_d2r("float")*msl_for_pre&lat))
   msl_anom_wgt = msl_anom*conform(msl_anom,clat,1)
      msl_anom_wgt!0 = "year"
      msl_anom_wgt!1 = "lat"
      msl_anom_wgt!2 = "lon"
      msl_anom_wgt&year = ispan(1979,2022,1)
      msl_anom_wgt&lat = msl_for_pre&lat
      msl_anom_wgt&lon = msl_for_pre&lon
      copy_VarAtts(msl_for_pre,msl_anom_wgt)
      delete(clat)
   
   msl_eof = eofunc_n_Wrap(msl_anom_wgt,1,True,0)
   msl_ets = eofunc_ts_n_Wrap(msl_anom_wgt,msl_eof,False,0)
   msl_ets := msl_ets(0,:)/stddev(msl_ets(0,:))
      msl_ets!0 = "year"
      msl_ets&year = ispan(1979,2022,1)
   
   ym = cd_calendar(sic_anom&time,-1)

   msl_sic_a = new((/dimsizes(tlg),dimsizes(sic_anom&lat),dimsizes(sic_anom&lon)/),"float")
   msl_sic_b = msl_sic_a
   
   do g = 0, dimsizes(tlg)-1
      iym_tlg = iym+g
      iym_tlg := iym_tlg(ind(iym_tlg.le.(dimsizes(ym)-1)))
      a = msl_ets					;  msl_ets      1d
      b = sic_anom(iym_tlg,:,:)		;  sic_anom_tlg 3d

      if (dimsizes(b&time).lt.dimsizes(a)) then
         a := a(:(dimsizes(b&time)-1))
      end if

      reg_tmp = regCoef_n(a,b,0,0)

      df = onedtond(reg_tmp@nptxy,dimsizes(reg_tmp))-1	;dof of auto linear regr = n-1
      tval = onedtond(reg_tmp@tval,dimsizes(reg_tmp))
      bb = new(dimsizes(reg_tmp),"float")
      bb = 0.5
      sig_tmp = betainc(df/(df+tval^2),df/2.,bb)
         copy_VarMeta(reg_tmp,sig_tmp)

      yintercept = onedtond(reg_tmp@yintercept,dimsizes(reg_tmp))
         copy_VarMeta(reg_tmp,yintercept)
      sig_tmp = where(sig_tmp.eq.0..or.sig_tmp.eq.1e+30,getFillValue(sig_tmp),sig_tmp)
      msl_sic_a(g,:,:) = where(sig_tmp.gt.0.05,0.,reg_tmp)
      msl_sic_b(g,:,:) = where(sig_tmp.gt.0.05,0.,yintercept)
         delete([/iym_tlg,a,b,reg_tmp,sig_tmp,df,tval,bb,yintercept/])
   end do
   delete([/ym,iym/])
end if

sic_eof := sic_eof(mon|:,evn|:,lat|:,lon|:)
sic_eof_apd = array_append_record(sic_eof,sic_eof,0)
sic_eof_apd := array_append_record(sic_eof_apd,sic_eof_apd,0)

ym = cd_calendar(sic&time,-1)
iym = ind(ym.eq.tym)
sic_anom_input = sic_anom(iym-11:iym,:,:)
sic_anom_obs = sic(iym:iym+12,:,:)-(/sic_clim_apd(tmm-1:tmm+11,:,:)/)

clat = sqrt(cos(get_d2r("float")*sic_anom_input&lat))
sic_anom_input_wgt = sic_anom_input*conform(sic_anom_input,clat^2,1)

sic_anom_pre = sic_anom_obs
sic_anom_pre = 0.

do g = 0, dimsizes(tlg)-1
do i = 0, neof-1
   pi = sum(sic_anom_input_wgt*conform(sic_anom_input_wgt,sie_cor(g,:),0)*sic_eof_apd(tmm:tmm+11,i,:,:))
   sic_anom_pre(g,:,:) = sic_anom_pre(g,:,:)+sic_eof_apd(tmm-1+g,i,:,:)*pi
end do
end do

if (optMSLP) then
   sic_anom_pre = sic_anom_pre+(msl_sic_a*msl_ets({tyr})+msl_sic_b)	;y=ax+b
end if

sic_anom_pre = where(ismissing(sic_anom_pre).and.sic_clim_apd(tmm-1:tmm+11,:,:).eq.0.,0.,sic_anom_pre)
sic_anom_pre = where(sic_anom_pre.gt. 1., 1.,sic_anom_pre)
sic_anom_pre = where(sic_anom_pre.lt.-1.,-1.,sic_anom_pre)

sic_obs = sic(iym:iym+12,:,:)
sic_pre = sic_anom_pre+sic_clim_apd(tmm-1:tmm+11,:,:)
sic_pre = where(sic_pre.gt.1.,1.,sic_pre)
sic_pre = where(sic_pre.lt.0.,0.,sic_pre)
   copy_VarMeta(sic_obs,sic_pre)
   copy_VarMeta(sic_obs,sic_anom_obs)
   copy_VarMeta(sic_obs,sic_anom_pre)

sie_obs = siecalcTLL(sic_obs)
sie_pre = siecalcTLL(sic_pre)

if (optMSLP) then
   fpath = "./output_msl/prediction_"+tym+".nc"
else 
   fpath = "./output/prediction_"+tym+".nc"
end if
system("echo "+fpath)
system("rm -f "+fpath)
out = addfile(fpath,"c")
out->sic_anom_obs = sic_anom_obs
out->sic_anom_pre = sic_anom_pre
out->sic_obs = sic_obs
out->sic_pre = sic_pre
out->sie_obs = sie_obs
out->sie_pre = sie_pre
