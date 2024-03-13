pro make_better_wee

  ; The osxp values from ell2 most "correct" as the engs match over files
  ; And trust the volume from the mf_fin_newel file so use them,
  ; and hence calc a new thermal energy for the files?
  ;
  ;
  ;
  ; 28-Feb-2024 IGH
  ; 13-Mar-2024 Add in vol error calc, but some len and wid err missing

  ; ospex correct in this one
  restgen,file='~/Desktop/wee_2007_bkhd/mf_ell2.genx',res
  ; vf_fit and vol correct in this one ????
  restgen,file='mf_fin_newel',rv

  ; GOES temps results
  restgen,file='goes_tem_all',tems

  ; Other data sources for filtering events
  ; this was laready calculated from mf_ell2.genx
  restore, file='modflux_48.dat'
  rat=nn/th
  lf=linfit(alog10([1d44,1d49]),alog10([1.5d45,3d49]))
  gfit=10d^(lf[0]+lf[1]*alog10(res.osx_p[0]*1d49))
  restgen,file='sig2back',bs

  ; Filter for the "good" events
  gdth=where((rv.vf_qflag eq 0.  or  rv.vf_qflag eq 8.) and rat lt 1 and $
    res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e2 and bs.sb48 gt 3.0 and $
    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and rv.vf_lpwid gt 2.3 and $
    res.bk_bf_flag eq 1,ngdth)
  gdnn=where((rv.vf_qflag eq 0.  or  rv.vf_qflag eq 8.) and rat lt 1 $
    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e-2 and $
    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and $
    res.osx_p[3] ne 1e-5 and bs.sb612 gt 3.0 and $
    res.osx_p[5] ne 20 and res.osx_p[5] gt 7.2 and abs(res.osx_p[5]-7.16) gt 0.1 and res.osx_p[5] ne 9. and $
    res.osx_p[6] ne 2 and res.osx_p[6] ne 12 and rv.vf_lpwid gt 2.3 and $
    res.bk_bf_flag eq 1,ngdnn)

  ;  Added in GOES restrictiond of "good" events
  gsgdth=where((rv.vf_qflag eq 0.  or  rv.vf_qflag eq 8.) and rat lt 1 $
    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e2 and $
    res.osx_p[1] gt 5. and res.osx_p[1] ne 11. and rv.vf_lpwid gt 2.3 $
    and res.bk_bf_flag eq 1 and res.gflx_bs gt 0. and bs.sb48 gt 3. and $
    tems.em ne tems[1].em and tems.tmk gt 4.0 and tems.tmk lt 20. and tems.bsub eq 1$
    and tems.em lt gfit ,ngd)
  gsgdnn=where((rv.vf_qflag eq 0.  or  rv.vf_qflag eq 8.) and rat lt 1 $
    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e-2 and $
    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and rv.vf_lpwid gt 2.3 $
    and res.gflx_bs gt 0. and bs.sb48 gt 3. and $
    tems.em ne tems[1].em and tems.tmk gt 4.0 and tems.tmk lt 20. and tems.bsub eq 1$
    and tems.em lt gfit and res.osx_p[3] ne 1e-5 and bs.sb612 gt 3.0 and $
    res.osx_p[5] ne 20 and res.osx_p[5] gt 7.2 and abs(res.osx_p[5]-7.16) gt 0.1 and res.osx_p[5] ne 9. and $
    res.osx_p[6] ne 2 and res.osx_p[6] ne 12,ngdnn)

  ec=eb2ec(res.osx_p[5],res.osx_p[6])

  ; Calc new thermal energy from ell2 osxp and vol from fin_newl
  kb=1.3806503d-23
  p=res.osx_p
  em=p[0]*1d49
  tk=p[1]*1d6

  rad=rv.vf_lpwid*0.5
  len=rv.vf_lparc
  ; this matches vf_vol in fin_newel file and fills in the 0 values.
  new_vol=!PI*rad^2*len*(725d5)^3
  new_ength=3*sqrt(em*new_vol)*kb*tk*1d7

  ; Do I actually trust these.... probably not
  raderr=rv.vf_lpwiderr*0.5
  lenerr=rv.vf_lparcerr
  rpart=(2*rad*len*raderr)^2
  lpart=(rad^2*lenerr)^2
  new_volerr=!PI*(rpart+lpart)^0.5*(725d5)^3

  ; Check fin_newel non-therm pow is same as save
  ; and also calc using "Ec" instead of eps_B
  nf=n_elements(res.fpeak)
  ee_nn=dblarr(nf)
  ee_nn_ec=dblarr(nf)
  f_1out=dblarr(nf)
  for i=0, nf-1 do begin
    p=res[i].osx_p
    f_50=p[3]
    eb=p[5]
    g1=p[4]
    g2=p[6]
    f_eb=f_50*(eb/50.)^(-1*g1)
    f_1=f_eb*eb^g2
    f_1out[i]=f_1
    ee_nn[i]=16.*$
      9.5d24*g2^2*(g2-1)*beta(g2-.5,1.5)*f_1*eb^(1-g2)
    ee_nn_ec[i]=16.*$
      9.5d24*g2^2*(g2-1)*beta(g2-.5,1.5)*f_1*ec[i]^(1-g2)
  endfor
  ; So caclulated version matches saved!
  ;  plot_oo,res.eng_nn,ee_nn,psym=1


  ; Check with histogram plots
  ; loop and volumes
  !p.multi=[0,2,1]
  loadct,0,/silent
  lhist=histogram(alog10(rv.vf_lparc),min=0,max=3,nbins=50,locations=lbin)
  whist=histogram(alog10(rv.vf_lpwid),min=0,max=3,nbins=50,locations=lbin)

  ; Basically same as paper
  window, 1, xsize=800,ysize=400
  plot,10^lbin,lhist,psym=10,/xlog,xtitle='Loop Size [arcsec]'
  oplot,10^lbin,whist,psym=10,color=150
  vhist=histogram(alog10(new_vol),min=24,max=29,nbins=50,locations=vbin)
  plot,10^vbin,vhist,psym=10,/xlog,xtit='Loop Volume [cm^3]'

  thist=histogram(res[gdth].osx_p[1],min=2,max=22,nbins=50,locations=tbin)
  emhist=histogram(alog10(res[gdth].osx_p[0]*1d49),min=44,max=48,nbins=50,locations=embin)

  ; Basically same as paper
  window, 2, xsize=800,ysize=400
  plot,tbin,thist,psym=10,xtitle='Temperature [MK]'
  plot,10^embin,emhist,psym=10,/xlog,xtit='Emission Measure [cm^-3]'

  ghist=histogram(res[gdnn].osx_p[6],min=2,max=12,nbins=50,locations=gbin)
  ebhist=histogram(res[gdnn].osx_p[5],min=5,max=16,nbins=50,locations=ebbin)
  eb2hist=histogram(rv[gdnn].osx_p[5],min=5,max=16,nbins=50,locations=ebbin)
  echist=histogram(ec[gdnn],min=5,max=20,nbins=50,locations=ecbin)

  !p.multi=[0,3,1]
  ; Not quite same as paper -> something in epsilon_b?
  window, 3, xsize=1200,ysize=400
  plot,gbin,ghist,psym=10,xtitle='Gamma',title=mean(res[gdnn].osx_p[6])
  plot,ebbin,ebhist,psym=10,xtit='Break Energy [keV]',title=mean(res[gdnn].osx_p[5])
  ;  oplot,ebbin,eb2hist,psym=10,color=150
  plot,ecbin,echist,psym=10,xtit='Low Energy Cutoff* [keV]',title=mean(ec[gdnn])

  maxa=31
  mina=24
  nbs=70.
  binsize=(maxa-mina)/nbs

  hth=histogram(alog10(new_ength[gdth]),min=mina,max=maxa,binsize=binsize,locations=xa)
  heb=histogram(alog10(ee_nn[gdnn]),min=mina,max=maxa,binsize=binsize,locations=xa)
  hec=histogram(alog10(ee_nn_ec[gdnn]),min=mina,max=maxa,binsize=binsize,locations=xa)

  eda=10^[xa,maxa]
  midsa=get_edges(eda,/mean)
  ewid=get_edges(eda,/width)
  tdiff=anytim('1-Mar-2007')-anytim('1-Mar-2002')
  restgen,file='bimon_frac',frac
  mnfrac=mean(frac)
  tottime=tdiff*mnfrac
  factor=2*!PI*(6.955e10^2)*(tottime)/1d50

  hth=hth/(ewid*factor)
  heb=heb/(ewid*factor)
  hec=hec/(ewid*factor)
  !p.multi=0
  ; Looks reasonable compared to the paper
  window, 4,xsize=600,ysize=400
  plot,midsa,hth,psym=10,xrange=[5e25,5e30],/xlog,$
    ystyle=17,yrange=[1e-10,1e-4],$
    ytit='Flare Frequency [10!U-50!N s!U-1!N cm!U-2!N ergs!U-1!N]',$
    xtit='Energy [ergs]',$
    xtickf='exp1',/ylog,ytickf='exp1',xstyle=17
  oplot,midsa,hec,psym=10,color=150

  nf=n_elements(res.osx_p[3])
  ph4_8=dblarr(nf)
  ph12=dblarr(nf)
  for i=0, nf-1 do begin
    ph4_8[i]=f_vth([4,8],[res[i].osx_p[0],res[i].osx_p[1]*0.086164])
    ph12[i]=f_bpow([12.],res[i].osx_p[3:6])
  endfor

  rr={fstart:anytim(res.fstart,/ccsds),fend:anytim(res.fend,/ccsds),fpeak:anytim(res.fpeak,/ccsds),$
    fpeak_tr:anytim(res.fpeak_tr,/ccsds),bk_bf_tr:anytim(res.bk_bf_tr,/ccsds),$
    tmk:res.osx_p[1],em:res.osx_p[0]*1d49,$
    norm:res.osx_p[3],f_1:f_1out,g1:1.5,eb:res.osx_p[5],g2:res.osx_p[6],ec:ec,$
    ph4_8:ph4_8,ph12:ph12,$
    vol:new_vol,volerr:new_volerr,lp_len:rv.vf_lparc, lp_wid:rv.vf_lpwid, lp_qflag:rv.vf_qflag,vx:rv.vf_fit.srcx,vy:rv.vf_fit.srcy,$
    gflx:res.gflx,gflx_bs:res.gflx_bs,gtmk:tems.tmk,gem:tems.em,$
    eng_th:new_ength,eng_nn_eb:res.eng_nn, eng_nn_ec:ee_nn_ec,$
    idgdth:gdth,idgdnn:gdnn,idgsgdth:gsgdth,idgsgdnn:gsgdnn}

  savegen,file='wee_all_v2.genx',rr

  stop
end