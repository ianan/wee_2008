pro rhessi_weeflares_results
  ;
  ;  ;   Taking old results files and making simiplier output file for sharing
  ;
  ;  ;    28-Apr-2023 IGH
  ;  ;    04-May-2023 Added more info to output genx
  ;  ;
  ;  ;
  ;;; *******************************************
  ;;; *******************************************
  ;;; *******************************************
  ;;; *******************************************
  
  ;  ; GOES temps results
  ;  restgen,file='goes_tem_all',tems
  ;
  ;  ; RHESSI results
  ;  restgen,file='mf_fin_newel',res
  ;
  ;  ; Other data sources for filtering events
  ;  restore, file='modflux_48.dat'
  ;  rat=nn/th
  ;  lf=linfit(alog10([1d44,1d49]),alog10([1.5d45,3d49]))
  ;  gfit=10d^(lf[0]+lf[1]*alog10(res.osx_p[0]*1d49))
  ;  restgen,file='sig2back',bs
  ;
  ;  ; Filter for the "good" events
  ;  gdth=where((res.vf_qflag eq 0.  or  res.vf_qflag eq 8.) and rat lt 1 and $
  ;    res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e2 and bs.sb48 gt 3.0 and $
  ;    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and res.vf_lpwid gt 2.3 and $
  ;    res.bk_bf_flag eq 1,ngdth)
  ;  gdnn=where((res.vf_qflag eq 0.  or  res.vf_qflag eq 8.) and rat lt 1 $
  ;    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e-2 and $
  ;    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and $
  ;    res.osx_p[3] ne 1e-5 and bs.sb612 gt 3.0 and $
  ;    res.osx_p[5] ne 20 and res.osx_p[5] gt 7.2 and abs(res.osx_p[5]-7.16) gt 0.1 and res.osx_p[5] ne 9. and $
  ;    res.osx_p[6] ne 2 and res.osx_p[6] ne 12 and res.vf_lpwid gt 2.3 and $
  ;    res.bk_bf_flag eq 1,ngdnn)
  ;
  ;  ;  Added in GOES restrictiond of "good" events
  ;  gsgdth=where((res.vf_qflag eq 0.  or  res.vf_qflag eq 8.) and rat lt 1 $
  ;    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e2 and $
  ;    res.osx_p[1] gt 5. and res.osx_p[1] ne 11. and res.vf_lpwid gt 2.3 $
  ;    and res.bk_bf_flag eq 1 and res.gflx_bs gt 0. and bs.sb48 gt 3. and $
  ;    tems.em ne tems[1].em and tems.tmk gt 4.0 and tems.tmk lt 20. and tems.bsub eq 1$
  ;    and tems.em lt gfit ,ngd)
  ;  gsgdnn=where((res.vf_qflag eq 0.  or  res.vf_qflag eq 8.) and rat lt 1 $
  ;    and res.osx_p[0] ne 1e-5 and res.osx_p[0] ne 1e-2 and $
  ;    res.osx_p[1] ne 5. and res.osx_p[1] ne 11. and res.vf_lpwid gt 2.3 $
  ;    and res.gflx_bs gt 0. and bs.sb48 gt 3. and $
  ;    tems.em ne tems[1].em and tems.tmk gt 4.0 and tems.tmk lt 20. and tems.bsub eq 1$
  ;    and tems.em lt gfit and res.osx_p[3] ne 1e-5 and bs.sb612 gt 3.0 and $
  ;    res.osx_p[5] ne 20 and res.osx_p[5] gt 7.2 and abs(res.osx_p[5]-7.16) gt 0.1 and res.osx_p[5] ne 9. and $
  ;    res.osx_p[6] ne 2 and res.osx_p[6] ne 12,ngdnn)
  ;
  ;  ec=eb2ec(res.osx_p[5],res.osx_p[6])
  ;
  ;  nf=n_elements(res.osx_p[3])
  ;  ph4_8=dblarr(nf)
  ;  ph12=dblarr(nf)
  ;  for i=0, nf-1 do begin
  ;    ph4_8[i]=f_vth([4,8],[res[i].osx_p[0],res[i].osx_p[1]*0.086164])
  ;    ph12[i]=f_bpow([12.],res[i].osx_p[3:6])
  ;  endfor
  ;
  ;  rr={tmk:res.osx_p[1],em:res.osx_p[0]*1d49,$
  ;    norm:res.osx_p[3],g1:1.5,eb:res.osx_p[5],g2:res.osx_p[6],$
  ;    ec:ec,ph4_8:ph4_8,ph12:ph12,vol:res.vf_vol,vflx4_8:res.vf_fit.srcflux,vx:res.vf_fit.srcx,vy:res.vf_fit.srcy,$
  ;    gflx_bs:res.gflx_bs,gtmk:tems.tmk,gem:tems.em,$
  ;    eng_th:res.eng_th,eng_nn:res.eng_nn,$
  ;    idgdth:gdth,idgdnn:gdnn,idgsgdth:gsgdth,idgsgdnn:gsgdnn}
  ;
  ;  savegen,file='wee_all.genx',rr

  ;;; *******************************************
  ;;; *******************************************
  ;;; *******************************************
  ;;; *******************************************

  restgen,file='wee_all.genx',rr

  ;  help,rr,/str

  ; The "good" thermal fits
  id=rr.idgdth
  ; The "good" thermal + non-thermal fits
  idnn=rr.idgdnn

  ;  ;---------------------------------
  ;  ;   Plot of the T vs EM
  ;  window,0,xsize=600,ysize=600
  ;  !p.multi=0
  ;  !p.font=0
  ;  !p.charsize=2
  ;  plot,rr.tmk[id],alog10(rr.em[id]),psym=3,$
  ;    yrange=[44,49],xrange=[2,25],xtit='T [MK]',ytit='log!D10!N EM [cm!U-3!N]'


  ;---------------------------------
  ; Plot RHESSI TMK, EM vs Non-thermal model flux at 12keV
  window,1,xsize=1200,ysize=600
  !p.multi=[0,2,1]
  plot,alog10(rr.ph12[idnn]),rr.tmk[idnn],psym=3,xstyle=17,ystyle=17,$
    yrange=[5,24.5],xrange=[-2.3,1.9],ytit='T [MK]',xtit='log!D10!N Non-thermal 12 keV'

  plot,alog10(rr.ph12[idnn]),alog10(rr.em[idnn]),psym=3,xstyle=17,ystyle=17,$
    yrange=[44,48.6],xrange=[-2.3,1.9],$
    ytit='log!D10!N EM [cm!U-3!N]',xtit='log!D10!N Non-thermal 12 keV'

  ;----------------------------------
  ; Plot model thermal
  ;    window,3,xsize=600,ysize=600
  ;    !p.multi=0
  ;    !p.font=0
  ;    !p.charsize=2
  ;    plot,alog10(rr.ph4_8[idnn]),alog10(rr.ph12[idnn]),psym=3,xstyle=17,ystyle=17,$
  ;      yrange=[-2.3,1.9],xrange=[0,3.9],xtit='log!D10!N Thermal 4-8 keV',ytit='log!D10!N Non-thermal 12 keV'

  ;  ;---------------------------------
  ;  ; Plot GOES vs RHESSI TMK and EM
  ;  ; The "good" thermal fits and "good" GOES
  ;  id2=rr.idgsgdth
  ;
  ;  window,4,xsize=1200,ysize=600
  ;  !p.multi=[0,2,1]
  ;  plot,rr.tmk[id2],rr.gtmk[id2],psym=3,$
  ;    yrange=[2,25],xrange=[2,25],xtit='RHESSI T [MK]',ytit='GOES T [MK]'
  ;  oplot,[2,25],[2,25],lines=2
  ;
  ;  plot,alog10(rr.em[id2]),alog10(rr.gem[id2]),psym=3,$
  ;    yrange=[44,49],xrange=[44,49],$
  ;    xtit='RHESSI log!D10!N EM [cm!U-3!N]',ytit='GOES log!D10!N EM [cm!U-3!N]'
  ;  oplot,[44,49],[44,49],lines=2


  stop
end