pro help_me
print,"calling sequence: sync,tag=tag,n0=n0,b0=b0,t0=t0,zz=zz,mach1=mach1,mach2=mach2,time_tot=time_tot,time_reacc=time_reacc,lsize=lsize,dist=dist,freq=freq,vshock1=vshock1,vshock2=vshock2,fi1=fi1,fi2=fi2,target_flux=target_flux,dt=dt"
print,"tag= ID of this run"
print,"n0=ambient density in part/cm^3"
print,"b0=mag.field in muG"
print,"t0=T in K"
print,"prepost = pre/post : are the density and Temperature pre or post-shock?"
print,"z=redshift"
print,"mach1= Mach number of first (initial) shock"
print,"mach2= Mach number of second  shock. If Mach2=0, there is no second shock"
print,"time_tot=total time of the run [Gyr]"
print,"time_reacc=epoch of second shock (if any) [Gyr]"
print,"lsize= size of a cube giving the desired volume [kpc]
print,"dist=luminosity distance [Mpc]"
print,"freq=radio frequency [MHz]"
print,"vshock1,vshock2=shock velocity of 1st and 2nd shock [km/s]"
print,"fi1,fi2=accel.efficiency of CR electrons of 1st and 2nd shock"
print,"target_flux= desired flux (of observation) to match at freq [W/Hz]"
print,"dt=timestep [Gyr]"
print,"calling sequence: sync,tag=tag,n0=n0,b0=b0,t0=t0,zz=zz,mach1=mach1,mach2=mach2,time_tot=time_tot,time_reacc=time_reacc,lsize=lsize,dist=dist,freq=freq,vshock1=vshock1,vshock2=vshock2,fi1=fi1,fi2=fi2,target_flux=target_flux,dt=dt"
print,"see source code for default values"
end

;....by F.Vazza, 2020 
pro fp;,tag=tag,n0=n0,b0=b0,t0=t0,prepost=prepost,zz=zz,mach1=mach1,mach2=mach2,time_tot=time_tot,time_reacc=time_reacc,lsize=lsize,dist=dist,freq=freq,vshock1=vshock1,vshock2=vshock2,fi1=fi1,fi2=fi2,target_flux=target_flux,dt=dt
time0=systime(1)

    
 main="/Users/francovazza/Dropbox/fokker_planck_solver/GIAN/"
 tag="test0" ;...ID to attach to any new output
 
 input=main+"input_fam2.dat.txt"   ;input data in column format 
 readcol,input,ti,shocked,machi,bi,tempi,ni,zi

    ns=size(ti)
    ntime=ns(1)   ;...number of timesteps
    time_tot=max(ti)-min(ti)     ;...total time sampled by this run (in Gyr)
    dt=time_tot/float(ntime)     ;...timestep in between outputs 
    tmax=time_tot ;...in Gyr
    epoch=0
 
    lsize=100  ;...[kpc] initial 1D size of each tracers (assuming cubic cell geometry)
    dist=100   ;..Mpc    assumed distance of the source 
    
 
   ;...constants and useful numerical conversions
    vol_trac=lsize ;...tracer size in kpc (1D)
    mp=1.67e-24
    me=0.1093897e-28   ;g
    vc=2.99792458e10   ;cm s-
    kpctocm=3.086e21  ;....i
    kb=1.380658e-16    ;erg k-1
    volume=3*alog10(vol_trac*kpctocm) ;;;vol im cm^3case of
    ntrac=1    ;....in principle, the routine can also run over ntrac>1 different tracers at the same time
    dist*=(1e3*kpctocm)   ;..in cm
   
     
   ;....output parameters for computing radio spectra
    n_epochs=5   
    dt_epochs=time_tot/float(n_epochs)
    epochs=dt_epochs*indgen(n_epochs)
    
    ;...output plot
    set_plot,'ps'
    device,filename=main+"run_test_"+tag+"3.eps",/color,xsize=16,ysize=7,/encapsul,font_size=14
    !p.multi=[0,2,0]
    !p.font=1
    loadct,3
       
  ;....

  ;.....initialisation of tracer spectrum parameters
    g_max=500000  ;...maximum gamma factor
    g_min=10        ;...minumum gamma factor
    dg=50           ;...resolution in gamma factor
    part=2          ;...1=proton  2=electron
    ngamma=uint(((g_max)-(g_min))/double(dg))
   
    a=initialise_spectrum(g_max,g_min,dg,ntrac,nn_trac,gammaval) ;...all spectra are created=0

  
    spectrum_total=fltarr(ngamma)
    spectrum_total(*)=1e-30
    spectrum_final=fltarr(ngamma,n_epochs)
    spectrum_final(*)=1e-30
  
 loadct,13
 
 dtmax=dt
 dti=fltarr(ntime) 
 dti(*)=dt
 time_run=0.

   
   
  for t=0,ntime-1 do begin 
  if t eq 0 then shock=1   ;..we assume the electron is shocked at t=0
   time_run+=dti(t)
   print,"time=",time_run,"Gyr, redshift=",zi(t)

   n0=ni(t)
   t0=tempi(t)
   b0=bi(t)
   zz=zi(t)
   mach1=machi(t)

     ;....pre and post-shock values  
    n0/=float(0.6*mp)   ;...densities are turned in to part/cm^3
    npre=n0*(mach1^2.+3.)/float(4.*mach1^2.)
    tpre=t0*(16*mach1^2.)/float((mach1^2.+3)*(5.*mach1^2.-1.))
    npost=n0
    tpost=t0
         
      ;....quantities related to DSA injection 
    cspre=1e-5*sqrt(1.666*tpre*kb/double(mp*0.6))    ;...preshock
    vshock1=mach1*cspre   
    mlow=1.3 ;...minimum mach number to inject shocks
    mach_inj=mach1
    vpre=vshock1 ;cm/s preshock velocity
    f=(4.*mach1*mach1)/float(mach1*mach1+3.)
    vpost=vpre/float(f)
    ;....shock kinetic energy flux
    Ecr_inj=10+alog10(0.5*npre*vshock1^3.0/float(vpost))+alog10(0.6*mp)+volume  ;...the factor 1e10 in front is to compensate for vshock1^3/vpost, in which both are in km/s instead of cm/s
   
    v=vshock1

     if t eq 0 then plot,gammaval(*),30+alog10(nn_trac(*)),/xlog,xrange=[g_min,g_max],yrange=[50,63],/ystyle,thick=1,title="electron spectra",xtitle='!9g!6',ytitle='log[!9g!6 N(!9g!6)]',/xstyle,/nodata

    spectrum_plot=fltarr(ngamma)
    spectrum_plot(*)=1e-30
   
     norm=30  ;...log10(erg/s) - renormalisation factor to avoid NaN
     
    for p=0,ntrac-1 do begin  ;....loop over ntrac tracers 
     
      nth=n0 
      shock=0
      
      if t eq 0 then begin
         ecr=ecr_inj ;...log10(ecr) erg/s)
         shock=1   ;...shock=1 just inject new electrons via DSA
         m=mach1
         ecr=ecr-norm
         vshock=vshock1
         print,'injection by a M=',m,'shock'
      endif
      
           
      if t ge 1 and mach1 gt mlow then begin
         shock=2   ;...shock=2 includes both DSA acceleration and shock re-acceleration on fossil electrons 
         m=mach1
         ecr=ecr-norm
         vshock=vshock1 
         print,'re-acceleration by a M=',m,'shock'
      endif

      t2=dti(t)*1e9*3.154e7 ;..in seconds
      t1=0
     
      zzc=zz
       nn=nn_trac(*,p)   ;....for consistency with other routines running on many tracers
       delta_t=t2-t1
    
      ca=evolve_spectrum(zzc,v,t2,t1,nth,m,b0,tpost,ecr,nn,shock,delta_t,t,gend,volume,g_max,g_min,dg,ngamma,part)
        
      spectrum_plot=nn;(io)+=nn(io)
      spectrum_total+=nn
      
      oplot,gammaval,alog10(gammaval)+alog10(nn)+norm,linestyle=0,col=10+10*t,thick=1
     nn_trac(*,p)=nn

  endfor


endfor

   nn=spectrum_total
   

 ;.........synchrotron radio emission (only for the final step here)
  
  loadct,13,/silent
  nfr=4
  freqf=[140,610,1400,5000]*1e6  ;...Hz
  pradio=fltarr(nfr,2)

   for ff=0,nfr-1 do begin   ;...loop over observing frequencies 
     ptot=compute_sync_full(gammaval,nn,b0,freqf(ff),dg)  ;...in erg/s/Hz
     ptot=norm+23+alog10(ptot)-alog10(4*!pi)-2*alog10(dist)   ;...from erg/s/Hz to Jy 
      pradio(ff)=1e3*10.^ptot

   endfor 
   plot,1e-6*freqf(0:nfr-1),pradio(0:nfr-1),ytitle='mJy',xtitle='freq[MHz]',title='radio spectra',/xlog,/ylog,/nodata,/xstyle,/ystyle
   loadct,0


   oplot,1e-6*freqf(0:nfr-1),pradio(0:nfr-1),noclip=0,linestyle=0,col=90,thick=3
   print,pradio(1:nfr-1,0)
   print,1e-6*freqf(1:nfr-1)
 
  
device,/close
device,/encapsul

end


;.....computation of the synchrotron radio emission

function compute_sync_full,gammaval,nn,b0,freq,dg
  h=6.626d-27
  vc=2.99d10
  kb=1.3806d-16
  me=9.109d-28
  re=2.81d-13  ;...electron radius
  
  keverg=1.602e-9
  e=4.8d-10 ;....esu
  norm=1e30

  ng=size(gammaval)
  dp=dblarr(ng(1))
  psync=dblarr(ng(1))
  g2=ng(1)-2
  g1=4
  
  c1=3*e*0.5*sqrt(2)/double(4*!pi*me*vc)
  c2=1.732*e^3./double(me*vc^2.)
  bb=b0*1e-6
  ;.,,,,compute syncrothron emissivity
  psync=0
  dth=0.1
  to=0.5*!pi*dth*(indgen(10))  ;..angles
  
    integ=0
    for gg=g1,g2 do begin  ;...gamma factors
    freqc=c1*bb*(gammaval(gg))^2.
     vv=freq/double(freqc)
       for ang=0,9 do begin  ;...angles
        thi=to(ang)
        gu=float(gammaval(gg))
        sinthi=sin(thi)
          totf=0
      
        for y=0L,12 do begin  ;....K53 modified bessel function
          bbes=beselk((1+y*0.25)*vv,1.666)   
          totf=totf+bbes*0.25
        endfor

        ffn=vv*totf
        integ=integ+nn(gg)*ffn*sinthi^2*dth
      endfor
    psync+=integ*(bb*c2)
    endfor


  return,total(psync)

end


;...routine for the Fokker Planck evolution of particle spectra 
;...we assumed injection paramters from shocks as in Pinzke et al. 2013 formalism (derived from Kang & Ryu 2011)
;...radiative cooling according to Jaffe & Perola 1973

function evolve_spectrum,zzc,v,t2,t1,nth,m,b0,tpost,ecr,nn,shock,tstep,t,gend,volume,g_max,g_min,dg,ngamma,part
  gend=ngamma

  ;nn=fltarr(ngamma)
  ;nn=nn2
  ;zzc is redshift
  ;..v in km/s
  ;..t2, t1 in s
  ;...nth in gr/cm^3
  ;...m mach
  ;...in muG
  ;...ecr is log(Ecr) in erg/s
  ;...nn is spectrum (size ngamma)
  ;...tstep in s
  ;...t is time step
  ;...gend last nonzero gamma factor
  ;...assumed shock velocity

  ;define constants
  hp=6.6260755e-27; erg s
  pi=acos(-1)
  pmass=1.672e-24 ;g
  emass=0.1093897e-28   ;g
  eradius=2.8179e-13
  vc=2.99792458e10   ;cm s-
  kb=1.380658e-16    ;erg k-1
  h0=71*3.24e-20   ;..km/s/Mpc
  cmtoMpc=3.086e24
  sy=31556926 ;
  erest=0.510988946 ;MeV electron rest energy
  evtoerg=1.60218e-12
  ;.....constants for particle cooling
  b1=1.37e-20 ;...1/s 
  b2=1.9e-9  ;...1/s
  b3=7.2e-12 ;...1/s
  cou=nth*1.2e-12 ;..1/s

  gam=dblarr(ngamma) ;...array of Lorentz factors 
  gam=double(dg)*indgen(ngamma)+(g_min)

  vpre=1e5*v ;cm/s preshock velocity
  f=(4*m*m)/float(m*m+3)   ;...compression factor 
  vpost=vpre/float(f)
  delta=2*(m*m+1)/float(m*m-1)
  lcorr=10 ;...assumed correlation scale of B-field - only to compute gamma_mcut 
  erest=emass*vc^2.   ;...electron rest energy 
  gammag=dblarr(ngamma)
  gammag=gam

  n_inj=dblarr(ngamma)  ;...spectrum of injected energy
  
  ;...compute non-stationary solutions
  dt=tstep
 
  idt=1/double(dt)
  idg=1/double(dg)
  icou=1/double(cou)
  inth=1/double(nth)

  tend=uint((t2-t1)/double(tstep))


  q_inj=dblarr(ngamma)   ;...array of new injected electrons by the shock (if any) 
  q_inj(*)=0.

  for ts=1,1 do begin

    n_re=dblarr(ngamma)
    n_re(*)=0.
    ta=ts

    if ts eq 1 and shock ge 1  then begin
      ic_lose=dblarr(ngamma)

      ic_lose=b1*((gam*gam))*(1+zzc)^4

      f=(4*m*m)/double(m*m+3)
      delta=2*(m*m+1)/double(m*m-1)

   ;...here we compute the maximum gamma that can be reached at the shock -> has a very neligible effect on anything we are doing here as gamma_cut >=gamma_max, typically
   ;...diffusion coefficient - assuming Bohm diffusion
      diff=dblarr(ngamma)
      diff=2.3e29*(((gam)-1)*(0.511e-3))^(0.3333)*b0^(-0.3333)*(lcorr*0.05)^(0.66666) ; in cm^2/s
   ;....shock acceleation timescale 
      sh_gain=dblarr(ngamma)
      sh_gain=(gam)*1/double(f)*(vpre)^2.*((f-1)/double(f+1))*1/double(3*diff)

      ;...thermal leakage MODEL after Kang & Ryu 2011 and Pinzke+13
      pth=sqrt(2*tpost*kb*pmass)   ;...thermal momentum of protons
      pthe=sqrt(2*tpost*kb*emass)  ;...thermal momentum of electrons 
      eB=0.23
      xinj=(1.17*pmass*vpost/float(pth))*(1.+1.07/float(eB))*(m/float(3.))^(0.1)  ;...Eq.7 in Pizke+13 -
      
      pinje=pthe*xinj              ;....injection momentum for electrons 
      gamma_inj=sqrt(1.+pinje^2./double(emass*vc)^2.)
    
      ;...find gamma max equating losses and acceleration
     
      gi=where(ic_lose gt sh_gain)

      g=gi(0)

      gamo=gam(g) ;ngamma-1) ;g)

      q_inj(0:ngamma-1)=0 
      n_re(0:ngamma-1)=0

      gam_c=gamo
      
      gam=double(dg)*indgen(ngamma)+(g_min)
      gammag=gam
      n_inj=dblarr(ngamma)  ;...spectrum of injected energy
      n_inj(*)=-60
      
      ;............CR-proton efficiency from Kang & Jones 2007 - can be replaced with newer versions
      m4=m^(-4)
      if m lt 2  then eta=1.96e-3*(m*m-1.)
      if m ge 2 then  eta=5.46*(m4)-9.78*(m-1.)*(m4)+4.17*(m-1.)^2.*(m4)-0.334*(m-1.)^3.*(m4)+0.57*(m-1.)^4.*(m4);
      q=4.*m^2/float(m^2.-1)
      Kpe=(pmass/double(emass))^(0.5*(q-3.))   ;...proton to electron ratio From Kang+11, eq.10
      Kep=1./float(Kpe)                        ;...electron to proton ratio to be used below
     
     
          inte=emass*vc^2*(gamma_inj^(-delta+2))/float(delta-2)
          Ke=Kep*eta*(10.^(Ecr-alog10(inte)))
          
;      
      n_inj=alog10(Ke)-delta*alog10(gammag);+alog10((1.-gammag/double(gam_c ))^(delta-2))  ;...final injection spectrum of electron

      n_inj(where(finite(n_inj) ne 1))=-60

      q_inj=(10.^(n_inj))
  
    endif else begin


      q_inj(*)=0 
      n_re(*)=0  
     

    endelse

    if shock eq 2 then begin
      print,"reacceleration"
    ;....re-acceleration by shocks as in Kang & Ryu 2011 
      n_re(0)=nn(0)
      gmin_re=1

      for gag=gmin_re,ngamma-2 do begin
        cutoff=(1.-(gammag(gag)/double(gammag(ngamma-1))))
        n_re(gag)=(delta+2.)*(gammag(gag))^(-delta)*(total(nn(gmin_re:gag)*dg*(gammag(gmin_re:gag))^(delta-1.))); as in Eq.6 of Kang & Ryu 2011
        if finite(n_re(gag)) ne 1  then print,'problem',n_re(gag),gammag(gag),nn(gag)
      endfor

      nn(*)=n_re(*)

    endif

    gend=ngamma-1



  ;...loop over all gammas to evolve the spectrum 
    for gga=1,gend do begin

      gg=uint(gend-gga)

      ga=gammag(gg)

      iga=1/double(ga)
      g1=ga
      g2=ga+dg
  ;....cooling constant 
      aa1=cou*(1+(alog(g1/double(nth)))/double(75)) ;...coulomb loss coefficient for the two extreme of the en.bin
      aa2=cou*(1+(alog(g2/double(nth)))/double(75)) 
      bb=b2*(0.666*b0*b0*1e-12+b3*(1+zzc)^4.)       ;....Jaffe & Perola 73, isotropic pitch angle scattering, also including Inverse Compton Losses
      if part eq 1 then bb=0
      if gg ge ngamma-1 then gg = ngamma-1
      if gg lt 0 then gg = 0
      gg2=gg+1
      if gg2 ge ngamma-1 then gg2=ngamma-1

      Nn(gg)=(1/double(idt+idg*(aa1+bb*g1*g1)))*(idt*q_inj(gg)+nn(gg)*idt+nn(gg2)*idg*(aa2+bb*g2*g2))


    endfor


  endfor

end

function initialise_spectrum,g_max,g_min,dg,ntrac,nn_trac,gammaval
 
  cmtokpc=3e21

  ngamma=uint(((g_max)-(g_min))/float(dg))

  gam=dblarr(ngamma) ;...logarithm of gamma factor
  gam=dg*indgen(ngamma,/long)+(g_min)
  gammaval=dblarr(ngamma)
  gammaval=gam


  nn_trac=dblarr(ngamma,ntrac)
  nn_trac(*,*)=0.
  
  end
  
