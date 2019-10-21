import numpy as np
import datetime
from flask import Flask, render_template, request, redirect, url_for, session
from etc_kcwi import do_calc as kcwi_docalc
import json
import html

def instr_calc(allowed):
    instr = 'KCWI'
    slicer = request.form.get("slicer",'Small')
    grating = request.form.get("grating",'BL')
    gratwave = request.form.get("gratwave",4500.)
    seeing = request.form.get("seeing",0.75)
    exptime = request.form.get("exptime",3600.)
    ccdbin = request.form.get("ccdbin",'1x1')
    magAB = request.form.get("magAB",23.)
    flux = request.form.get("flux",None)
    nframes = request.form.get("nframes",None)
    ccdspeed = request.form.get("ccdspeed",None)
    spat_binx = request.form.get("spat_binx",None)
    spat_biny = request.form.get("spat_biny",None)
    print(spat_binx)
    spat_bin = [spat_binx,spat_biny]
    spec_bin = request.form.get("spec_bin",None)
    surf_bright = request.form.get("surf_bright",None)
    emline_w = request.form.get("emline_w",None)
    calculate = request.form.get("calculate")
    print('FLUX',flux)
    if calculate:
        gratwave = float(gratwave)
        seeing = float(seeing)
        exptime = float(exptime)
        if magAB != 'None':
          magAB = float(magAB)
        print(nframes,type(nframes))
        if nframes != 'None':
          nframes = float(nframes)
        try:
          spat_binx = float(spat_binx)
          spat_biny = float(spat_biny)
          print('floating spatbin')
        except:
         	pass
        spat_bin = [spat_binx,spat_biny]
        print(spat_bin,'a')
        # if spat_binx == '':
        # 	spat_binx = 'None'
        # if spat_biny == '':
        # 	spat_biny = 'None'
        if flux != 'None':
          flux = float(flux)
        if spec_bin != 'None':
          spec_bin = float(spec_bin)
        if surf_bright != 'None':
          surf_bright = float(surf_bright)
        if emline_w != 'None' and emline_w != None:
          emline_w = float(emline_w)
        ###### FILL WITH ESI DOCALC ######
        script_snr,div_snr,script_obj,div_obj,script_sky,div_sky,script_rn,div_rn,script_os,div_os,script_flux,div_flux = kcwi_docalc(slicer,grating,gratwave,seeing,exptime,ccdbin,
                          magAB,flux,nframes,ccdspeed,spat_bin,spec_bin,
                          surf_bright,emline_w)
        # div_snr = html.escape(div_snr,quote=True)

        return render_template('etc_kcwi.html',instrument=instr,slicer=slicer,grating=grating,
                               gratwave=gratwave,seeing=seeing,exptime=exptime,ccdbin=ccdbin,
                               magAB=magAB,flux=flux,nframes=nframes,ccdspeed=ccdspeed,
                               spat_bin=spat_bin,spat_binx=spat_binx,spat_biny=spat_biny,spec_bin=spec_bin,surf_bright=surf_bright,
                               emline_w=emline_w,calculate=calculate,script_snr=script_snr,
                               div_snr=div_snr,script_obj=script_obj,div_obj=div_obj,
                               script_sky=script_sky,div_sky=div_sky,script_rn=script_rn,
                               div_rn=div_rn,script_os=script_os,div_os=div_os,
                               script_flux=script_flux,div_flux=div_flux,
                               allowed=allowed)

    return render_template('etc_kcwi.html',instrument=instr,slicer=slicer,grating=grating,
                           gratwave=gratwave,seeing=seeing,exptime=exptime,ccdbin=ccdbin,
                           magAB=magAB,flux=flux,nframes=nframes,ccdspeed=ccdspeed,
                           spat_binx=spat_binx,spat_biny=spat_biny,spec_bin=spec_bin,surf_bright=surf_bright,
                           emline_w=emline_w,calculate=calculate,allowed=allowed)