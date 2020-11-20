import numpy as np
import datetime
from flask import Flask, render_template, request, redirect, url_for, session
#import nirc2
#import nires
import json
import importlib

def etc_gui():
    instr = ''
    #get instrument if selected
    if request.form: instr = request.form.get('instrument')
    elif request.args:
        instr = request.args.get('instrument')
    #current instruments with exposure time calculators
    allowed = ['','HIRES','KCWI','NIRC2', 'NIRES']
    #if instrument chosen, load exposure time calculator
    print(instr)
    if instr in allowed and instr != '':
        m = importlib.__import__(instr.lower())
        return m.instr_calc(allowed)
    print('none')
    return render_template('etc_noinstr.html',instrument=instr, allowed=allowed)
