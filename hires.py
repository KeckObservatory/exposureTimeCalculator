import numpy as np
import datetime
from flask import Flask, render_template, request, redirect, url_for, session
from etc_nirc2 import do_calc as nirc2_docalc
import json

def instr_calc(allowed):
	instr = 'HIRES'
