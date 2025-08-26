print('--- Start Final2DFit ---')

import os, sys

# --- make a directory having the script become a working directory ---
script_path = os.path.abspath(sys.argv[0])
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

from ROOT import gSystem, gROOT
gROOT.SetBatch(True) # no GUI
gSystem.Load('Analysis.so') # load compiled 2D fit codes

from ROOT import Final2DFit

# make an instance
fit = Final2DFit(
    ptLow=6.5, ptHigh=50,
    yLow=0, yHigh=2.4,
    cLow=0, cHigh=180,
    cosLow=-1.0, cosHigh=1.0,
    PR=1,
    PRw=1,
    fEffW=False,
    fAccW=False,
    isPtW=False,
    isTnP=False,
    DATE = '5p36TeV_OO'
)

# set user custom values if you want

# fit.inputFilePath = ''


# ===== prepare and perform the fit =====
fit.init() # up to bulding models

# --- set parameters ---
# fit.initVar('b_Jpsi', 0.20, 0, 0.8) # by experience, < 0.80


fit.run() # fit and save

print('--- Finished Final2DFit ---')


# -----------------------
# --- free user notes ---
# -----------------------