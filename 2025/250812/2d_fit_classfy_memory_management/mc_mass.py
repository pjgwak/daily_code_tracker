print('--- Start McMassFit ---')

import os, sys

# --- make a directory having the script become a working directory ---
script_path = os.path.abspath(sys.argv[0])
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

from ROOT import gSystem, gROOT
gROOT.SetBatch(True) # no GUI
gSystem.Load('Analysis.so') # load compiled 2D fit codes

from ROOT import McMassFit

# make an instance
fit = McMassFit(
    ptLow=6.5, ptHigh=50,
    yLow=0, yHigh=2.4,
    cLow=0, cHigh=180,
    cosLow=-1.0, cosHigh=1.0,
    PR=2,
    PRw=1,
    fEffW=False,
    fAccW=False,
    isPtW=False,
    isTnP=False,
    DATE = '5p36TeV_OO'
)

# === set user custom values if you want ===
# fit.inputFilePath = ''

# ===== prepare and perform the fit =====
fit.init() # up to bulding models

# --- set parameters ---
fit.fitLimit = 3.23 # max of mass for fit
fit.initVar('N_Jpsi', 2000000, 1000000, 5000000)
fit.initVar('mean', 3.096, 3.086, 3.106)
fit.initVar('x_A', 1.1,1,3)
fit.initVar('sigma_1_A', 0.01,0.001,0.1)
fit.initVar('alpha_1_A', 1.5,0.8,5)
fit.initVar('n_1_A', 1.5,0.8,5)
fit.initVar('f', 0.6,0.05,0.95)


fit.run() # fit and save

print('--- Finished McMassFit ---')


# -----------------------
# --- free user notes ---
# -----------------------