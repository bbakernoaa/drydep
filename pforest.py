#!/usr/bin/env python
#
# pforest.py
#
import os
import sys
import matplotlib.pylab as plt
from numpy import array, zeros, linspace
import seaborn as sns

# measurements
def get_data_Vong2010():
    vd = array([0.18, 0.29, 0.46, 0.63, 0.61, 0.38])
    dp = array([0.26, 0.30, 0.33, 0.40, 0.49, 0.78])

    vd_err = zeros(6)
    dp_err = zeros(6)

    return vd, vd_err, dp, dp_err

def get_data_Gordon2011():
    vd = array([0.40, 0.58, 0.39, 0.26, 0.16, 0.12, 0.11, 0.105, 0.10, 0.095, 0.09, 0.09, 0.08, 0.58, 0.46])
    dp = array([0.017, 0.024, 0.029, 0.035, 0.043, 0.052, 0.062, 0.076, 0.091, 0.11, 0.133, 0.16, 0.194, 0.282, 0.34])

    vd_err = zeros(15)
    dp_err = zeros(15)

    return vd, vd_err, dp, dp_err

def get_data_Mammarella2011():
    vd = array([0.4, 0.34, 0.29, 0.26, 0.21, 0.17, 0.14, 0.12, 0.08, 0.07, 0.08, 0.084, 0.08, 0.11, 0.17, 0.20])
    dp = array([0.02, 0.023, 0.028, 0.032, 0.039, 0.047, 0.054, 0.065, 0.075, 0.091, 0.11, 0.13, 0.17, 0.19, 0.21, 0.24])

    vd_err = zeros(16)
    dp_err = zeros(16)

    return vd, vd_err, dp, dp_err

def get_data_Lavi2013():
    vd = array([0.38, 0.48, 0.57, 0.70, 0.68, 0.64, 0.52, 0.51])
    dp = array([0.27, 0.29, 0.33, 0.38, 0.43, 0.48, 0.54, 0.57])

    vd_err = zeros(8)
    dp_err = zeros(8)

    return vd, vd_err, dp, dp_err

def get_data_Pryor2006():
    vd = array([0.47,  0.31,  0.30,  0.26,  0.15])
    dp = array([0.025, 0.035, 0.045, 0.055, 0.065])

    vd_err = array([ [0.25, 0.20, 0.18, 0.16, 0.03] , [0.29, 0.24, 0.25, 0.27, 0.38] ])
    dp_err = zeros(5)
    
    return vd, vd_err, dp, dp_err

def get_data_Pryor2009():
    vd = array([0.320, 0.220, 0.300, 0.200, 0.196, 0.120, 0.200, 0.220, 0.160, 0.100, 0.100, 0.080, 0.072, 0.060, 0.060, 0.060])
    dp = array([0.011, 0.012, 0.015, 0.017, 0.019, 0.022, 0.025, 0.029, 0.033, 0.039, 0.043, 0.051, 0.060, 0.070, 0.080, 0.092])

    vd_err = zeros(16)
    dp_err = zeros(16)
 
    return vd, vd_err, dp, dp_err

def get_data_Pryor2007():
    vd = array([0.420, 0.300, 0.210, 0.190, 0.170, 0.150, 0.140, 0.130])
    dp = array([0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085])

    vd_err = array([ [0.23, 0.17, 0.12, 0.11, 0.09, 0.08, 0.07, 0.06] , [0.38, 0.25, 0.21, 0.19, 0.16, 0.17, 0.16, 0.16] ])
    dp_err = zeros(8)
    
    return vd, vd_err, dp, dp_err

def get_data_Hofken1982():
    vd = array([1.80, 1.10, 1.00, 0.90, 0.70, 1.00, 1.00, 1.30])
    dp = array([0.26, 0.63, 0.70, 0.78, 1.10, 1.50, 2.20, 2.40])

    vd_err = zeros(8)
    dp_err = zeros(8)
 
    return vd, vd_err, dp, dp_err

def get_data_Gronholm2007():
    vd = array([1.770, 1.490, 1.180, 1.230, 0.900, 0.940, 0.740, 0.620, 0.660, 0.990, 0.780, 0.690, 0.250])
    dp = array([0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.100, 0.100, 0.150])

    vd_err = array([ [1.76, 1.12, 1.15, 0.92, 0.57, 0.67, 0.40, 0.47, 0.60, 0.87, 0.60, 0.51, 0.24] ,
                     [2.64, 0.97, 1.03, 0.85, 0.90, 0.59, 0.56, 0.49, 0.48, 0.73, 0.55, 0.37, 0.59] ])
    dp_err = zeros(13)
 
    return vd, vd_err, dp, dp_err

def get_data_Beswick1991():
    vd = array([1.15, 2.51, 1.93, 1.73, 1.83, 1.83, 1.31, 3.13, 3.65, 3.81, 3.23, 5.99, 9.40])
    dp = array([3.00, 5.00, 7.00, 9.00, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0])

    vd_err = array([ [0.86, 0.91, 0.67, 0.52, 0.54, 0.44, 0.48, 0.74, 0.86, 1.58, 1.84, 4.24, 8.54] ,
                     [1.24, 0.81, 0.73, 0.52, 0.54, 0.44, 0.48, 0.74, 0.86, 1.32, 1.76, 5.27, 5.80] ])
    dp_err = zeros(13)
 
    return vd, vd_err, dp, dp_err

def get_data_Gronholm2009():
    vd = array([0.470, 0.160, 0.060, 0.050, 0.100, 0.085])
    dp = array([0.010, 0.020, 0.030, 0.040, 0.050, 0.060])

    vd_err = array([ [0.24, 0.10, 0.03, 0.02, 0.015, 0.06] , [0.25, 0.18, 0.09, 0.08, 0.05, 0.04] ])
    dp_err = zeros(6)
 
    return vd, vd_err, dp, dp_err

def get_data_Gallagher1997():
    vd = array([0.150, 0.030, 0.130, 0.180, 0.120, 0.180, 0.300, 0.850, 0.850, 1.800, 2.300, 0.600, 7.000, 1.100, 1.200, 1.700, 0.700, 11.00, 1.300, 0.800])
    dp = array([0.100, 0.110, 0.130, 0.140, 0.160, 0.180, 0.210, 0.220, 0.250, 0.290, 0.320, 0.380, 0.420, 0.480, 0.530, 0.630, 0.830, 0.900, 1.200, 2.100])

    vd_err = zeros(20)
    dp_err = zeros(20)
 
    return vd, vd_err, dp, dp_err

def get_data_Waraghai1989():
    vd = array([1.500, 1.100, 0.860, 1.100, 0.93, 1.64])
    dp = array([0.750, 1.500, 3.500, 6.500, 9.00, 10.00])

    vd_err = array([0.34, 0.36, 0.39, 0.41, 0.70, 0.59])
    dp_err = array([0.25, 0.50, 1.50, 1.50, 1.00, 0.00])

    return vd, vd_err, dp, dp_err

def get_data_Grosch1988():
    vd = array([1.85, 1.30, 1.70])
    dp = array([0.4, 1.1, 4.4])

    vd_err = array([0.75, 0.40, 0.5])
    dp_err = zeros(3)

    return vd, vd_err, dp, dp_err

def get_data_Lorenz1989():
    vd = array([0.34, 0.78, 0.92])
    dp = array([0.75, 1.5, 3.5])

    vd_err = array([0.35, 0.64, 0.76])
    dp_err = array([0.25, 0.5, 1.5])

    return vd, vd_err, dp, dp_err

def get_data_Gaman2004():
    vd = array([0.05])
    dp = array([0.43])

    vd_err = array([0.06])
    dp_err = zeros(1)

    return vd, vd_err, dp, dp_err

def get_data_Zhang2014():
    vd = array([0.24, 0.61, 0.78, 0.91, 2.00, 4.13, 4.37, 5.99, 10.7, 0.92, 1.58, 4.30, 3.60, 3.50, 5.50, 3.80, 4.30, 13.0])
    dp = array([1.00, 2.25, 4.00, 7.50, 12.5, 17.5, 22.5, 27.5, 40.0, 1.00, 2.25, 4.00, 7.50, 12.5, 17.5, 22.5, 27.5, 40.0])

    vd_err = array([1.93, 1.61, 0.00, 1.35, 0.00, 4.50, 0.00, 0.00, 1.40, 3.90, 4.10, 3.10, 4.30, 0.00, 13.0, 11.5, 7.5, 11.4])
    dp_err = array([0.50, 0.75, 1.00, 2.50, 2.50, 2.50, 2.50, 2.50, 10.0, 0.50, 0.75, 1.00, 2.50, 2.50, 2.50, 2.50, 2.50, 10.0])

    return vd, vd_err, dp, dp_err

def plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl):
    if (showebars == True):
       plt.errorbar(dp, vd, yerr=vd_err, xerr=dp_err, fmt=sym, color=clr,
                    ms=msz, label=lbl, elinewidth=1.0)
    else:
       plt.plot(dp, vd, marker=sym, color=clr, ms=msz, ls='None', label=lbl)

    return

# main
def main(argv=None):
    import slinn
    import cmaq
    import zhang01
    import pz2010
    import proto
    import wesely

    if argv is None:
        argv = sys.argv

    # enforce proper number of arguments passed
    if len(argv) != 3:
        print "usage: %s USTARVALUE OUTTYPE" % os.path.basename(sys.argv[0])
        return 2

    # extract argument
    ustar = float(argv[1])
    otype = argv[2]

    # define data dictionary
    get_data = { "Pryor2006"      : get_data_Pryor2006,
                 "Pryor2009"      : get_data_Pryor2009,
                 "Pryor2007"      : get_data_Pryor2007,
                 "Hofken1982"     : get_data_Hofken1982,
                 "Gronholm2007"   : get_data_Gronholm2007,
                 "Gronholm2009"   : get_data_Gronholm2009,
                 "Gallagher1997"  : get_data_Gallagher1997,
                 "Waraghai1989"   : get_data_Waraghai1989,
                 "Grosch1988"     : get_data_Grosch1988,
                 "Lorenz1989"     : get_data_Lorenz1989,
                 "Beswick1991"    : get_data_Beswick1991,
                 "Zhang2014"      : get_data_Zhang2014,
                 "Vong2010"       : get_data_Vong2010,
                 "Gordon2011"     : get_data_Gordon2011,
                 "Mammarella2011" : get_data_Mammarella2011,
                 "Lavi2013"       : get_data_Lavi2013,
                 "Gaman2004"      : get_data_Gaman2004 } 

    # extra settings from seaborn
    sns.set_context("talk")
    sns.set_style("ticks")

    # formatting parameters
    tfsize = 18     # plot title font size
    tyloc  = 1.03   # title y location
    lfsize = 9      # legend font size
    yfsize = 18     # y-axis title font size
    xfsize = 18     # x-axis title font size
    tlmaj  = 6      # major tick length
    tlmin  = 4      # minor tick length
    tlbsize = 16    # tick label size
    msize  = 7      # marker size

    f, ax = plt.subplots()

#   dp = array([0.01, 0.1, 1., 10.0])
    dp = linspace(0.01, 100., 100000)
#   dp = linspace(0.01, 100., 5)
    dp = dp*1.0e-06

    # parameters
    rhop = 1000.0
    t    = 300.0
    luc  = 3

    # track number of model line plotted
    nmodels=0

    # Slinn
    uref = 2.0
    zol = 0.0
    zrefm = 10.0
    vd, eb, ein, eim, rfac, vg, rc = slinn.model(dp, rhop, t, ustar, zol, uref, zrefm, luc)
    plt.loglog(dp*1.0e+06, vd*100., '-', color='violet', label='Slinn')
    nmodels+=1

    # CMAQ
    wstar = 0.0
    zol   = 0.0
    zrefm = 10.0 
    vd, eb, ein, eim, vg, rc = cmaq.model(dp, rhop, t, ustar, wstar, zol, zrefm, luc)
    plt.loglog(dp*1.0e+06, vd*100., '-', color='orange', label='CMAQ')
    nmodels+=1

    # Zhang et al. (2001)
    zol   = 0.0
    zrefm = 10.0 
    vd, eb, ein, eim, rfac, vg, rc = zhang01.model(dp, rhop, t, ustar, zol, zrefm, luc)
    plt.loglog(dp*1.0e+06, vd*100., '-', color='black', label='Zhang01')
    nmodels+=1

    # Petroff and Zhang (2010)
    zol   = 0.0
    zrefm = 10.0 
    vd, eb, ebin, eim, eit, egb, egt, qfac, vg, rc = pz2010.model(dp, rhop, t, ustar, zol, zrefm, luc)
    plt.loglog(dp*1.0e+06, vd*100., '-', color='purple', label='PZ2010')
    nmodels+=1

    # proto
    zol = 0.0
    zrefm = 10.0 
    uref = 1.0
    vd, eb, ein, eim, euk, rfac, vg, rc = proto.model(dp, rhop, t, ustar, zol, uref, zrefm, luc)
    plt.loglog(dp*1.0e+06, vd*100., '-', color='green', label='proto')
    nmodels+=1

    # Wesely
#   c = 0.63
#   zol = 0.0
#   vd = wesely.model(c, dp*1.0e+06, ustar, zol)
#   plt.loglog(dp*1.0e+06, vd, '-', color='crimson', label='Wesely')
#   nmodels+=1

    # plot data
    vd, vd_err, dp, dp_err = get_data.get("Beswick1991", lambda: None)()
    sym = 'o'
    clr = 'peru'
    msz = msize
    lbl = 'Beswick et al. (1991)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Pryor2007", lambda: None)()
    sym = 's'
    clr = 'green'
    msz = msize
    lbl = 'Pryor et al. (2007)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Hofken1982", lambda: None)()
    sym = '^'
    clr = 'red'
    msz = msize
    lbl = 'Hofken et al. (1982)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Gronholm2007", lambda: None)()
    sym = 'h'
    clr = 'orange'
    msz = msize
    lbl = 'Gronholm et al. (2007)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Gronholm2009", lambda: None)()
    sym = 'v'
    clr = 'black'
    msz = msize
    lbl = 'Gronholm et al. (2009)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Gallagher1997", lambda: None)()
    sym = 'o'
    clr = 'magenta'
    msz = msize
    lbl = 'Gallagher et al. (1997)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Waraghai1989", lambda: None)()
    sym = 'D'
    clr = 'blue'
    msz = msize-1.0
    lbl = 'Waraghai/Gravenhorst (1989)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Grosch1988", lambda: None)()
    sym = 'v'
    clr = 'brown'
    msz = msize
    lbl = 'Grosch/Schmitt (1988)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Lorenz1989", lambda: None)()
    sym = 's'
    clr = 'cyan'
    msz = msize
    lbl = 'Lorenz/Murphy (1989)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Gaman2004", lambda: None)()
    sym = 'h'
    clr = 'greenyellow'
    msz = msize
    lbl = 'Gaman et al. (2004)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Pryor2006", lambda: None)()
    sym = '^'
    clr = '0.50'
    msz = msize
    lbl = 'Pryor et al. (2006)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)


    vd, vd_err, dp, dp_err = get_data.get("Pryor2009", lambda: None)()
    sym = 'o'
    clr = '0.75'
    msz = msize
    lbl = 'Pryor et al. (2009)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Zhang2014", lambda: None)()
    sym = 'D'
    clr = 'goldenrod'
    msz = msize-1.0
    lbl = 'Zhang et al. (2014)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Vong2010", lambda: None)()
    sym = '*'
    clr = 'darkgreen'
    msz = msize+1.0
    lbl = 'Vong et al. (2010)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Gordon2011", lambda: None)()
    sym = 'v'
    clr = 'mediumpurple'
    msz = msize
    lbl = 'Gordon et al. (2011)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Mammarella2011", lambda: None)()
    sym = 'o'
    clr = 'saddlebrown'
    msz = msize
    lbl = 'Mammarella et al. (2011)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    vd, vd_err, dp, dp_err = get_data.get("Lavi2013", lambda: None)()
    sym = 'h'
    clr = 'sage'
    msz = msize
    lbl = 'Lavi et al. (2013)'
    showebars=True
    plotdata(plt, showebars, dp, vd, vd_err, dp_err, sym, clr, msz, lbl)

    # make everything tidy
    plt.xticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])
    plt.yticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])

    ax.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=0.5)
    ax.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=0.25)

    # hack from Stack Overflow to only plot symbols (not errorbars) on legend
    handles, labels = ax.get_legend_handles_labels()
    hnew = handles[nmodels:]
    i=nmodels-1
    for h in hnew:
       i+=1
       handles[i]=h[0]
    plt.legend(handles, labels, loc=2, fontsize=lfsize, ncol=2)
#   plt.legend(handles, labels, loc=4, fontsize=lfsize, bbox_to_anchor=(1.10, 0.0))

    # grid and ticks
    plt.ylim(0.001, 1000.0)
    plt.grid(b=True, which="major", color="gray", linewidth=1.0, alpha=0.5)
    plt.grid(b=True, which="minor", color="gray", linewidth=0.5, alpha=0.25)

    ax.tick_params(which="both", direction="out")
    ax.tick_params(which="major", length=tlmaj)
    ax.tick_params(which="minor", length=tlmin)
    ax.tick_params(axis="both", which="major", pad=2)
    ax.tick_params(which="both", labelsize=tlbsize)

    # labels and title
    plt.xlabel('D$_p$ ($\mu$m)', fontsize=xfsize, labelpad=10)
    plt.ylabel('v$_d$ (cm/s)', fontsize=yfsize, labelpad=-10)
    plt.title('Forest', fontsize=tfsize, y=tyloc)

    if (otype == "pdf"):
       plt.savefig('./pforestfig.pdf')
    elif (otype == "png"):
       plt.savefig('./pforestfig.png')
    else:
       plt.show()

    # all done!
    return 0

if __name__ == "__main__":
    sys.exit(main())
