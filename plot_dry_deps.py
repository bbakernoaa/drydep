from matplotlib.pylab import *
from numpy import linspace, array

# NOTE THAT THIS CODE HAS DATA ALREADY LOADED FROM THIS FUNCTION


def ice_data():
    from numpy import array
    # data from [ duan1988,duan1988,nilsson2001,nilsson2001,ibrahim1983,ibrahim1983,ibrahim1983 ]
    dp = array([.225, .75, .023, .049, .69, .78, 7.2])
    dperror = array([.075, .021, .001, .017, .001, .001, 2.4])
    vd = array([.034, .021, .15, .055, .035, .05, .14])
    vderror = array([.0001, .00001, .575, .251, .01, .03, .05])

    return dp, dperror, vd, vderror


def water_data():
    from numpy import array, zeros

    vdd = array(
        [.06, .04, .023, .023, .02, .018, .012, .02, .013, .018, .011, .018, .015, .02, 10.1, 3.0, 1.6, .14, .022, .018,
         .007, .0043, .0051, .006, .0055])
    vderror = zeros(vdd.shape[0])
    dp = array(
        [.031, .035, .038, .041, .045, .053, .09, .11, .17, .24, .35, .45, .8, 1.1, 37., 20, 17, 6.9, 1.8, 1.0, .55, .3,
         .17, .091, .051])
    dperror = zeros(vdd.shape[0])
    return dp, dperror, vdd, vderror


def coniferous_data():
    from numpy import array, zeros

    vdd = array([1.8, 1.4, 1.3, 1.45, 1.6, 1.4, 2.8, 3, 3.4, 3.3, 5, 6])
    dp = array([5, 7, 9, 10.1, 10.25, 10.4, 10.6, 10.8, 11, 11.1, 11.3, 11.6])
    vderror = array(
        [(.7 - 9.7) / 2, (3.8 - 11.9) / 2, ])
    dperror = array(
        [(3 - 27) / 2, (7 - 25) / 2, ])
    return dp, dperror, vdd, vderror


def grass_data():
    from numpy import array, zeros
    # 0:13 chamberlain 1967
    #13:18 nemitz 2002
    vdd = array([6.05, 3.3, .68, .09, .036, .0495, 4.1, .04, .017, .035, .063, .151])
    dp = array([23.6, 19, 5, 2, 1, .08, 23.6, .0255, .125, .155, .23, .425])
    vderror = array([3.15, .6, .52, 0, .017, .0265, 2.7, .02, .015, .031, .051, .119])
    dperror = array([0, 0, 0, 0, 0, 0, 0, .0145, .05, .05, .1, .025])
    return dp, dperror, vdd, vderror


def makeiceplot():
    from matplotlib.pylab import loglog, errorbar, xlabel, ylabel, title
    from prettyplotlib import legend
    from numpy import linspace, array
    # I HAVEN"T UPDATED THIS FUNCTION TO USE THE NEW DRY DEPOSITION SCHEMES
    #PLEASE SEE makewaterplot() and makegrassplot to see how to modify this function
    ust = .17  #cm/s
    mlu = 1  #ice
    z0 = .001  #surface roughness (m)
    vdd = []
    for i in linspace(.01e-6, 100.e-6, 10000):
        vdd.append(vd(ust, i, mlu, z0, mod=False)[0])
    vdzhang17 = vdd
    vdd = []
    for i in linspace(.01e-6, 100.e-6, 10000):
        vdd.append(vd(ust, i, mlu, z0, mod=True)[0])
    vdmod17 = vdd
    ust = .36  #cm/s
    vdd = []
    for i in linspace(.01e-6, 100.e-6, 10000):
        vdd.append(vd(ust, i, mlu, z0, mod=False)[0])
    vdzhang36 = vdd
    vdd = []
    for i in linspace(.01e-6, 100.e-6, 10000):
        vdd.append(vd(ust, i, mlu, z0, mod=True)[0])
    vdmod36 = vdd
    loglog(linspace(.01, 100, 10000), array(vdmod17) * 100, '-', label='Mod u*=17cm/s')
    loglog(linspace(.01, 100, 10000), array(vdzhang17) * 100, '--', label='Zhang01 u*=17cm/s')
    loglog(linspace(.01, 100, 10000), array(vdmod36) * 100, '-', label='Mod u*=36cm/s')
    loglog(linspace(.01, 100, 10000), array(vdzhang36) * 100, '--', label='Zhang01 u*=36cm/s')
    errorbar([.225, .75], [.034, .021], xerr=[.075, .25], fmt='o', label='Duan 1988 u*=12cm/s')
    errorbar([.023, .049], [.15, .055], xerr=[.001, .017], fmt='>', label='Nilsson 01 u*=18cm/s', color='green')
    errorbar([.69, .78, 7.2], [.035, .05, .14], xerr=[.001, .001, 2.4], fmt='<', label='Ibrahim1983 u*=16cm/s',
             color='red')
    legend(loc=2, fontsize=9)
    xlabel('dp ($\mu$m)')
    ylabel('V$_d$ (cm/s)')
    title('Ice')


def makewaterplot():
    from matplotlib.pylab import loglog, errorbar, xlabel, ylabel, title
    from dry_deposition import pleim, zhang01_class, petroff_zhang2010
    import seaborn as sns

    sns.set_context("talk")
    sns.set_style("whitegrid", {'xtick.major.size': 1, 'xtick.minor.size': .3})
    p = pleim.pleimdep()
    z = zhang01_class.zhangdep()
    pz = petroff_zhang2010.pzdep()

    p.do_all()
    z.do_all()
    pz.do_all()
    figure(figsize=(10, 7))
    f, ax = subplots(1, 1)
    loglog(p.dp, p.vd * 100, '-', label='CMAQ 11cm/s')
    loglog(pz.dp, pz.vd * 100, '--', label='Petroff 11cm/s')
    loglog(z.dp, z.vd * 100, '-.', label='Zhang01 11cm/s', linewidth=1)
    z.mod = True
    z.do_all()
    loglog(z.dp, z.vd * 100, '-.', label='MOD 11cm/s', linewidth=2)
    xticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])
    yticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])
    ax.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=.5)
    ax.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=.25)

    z.mod = False
    p.ustar = 1.17
    p.do_all()
    z.ustar = 1.17
    z.do_all()
    pz.ustar = 1.17
    pz.do_all()

    loglog(p.dp, p.vd * 100, '-', label='CMAQ 117cm/s')
    loglog(pz.dp, pz.vd * 100, '--', label='Petroff 117cm/s')
    loglog(z.dp, z.vd * 100, '-.', label='Zhang01 117cm/s', linewidth=1)
    z.mod = True
    z.do_all()
    loglog(z.dp, z.vd * 100, '-.', label='MOD 117cm/s', linewidth=2)

    dp, dperror, vddata, vderror = water_data()
    errorbar(dp[:14], vddata[:14], yerr=vderror[:14], xerr=dperror[:14], fmt='>', color='red',
             label='Moller70 u*=40cm/s')
    errorbar(dp[14:14 + 11], vddata[14:14 + 11], yerr=vderror[14:14 + 11], xerr=dperror[14:14 + 11], fmt='o',
             color='blue',
             label='Caffery98 u*=11-16cm/s', lolims=True, uplims=True)

    legend(loc=4, fontsize=9)
    xlabel('dp ($\mu$m)')
    ylabel('V$_d$ (cm/s)')
    title('water')


def makeshortgrassplot():
    from matplotlib.pylab import loglog, errorbar, xlabel, ylabel, title
    from dry_deposition import pleim, zhang01_class, petroff_zhang2010  # import dry deposition models
    import seaborn as sns  # this is an extra plot library not necessary if you only have matplotlib

    sns.set_context("talk")  # extra setting using seaborn package
    sns.set_style("whitegrid", {'xtick.major.size': 1, 'xtick.minor.size': .3})

    # load dry deposition modules
    p = pleim.pleimdep()
    z = zhang01_class.zhangdep()
    pz = petroff_zhang2010.pzdep()

    #change land use
    p.LUC, z.LUC, pz.LUC = 'LUC#13', 'LUC#13', 'LUC#13N'
    p.ustar, z.ustar, pz.ustar = .4, .4, .4
    p.do_all()
    z.do_all()
    pz.do_all()

    #plot all the functions
    figure(figsize=(10, 7))
    f, ax = subplots(1, 1)
    loglog(p.dp, p.vd * 100, '-', label='CMAQ 40cm/s')
    loglog(pz.dp, pz.vd * 100, '--', label='Petroff 40cm/s')
    loglog(z.dp, z.vd * 100, '-.', label='Zhang01 40cm/s', linewidth=1)

    #turn on the zhang01 modifications
    z.mod = True
    z.do_all()

    #plot zhang01 mod velocity deposition
    loglog(z.dp, z.vd * 100, '-.', label='MOD 40cm/s', linewidth=2)
    xticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])
    yticks([.01, .1, 1, 10, 100], ['.01', '.1', '1', '10', '100'])
    ax.grid(b=True, which='major', color='gray', linewidth=1.0, alpha=.5)
    ax.grid(b=True, which='minor', color='gray', linewidth=0.5, alpha=.25)

    #have to change back to turn off zhang01 modifications
    z.mod = False
    p.ustar, z.ustar, pz.ustar = .8, .8, .8
    p.do_all()
    z.do_all()
    pz.do_all()

    loglog(p.dp, p.vd * 100, '-', label='CMAQ 80cm/s')
    loglog(pz.dp, pz.vd * 100, '--', label='Petroff 80cm/s')
    loglog(z.dp, z.vd * 100, '-.', label='Zhang01 80cm/s', linewidth=1)
    z.mod = True
    z.do_all()
    loglog(z.dp, z.vd * 100, '-.', label='MOD 80cm/s', linewidth=2)

    #get the grass data
    dp, dperror, vddata, vderror = grass_data()

    #plot the grass data
    errorbar(dp[:7], vddata[:7], yerr=vderror[:7], xerr=dperror[:7], fmt='o', color='blue',
             label='Chamberlain67', lolims=True, uplims=True)
    errorbar(dp[7:7 + 5], vddata[7:7 + 5], yerr=vderror[7:7 + 5], xerr=dperror[7:7 + 5], fmt='*', color='green',
             label='Nemitz2002', lolims=True, uplims=True)


    #legend and axes label
    legend(loc=4, fontsize=9, bbox_to_anchor=(1.1, 0))
    xlabel('dp ($\mu$m)')
    ylabel('V$_d$ (cm/s)')
    title('Short Grass')




