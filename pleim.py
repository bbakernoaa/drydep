from numpy import log, sqrt, mean, exp, pi, size, where, around, linspace
import common_calculations as com


roarow, amfp = 1.19, 0.067e-6
kb = 1.3806044503487214E-23


class pleimdep():
    # =======================================================
    # Parameters easy to change at runtime
    # dp              Particle diameter
    # rhop            particle density
    # ustar           friction velocity
    # LUC             Land Use Category (see function get_land_parameters() for information on LUC)
    # H               Surface Heat Flux (not used in all schemes)
    # Temp            Average temperature in surface layer
    # mod             Activates the modifications to the Zhang01 scheme using
    #                 resistance terms from Petroff and Zhang 2010
    # =========================================================
    # !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # !LUC#1   water
        # !LUC#2   ice
        # !LUC#3   inland lake
        # !LUC#4   evergreen needleleaf
        # !LUC#5   evergreen broadleaf
        # !LUC#6   deciduous needleleaf
        # !LUC#7   deciduous broadleaf
        # !LUC#8   tropical broadleaf
        # !LUC#9   drought deciduous forest
        # !LUC#10  evergreen broadleaf shrubs
        # !LUC#11  deciduous shrubs
        # !LUC#12  thorn shrubs
        # !LUC#13N short grass and forbs (needle like shape)
        # !LUC#13L short grass and forbs (leaf like shape)
        # !LUC#14  long grass
        # !LUC#15  crops
        # !LUC#16  rice
        # !LUC#17  sugar
        # !LUC#18  maize
        # !LUC#19  cotton
        # !LUC#20  irrigated crops
        # !LUC#21N  urban (coniferous trees)
        # !LUC#21L  urban (broadleaf trees)
        # !LUC#22  tundra
        # !LUC#23  swamp
        # !LUC#24  desert
        # !LUC#25  mixed wood forest
    #
    #
    # To use this scheme do the following
    #
    # from pleim import pleimdep
    # vd = pleimdep(dp=linspace(.01, 100., 100000), rhop=1500, ustar=.11, LUC='LUC#1', H=-.00001, Temp=293.,mod=False)
    # vd.do_all()
    #
    # To change LUC from LUC 1 to 10 and rerun
    #
    # from pleim import pleimdep
    # vd = pleimdep()
    # vd.LUC='LUC#10'
    # vd.do_all()
    def __init__(self, dp=linspace(.01, 100., 100000), rhop=1500, ustar=.11, LUC='LUC#1', H=-.00001, Temp=293.,
                 mod=False):

        self.amu = com.dynamic_viscosity(com.mua)
        self.rd = 0
        self.mod = mod
        self.eim = 0
        self.ein = 0
        self.st = 0
        self.vdrift = None
        self.lai_min = None
        self.dp = dp
        self.dpm = dp * 1e-6
        self.g = com.g
        self.ustar = ustar
        self.wstar = 0
        self.ra = 0
        self.pllp = None
        self.lai_f = 0
        self.vd = None
        self.LUC = LUC
        self.rhop = rhop
        self.H = H
        self.Tpe = Temp
        self.aest = None
        self.lai_max = None
        self.cin = None
        self.ll = None
        self.gamma = None
        self.pllp1 = None
        self.pllp2 = None
        self.cb = None
        self.AltiRef = 10.
        self.zl = 0
        self.LO = 1e-10
        self.dair = 0.369 * 29. + 6.29
        self.dh2o = 0.369 * 18. + 6.29
        self.zrough = 1e-10
        self.sc = None
        self.nua = com.nua

    def get_land_parameters(self):
        # !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # !LUC#1   water
        # !LUC#2   ice
        # !LUC#3   inland lake
        # !LUC#4   evergreen needleleaf
        # !LUC#5   evergreen broadleaf
        # !LUC#6   deciduous needleleaf
        # !LUC#7   deciduous broadleaf
        # !LUC#8   tropical broadleaf
        # !LUC#9   drought deciduous forest
        # !LUC#10  evergreen broadleaf shrubs
        # !LUC#11  deciduous shrubs
        # !LUC#12  thorn shrubs
        # !LUC#13N short grass and forbs (needle like shape)
        # !LUC#13L short grass and forbs (leaf like shape)
        # !LUC#14  long grass
        # !LUC#15  crops
        # !LUC#16  rice
        # !LUC#17  sugar
        # !LUC#18  maize
        # !LUC#19  cotton
        # !LUC#20  irrigated crops
        # !LUC#21N  urban (coniferous trees)
        # !LUC#21L  urban (broadleaf trees)
        # !LUC#22  tundra
        # !LUC#23  swamp
        # !LUC#24  desert
        # !LUC#25  mixed wood forest

        if (self.LUC == 'LUC#1') | (self.LUC == 'LUC#3'):
            self.aest = 100.
            self.gamma = 0.5
            self.pllp1 = -0.9
            self.pllp2 = -0.9
            self.cb = 0.0
            self.lai_max = 0.0
            self.lai_min = 0.0
            self.cin = 0
            self.ll = 0
            self.zrough = 0.11 * self.nua / self.ustar + 0.011 * self.ustar ** 2 / self.g
            print 'Water surface'

        elif self.LUC == 'LUC#2':
            self.aest = 50.
            self.gamma = 0.54
            self.pllp1 = -0.9
            self.pllp2 = -0.9
            self.cb = 0.0
            self.lai_max = 0.0
            self.lai_min = 0.0
            self.cin = 0
            self.ll = 0
            self.zrough = 0.9
            print 'Ice'

        elif self.LUC == 'LUC#4':
            self.aest = 1.
            self.gamma = 0.56
            self.pllp1 = 2
            self.pllp2 = 2
            self.cb = 0.887
            self.lai_max = 5
            self.lai_min = 5
            self.cin = 0.81
            self.ll = 0.15
            self.zrough = 0.01
            print 'Needle forest'

        elif self.LUC == 'LUC#5':
            self.aest = 0.8
            self.gamma = 0.56
            self.pllp1 = 5
            self.pllp2 = 5
            self.cb = 0.887
            self.lai_max = 5
            self.lai_min = 5
            self.cin = .216
            self.ll = 0.4
            self.zrough = 2.
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#6':
            self.aest = 1.1
            self.gamma = 0.56
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = 0.887
            self.lai_max = 6
            self.lai_min = 6
            self.cin = 0.81
            self.ll = 0.15
            self.zrough = 0.9
            print 'Needle forest'

        elif self.LUC == 'LUC#7':
            self.aest = .8
            self.gamma = .56
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = 1.262
            self.lai_max = 5
            self.lai_min = .1
            self.cin = .216
            self.ll = 3
            self.zrough = 1.
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#8':
            self.aest = 0.6
            self.gamma = .58
            self.pllp1 = 5
            self.pllp2 = 5
            self.cb = 1.262
            self.lai_max = 6
            self.lai_min = 6
            self.cin = .216
            self.ll = 4
            self.zrough = 2.5
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#9':
            self.aest = 1
            self.gamma = .56
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = 1.262
            self.lai_max = 4
            self.lai_min = 4
            self.cin = .216
            self.ll = 3
            self.zrough = 0.6
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#10':
            self.aest = 1.1
            self.gamma = .55
            self.pllp1 = 5
            self.pllp2 = 5
            self.cb = .93
            self.lai_max = 3
            self.lai_min = 3
            self.cin = .14
            self.ll = 2
            self.zrough = 0.2
            print 'Broadleaf shrubs'

        elif self.LUC == 'LUC#11':
            self.aest = 1.1
            self.gamma = .55
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = 0.93
            self.lai_max = 3
            self.lai_min = 0.5
            self.cin = .14
            self.ll = 2
            self.zrough = 0.2
            print 'Broadleaf shrubs'

        elif self.LUC == 'LUC#12':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = 0.93
            self.lai_max = 3
            self.lai_min = 3
            self.cin = .14
            self.ll = 2
            self.zrough = 0.2
            print 'Broadleaf shrubs'

        elif self.LUC == 'LUC#13':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = .996
            self.lai_max = 1
            self.lai_min = 1
            self.cin = .7
            self.ll = 0.5
            self.zrough = 0.04
            print 'short grass'

        elif self.LUC == 'LUC#14':
            self.aest = 1.2
            self.gamma = .55
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = .996
            self.lai_max = 2
            self.lai_min = 0.5
            self.cin = 0.162
            self.ll = 1
            self.zrough = 0.04
            print 'Leaf long grass'

        elif self.LUC == 'LUC#15':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = .996
            self.lai_max = 4
            self.lai_min = .1
            self.cin = .162
            self.ll = 3
            self.zrough = 0.1
            print 'Crops'

        elif self.LUC == 'LUC#16':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = .996
            self.lai_max = 6
            self.lai_min = 1 / 10
            self.cin = .162
            self.ll = 2
            self.zrough = 0.1
            print 'Rice'

        elif self.LUC == 'LUC#17':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = .996
            self.lai_max = 5
            self.lai_min = .1
            self.cin = .162
            self.ll = 4
            self.zrough = 0.1
            print 'Sugar'

        elif self.LUC == 'LUC#18':
            self.aest = 1.1
            self.gamma = .55
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = .996
            self.lai_max = 4
            self.lai_min = .1
            self.cin = .162
            self.ll = 5
            self.zrough = 0.1
            print 'Maize'

        elif self.LUC == 'LUC#19':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 5
            self.pllp2 = 10
            self.cb = .996
            self.lai_max = 5
            self.lai_min = .1
            self.cin = .162
            self.ll = 7
            self.zrough = 0.2
            print 'Cotton'

        elif self.LUC == 'LUC#20':
            self.aest = 1.2
            self.gamma = .54
            self.pllp1 = 2
            self.pllp2 = 5
            self.cb = .996
            self.lai_max = 1
            self.lai_min = 1
            self.cin = .162
            self.ll = 3
            self.zrough = 0.05
            print  'Irrigated crops'

        elif self.LUC == 'LUC#21':
            self.aest = 1.5
            self.gamma = .56
            self.pllp1 = 10
            self.pllp2 = 10
            self.cb = 1.262
            self.lai_max = 1
            self.lai_min = .10
            self.cin = 0.81
            self.ll = 0.15
            self.zrough = 1.
            print 'Urban center'

        elif self.LUC == 'LUC#22':
            self.aest = 50
            self.gamma = .54
            self.pllp1 = -0.9
            self.pllp2 = 0.7
            self.cb = 2
            self.lai_max = 0.1
            self.lai_min = 0.7
            self.cin = .81
            self.ll = 1
            self.zrough = 0.03
            print 'Tundra'

        elif self.LUC == 'LUC#23':
            self.aest = 2
            self.gamma = .54
            self.pllp1 = 10
            self.pllp2 = 10
            self.cb = 0.996
            self.lai_max = 4
            self.lai_min = 4
            self.cin = 0.162
            self.ll = 3
            self.zrough = 0.1
            print 'Swamp'

        elif self.LUC == 'LUC#24':
            self.aest = 50
            self.gamma = .54
            self.pllp1 = -.9
            self.pllp2 = 0
            self.cb = 0
            self.lai_max = 0
            self.lai_min = 0
            self.cin = 0
            self.ll = 2
            self.zrough = 0.04
            print 'Desert'

        elif (self.LUC == 'LUC#25') | (self.LUC == 'LUC#26'):
            self.aest = 0.8
            self.gamma = .56
            self.pllp1 = 5
            self.pllp2 = 5
            self.cb = 1.26
            self.lai_max = 5
            self.lai_min = 3
            self.cin = 0.81
            self.ll = 5
            self.zrough = 0.9

    def calc_atmosphere(self):
        self.LO = com.calc_monin_obukhov(self.ustar, self.Tpe, self.H)
        self.zl = self.AltiRef / self.LO
        self.zrl = self.zrough / self.LO
        # self.zl=0

    def aerosol_characteristics(self):
        # Cunningham slip factor
        self.cu = com.cfac(self.dpm)
        # Brownian Diffusivity
        self.db = com.brownian_diffusion(self.dpm, self.Tpe)
        # Schmidt number
        self.sc = com.schmidt_number(self.dpm, self.Tpe)
        # relaction time
        self.trel = com.relaxation_time(self.dpm, self.rhop)
        #settling velocity
        self.vg = self.g * self.trel

    def get_aerodynamic_resistance(self):
        if self.zl >= 0:
            if self.zl >= 1:
                self.psih0 = 1 - 8.21 - self.zl
            else:
                self.psih0 = - 8.21 * self.zl

            if self.zrough / self.LO > 1:
                self.psih = self.psih0 - (1 - 8.21 - self.zrl)
            else:
                self.psih = self.psih0 + 8.21 * self.zrl
        else:
            self.psih = 2.0 * log((1 + sqrt(1 - 11.60 * self.zl)) / (1 + sqrt(1 - 11.6 * self.zrl)))

        self.ra = 0.95 * (log(self.AltiRef / self.zrough) - self.psih) / (com.kappa * self.ustar)

        if self.zl >= 0:
            self.ra = (.74 * log(self.AltiRef / self.zrough) + 4.7 * self.zl) / (0.4 * self.ustar)
        else:
            self.ra = .74 / com.kappa / self.ustar * (
                log(self.AltiRef / self.zrough) - 2 * log((1 + sqrt(1 - 9 * self.zl)) * .5))

        self.ra = max([self.ra, 5])

    def funct_PSIM(self, xx):
        from numpy import log, arctan

        if (xx >= -2.) & (xx < 0.):
            PSIM = 2. * log(0.5 * (1. + (1. - 16. * xx) ** 0.25 )) + log(
                0.5 * (1. + (1. - 16. * xx) ** 0.5)) - 2. * arctan((1. - 16. * xx) ** 0.25) + self.PI / 2.
        elif ( xx >= 0.) & (xx <= 1.):
            PSIM = -5. * xx
        else:
            print 'out of stability range for PSIM'
        return PSIM

    def funct_PSIH(self, xx):
        from numpy import log

        if (xx >= -2.) & (xx < 0.):
            PSIH = 2. * log(0.5 * (1. + (1. - 16. * xx) ** (1. / 2.)))
        elif ( xx >= 0.) & (xx <= 1.):
            PSIH = -5. * xx
        else:
            print 'out of stability range for PSIH'
        return PSIH

    def funct_PHIM(self, xx):
        if (xx >= -2.) & (xx < 0.):
            PHIM = (1. - 16. * xx) ** (-1. / 4.)
        elif ( xx >= 0.) & (xx <= 1.):
            PHIM = 1. + 5. * xx
        else:
            print 'out of stability range for Phim'
        return PHIM

    def funct_PHIH(self, xx):
        # print xx
        if (xx >= -2.) & (xx < 0.):
            PHIH = (1. - 16. * xx) ** (-1. / 2.)
        elif ( xx >= 0.) & (xx <= 1.):
            PHIH = 1. + 5. * xx
        else:
            print 'out of stability range for Phih'
        return PHIH

    def funct_PSIH(self, xx):
        from numpy import log

        if (xx >= -2.) & (xx < 0.):
            PSIH = 2. * log(0.5 * (1. + (1. - 16. * xx) ** (1. / 2.)))
        elif ( xx >= 0.) & (xx <= 1.):
            PSIH = -5. * xx
        else:
            print 'out of stability range for PSIH'
        return PSIH

    def funct_Ra(self, sctN, kappa, ustar, z1, z2, PSIHz1, PSIHz2):
        from numpy import log

        Ra = sctN / (kappa * ustar) * (log(z1 / z2) - PSIHz1 + PSIHz2)
        return Ra

    def deposition_model(self):
        self.st = self.trel * self.ustar ** 2 / self.nua
        self.eim = self.st ** 2 / (400. + self.st ** 2)
        #self.eim[self.eim > 1] = 1.
        self.ein = 0
        self.rb = ((1 + .24 * self.wstar ** 2 / self.ustar ** 2) * self.ustar * (
            self.sc ** (-2 / 3) + self.eim + self.ein)) ** -1
        # 1 / ((self.sc ** (-2 / 3) + 10 ** (-3 / self.st))* (self.ustar + 0.24 * self.wstar ** 2 / self.ustar))

        if self.zl >= 0:
            self.ra = (.74 * log(self.AltiRef / self.zrough) + 4.7 * self.zl) / (0.4 * self.ustar)
        else:
            self.ra = .74 / com.kappa / self.ustar * (
                log(self.AltiRef / self.zrough) - 2 * log((1 + sqrt(1 - 9 * self.zl)) * .5))

        self.ra = max([self.ra, 5])
        #self.get_aerodynamic_resistance()


        self.vd = self.vg / (1 - exp(-self.vg * (self.ra + self.rb)))


    def do_all(self):
        self.get_land_parameters()
        self.calc_atmosphere()
        self.get_aerodynamic_resistance()
        self.aerosol_characteristics()
        self.deposition_model()

