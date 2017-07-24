__author__ = 'barry'  #


class pzdep:
    from numpy import linspace, pi
    # =======================================================
    # Parameters easy to change at runtime
    # dp              Particle diameter
    # rhop            particle density
    # ustar           friction velocity
    # LUC             Land Use Category (see function get_land_parameters() for information on LUC)
    # H               Surface Heat Flux (not used in all schemes)
    # Temp            Average temperature in surface layer
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
    # from petroff_zhang2010 import pzdep
    # vd = pzdep(dp=linspace(.01, 100., 100000), rhop=1500, ustar=.11, LUC='LUC#1', H=-.00001, Temp=293.,mod=False)
    # vd.do_all()
    #
    # To change LUC from LUC 1 to 10 and rerun
    #
    # from petroff_zhang2010 import pzdep
    # vd = pzdep()
    # vd.LUC='LUC#10'
    # vd.do_all()


    def __init__(self, dp=linspace(.01, 100., 10000), rhop=1500, ustar=.11, LUC='LUC#1', H=-.1, Tpe=293.):
        self.dp = dp
        self.dpm = 0.000001 * self.dp
        self.H = H
        self.Tpe = Tpe
        self.LUC = LUC
        self.ustar = ustar
        self.rhop = rhop
        self.eb = 0
        self.ein = 0
        self.eit = 0
        self.eim = 0
        self.kappa = 0.4
        from numpy import pi

        self.PI = pi
        self.kb = 1.38e-23
        self.rhoa = 1.2
        self.cpa = 1000.
        self.egt = 1e-10
        self.egb = 1e-10
        self.lpm = 0.067E-6
        self.g = 9.81
        self.mua = 1.89E-5
        self.nua = 1.57E-5
        self.LO = None
        self.cu = None  # Cunningham slip factor
        self.Db = None  # Brownian Diffusivity (m^2 s^-1)
        self.sc = None  # Schmidt Number
        self.trel = None  # partical relaction time (s)
        self.ws = None  # sedimentation velocity (m s^-1)
        self.Uh, self.lm, self.alphaVeg = None, None, None
        self.xb, self.xin, self.xim, self.xit, self.lai, self.ObstSize = 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10
        self.ObstShape, self.hcanop, self.zrough, self.hdepl = 1e-10, 1e-10, 1e-10, 1e-10
        self.AltiRef, self.kxN, self.xbN, self.xinN, self.ximN, self.ObstSizeN = 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10
        self.kxL, self.xitN, self.xbL, self.xinL = 1e-10, 1e-10, 1e-10, 1e-10
        self.ximL, self.xitL, self.ObstSizeL, self.kx = 1e-10, 1e-10, 1e-10, 1e-10

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
            self.AltiRef = 10.
            self.hdepl = 0.  # ! reference height (m), displacement height (m)
            self.zrough = 0.11 * self.nua / self.ustar + 0.011 * self.ustar ** 2 / self.g  #! roughness length (m)
            self.xit = 1.
            print 'Water surface'

        elif self.LUC == 'LUC#2':
            self.AltiRef = 10.
            self.hdepl = 0.
            self.zrough = 0.01
            self.xit = 1.
            print 'Ice'

        elif self.LUC == 'LUC#4':
            self.ObstShape = 'needle'
            self.kx = 0.27  # !plagiophile
            self.xb = 1.9
            self.xin = 1.5
            self.xim = 0.6
            self.xit = 0.
            self.hcanop = 15.
            self.zrough = 0.9
            self.hdepl = 12.
            self.lai = 10.
            self.ObstSize = 1.5e-3
            self.AltiRef = 2. * self.hcanop  # !height(m)of the canopy top, roughness, length(m),! displacement height(m), total
            # leaf area index(-), characteristic size of the obstacle(m) reference height(m)
            print 'Needle forest'

        elif self.LUC == 'LUC#5':
            self.ObstShape = 'leaf'
            self.kx = 0.216  # !plagiophile
            self.xb = 1.9
            self.xin = 2.
            self.xim = 0.6
            self.xit = 0.4
            self.hcanop = 33.33
            self.zrough = 2.
            self.hdepl = 26.67
            self.lai = 12.
            self.ObstSize = 4e-2
            self.AltiRef = 2. * self.hcanop
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#6':
            self.ObstShape = 'needle'
            self.kx = 0.27  #!plagiophile
            self.xb = 1.9
            self.xin = 1.5
            self.xim = 0.6
            self.xit = 0.
            self.hcanop = 15.
            self.AltiRef = 2. * self.hcanop
            self.zrough = 0.9
            self.hdepl = 12.
            self.lai = 10.
            self.ObstSize = 1.5e-3  #! max roughness and self.lai
            print 'Needle forest'

        elif self.LUC == 'LUC#7':
            self.ObstShape = 'leaf'
            self.kx = 0.216  #!plagiophile
            self.xb = 1.9
            self.xin = 2.
            self.xim = 0.6
            self.xit = 0.4
            self.hcanop = 16.67
            self.AltiRef = 2. * self.hcanop  #
            self.zrough = 1.
            self.hdepl = 13.33
            self.lai = 10.
            self.ObstSize = 3e-2  #! maxroughness and self.lai
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#8':
            self.ObstShape = 'leaf'
            self.kx = 0.216  #!plagiophile
            self.xb = 1.9
            self.xin = 2.
            self.xim = 0.6
            self.xit = 0.4
            self.hcanop = 41.67
            self.zrough = 2.5
            self.hdepl = 33.33
            self.lai = 12.
            self.ObstSize = 4e-2
            self.AltiRef = 2. * self.hcanop
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#9':
            self.ObstShape = 'leaf'
            self.kx = 0.216  #!plagiophile
            self.xb = 1.9
            self.xin = 2.
            self.xim = 0.6
            self.xit = 0.4
            self.hcanop = 16.67
            self.zrough = 0.6
            self.hdepl = 13.33
            self.lai = 8.
            self.ObstSize = 3e-2
            self.AltiRef = 2. * self.hcanop
            print 'Broadleaf forest'

        elif self.LUC == 'LUC#10':
            self.ObstShape = 'leaf'
            self.kx = 0.216  # !plagiophile
            self.xb = 1.4
            self.xin = 1.3
            self.xim = 0.4
            self.xit = 0.1
            self.hcanop = 1.54
            self.zrough = 0.2
            self.hdepl = 0.98
            self.lai = 6.
            self.ObstSize = 2e-2
            self.AltiRef = 10.
            print 'Broadleaf shrubs'

        elif self.LUC == 'LUC#11':
            self.ObstShape = 'leaf'
            self.kx = 0.216  # !plagiophile
            self.xb = 1.4
            self.xin = 1.3
            self.xim = 0.4
            self.xit = 0.1
            self.AltiRef = 10.
            self.hcanop = 1.54
            self.zrough = 0.2
            self.hdepl = 0.98
            self.lai = 6.
            self.ObstSize = 2e-2  #!max roughness and self.lai
            print 'Broadleaf shrubs'

        elif self.LUC == 'LUC#12':
            self.ObstShape = 'leaf'
            self.kx = 0.216  #!plagiophile
            self.xb = 1.4
            self.xin = 1.3
            self.xim = 0.4
            self.xit = 0.1
            self.hcanop = 1.54
            self.zrough = 0.2
            self.hdepl = 0.98
            self.lai = 6.
            self.ObstSize = 2e-2
            self.AltiRef = 10.
            print 'Broadleaf shrubs'

        elif (self.LUC == 'LUC#13N'):
            self.ObstShape = 'needle'
            self.kx = 1. / self.PI  #!vertical needle
            self.xb = 1.5
            self.xin = 1.1
            self.xim = 0.6
            self.xit = 0.3
            self.hcanop = 0.31
            self.zrough = 0.04
            self.hdepl = 0.20
            self.lai = 2.
            self.ObstSize = 5e-3
            self.AltiRef = 10.
            print 'Needle short grass'

        elif (self.LUC == 'LUC#13L'):
            self.ObstShape = 'leaf'
            self.kx = 1. / self.PI  #! verticalleaf
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.6
            self.xit = 0.3
            self.hcanop = 0.31
            self.zrough = 0.04
            self.hdepl = 0.20
            self.lai = 2.
            self.ObstSize = 5e-3
            self.AltiRef = 10.
            print 'Leaf short grass'

        elif (self.LUC == 'LUC#14'):
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 4.
            self.ObstSize = 1e-2  #!max  self.hcanop, rough, self.hdepl and self.lai
            print 'Leaf long grass'

        elif (self.LUC == 'LUC#15'):
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 8.
            self.ObstSize = 3e-2  #!max  self.hcanop, rough, self.hdepl and self.lai
            print 'Crops'

        elif self.LUC == 'LUC#16':
            self.ObstShape = 'leaf'
            self.kx = 0.270  # leaf        erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 12.
            self.ObstSize = 2e-2  #!max        self.hcanop, rough, self.hdepl and self.lai
            print 'Rice'

        elif self.LUC == 'LUC#17':
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 10.
            self.ObstSize = 4e-2  #!max        self.hcanop, rough, self.hdepl and self.lai
            print 'Sugar'

        elif self.LUC == 'LUC#18':
            self.ObstShape = 'leaf'
            self.kx = 0.270  # !leaf erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 8.
            self.ObstSize = 5e-2  #!max self.hcanop, rough, self.hdepl and self.lai
            print 'Maize'

        elif (self.LUC == 'LUC#19'):
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf        erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 1.54
            self.zrough = 0.2
            self.hdepl = 0.98
            self.lai = 10.
            self.ObstSize = 7e-2  #!max        self.hcanop, rough, self.hdepl and self.lai
            print 'Cotton'

        elif self.LUC == 'LUC#20':
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf        erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.38
            self.zrough = 0.05
            self.hdepl = 0.25
            self.lai = 10.
            self.ObstSize = 3e-2  #!max        self.hcanop  rough  self.hdepl and self.lai
            print  'Irrigated crops'

        elif self.LUC == 'LUC#21N':
            self.ObstShape = 'needle'
            self.kx = 0.27  #!plagiophile
            self.xb = 1.9
            self.xin = 1.5
            self.xim = 0.6
            self.xit = 0.
            self.hcanop = 17.
            self.zrough = 1.
            self.hdepl = 11.9
            self.lai = 1.
            self.ObstSize = 1.5e-3
            self.AltiRef = 2. * self.hcanop
            print 'Urban center with coniferous trees'

        elif self.LUC == 'LUC#21L':
            self.ObstShape = 'leaf'
            self.kx = 0.216  #!plagiophile
            self.xb = 1.9
            self.xin = 2.
            self.xim = 0.6
            self.xit = 0.4
            self.hcanop = 17.
            self.zrough = 1.
            self.hdepl = 11.9
            self.lai = 1.
            self.ObstSize = 3.e-2
            self.AltiRef = 2. * self.hcanop
            print 'Urban center with broadleaf trees'

        elif self.LUC == 'LUC#22':
            self.ObstShape = 'needle'
            self.kx = 1. / self.PI  #!vertical        needle
            self.xb = 1.5
            self.xin = 1.1
            self.xim = 0.6
            self.xit = 0.3
            self.AltiRef = 10.
            self.hcanop = 0.23
            self.zrough = 0.03
            self.hdepl = 0.14
            self.lai = 4.
            self.ObstSize = 5e-3  #!max
            self.lai
            print 'Tundra'

        elif self.LUC == 'LUC#23':
            self.ObstShape = 'leaf'
            self.kx = 0.270  #!leaf        erectophile
            self.xb = 1.5
            self.xin = 1.2
            self.xim = 0.3
            self.xit = 0.4
            self.AltiRef = 10.
            self.hcanop = 0.77
            self.zrough = 0.1
            self.hdepl = 0.49
            self.lai = 8.
            self.ObstSize = 4e-2  #!max        obstacle        size
            print 'Swamp'

        elif self.LUC == 'LUC#24':
            self.AltiRef = 10.
            self.hdepl = 0.
            self.zrough = 0.04
            self.xit = 1.
            print 'Desert'

        elif (self.LUC == 'LUC#25') | (self.LUC == 'LUC#26'):
            self.ObstShape = 'mixed'
            self.hcanop = 15.
            self.zrough = 0.9
            self.hdepl = 12.
            self.lai = 10.
            self.AltiRef = 2. * self.hcanop
            self.kxN = 0.27  # #!plagiophile
            self.xbN = 1.9
            self.xinN = 1.5
            self.ximN = 0.6
            self.xitN = 0.
            self.ObstSizeN = 1.5e-3
            self.kxL = 0.216  # #!plagiophile
            self.xbL = 1.9
            self.xinL = 2.
            self.ximL = 0.6
            self.xitL = 0.4
            self.ObstSizeL = 3e-2
            self.mix_ratio = .5
            self.kx = self.mix_ratio * self.kxN + (1. - self.mix_ratio) * self.kxL

    def calc_atmosphere(self):  # calculation of the atmosperhic characteristics
        # self.Tpe is taken as an approximation of the virtual potential temperature of the atmosphere
        # self.LO is the monin obhukov length
        # self.kappa is the von Karmin number
        from numpy import log

        self.LO = -self.ustar ** 3. / self.kappa * self.Tpe / self.g * self.rhoa * self.cpa / self.H
        self.Uh = self.ustar / self.kappa * (
            log((self.hcanop - self.hdepl) / self.zrough) - self.funct_PSIM(
                (self.hcanop - self.hdepl) / self.LO) + self.funct_PSIM(self.zrough / self.LO))
        self.lm = self.kappa * (self.hcanop - self.hdepl) / self.funct_PHIH((self.hcanop - self.hdepl) / self.LO)
        #! extinction coefficient for aerodynamic characteristics within vegetative canopies
        self.alphaVeg = (self.kx * self.lai / (12. * self.kappa ** 2. * (1. - self.hdepl / self.hcanop) ** 2.)) ** (
            1. / 3.) * (self.funct_PHIM((self.hcanop - self.hdepl) / self.LO)) ** (2. / 3.)

    def aerosol_characteristics(self):
        from numpy import exp

        # Cunningham slip factor
        self.cu = 1. + 2. * self.lpm / self.dpm * (1.257 + 0.4 * exp(-1.1 * self.dpm / (2. * self.lpm)))
        # Brownian Diffusivity (m^2 s^-1)
        self.Db = (self.cu * self.kb * self.Tpe) / (3. * self.PI * self.mua * self.dpm)
        # Schmidt Number
        self.sc = self.nua / self.Db
        # partical relaction time (s)
        self.trel = self.rhop * (self.dpm ** 2.) * self.cu / (18. * self.mua)
        # sedimentation velocity (m s^-1)
        self.ws = self.g * self.trel

    def calc_ground_efficiencies(self):
        from numpy import sqrt, log, exp, max, arctan, size

        # ground deposition
        self.egb = self.sc ** (-2. / 3.) / (14.5 * (self.PI / (6. * sqrt(3.)) + 1. / sqrt(3.) * arctan(
            (2. * self.sc ** (1. / 3.) / 2.9 - 1.) / sqrt(3.)) + 1. / 6. * log(
            (1. + (self.sc ** (1. / 3.) / 2.9)) ** 2. / (
                1. - self.sc ** (1. / 3.) / 2.9 + (self.sc ** (1. / 3.) / 2.9) ** 2.))))

        self.egt = (3.5E-4 / (self.nua ** 2.)) * self.trel ** 2. * (self.ustar * exp(-self.alphaVeg)) ** 4.
        n = range(size(self.trel))
        index = (self.trel * (self.ustar * exp(-self.alphaVeg)) ** 2. / self.nua) < 20.
        self.egt[index] = (3.5E-4 / (self.nua ** 2.)) * self.trel[index] ** 2. * (self.ustar * exp(-self.alphaVeg)) ** 4
        index = (self.trel * (self.ustar * exp(-self.alphaVeg)) ** 2. / self.nua) >= 20.
        self.egt[index] = 0.14
        # for i in n:
        #     if (self.trel * (self.ustar * exp(-self.alphaVeg)) ** 2. / self.nua)[i] < 20.:
        #         self.egt[i] = (3.5E-4 / (self.nua ** 2.)) * self.trel[i] ** 2. * (
        #                                                                              self.ustar * exp(-self.alphaVeg)) ** 4.
        #     else:
        #         self.egt[i] = 0.14

    def calc_veg_efficiencies(self):
        """

        @type self: object
        """
        from numpy import log, size

        n = size(self.trel)
        # turbulent impaction
        if self.ObstShape == 'mixed':
            self.eit = range(n)
            for i in range(n):
                if (self.trel[i] * self.ustar ** 2. / self.nua) < 20.:
                    self.eit[i] = (self.mix_ratio * self.xitN + (1. - self.mix_ratio) * self.xitL) * (
                        3.5E-4 / (self.nua ** 2.)) * self.trel ** 2. * self.ustar ** 4.
                else:
                    self.eit[i] = (self.mix_ratio * self.xitN + (1. - self.mix_ratio) * self.xitL) * 0.14
        else:
            self.eit = range(n)
            for i in range(n):
                if (self.trel[i] * self.ustar ** 2. / self.nua) < 20.:
                    self.eit[i] = self.xit * (3.5E-4 / (self.nua ** 2.)) * self.trel[i] ** 2. * self.ustar ** 4.
                else:
                    self.eit[i] = self.xit * 0.14

        if self.ObstShape == 'leaf':
            # brownian Diffusion efficiency
            self.eb = self.xb * 0.664 * (self.nua / self.Db) ** (-2. / 3.) * (self.Uh * self.ObstSize / self.nua) ** (
                -1. / 2.)
            # interception efficiency
            self.ein = self.xin * self.kx / 2. * self.dp * 1.e-6 / self.ObstSize * (
                2. + log(4. * self.ObstSize / (self.dp * 1.e-6)))
            # Inertial Impaction
            self.eim = self.xim * self.kx / (1. + 0.47 * self.ObstSize / (self.trel * self.Uh) ) ** 2.
        elif self.ObstShape == 'needle':
            self.eb = self.xb * 0.467 * (self.nua / self.Db) ** (-2. / 3.) * (self.Uh * self.ObstSize / self.nua) ** (
                -1. / 2.)
            self.ein = self.xin * 2. * self.kx * self.dp * 1.e-6 / self.ObstSize
            self.eim = self.xim * self.kx / (1. + 0.6 * self.ObstSize / (self.trel * self.Uh) ) ** 2.
        elif (self.ObstShape == 'mixed'):
            self.eb = self.mix_ratio * self.xbN * 0.467 * (self.nua / self.Db) ** (-2. / 3.) * (
                                                                                                   self.Uh * self.ObstSizeN / self.nua) ** (
                                                                                                   -1. / 2.) + (
                                                                                                                   1. - self.mix_ratio) * self.xbL * 0.664 * (
                                                                                                                                                                 self.nua / self.Db) ** (
                                                                                                                                                                 -2. / 3.) * (
                                                                                                                                                                                 self.Uh * self.ObstSizeL / self.nua) ** (
                                                                                                                                                                                 -1. / 2.)
            self.ein = self.mix_ratio * self.xinN * 2. * self.kxN * self.dp * 1.e-6 / self.ObstSizeN + (
                                                                                                           1. - self.mix_ratio) * self.xinL * self.kxL / 2. * self.dp * 1.e-6 / self.ObstSizeL * (
                                                                                                           2. + log(
                                                                                                               4. * self.ObstSizeL / (
                                                                                                                   self.dp * 1.e-6)))
            self.eim = self.mix_ratio * self.ximN * self.kxN / (1. + 0.6 * self.ObstSizeN / (
                self.trel * self.Uh)) ** 2. + (1. - self.mix_ratio) * self.ximL * self.kxL / (1. + 0.47 * self.ObstSizeL / (
                self.trel * self.Uh)) ** 2.

    def deposition_model(self):
        from numpy import tanh, sqrt

        self.Qveg = self.hcanop / self.lm * self.lai * (
            self.Uh / self.ustar * (self.eb + self.ein + self.eim) + self.eit)
        self.Qsol = self.hcanop / self.lm * (self.egb + self.egt)
        self.eta = sqrt(1. + 4. * self.Qveg / self.alphaVeg ** 2.)
        self.vds = self.ustar * ( self.egb + self.egt ) * (
            1. + (self.Qveg / self.Qsol - self.alphaVeg / 2.) * 2. / (self.alphaVeg * self.eta) * tanh(
                self.alphaVeg * self.eta / 2.)) / (
                       1. + (self.Qsol + self.alphaVeg / 2.) * 2. / (self.alphaVeg * self.eta) * tanh(
                           self.alphaVeg * self.eta / 2.))

        #!!!!!! Drift velocity accounting for the effects of gravity and phoretic effects.
        if (self.LUC == 'LUC#1') | (self.LUC == 'LUC#3') | (self.LUC == 'LUC#23'):
            self.vdrift = self.ws + 5.e-5
        elif self.LUC == 'LUC#2':
            self.vdrift = self.ws + 2.e-4
        else:
            self.vdrift = self.ws

        #!!!!!!!!!!! Formulation of vd(m.s - 1) without vegetation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (self.LUC == 'LUC#1') | (self.LUC == 'LUC#2') | (self.LUC == 'LUC#3') | (self.LUC == 'LUC#24'):
            z1 = self.AltiRef - self.hdepl
            z2 = self.zrough
            PSIHz1 = self.funct_PSIH((self.AltiRef - self.hdepl) / self.LO)
            PSIHz2 = self.funct_PSIH(self.zrough / self.LO)
            sctN = 1
            ra = self.funct_Ra(sctN, self.kappa, self.ustar, z1, z2, PSIHz1, PSIHz2)
            self.vd = self.vdrift + 1. / (ra + 1. / (self.ustar * (self.eit + self.egb)))
        #!!!!!!!!!!! Formulation of vd(m.s - 1) with vegetation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            z1 = self.AltiRef - self.hdepl
            z2 = self.hcanop - self.hdepl
            PSIHz1 = self.funct_PSIH((self.AltiRef - self.hdepl) / self.LO)
            PSIHz2 = self.funct_PSIH((self.hcanop - self.hdepl) / self.LO)
            sctN = 1
            ra = self.funct_Ra(sctN, self.kappa, self.ustar, z1, z2, PSIHz1, PSIHz2)
            self.vd = self.vdrift + 1. / (ra + 1. / self.vds)

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

    def do_all(self):
        self.get_land_parameters()
        self.calc_atmosphere()
        self.aerosol_characteristics()
        self.calc_ground_efficiencies()
        self.calc_veg_efficiencies()
        self.deposition_model()

