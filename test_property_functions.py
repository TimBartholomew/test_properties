# test property functions
class Prop:
    def __init__(self,mass_frac,temperature):
        self.mass_frac = mass_frac  # mass fraction [-]
        self.temperature = temperature + 273.15  # temperature [K]

    def calculate(self):
        self.dens_mass = self._dens_mass()  # density [g/m3]
        self.dens_mass_comp = self._dens_mass_comp()  # concentration [g/m3]
        self.pressure_osm = self._pressure_osm()  # osmotic pressure [Pa]

    def result(self):
        self.calculate()
        print('Density = %.1f kg/m3' % (self.dens_mass / 1000))
        print('Concentration = %.1f g/L' % (self.dens_mass_comp / 1000))
        print('Osmotic Pressure = %.1f bar' % (self.pressure_osm / 1e5))

    def _dens_mass(self): # density [g/m3]
        return 1000 * (995 + 756 * self.mass_frac)

    def _dens_mass_comp(self):  # concentration [g/m3]
        return self.dens_mass * self.mass_frac

    def _osm_coeff(self):  # osmotic coefficient [-]
        c = self.dens_mass_comp / 1000  # concentration [g/L]
        return 3.14e-6 * c ** 2 + 2.13e-4 * c + 0.917

    def _pressure_osm(self):  # osmotic pressure [Pa]
        osm_coeff = self._osm_coeff()
        i = 2  # number of ionic species
        R = 8.314  # gas constant [J/mol-K]
        MW = 58.44 # molecular weight [g/mol]
        return i * osm_coeff * self.dens_mass_comp / MW * R * self.temperature

prop1 = Prop(0.035,25)
prop1.result()




# def osm_coef(X,T):
#     s = X * 1000
#     t = T - 273.15
#     osm_coeff_data = {'a1': 8.9453e-1,
#                       'a2': 4.1561e-4,
#                       'a3': -4.6262e-6,
#                       'a4': 2.2211e-11,
#                       'a5': -1.1445e-1,
#                       'a6': -1.4783e-3,
#                       'a7': -1.3526e-8,
#                       'a8': 7.0132,
#                       'a9': 5.696e-2,
#                       'a10': -2.8624e-4}
#     d = osm_coeff_data
#     osm_coeff = (d['a1'] + d['a2'] * t + d['a3'] * t ** 2
#                  + d['a4'] * t ** 4 + d['a5'] * s + d['a6'] * s * t
#                  + d['a7'] * s * t ** 3 + d['a8'] * s ** 2
#                  + d['a9'] * s ** 2 * t
#                  + d['a10'] * s ** 2 * t ** 2)

