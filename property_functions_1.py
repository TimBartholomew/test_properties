# test property functions
import numpy as np
import matplotlib.pyplot as plt


class Prop:
    """
    Properties for NaCl solutions with no temperature dependence.
    Bartholomew and Mauter. Journal of Membrane Science. 2019, 573, 682-693
    """

    def __init__(self, mass_frac=0.035, temperature=298):
        self.mass_frac = mass_frac  # mass fraction [-]
        self.temperature = temperature  # temperature [K]

    def calculate(self):
        self.dens_mass = self._dens_mass()  # density [g/m3]
        self.viscosity = self._viscosity()  # viscosity [Pa-s]
        self.diffusivity = self._diffusivity()  # diffusivity [m2/s]
        self.dens_mass_comp = self._dens_mass_comp()  # concentration [g/m3]
        self.osm_coeff = self._osm_coeff()  # osmotic coefficient [-]
        self.pressure_osm = self._pressure_osm()  # osmotic pressure [Pa]

    def show(self):
        self.calculate()
        print('Mass fraction = %.3f' % self.mass_frac)
        print('Temperature = %.1f C' % (self.temperature - 273))
        print('Density = %.1f kg/m3' % (self.dens_mass / 1000))
        print('Viscosity = %.1f mPa-s' % (self.viscosity * 1000))
        print('Diffusivity = %.2f 1E-9 m2/s' % (self.diffusivity / 1e-9))
        print('Concentration = %.1f g/L' % (self.dens_mass_comp / 1000))
        print('Osmotic Pressure = %.1f bar' % (self.pressure_osm / 1e5))

    def _dens_mass(self):  # density [g/m3]
        return 1000 * (995 + 756 * self.mass_frac)

    def _viscosity(self):
        return 2.15e-3 * self.mass_frac + 9.80e-4

    def _diffusivity(self):
        return ((153 * self.mass_frac ** 4 - 122 * self.mass_frac ** 3
                 + 30.1 * self.mass_frac ** 2 - 2.00 * self.mass_frac + 1.51)
                * 1e-9)

    def _dens_mass_comp(self):  # concentration [g/m3]
        return self.dens_mass * self.mass_frac

    def _osm_coeff(self):  # osmotic coefficient [-]
        c = self.dens_mass_comp / 1000  # concentration [g/L]
        return 3.14e-6 * c ** 2 + 2.13e-4 * c + 0.917

    def _pressure_osm(self):  # osmotic pressure [Pa]
        i = 2  # number of ionic species
        R = 8.314  # gas constant [J/mol-K]
        MW = 58.44  # molecular weight [g/mol]
        return (i * self.osm_coeff * self.dens_mass_comp / MW
                * R * self.temperature)

    def visualize_2D(self, var_x=None, var_case=None, var_sim=None):
        # set up simulation
        x_arr_scaled = np.linspace(var_x['range'][0],
                                   var_x['range'][1],
                                   var_x['N'])  # var array in scaled units
        x_arr = scale_to_si(
            x_arr_scaled, var_x['scale'])  # var array in si
        case_arr_scaled = np.linspace(var_case['range'][0],
                                      var_case['range'][1],
                                      var_case[
                                          'N'])  # var array in scaled units
        case_arr = scale_to_si(
            case_arr_scaled, var_case['scale'])  # var array in si
        # allocate results
        results = {}
        for var in var_sim.keys():
            results[var] = np.empty((var_x['N'], var_case['N']))
            results[var].fill(np.nan)

        # simulate results
        for j in range(len(case_arr)):
            setattr(self, var_case['name'], case_arr[j])
            for i in range(len(x_arr)):
                setattr(self, var_x['name'], x_arr[i])
                self.calculate()
                for var in var_sim.keys():
                    results[var][i][j] = getattr(self, var)

        # plot results
        for var in var_sim.keys():
            result = results[var]
            result_scaled = scale_from_si(result, var_sim[var]['scale'])
            plt.plot(x_arr_scaled, result_scaled)
            plt.xlabel(var_x['label'])
            plt.ylabel(var_sim[var]['label'])
            plt.legend(case_arr_scaled, title=var_case['label'])
            plt.show()


def scale_from_si(var, var_scale):
    return (var + var_scale[0]) / var_scale[1]


def scale_to_si(var, var_scale):
    return var * var_scale[1] - var_scale[0]


def visualize_properties(x_name='mass_frac', x_N=50,
                         case_name='temperature', case_N=3,
                         var_lst=['dens_mass', 'viscosity', 'diffusivity',
                                  'dens_mass_comp','osm_coeff',
                                  'pressure_osm']):
    # data
    data_var_indep = {
        'mass_frac':
            {'range': (0, 26),
             'scale': (0, 1 / 100),
             'label': 'Salinity [wt %]'},
        'temperature':
            {'range': (0, 100),
             'scale': (-273, 1),
             'label': 'Temperature [C]'}
    }
    data_var_dep = {
        'dens_mass':
            {'scale': (0, 1000),
             'label': 'Density [kg/m3]'},
        'viscosity':
            {'scale': (0, 1e-3),
             'label': 'Viscosity [mPa-s]'},
        'diffusivity':
            {'scale': (0, 1e-9),
             'label': 'Diffusivity [1E-9 m2/s]'},
        'dens_mass_comp':
            {'scale': (0, 1000),
             'label': 'Concentration [g/L]'},
        'osm_coeff':
            {'scale': (0, 1),
             'label': 'Osmotic coefficient [-]'},
        'pressure_osm':
            {'scale': (0, 1e5),
             'label': 'Osmotic pressure [bar]'}
    }
    # assigning variables
    var_x = data_var_indep[x_name]
    var_x['name'] = x_name
    var_x['N'] = x_N
    var_case = data_var_indep[case_name]
    var_case['name'] = case_name
    var_case['N'] = case_N
    var_sim = {}
    for var in var_lst:
        var_sim[var] = data_var_dep[var]
    prop = Prop(0.035, 25 + 273)
    prop.visualize_2D(var_x=var_x, var_case=var_case, var_sim=var_sim)


prop1 = Prop(mass_frac=0.1)
prop1.show()
visualize_properties()
