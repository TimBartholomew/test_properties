# test property functions
import numpy as np
import matplotlib.pyplot as plt


class Prop:
    """
    Properties for NaCl solutions with temperature dependence.
    Sharqawy et al. Desalination and Water Treatment, 2010, 16, 354-380
    """

    def __init__(self, mass_frac=0.035, temperature=298):
        self.mass_frac = mass_frac  # mass fraction [-]
        self.temperature = temperature  # temperature [K]

    # data
    data_var_indep = {
        'mass_frac':
            {'label': 'Salinity',
             'unit': 'wt %',
             'scale': (0, 1 / 100),
             'range': (0, 18)},
        'temperature':
            {'label': 'Temperature',
             'unit': 'C',
             'scale': (-273.15, 1),
             'range': (0, 100)}
    }
    data_var_dep = {
        'dens_mass':
            {'label': 'Density',
             'unit': 'kg/m3',
             'scale': (0, 1),
             'fmt': '%.f'},
        'viscosity':
            {'label': 'Viscosity',
             'unit': 'mPa-s',
             'scale': (0, 1e-3),
             'fmt': '%.2f'},
        'diffusivity':
            {'label': 'Diffusivity',
             'unit': '1E-9 m2/s',
             'scale': (0, 1e-9),
             'fmt': '%.2f'},
        'dens_mass_comp':
            {'label': 'Concentration',
             'unit': 'g/L',
             'scale': (0, 1),
             'fmt': '%.1f'},
        'osm_coeff':
            {'label': 'Osmotic coefficient',
             'unit': 'unitless',
             'scale': (0, 1),
             'fmt': '%.2f'},
        'pressure_osm':
            {'label': 'Osmotic pressure',
             'unit': 'bar',
             'scale': (0, 1e5),
             'fmt': '%.1f'},
        'pressure_vap':
            {'label': 'Vapor pressure',
             'unit': 'kPa',
             'scale': (0, 1e3),
             'fmt': '%.1f'},
        'cp':
            {'label': 'Specific heat',
             'unit': 'kJ/kg-K',
             'scale': (0, 1e3),
             'fmt': '%.2f'},
        'therm_cond':
            {'label': 'Thermal conductivity',
             'unit': 'W/m-K',
             'scale': (0, 1),
             'fmt': '%.3f'},
        'hvap_mass':
            {'label': 'Specific heat of vaporization',
             'unit': 'kJ/kg',
             'scale': (0, 1e3),
             'fmt': '%.f'},
        'enth_mass_liq':
            {'label': 'Specific enthalpy',
             'unit': 'kJ/kg',
             'scale': (0, 1e3),
             'fmt': '%.f'}
    }

    def calculate(self):
        for k in Prop.data_var_dep.keys():  # loops through dependent variables
            v = getattr(self,'_'+k)()  # calculates dependent variable in SI
            setattr(self,k,v)  # sets dependent variable

    def show(self):
        self.calculate()
        for k in Prop.data_var_dep.keys():  # loops through dependent variables
            label = Prop.data_var_dep[k]['label']
            unit = Prop.data_var_dep[k]['unit']
            scale = Prop.data_var_dep[k]['scale']
            fmt = Prop.data_var_dep[k]['fmt']
            v = getattr(self,k)
            v_scaled = scale_from_si(v,scale)
            print((label + ' = ' + fmt + ' ' + unit) % v_scaled)

    def _dens_mass(self):  # density [kg/m3], Eq. 5
        t = self.temperature - 273.15
        S = self.mass_frac * 1000
        A = (2 * t - 200) / 160
        B = (2 * S - 150) / 150
        F1 = 0.5
        F2 = A
        F3 = 2 * A ** 2 - 1
        F4 = 4 * A ** 3 - 3 * A
        G1 = 0.5
        G2 = B
        G3 = 2 * B ** 2 - 1
        A1 = 4.032 * G1 + 0.115 * G2 + 3.26E-4 * G3
        A2 = -0.108 * G1 + 1.571e-3 * G2 - 4.23e-4 * G3
        A3 = -0.012 * G1 + 1.74e-3 * G2 - 9e-6 * G3
        A4 = 6.92e-4 * G1 - 8.7e-5 * G2 - 5.3e-5 * G3
        return 1e3 * (A1 * F1 + A2 * F2 + A3 * F3 + A4 * F4)

    def _viscosity(self):  # dynamic viscosity [Pa-s], Eq. 22
        t = self.temperature - 273.15
        s = self.mass_frac
        mu_w = 4.2844e-5 + (0.157 * (t + 64.993) ** 2 - 91.296) ** -1
        A = 1.541 + 1.998e-2 * t - 9.52e-5 * t ** 2
        B = 7.974 - 7.561e-2 * t + 4.724e-4 * t ** 2
        return mu_w * (1 + A * s + B * s ** 2)

    def _diffusivity(self):  # diffusivity [m2/s]
        # regression from Bartholomew et al. Cost optimization of MD
        # data from Fell and Hutchison, Journal of Chemical Engineering Data
        # 1971, 16 (4)
        A = 3.847e-4
        B = -0.1984
        C = 26.54
        return 1e-9 * (A * self.temperature ** 2 + B * self.temperature + C)

    def _dens_mass_comp(self):  # concentration [kg/m3]
        return self.dens_mass * self.mass_frac

    def _osm_coeff(self):  # osmotic coefficient [-], eq. 49
        s = self.mass_frac  # typo in Sharqawy, s is just mass_frac
        t = self.temperature - 273.15
        d = osm_coeff_data = {'a1': 8.9453e-1,
                              'a2': 4.1561e-4,
                              'a3': -4.6262e-6,
                              'a4': 2.2211e-11,
                              'a5': -1.1445e-1,
                              'a6': -1.4783e-3,
                              'a7': -1.3526e-8,
                              'a8': 7.0132,
                              'a9': 5.696e-2,
                              'a10': -2.8624e-4}
        osm_coeff = (d['a1'] + d['a2'] * t + d['a3'] * t ** 2
                     + d['a4'] * t ** 4 + d['a5'] * s + d['a6'] * s * t
                     + d['a7'] * s * t ** 3 + d['a8'] * s ** 2
                     + d['a9'] * s ** 2 * t
                     + d['a10'] * s ** 2 * t ** 2)
        return osm_coeff

    def _pressure_osm(self):  # osmotic pressure [Pa]
        i = 2  # number of ionic species
        R = 8.314  # gas constant [J/mol-K]
        MW = 58.44  # molecular weight [g/mol]
        return (i * self.osm_coeff * self.dens_mass_comp / MW
                * R * self.temperature)

    def _pressure_vap(self):  # vapor pressure [Pa], Eq. 29 and 53
        T = self.temperature
        S = self.mass_frac * 1000
        d = vap_pressure_data = {'a1': -5800,
                                 'a2': 1.391,
                                 'a3': -4.846e-2,
                                 'a4': 4.176e-5,
                                 'a5': -1.445e-8,
                                 'a6': 6.545}
        log_pv_w = (d['a1'] / T + d['a2'] + d['a3'] * T + d['a4'] * T ** 2
                   + d['a5'] * T ** 3 + d['a6'] * np.log(T))
        pv_w = np.exp(log_pv_w)
        pv_sw = pv_w / (1 + 0.57357 * (S / (1000 - S)))
        return pv_sw

    def _cp(self):  # specific heat [J/kg-K], Eq. 9
        T = self.temperature
        S = self.mass_frac * 1000
        A = 5.328 - 9.76e-2 * S + 4.04e-4 * S ** 2
        B = -6.913e-3 + 7.351e-4 * S - 3.15e-6 * S ** 2
        C = 9.6e-6 - 1.927e-6 * S + 8.23e-9 * S ** 2
        D = 2.5e-9 + 1.666e-9 * S - 7.125e-12 * S ** 2
        return 1000 * (A + B * T + C * T ** 2 + D * T ** 3)

    def _therm_cond(self):  # thermal conductivity [W/m-K], Eq. 13
        t = self.temperature - 273.15
        S = self.mass_frac * 1000
        log10_k = (np.log10(240 + 0.0002 * S)
                   + 0.434 * (2.3 - (343.5 + 0.037 * S) / (t + 273.15))
                   * (1 - (t + 273.15)/(647 + 0.03 * S)) ** (1/3))
        return 1e-3 * 10 ** log10_k

    def _hvap_mass(self):  # specific heat of vaporization [J/kg], Eq. 37 and 54
        t = self.temperature - 273.15
        s = self.mass_frac
        hvap_w = (2.501e6 - 2.369e3 * t + 2.678e-1 * t ** 2 - 8.103e-3 * t ** 3
                  - 2.079e-5 * t ** 4)
        return hvap_w * (1 - s)

    def _enth_mass_liq(self):  # specific enthalpy [J/kg], Eq. 43 and 55
        t = self.temperature - 273.15
        S = self.mass_frac
        h_w = 124.790 + 4203.075 * t - 0.552 * t ** 2 + 0.004 * t ** 3
        h_sw = (h_w - (S * (27062.623 + S) + S * (4835.675 + S) * t))
        return h_sw

    def visualize_2D(self, var_x=None, var_case=None, var_sim=None):
        # set up simulation
        x_arr_scaled = np.linspace(var_x['range'][0],
                                   var_x['range'][1],
                                   var_x['N'])  # var array scaled
        x_arr = scale_to_si(
            x_arr_scaled, var_x['scale'])  # var array in si
        case_arr_scaled = np.linspace(var_case['range'][0],
                                      var_case['range'][1],
                                      var_case['N'])  # var array scaled
        case_arr = scale_to_si(
            case_arr_scaled, var_case['scale'])  # var array in si
        # allocate results
        results = {}
        for var in var_sim.keys():
            results[var] = np.empty((var_x['N'], var_case['N']))

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
            plt.xlabel((var_x['label'] + ' [%s]') % var_x['unit'])
            plt.ylabel((var_sim[var]['label'] + ' [%s]') % var_sim[var]['unit'])
            plt.legend(case_arr_scaled,
                       title=(var_case['label'] + ' [%s]') % var_case['unit'])
            plt.show()


def scale_from_si(var, var_scale):
    return (var + var_scale[0]) / var_scale[1]


def scale_to_si(var, var_scale):
    return var * var_scale[1] - var_scale[0]


def visualize_properties(x_name='mass_frac', x_N=50,
                         case_name='temperature', case_N=5,
                         var_lst=[k for k in Prop.data_var_dep.keys()]):
    # assigning variables
    var_x = Prop.data_var_indep[x_name]
    var_x['name'] = x_name
    var_x['N'] = x_N
    var_case = Prop.data_var_indep[case_name]
    var_case['name'] = case_name
    var_case['N'] = case_N
    var_sim = {}
    for var in var_lst:
        var_sim[var] = Prop.data_var_dep[var]
    prop = Prop(0.035, 25 + 273.15)
    prop.visualize_2D(var_x=var_x, var_case=var_case, var_sim=var_sim)


prop1 = Prop(mass_frac=0)
# prop1.calculate()
prop1.show()
# print(prop1.cp)

# visualize_properties()

# plotting
var_lst = ['enth_mass_liq']
# salinity x axis
visualize_properties(var_lst=var_lst)
# temperature x axis
visualize_properties(x_name='temperature',
                     case_name='mass_frac',
                     var_lst=var_lst)


# prop1 = Prop(mass_frac=0.1)