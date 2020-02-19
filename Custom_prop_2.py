"""
Initial property package for H2O-NaCl system
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, log, NonNegativeReals,\
    Var, Set, Param, sqrt, log10, TerminationCondition
from pyomo.opt import SolverFactory

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference, extract_data
from idaes.core.util.model_statistics import degrees_of_freedom, \
                                             number_unfixed_variables

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("NaClParameterBlock")
class NaClParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(NaClParameterData, self).build()

        self.state_block_class = IdealStateBlock

        self.component_list = Set(initialize=['H2O','NaCl'])

        self.phase_list = Set(initialize=['Liq'],ordered=True)

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list}

        # Source: google
        mw_comp_data = {'H2O': 18.0E-3,'NaCl': 58.4E-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mass': {'method': None, 'units': 'g/s'},
             'mass_frac': {'method': None, 'units': 'none'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'dens_mass': {'method': '_dens_mass', 'units': 'kg/m3'},
             'viscosity': {'method': '_viscosity', 'units': 'Pa-s'},
             'dens_mass_comp': {'method': '_dens_mass_comp', 'units': 'kg/m3'},
             'osm_coeff': {'method': '_osm_coeff', 'units': 'none'},
             'pressure_osm': {'method': '_pressure_osm', 'units': 'Pa'},
             'enth_mass_liq': {'method': '_enth_mass_liq', 'units': 'J/kg'}
             })

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'g'})
    # parameters
    R = 8.314  # gas constant [J/mol-K]
    MW = 58.44  # molecular weight [g/mol]
    osm_coeff_data = {'a1': 8.9453e-1,
                      'a2': 4.1561e-4,
                      'a3': -4.6262e-6,
                      'a4': 2.2211e-11,
                      'a5': -1.1445e-1,
                      'a6': -1.4783e-3,
                      'a7': -1.3526e-8,
                      'a8': 7.0132,
                      'a9': 5.696e-2,
                      'a10': -2.8624e-4}


class _IdealStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=1,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provied at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mol_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
            optarg : solver options dictionary object (default=None)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """

        _log.info('Starting {} initialization'.format(blk.name))

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")
        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        if optarg is None:
            sopt = {'tol': 1e-8}
        else:
            sopt = optarg

        opt = SolverFactory('ipopt')
        opt.options = sopt
# ---------------------------------------------------------------------
        # Initialize flow rates and compositions

        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])
        if free_vars > 0:
            try:
                results = solve_indexed_blocks(opt, [blk], tee=stee)
            except:
                results = None
        else:
            results = None

        if outlvl > 0:
            if results is None or results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Property initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Property initialization for "
                             "{} failed".format(blk.name))

        # ---------------------------------------------------------------------
        # Return state to initial conditions
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            _log.info("Initialization completed for {}".format(blk.name))

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} states released.'.format(blk.name))

@declare_process_block_class("IdealStateBlock",
                             block_class=_IdealStateBlock)
class IdealStateBlockData(StateBlockData):
    """An example property package for ideal VLE."""

    def build(self):
        """Callable method for Block construction."""
        super(IdealStateBlockData, self).build()

        # Add state variables
        self.flow_mass = Var(
            initialize=0.5,
            bounds=(1e-8, 100),
            doc='mass flow rate [g/s]')

        self.mass_frac = Var(
            initialize=0.1,
            bounds=(1e-8, 1),
            doc='mass fraction [unitless]')

        self.pressure = Var(
            initialize=101325,
            bounds=(101325, 400000),
            domain=NonNegativeReals,
            doc='State pressure [Pa]')

        self.temperature = Var(
            initialize=298.15,
            bounds=(298.15, 1000),
            domain=NonNegativeReals,
            doc='State temperature [K]')

        # Add supporting variables
        # def flow_mass_comp(b):
        #     return b.flow_mass * b.mass_frac
        # self.flow_mass_comp = Expression(rule=flow_mass_comp,
        #                                  doc='mass fraction [unitless]')

        #

# -----------------------------------------------------------------------------
# Property Methods
    def _dens_mass(self):
        self.dens_mass = Var(
            initialize=1e3,
            bounds=(1e-6, 1e6),
            doc="Mass density [kg/m3]")

        def rule_dens_mass(b):  # density [kg/m3]
            t = b.temperature - 273.15
            S = b.mass_frac * 1000
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
            return b.dens_mass == \
                   1e3 * (A1 * F1 + A2 * F2 + A3 * F3 + A4 * F4)

        self.eq_dens_mass = Constraint(rule=rule_dens_mass)

    def _viscosity(self):
        self.viscosity = Var(
            initialize=1e-3,
            bounds=(1e-8, 1),
            doc="Viscosity [Pa-s]")

        def rule_viscosity(b):  # dynamic viscosity [Pa-s]
            t = b.temperature - 273.15
            s = b.mass_frac
            mu_w = 4.2844e-5 + (0.157 * (t + 64.993) ** 2 - 91.296) ** -1
            A = 1.541 + 1.998e-2 * t - 9.52e-5 * t ** 2
            B = 7.974 - 7.561e-2 * t + 4.724e-4 * t ** 2
            return b.viscosity == mu_w * (1 + A * s + B * s ** 2)
        self.eq_viscosity = Constraint(rule=rule_viscosity)

    # TODO: scale diffusivity
    # def _diffusivity(self):
    #     self.diffusivist = Var(
    #         initialize=1e6,
    #         bounds=(1e-8, 1e2),
    #         doc="Diffiusivity [m2/s]")
    #
    #     def rule_diffusivity(b):  # diffusivity [m2/s]
    #         A = 3.847e-4
    #         B = -0.1984
    #         C = 26.54
    #         return 1e-9 * (A * b.temperature ** 2 + B * b.temperature + C)
    #
    #     self.eq_diffusivity = Constraint(rule=rule_diffusivity)

    def _dens_mass_comp(self):
        self.dens_mass_comp = Var(
                        initialize=1e2,
                        bounds=(1e-6, 1e6),
                        doc="Mass concentration [kg/m3]")

        self.eq_dens_mass_comp = Constraint(
            expr=self.dens_mass_comp == self.dens_mass * self.mass_frac)

    def _osm_coeff(self):
        self.osm_coeff = Var(
            initialize=1,
            bounds=(1e-8, 10),
            doc="Osmotic coefficient [unitless]")

        def rule_osm_coeff(b):  # osmotic coefficient [-], eq. 49
            s = b.mass_frac  # typo in Sharqawy, s is just mass_frac
            t = b.temperature - 273.15
            d = NaClParameterData.osm_coeff_data
            osm_coeff = (d['a1'] + d['a2'] * t + d['a3'] * t ** 2
                         + d['a4'] * t ** 4 + d['a5'] * s + d['a6'] * s * t
                         + d['a7'] * s * t ** 3 + d['a8'] * s ** 2
                         + d['a9'] * s ** 2 * t
                         + d['a10'] * s ** 2 * t ** 2)
            return b.osm_coeff == osm_coeff
        self.eq_osm_coeff = Constraint(rule=rule_osm_coeff)

    def _pressure_osm(self):
        self.pressure_osm = Var(
            initialize=1e6,
            bounds=(1, 1e8),
            doc="Osmotic pressure [Pa]")

        def rule_pressure_osm(b):  # osmotic pressure [Pa]
            i = 2  # number of ionic species
            R = NaClParameterData.R
            MW = NaClParameterData.MW
            return b.pressure_osm == \
                   (i * b.osm_coeff * b.dens_mass_comp * 1000 / MW
                    * R * b.temperature)
        self.eq_pressure_osm = Constraint(rule=rule_pressure_osm)

    # TODO: add vapor pressure, specific heat, thermal conductivitiy,
    #   and heat of vaporization

    def _enth_mass_liq(self):
        self.enth_mass_liq = Var(
            initialize=1e6,
            bounds=(1, 1e9),
            doc="Specific enthalpy [J/kg]")

        def rule_enth_mass_liq(b):  # specific enthalpy [J/kg]
            t = b.temperature - 273.15
            S = b.mass_frac
            h_w = 124.790 + 4203.075 * t - 0.552 * t ** 2 + 0.004 * t ** 3
            h_sw = (h_w - (S * (27062.623 + S) + S * (4835.675 + S) * t))
            return b.enth_mass_liq == h_sw
        self.eq_enth_mass_liq = Constraint(rule=rule_enth_mass_liq)


# -----------------------------------------------------------------------------
# General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if j == 'NaCl':
            return self.flow_mass * self.mass_frac
        elif j == 'H2O':
            return self.flow_mass * (1 - self.mass_frac)
        else:
            raise Exception("Property package only supports solutions"
                            "with NaCl and H2O")

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.flow_mass * self.enth_mass_liq

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        # if j in self._params.component_list:
        #     return self.dens_mol_phase[p] * self.mole_frac_phase_comp[p, j]
        # else:
        #     return 0
        pass

    def get_enthalpy_density_terms(self, p):
        """Create enthalpy density terms."""
        # return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]
        pass

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mass": self.flow_mass,
                "mass_frac": self.mass_frac,
                "temperature": self.temperature,
                "pressure": self.pressure}
