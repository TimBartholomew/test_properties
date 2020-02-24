from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Expression, Var
from pyomo.network import Port
from idaes.core import (FlowsheetBlock)
import Custom_prop_2 as props
# import Mod2_hda_ideal_VLE as props
from idaes.unit_models import Mixer
from idaes.unit_models.mixer import MixingType, MomentumMixingType
from idaes.core.util.model_statistics import degrees_of_freedom

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = props.NaClParameterBlock()
# m.fs.properties = props.HDAParameterBlock()
m.fs.Mixer = Mixer(default={"property_package": m.fs.properties,
                            "inlet_list": ["feed_1", "feed_2"],
                            "energy_mixing_type": MixingType.extensive,
                            "momentum_mixing_type": MomentumMixingType.none})


m.fs.Mixer.feed_1.flow_mass[0].fix(0.5)
m.fs.Mixer.feed_1.mass_frac[0].fix(0.1)
m.fs.Mixer.feed_1.temperature[0].fix(273.15 + 50)
m.fs.Mixer.feed_2.flow_mass[0].fix(0.5)
m.fs.Mixer.feed_2.mass_frac[0].fix(0.035)
m.fs.Mixer.feed_2.temperature[0].fix(273.15 + 25)
# m.fs.Mixer.outlet.temperature[0].fix(273.15 + 40)
# m.fs.Mixer.mixed_state[0].dens_mass
# m.fs.Mixer.mixed_state[0].viscosity
# m.fs.Mixer.mixed_state[0].dens_mass_comp
m.fs.Mixer.mixed_state[0].pressure_osm
# m.fs.Mixer.mixed_state[0].osm_coeff
# m.fs.Mixer.mixed_state[0].enth_mass_liq
m.fs.Mixer.initialize(outlvl=0)
print("Degrees of Freedom =", degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-6, 'max_iter': 5000}
results = solver.solve(m, tee=False)
assert results.solver.termination_condition == TerminationCondition.optimal
# m.display()
#
m.fs.Mixer.display()
for p in m.fs.Mixer.component_objects(Port, descend_into=True):
    p.display()
# m.fs.Mixer.pprint()