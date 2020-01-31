from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value
from idaes.core import FlowsheetBlock
import Mod2_hda_ideal_VLE as thermo_props
from idaes.unit_models import Mixer
from idaes.core.util.model_statistics import degrees_of_freedom

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = thermo_props.HDAParameterBlock()
m.fs.Mixer = Mixer(default={"property_package": m.fs.properties,
                            "inlet_list": ["feed_1", "feed_2"]})
#
for comp in m.fs.properties.component_list:
    for phase in m.fs.properties.phase_list:
        m.fs.Mixer.feed_1.flow_mol_phase_comp[0,phase,comp].fix(1e-5)
        m.fs.Mixer.feed_2.flow_mol_phase_comp[0, phase, comp].fix(1e-5)


print("Degrees of Freedom =", degrees_of_freedom(m))
# m.fs.Mixer.report()