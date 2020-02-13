from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
import Custom_prop_1 as props

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = props.NaClParameterBlock()
m.fs.stream = m.fs.properties.state_block_class([0],default={"parameters":m.fs.properties})


m.display()
m.fs.stream[0].flow_mass.fix(1)
x_NaCl = 0.035
m.fs.stream[0].mass_frac_comp['NaCl'].fix(x_NaCl)
m.fs.stream[0].mass_frac_comp['H2O'].fix(1 - x_NaCl)

m.fs.stream[0].dens_mass
m.fs.stream[0].dens_mass_comp
m.fs.stream[0].pressure_osm
m.display()

print(degrees_of_freedom(m))

solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-6, 'max_iter': 5000}
results = solver.solve(m, tee=False)
assert results.solver.termination_condition == TerminationCondition.optimal
m.display()
