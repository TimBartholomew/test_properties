"""
Initial property package for H2O-NaCl system
"""

# Import Python libraries
import logging

# Import Pyomo libraries
# from pyomo.environ import Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
# from idaes.core.util.misc import extract_data

# Some more information about this module
__author__ = "Tim Bartholomew"
__version__ = "0.0.1"

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("HDAParameterBlock")
class HDAParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(HDAParameterData, self).build()

        self.state_block_class = IdealStateBlock

        self.component_list = Set(initialize=['H2O',
                                              'NaCl'])

        self.phase_list = Set(initialize=['Liq', 'Vap'],
                              ordered=True)

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # Source: google
        mw_comp_data = {'H2O': 18.0E-3,
                        'NaCl': 58.4E-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")