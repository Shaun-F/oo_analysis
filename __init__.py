"""
__init__.py: main init file for analysis

Created by: Erik Lentz
Creation Date: 10/26/18
"""
# undertake analysis calculations
from control import core

analysis = core.core_analysis()
analysis.execute()

# more on analytics?
