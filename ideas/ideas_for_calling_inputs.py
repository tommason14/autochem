"""
mol
  - opt
  - spec
  - hess
"""

for path, dirs, files in os.walk('.'):
    if path.split('.')[-1].endswith('opt'):
        Gamess_input(..., run="optimize", ...)
    if path.split('.')[-1].endswith('spec'):
        Gamess_input(..., run="single point", ...)
    if path.split('.')[-1].endswith('hess'):
        Gamess_input(..., run="hessian", ...)
