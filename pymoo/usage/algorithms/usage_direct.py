from pymoo.algorithms.soo.nonconvex.direct import DIRECT
from pymoo.optimize import minimize
from pymoo.problems.single import Himmelblau

problem = Himmelblau()

algorithm = DIRECT()

ret = minimize(problem,
               algorithm,
               verbose=True)

