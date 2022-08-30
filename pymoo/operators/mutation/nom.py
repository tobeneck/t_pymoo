from pymoo.core.mutation import Mutation
from pymoo.core.sampling import TracingTypes


class NoMutation(Mutation):

    def _do(self, problem, X, **kwargs):
        return X
