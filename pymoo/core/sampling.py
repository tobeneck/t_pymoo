import numpy as np

from abc import abstractmethod

from pymoo.core.population import Population

from pymoo.core.individual import TraceTuple, TraceList


class Sampling:

    def __init__(self) -> None:
        """
        This abstract class represents any sampling strategy that can be used to create an initial population or
        an initial search point.
        """
        super().__init__()

    def do(self, problem, n_samples, pop=Population(), **kwargs):
        """
        Sample new points with problem information if necessary.

        Parameters
        ----------

        problem : :class:`~pymoo.core.problem.Problem`
            The problem to which points should be sampled. (lower and upper bounds, discrete, binary, ...)

        n_samples : int
            Number of samples

        pop : :class:`~pymoo.core.population.Population`
            The sampling results are stored in a population. The template of the population can be changed.
            If 'none' simply a numpy array is returned.

        Returns
        -------
        X : numpy.array
            Samples points in a two dimensional array

        """

        val = self._do(problem, n_samples, **kwargs)

        if pop is None:
            return val
        
        #create the offspring flags
        IsOff = np.zeros(len(val), dtype=bool) #in the initial population, all individuals are marked to be non offspring (or parents)

        #create the trace vectors
        T = []
        for indIndex in range(0, n_samples): #trace lists for all individuals
            curr_T = []
            for genomeIndex in range(0, problem.n_var): #trace list for each gene
                curr_trace_list = np.array( [ TraceTuple(indIndex, 1.0) ], dtype=TraceTuple )
                curr_T.append( TraceList(curr_trace_list) )
            currT = np.array(curr_T, dtype=TraceList)
            T.append(curr_T)

        T = np.array(T, dtype=TraceTuple)

        return Population.new("X", val, "IsOff", IsOff, "T", T)

    @abstractmethod
    def _do(self, problem, n_samples, **kwargs):
        pass



