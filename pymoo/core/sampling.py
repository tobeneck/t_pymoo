import numpy as np

from abc import abstractmethod

from pymoo.core.population import Population

from pymoo.core.individual import TraceTuple, TraceList

from enum import Enum

class TracingTypes(Enum):
    NO_TRACING = 0
    TRACE_ID = 1
    TRACE_LIST = 2
    TRACE_VECTOR = 3

class Sampling:

    def __init__(self) -> None:
        """
        This abstract class represents any sampling strategy that can be used to create an initial population or
        an initial search point.
        """


    def do(self, problem, n_samples, pop=Population(), tracing_type=TracingTypes.NO_TRACING, **kwargs):
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

        ParentIDs = [] # the IDs of the parent from the previous generation #TODO: implement

        Birthday = 0 # the generation the individual was generated in #TODO: implement



        # build the tracing information based on the trace types
        if tracing_type == TracingTypes.NO_TRACING:
            return Population.new("X", val, "IsOff", IsOff)
        
        if tracing_type == TracingTypes.TRACE_ID:
            T = np.zeros( (n_samples, problem.n_var) )
            for i in range(0, n_samples):
                T[i] = np.zeros(problem.n_var) + i + 1
            return Population.new("X", val, "IsOff", IsOff, "T", T, "IsOff", IsOff)
        
        if tracing_type == TracingTypes.TRACE_LIST:
            T = []
            for indIndex in range(0, n_samples): #trace lists for all individuals
                curr_T = []
                for genomeIndex in range(0, problem.n_var): #trace list for each gene
                    curr_trace_list = np.array( [ TraceTuple(indIndex + 1, 1.0) ], dtype=TraceTuple )
                    curr_T.append( TraceList(curr_trace_list) )
                currT = np.array(curr_T, dtype=TraceList)
                T.append(curr_T)
            T = np.array(T, dtype=TraceTuple)
            return Population.new("X", val, "IsOff", IsOff, "T", T, "IsOff", IsOff)
        
        if tracing_type == TracingTypes.TRACE_VECTOR:
            T = np.zeros( (n_samples, problem.n_var, n_samples + 1) ) #TODO: this is wrong!
            for i in range(0, n_samples):
                T[i, :, i] = 1.0
            return Population.new("X", val, "IsOff", IsOff, "T", T, "IsOff", IsOff)



    @abstractmethod
    def _do(self, problem, n_samples, **kwargs):
        pass



