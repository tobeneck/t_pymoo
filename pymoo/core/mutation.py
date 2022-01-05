import numpy as np

from pymoo.core.individual import TraceTuple, TraceList

class Mutation:

    def __init__(self, tracing_activated=True, accumulate_mutations=True) -> None: #TODO: tracing_activated to false
        super().__init__()
        self.algorithm = None
        self.problem = None
        self.tracing_activated=tracing_activated
        self.accumulate_mutations=accumulate_mutations
        self.mutation_counter = -1


    def mutateTraceList(self, p_x, new_c_x, p_tl):
        """
        calculates the trace list of the offspring gene based on the gene values of the offspring and the parent

        Parameters
        ----------
        p_x: Gene Value
            the gene value of the parent gene
        new_c_x: Gene Value
            the gene value of the offspring
        p_tl: numpy.array
            p_tl contains the trace lists of the parents
            

        Returns
        -------
        new_traceList : 
            the trace list of the child of parent 1 and parent 2, based on the distance between the parents and the child value        
        """

        #inclrease the mutation counter if we do not accumulate the mutations
        if not self.accumulate_mutations:
            self.mutation_counter += 1

        if p_x == new_c_x : #return the parents trace list if the value was not altered
            return p_tl

        influence_factor_mut = abs(p_x - new_c_x) / ( abs(p_x) + abs(p_x - new_c_x) )
        influence_factor_old = 1 - influence_factor_mut

        new_traceList = TraceList()

        i = 0 #index for this traceList

        #check if we accumulate and the first element is already the mutation
        if self.accumulate_mutations and p_tl.get(0).traceID == self.mutation_counter: #if we accumulate mutation and the gene is influenced by mutation, accumulate that
            oldScaledMutImpact = p_tl.get(0).influenceFactor * influence_factor_old
            new_traceList.append(TraceTuple(self.mutation_counter, influence_factor_mut + oldScaledMutImpact))
            i = i + 1 #increment i, as the first element of the traceList from p_tl is already added
        else: #add the new mutation ID if we don't accumulate or there is no mutation present bevore
            new_traceList.append(TraceTuple(self.mutation_counter, influence_factor_mut))

        while i < p_tl.len(): #this iterates over the traceList of this individual
            currentAID = p_tl.get(i).traceID
            currentAImpact = p_tl.get(i).influenceFactor

            new_traceList.append(TraceTuple(currentAID, influence_factor_old * currentAImpact))
            i = i + 1

        return new_traceList


    #calculates the new trace lists from the parents values as well as the new child  values
    def calculateOffspringTraceLists(self, parents_X, children_X, parents_T):
        """
        this methid calculates the new trace lists for the offsprings of the mutation operation based on the parent and offspring gene value.

        Parameters
        ----------
        parents_X: numpy.array
            parents_X contains the genome values of the parents
        children_X: numpy.array
            children_X contains the genome values of the offspring
        parents_T: numpy.array
            parents_T contains the trace lists of the parents
            

        Returns
        -------
        offspring_TL : 
            the trace list of the child of parent 1 and parent 2, based on the distance between the parents and the child value        
        """

        offspring_T = []

        for indIndex in range(0, len(parents_X)):
            parent_X = parents_X[indIndex]
            parent_T = parents_T[indIndex]
            child_X = children_X[indIndex]

            child_T = []

            for genomeIndex in range(0, len(parent_X)):
                parent_current_gene_value = parent_X[genomeIndex]
                parent_current_TL = parent_T[genomeIndex]
                child_current_gene_value = child_X[genomeIndex]

                child_T.append( self.mutateTraceList(parent_current_gene_value, child_current_gene_value, parent_current_TL) )
            
            offspring_T.append(child_T)

        return np.array(offspring_T, dtype=TraceList)

    def do(self, problem, pop, **kwargs):
        """

        Mutate variables in a genetic way.

        Parameters
        ----------
        problem : class
            The problem instance - specific information such as variable bounds might be needed.
        pop : Population
            A population object

        Returns
        -------
        X_m : Population
            The mutated population.

        """
        X_p = pop.get("X") #the parent genome values
        X_m = self._do(problem, X_p, **kwargs) #the mutated genome values

        #calculate the new trace lists if tracing is activated
        if self.tracing_activated:
            T_m = self.calculateOffspringTraceLists(X_p, X_m, pop.get("T")) #the trace lists of the mutated offspring
        

        # create and return a population object depending on if tracing is activated
        if self.tracing_activated:
            return pop.new("X", X_m, "T", T_m) #TODO: offspring flaggs
        else:
            return pop.new("X", X_m) #TODO: offspring flaggs

    def _do(self, problem, X, **kwargs):
        pass
