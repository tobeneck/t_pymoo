import numpy as np

from pymoo.core.population import Population

from pymoo.core.individual import TraceTuple, TraceList


class Crossover:
    """
    The crossover combines parents to offsprings. Some crossover are problem specific and use additional information.
    This class must be inherited from to provide a crossover method to an algorithm.
    """

    def __init__(self, n_parents, n_offsprings, prob=0.9, tracing_activated=True):
        self.prob = prob
        self.n_parents = n_parents
        self.n_offsprings = n_offsprings
        self.tracing_activated = tracing_activated




    def recombineTraceLists(self, p1_x, p2_x, p1_tl, p2_tl, new_c_x):
        """
        this method calculates the new trace list for for the gene values a and b based on the new_c_x.

        Parameters
        ----------
        p1_x: Genome Value
            the genome value of parent 1
        p1_tl: TraceList
            the trace vector of parent 1
        p2_x: Genome Value
            the genome value of parent 2
        p2_tl: TraceList
            the trace vector of parent 2
        new_c_x: Genome Value
            the genome value of the child of parent 1 and parent 2
            

        Returns
        -------
        new_c_tl : TraceList
            the trace vector of the child of parent 1 and parent 2, based on the distance between the parents and the child value
        """

        #if aVal = bVal = newVal, both genes should have 50% impact
        influence_factor_p1 = 0.5
        influence_factor_p2 = 0.5

        if(p1_x == new_c_x and p2_x != new_c_x): #return a if it has 100% influence. Faster and avoids 0.0 impact traceTuples
            return p1_tl
        if(p2_x == new_c_x and p1_x != new_c_x): #return b if it has 100% influence. Faster and avoids 0.0 impact traceTuples
            return p2_tl
        if(p1_x != new_c_x and p2_x != new_c_x): #compute the new values if non of them are equal
            #if you don't cast here, the resulting values are integers and will be roundet!
            influence_factor_p1 = 1.0 - ( abs(p1_x - new_c_x) / ( abs(p1_x - new_c_x) + abs(p2_x - new_c_x) ) )
            influence_factor_p2 = 1.0 - ( abs(p2_x - new_c_x) / ( abs(p1_x - new_c_x) + abs(p2_x - new_c_x) ) )

        i = 0 #index for trace vector a
        j = 0 #index for trace vector b

        new_c_tl = TraceList()

        while True: #this iterates over the traceList of this individual
            if i >= p1_tl.len() and j >= p2_tl.len(): #stop if both vectors are empty
                break
            elif i >= p1_tl.len() and not (j >= p2_tl.len()):#append if the a vector is empty and b vector is not.
                currentP2ID = p2_tl.get(j).traceID
                currentP2Influence = p2_tl.get(j).influenceFactor
                new_c_tl.append(TraceTuple(currentP2ID, influence_factor_p2 * currentBInfluence))
                j = j + 1
            elif not (i >= p1_tl.len()) and j >= p2_tl.len(): #append if the b vector is empty and a vector is not.
                currentP1ID = p1_tl.get(i).traceID
                currentAInfluence = p1_tl.get(i).influenceFactor
                new_c_tl.append(TraceTuple(currentP1ID, influence_factor_p1 * currentAInfluence))
                i = i + 1
            else: #if both arrays are not empty, append the next traceID:
                currentP1ID = p1_tl.get(i).traceID
                currentP2ID = p2_tl.get(j).traceID

                currentAInfluence = p1_tl.get(i).influenceFactor
                currentBInfluence = p2_tl.get(j).influenceFactor

                if currentP1ID == currentP2ID: #combine the two if equal
                    new_c_tl.append(TraceTuple(currentP1ID, influence_factor_p1 * currentAInfluence + influence_factor_p2 * currentBInfluence))
                    i = i + 1
                    j = j + 1

                if currentP1ID < currentP2ID: #add the traceID of a if its smaller than the traceID of b
                    new_c_tl.append(TraceTuple(currentP1ID, influence_factor_p1 * currentAInfluence))
                    i = i + 1

                if currentP2ID < currentP1ID: #add the traceID of b if its smaller than the traceID of a
                    new_c_tl.append(TraceTuple(currentP2ID, influence_factor_p2 * currentBInfluence))
                    j = j + 1

        return new_c_tl



    def calculateOffspringTraceLists(self, parents_X, offspring_X, parents_T):
        """
        this methid calculates the new trace lists for the offsprings of the crossover operation based on the parents and offspring gene values.

        Parameters
        ----------
        parents_X: numpy.array
            parents_X[0] contains the genome values for the first half of parents and parents_X[1] for the second half 
        offspring_X: numpy.array
            offspring_X[0] contains the genome values for the first half of offspring and offspring_X[1] for the second half 
        parents_T: numpy.array
            parents_T[0] contains the trace lists for the first half of parents and parents_T[1] for the second half 
            

        Returns
        -------
        offspring_T : 
            the trace list of the child of parent 1 and parent 2, based on the distance between the parents and the child value        
        """

        parents1_X = parents_X[0]
        parents2_X = parents_X[1]

        parents1_T = parents_T[0]
        parents2_T = parents_T[1]

        children1_X = np.copy(offspring_X[0])
        children2_X = np.copy(offspring_X[1])

        children1_T = []
        children2_T = []

        for indIndex in range(0, len(parents1_X)): #iterate over the individuals
            parent1_X = parents1_X[indIndex]
            parent2_X = parents2_X[indIndex]
            parent1_T = parents1_T[indIndex]
            parent2_T = parents2_T[indIndex]
            child1_X  = children1_X[indIndex]
            child2_X  = children2_X[indIndex]

            child1_T = []
            child2_T = []

            for genomeIndex in range(0, len(parent1_X)): #iterate over the genome vector
                #get the genome values and trace lists of the current parents and the genome values of the offspring
                parent1_current_gene_x = parent1_X[genomeIndex]
                parent2_current_gene_x = parent2_X[genomeIndex]
                parent1_current_gene_TL = parent1_T[genomeIndex]
                parent2_current_gene_TL = parent2_T[genomeIndex]
                child1_current_gene_x  = child1_X[genomeIndex]
                child2_current_gene_x  = child2_X[genomeIndex]

                #append the new trace list to the trace lists vector of the two offspring individuals
                child1_T.append( self.recombineTraceLists(parent1_current_gene_x, parent2_current_gene_x, parent1_current_gene_TL, parent2_current_gene_TL ,child1_current_gene_x) )
                child2_T.append( self.recombineTraceLists(parent1_current_gene_x, parent2_current_gene_x, parent1_current_gene_TL, parent2_current_gene_TL ,child2_current_gene_x) )
            
            #append the child trace lists to the children trace lists vector
            children1_T.append(child1_T)
            children2_T.append(child2_T)
        
        return np.array([children1_T, children2_T], dtype=TraceList)

    def do(self, problem, pop, parents, **kwargs):
        """

        This method executes the crossover on the parents. This class wraps the implementation of the class
        that implements the crossover.

        Parameters
        ----------
        problem: class
            The problem to be solved. Provides information such as lower and upper bounds or feasibility
            conditions for custom crossovers.

        pop : Population
            The population as an object

        parents: numpy.array
            The select parents of the population for the crossover

        kwargs : dict
            Any additional data that might be necessary to perform the crossover. E.g. constants of an algorithm.

        Returns
        -------
        offsprings : Population
            The off as a matrix. n_children rows and the number of columns is equal to the variable
            length of the problem.

        """

        if self.n_parents != parents.shape[1]:
            raise ValueError('Exception during crossover: Number of parents differs from defined at crossover.')

        # get the design space matrix as well as the traceLists form the population and parents
        parents_X  = pop.get("X")[parents.T].copy()
        if self.tracing_activated:
            parents_T = pop.get("T")[parents.T].copy()

        # now apply the crossover probability
        do_crossover = np.random.random(len(parents)) < self.prob

        # execute the crossover
        offspring_X = self._do(problem, parents_X, **kwargs)

        offspring_X[:, do_crossover, :] = offspring_X[:, do_crossover, :] #here the old parent values are replaced with the new child values where _do_crossover is true

        #create the trace lists of the offspring if tracing is enabled
        if self.tracing_activated:
            offspring_T = self.calculateOffspringTraceLists(parents_X, offspring_X, parents_T)

        # flatten the genome array to become a 2d-array
        offspring_X = offspring_X.reshape(-1, offspring_X.shape[-1])

        # also flatten the trace list array if tracing is enabled
        if self.tracing_activated:
            offspring_T = offspring_T.reshape(-1, offspring_T.shape[-1])

        # create and return a population object depending on if tracing is activated
        if self.tracing_activated:
            return Population.new("X", offspring_X, "T", offspring_T)
        else:
            return Population.new("X", offspring_X)


class ElementwiseCrossover(Crossover):
    """
    The purpose of this crossover meta-crossover is the convenience of applying the operation an individual level
    (elementwise). This wrapper transform all crossovers in pymoo to such an operation.
    """

    def __init__(self, crossover):
        self.crossover = crossover

    def do(self, problem, *args, **kwargs):
        pop = Population.create(*args)
        parents = np.arange(len(args))[None, :]
        return self.crossover.do(problem, pop, parents, **kwargs)
