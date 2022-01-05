# Changes
This File describes the chances made to the original pymoo framework. The two features wich were added are gene heritage tracking capabilitys with the Traceable Evolutionary Algorithm (T-EA), and tracking of the amount of offspring surviving to the next generation.

## Gene Tracing (T-EA)
### individual.py
To store the tracing information, this file now also includes the classes "traceList" and "traceTuple". The traceTuple is the object which stores the traceID and influence factor of the traceID. The traceList is basically a nest for a list of trace tuples. This is needed, as numpy can not cope with variable dimensions in their ndarray.

In the class "individual", a new member T is added as a paralel structure to the genome values (X), storing the heritage information for each genome. T is a numpy array with the same size as X containing the trace lists of the genes. TLs[0] is the trace list of X[0], TLs[1] the trace list of X[1] and so on. 
Each trace list is a variably sized list containing traceTuples.

The TLs need to be updated externally within the crossover and mutation operation if they are used.

### crossover.py
To activate gene tracing, the "crossover" class now has a "tracing_activated" flag. If the flag is True, the gene tracing steps will be executed. If the flag is False, the steps are skipped. If gene tracing is not required, deactivating it is recommended for performance reasons.

To calculate the trace vectors for the offspring of the crossover operation, two new methods are added to crossover.py. In "calculateOffspringTraceLists" the traceLists of the offspring individuals are computed and returned. The "recombineTraciLists" is a helper method calculating the trace list of one single genome.

The "do" method is also altered, now calling the "calculateOffspringTraceLists" depending on the "tracing_activated" flag.

### mutation.py
Similar to the crossover.py, the mutation also has a new flag called "tracing_activated" which can acitvate or deactivate the trace list calculations for performance reasons.

Again, similar to the crossover, there are two new methods to calculate the new trace list.  In "calculateOffspringTraceLists" the traceLists of the offspring individuals are computed and returned. The "mutateTraceList" is a helper method calculating the trace list of the newly mutated genome.

If the tracing is activated, there is a second "accumulate_mutations" flag. If accumulating the mutations is activated, the new traceID for all mutations will be -1. Otherwise, every new mutation will get a new negative number counting down from -1. If we wand to distinglish between the impact of different mutations, this flag should be deactivated. However, as the trace lists are expanding larger and larger for every new mutation, the runtime can be significantly longer.

### sampling.py
The "do()" method in sampling.py needed to be altered to initialize the trace lists for the T-EA.

### pm.py
In the file of the polynomial mutation the "kwargs" input in the "__init__" was missing. This needed to be added for the "tracing_activated" flag to be passed to the parent.

## Offspring Survivors
### individual.py
The flag "isOff" is added. By default, each newly initialized individual is marked as offspring, except when specified otherwise. Similar to the gene tracing, the "IsOff" flag is updated externally.

### sampling.py
As all newly created individuals would be flagged as offspring by default, the "do()" method was altered to create non-offspring individuals.

### algorithm.py
After calling the callback for every generation, the "IsOff" flag needs to be changed to "False", as in the next generation all individuals will be parents for the new offspring. For this reason, these flags are all reset at the end of the "\_post_advance()" method, importantly after the callback call.