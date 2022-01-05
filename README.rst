
.. |python| image:: https://img.shields.io/badge/python-3.9-blue.svg
   :alt: python 3.9

.. |license| image:: https://img.shields.io/badge/license-apache-orange.svg
   :alt: license apache
   :target: https://www.apache.org/licenses/LICENSE-2.0




|python| |license|




t_pymoo: The Traceable Evolutionary Algorithm in pymoo
====================================================================

This repository is a fork of the popular `pymoo <https://github.com/anyoptimization/pymoo>`_ framework, currently of the release 5.0. It implements heritage tracking capabilitys with the `traceable evolutionary algorithm (T-EA) <https://ci.ovgu.de/Publications/CEC_RBM_2020-p-870.html>`_. With this method, heritage information of the initial population can be tracked throughout the run of an EA.

A list with all changes to the original framework can be found `here <Changes_to_pymoo.md>`_.


How to Install
********************************************************************************

The installation of this framework is basically the same as manually installing pymoo, a guide of which can be found in its `official documentation <https://pymoo.org/installation.html>`_. Basically, just clone this repository and install it manually with pip.

.. code:: bash

    git clone https://github.com/tobeneck/t_pymoo
    cd pymoo
    pip install .




How to Use
********************************************************************************

The framework still works like pymoo, if you are unfamiliar with it look into its fantastic `documentation <https://pymoo.org/index.html>`_. The tracing information for the individuals are found in the variable T and can be retrieved similarly to the genome values X or fitness information F, as seen in the example below. Here the traceIDs influencing the first gene of the first individual of the final population are printed with their respective influence factor on the gene.

.. code:: python

   from pymoo.algorithms.moo.nsga2 import NSGA2
   from pymoo.factory import get_problem, get_crossover, get_mutation, get_sampling
   from pymoo.optimize import minimize
   from pymoo.visualization.scatter import Scatter

   from pymoo.operators.sampling.rnd import FloatRandomSampling

   from my_callback import MyCallback

   problem = get_problem("dtlz1")

   algorithm = NSGA2(
       pop_size=100,
       sampling=get_sampling("real_random"),
       crossover=get_crossover("real_sbx", prob=0.9, eta=20, tracing_activated=True),
       mutation=get_mutation("real_pm", prob=1/problem.n_var, eta=20, tracing_activated=True)
       )

   res = minimize(problem,
                  algorithm,
                  ('n_gen', 100),
                  seed=1)

   trace_list_ind1_gene1 = res.pop.get("T")[0, 0].get_all() #trace list of the first gene from the first individual of the population
   for trace_tuple in trace_list_ind1_gene1:
       print("traceID:", trace_tuple.traceID, "; influence factor:", trace_tuple.influenceFactor)



You can deactivate gene tracing, if you set "tracing_activated" flag in the crossover and mutation operator to False. By default they are set to True. Furthermore, for the mutation operator, you can also set the flag "accumulate_mutations", which by default is True. If you set this flag to False, each new mutation will generate a new traceID, enabeling to distinglish between different kinds of mutations. However, the runtime and memory usage of the algorithm might be significantly increased, as the trace lists of each gene will get bigger and bigger in each generation.

Currently, only float and integer representations are supported in this framework. If you use other datatypes you will probably run into errors.


