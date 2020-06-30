.. _x_workshops:

IDAES Workshop Materials
========================
Examples presented at technical workshops.

- :doc:`Python and Pyomo Introduction (short) <Introduction_to_Python_and_Pyomo/Introduction_Short_Solution_doc>`

  - In this module we introduce the fundamentals of working with the IDAES process modeling toolset,
    and we will demonstrate how these tools can be applied for optimization applications.

- :doc:`Python and Pyomo Introduction (full) <Introduction_to_Python_and_Pyomo/Introduction_Short_Solution_doc>`

  - In this module we introduce the fundamentals of working with the IDAES process modeling toolset,
    and we will demonstrate how these tools can be applied for optimization applications. Contains some additional
    Python examples, compared to the "short" version above.

- :doc:`Flash Unit Model <Flash_Unit_Model/Flash_Unit_Solution_doc>`

  - In this module, we will familiarize ourselves with the IDAES framework by creating and working with a flowsheet
    that contains a single flash tank. The flash tank will be used to perform separation of Benzene and Toluene.

- :doc:`HDA Flowsheet Simulation and Optimization <HDA_Flowsheet_Optimization/HDA_Flowsheet_Solution_doc>`

  - This workshop includes:

    - Constructing a steady-state flowsheet using the IDAES unit model library
    - Connecting unit models in a flowsheet using Arcs
    - Using the SequentialDecomposition tool to initialize a flowsheet with recycle
    - Formulating and solving an optimization problem

- Parameter estimation

  - In this module, we use Pyomo's `parmest` tool in conjunction with IDAES models for parameter estimation.
    We demonstrate these tools by estimating the parameters associated with the NRTL property model for a
    benzene-toluene mixture. We demonstrate two different approaches to solving this problem:

    - :doc:`Using the NRTL State Block <Parameter_Estimation/Parameter_estimation_NRTL_using_state_block_solution_testing_doc>`
    - :doc:`Using unit models <Parameter_Estimation/Parameter_estimation_NRTL_using_unit_model_solution_testing_doc>`

.. toctree::
    :hidden:
    :maxdepth: 1

    Introduction_to_Python_and_Pyomo/Introduction_Short_Solution_doc
    Introduction_to_Python_and_Pyomo/Introduction_Solution_doc
    Flash_Unit_Model/Flash_Unit_Solution_doc
    HDA_Flowsheet_Optimization/HDA_Flowsheet_Solution_doc
    Parameter_Estimation/Parameter_estimation_NRTL_using_state_block_solution_testing_doc
    Parameter_Estimation/Parameter_estimation_NRTL_using_unit_model_solution_testing_doc

