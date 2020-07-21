Data Parameter Estimation for a Single Unit - Economizer
========================================================

This notebook demonstrates parameter estimation continuing from the data
reconciliation results in ``econ_recon.ipynb`` or
``boiler_recon.ipynb``.

1. Read Data
------------

The data used here is produced by the data reconciliation step. The data
tags are mostly systematically generated from stream names to minimize
effort in mapping data to the model. The same model is used for
parameter reconciliation and data reconciliation, so the stream names
are consistent and any data mapping can be reused. The bin information
columns where included in the data reconciliation output, so there is no
need to bin the data again here.

.. code:: ipython3

    import pandas as pd
    import idaes.dmf.model_data as da
    from idaes.logger import getLogger
    import logging
    getLogger('idaes.core').setLevel(logging.ERROR)

Either economizer only or full boiler data reconciliation results can be
used here. You can select below.

.. code:: ipython3

    recon_data = "econ_recon.csv"
    
    #Uncomment the next line to use boiler recon data.
    #redon_data = "boiler_recon.csv"

.. code:: ipython3

    # Since the data is already tagged to match the model and in the correct units, we directly read the data
    # into a Pandas data frame, and there is no need to use the data processing functions that were used in the
    # data reconciliation notebook (although they could be used here). 
    df = pd.read_csv(recon_data)
    
    # Calculate the standard deviations of the binned data 
    bin_stdev = da.bin_stdev(df, bin_no="bin_no")

Create a function to set the model data. In this case, we take the data
reconciliation results for model input to be correct and set those in
the model. Another approach would be to also estimate model inputs given
that there is some uncertainty in their measurements.

.. code:: ipython3

    # The 'set_data' function below takes data from the DataFrame and updates the model data parameters.
    def set_data(m, df, data_tags, index=None, indexindex=None):
        if index is None:
            index = df.index[indexindex]
        m.bin_no = df.iloc[index]["bin_no"]
        for t in data_tags:
            m.data[t] = df.iloc[index][t]
            m.data_stdev[t] = bin_stdev[m.bin_no][t]
        # Set the inlet streams from the data.
        input_tags = [
            "BFW_h",
            "BFW_P",
            "BFW_F",
            "FG_2_ECON_T",
            "FG_2_ECON_P",
            "FG_2_ECON_F[O2]",
            "FG_2_ECON_F[NO]",
            "FG_2_ECON_F[N2]",
            "FG_2_ECON_F[SO2]",
            "FG_2_ECON_F[CO2]",
            "FG_2_ECON_F[H2O]",
        ]
        for t in input_tags:
            m.data_tags[t].value = df.iloc[index][t]

2. Create Model Generator Function
----------------------------------

We use the Parmest tool from Pyomo to do the parameter estimation here,
which requires a function to generate a model for each case. The cases
are put together by Parmest to set up a parameter estimation problem.

.. code:: ipython3

    # Import models
    import os
    import pyomo.environ as pyo
    from idaes.core.util import model_serializer as ms
    from idaes.core import FlowsheetBlock
    from idaes.power_generation.properties.IdealProp_FlueGas import FlueGasParameterBlock
    from idaes.power_generation.unit_models.boiler_heat_exchanger import (
        BoilerHeatExchanger, 
        TubeArrangement, 
        DeltaTMethod
    )
    from idaes.generic_models.properties import iapws95
    import idaes.core.util.tables as ta

.. code:: ipython3

    # Add a function to get an instance of the economizer model.
    
    solver = pyo.SolverFactory('ipopt')
    def get_model(data=0):
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False, "time_units":pyo.units.s})
        m.fs.prop_water = iapws95.Iapws95ParameterBlock()
        m.fs.prop_fluegas = FlueGasParameterBlock()
    
        m.fs.econ = BoilerHeatExchanger(default={
                "side_1_property_package": m.fs.prop_water,
                "side_2_property_package": m.fs.prop_fluegas,
                "has_pressure_change": True,
                "has_holdup": False,
                "delta_T_method": DeltaTMethod.counterCurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False
            }
        )
        # Set inputs and initialize.  Since the initialization is repeated each time a
        # model is created, we'll save the results and reload them.
        if os.path.isfile("econ_init.json.gz"):
            ms.from_json(m, fname="econ_init.json.gz")
        else:
            h = iapws95.htpx(563.706, 2.5449e7)
            m.fs.econ.side_1_inlet.flow_mol[0].fix(24678.26) # mol/s
            m.fs.econ.side_1_inlet.enth_mol[0].fix(h) #J/mol         
            m.fs.econ.side_1_inlet.pressure[0].fix(2.5449e7) # Pa
    
            # Set the flue gas flow and composition
            fg_rate = 28.3876e3  # mol/s equivalent of ~1930.08 klb/hr
            fg_comp = { # mol fraction of flue gas components
                "H2O":8.69/100,
                "CO2":14.49/100,
                "O2":2.47/100,
                "NO":0.0006,
                "SO2":0.002,
            }
            # The rest is N2
            fg_comp["N2"] = 1 - sum(fg_comp[i] for i in fg_comp)
    
            # Set economizer inlets
            for c in fg_comp:
                m.fs.econ.side_2_inlet.flow_component[0, c].fix(fg_rate*fg_comp[c])    
            m.fs.econ.side_2_inlet.temperature[0].fix(682.335)  # K
            m.fs.econ.side_2_inlet.pressure[0].fix(100145)  # Pa
    
            # Set economizer design variables and parameters
            ITM = 0.0254  # inch to meter conversion
            # Based on NETL Baseline Report Rev4
            m.fs.econ.tube_thickness.fix(0.188*ITM)  # tube thickness
            m.fs.econ.tube_di.fix((2.0 - 2.0 * 0.188)*ITM) # calc inner diameter
            m.fs.econ.pitch_x.fix(3.5*ITM)
            m.fs.econ.pitch_y.fix(5.03*ITM)
            m.fs.econ.tube_length.fix(53.41*12*ITM)  # use tube length (53.41 ft)
            m.fs.econ.tube_nrow.fix(36*2.5) # use to match baseline performance
            m.fs.econ.tube_ncol.fix(130) # 130 from thermoflow
            m.fs.econ.nrow_inlet.fix(2)
            m.fs.econ.delta_elevation.fix(50)
            m.fs.econ.tube_r_fouling = 0.000176
            m.fs.econ.shell_r_fouling = 0.00088
            m.fs.econ.fcorrection_htc.fix(1.5)
            m.fs.econ.fcorrection_dp_tube.fix(1.0)
            m.fs.econ.fcorrection_dp_shell.fix(1.0)
            m.fs.econ.initialize(
                state_args_1={
                    "flow_mol": m.fs.econ.side_1_inlet.flow_mol[0].value,
                    "pressure": m.fs.econ.side_1_inlet.pressure[0].value,
                    "enth_mol": m.fs.econ.side_1_inlet.enth_mol[0].value,
                },
                state_args_2={
                    "flow_component":{
                        "H2O": m.fs.econ.side_2_inlet.flow_component[0, "H2O"].value,
                        "CO2": m.fs.econ.side_2_inlet.flow_component[0, "CO2"].value,
                        "N2": m.fs.econ.side_2_inlet.flow_component[0, "N2"].value,
                        "O2": m.fs.econ.side_2_inlet.flow_component[0, "O2"].value,
                        "NO": m.fs.econ.side_2_inlet.flow_component[0, "NO"].value,
                        "SO2": m.fs.econ.side_2_inlet.flow_component[0, "SO2"].value,
                    },
                    "temperature": m.fs.econ.side_2_inlet.temperature[0].value,
                    "pressure": m.fs.econ.side_2_inlet.pressure[0].value,
                }
            )
            ms.to_json(m, fname="econ_init.json.gz")
        
        #Add tags and data parameters
        stream_dict = ta.arcs_to_stream_dict(
            m, 
            additional={
                "BFW": m.fs.econ.side_1_inlet,
                "ECON_OUT": m.fs.econ.side_1_outlet,
                "FG_2_ECON": m.fs.econ.side_2_inlet,
                "FG_2_AIRPH": m.fs.econ.side_2_outlet,
            },
            sort=True,
        )
        state_dict = ta.stream_states_dict(stream_dict, time_point=0)
        m.data_tags = ta.tag_state_quantities(
            blocks=state_dict, 
            attributes=(
                "flow_mass", 
                "flow_mol", 
                "enth_mol", 
                "temperature", 
                "pressure", 
                ("flow_component", "O2"),
                ("flow_component", "NO"),
                ("flow_component", "N2"),
                ("flow_component", "SO2"),
                ("flow_component", "CO2"),
                ("flow_component", "H2O"),
            ), 
            labels=("_Fm", "_F", "_h", "_T", "_P", "_F[O2]", "_F[NO]", "_F[N2]", "_F[SO2]", "_F[CO2]", "_F[H2O]"),
        )
        m.data_tags["ECON_Q"] = m.fs.econ.heat_duty[0]
        
        m.data = pyo.Param(m.data_tags, mutable=True, doc="Process data for a specific point in time.")
        m.data_stdev = pyo.Param(m.data_tags, mutable=True, doc="Process data standard deviation.")
        @m.Expression(m.data_tags)
        def err(m, i):
            return (m.data[i] - m.data_tags[i])/m.data_stdev[i]    
    
        # Set the data
        set_data(m, df, data_tags=m.data_tags, index=data)
        solver.solve(m)
    
        return m

.. code:: ipython3

    # Try the get model function
    solver = pyo.SolverFactory('ipopt')
    print(df.index)
    m = get_model(0)
    
    # Solve the model at the first data point
    solver.solve(m)


.. parsed-literal::

    RangeIndex(start=0, stop=250, step=1)
    2020-07-21 06:19:19 [INFO] idaes.init.fs.econ.side_1: Initialization Complete
    2020-07-21 06:19:19 [INFO] idaes.init.fs.econ.side_2: Initialization Complete




.. parsed-literal::

    {'Problem': [{'Lower bound': -inf, 'Upper bound': inf, 'Number of objectives': 1, 'Number of constraints': 118, 'Number of variables': 118, 'Sense': 'unknown'}], 'Solver': [{'Status': 'ok', 'Message': 'Ipopt 3.13.2\\x3a Optimal Solution Found', 'Termination condition': 'optimal', 'Id': 0, 'Error rc': 0, 'Time': 0.14017415046691895}], 'Solution': [OrderedDict([('number of solutions', 0), ('number of solutions displayed', 0)])]}



.. code:: ipython3

    # Show the model result at the first data point
    from idaes.core.util.misc import svg_tag  # utility to place numbers/text in an SVG
    from IPython.display import SVG, display
    
    with open("econ.svg", "r") as f:
        s = svg_tag(svg=f, tags={"subtitle":"Initialized Model"})
        s = svg_tag(svg=s, tags=m.data_tags, outfile="econ_init.svg")
    display(SVG(s))



.. image:: output_12_0.svg


3. Set Up Parameter Estimation
------------------------------

Here we use the Parmest tool to solve the parameter estimation problem.
The theta\_names list is a list of parameters to estimate. The theta
names strings are the location of the parameters in the model. A
function ``sse()`` is also defined that creates the objective function
for each model instance. The objective from the individual cases is
summed to produce the overall parameter estimation objective.

.. code:: ipython3

    # List of parameters to estimate
    theta_names = [
        "fs.econ.fcorrection_htc",
        "fs.econ.fcorrection_dp_tube",
        "fs.econ.fcorrection_dp_shell",
    ]
    
    # Tags to include in the objective
    objective_tags = {
        "ECON_OUT_P", 
        "ECON_OUT_T", 
        "FG_2_AIRPH_T", 
        "FG_2_AIRPH_P",
    }
    
    # Return expressions for the objective
    def sse(model, data):
        return sum((model.err[i])**2 for i in objective_tags)

Run Parmest and record the results. Here we group the data by bin. Each
parameter in theta\_names will be estimated based on all the points in a
bin. This will allow us to examine whether the parameters have a
dependence on load.

.. code:: ipython3

    import pyomo.contrib.parmest.parmest as parmest
    import numpy as np
    
    parmest_results={}
    # run parmest tool for each power bin
    for i, group in df.groupby("bin_no"):
        pest = parmest.Estimator(get_model, list(group.index), theta_names, sse)
        obj, theta = pest.theta_est()
        print(f"Bin number: {i},  objective: {obj}")
        for k in theta:
            print(f"  {k}: {theta[k]}")
        parmest_results[i] = {'obj':obj, 'theta': theta}


.. parsed-literal::

    Bin number: 0.0,  objective: 56.836778128558095
      fs.econ.fcorrection_dp_shell: 0.041370074472566756
      fs.econ.fcorrection_dp_tube: 1.0130021553375035
      fs.econ.fcorrection_htc: 1.5352438283433547
    Bin number: 1.0,  objective: 21.514318331086947
      fs.econ.fcorrection_dp_shell: 5.095463591168366
      fs.econ.fcorrection_dp_tube: 1.023389846737705
      fs.econ.fcorrection_htc: 1.6233449672936595
    Bin number: 2.0,  objective: 19.54952337420911
      fs.econ.fcorrection_dp_shell: 4.023425622978943
      fs.econ.fcorrection_dp_tube: 1.0073324825371555
      fs.econ.fcorrection_htc: 1.4062962992653394
    Bin number: 3.0,  objective: 31.337351018147647
      fs.econ.fcorrection_dp_shell: -1.1369033273232936
      fs.econ.fcorrection_dp_tube: 1.0247444882292522
      fs.econ.fcorrection_htc: 1.491617859681654
    Bin number: 4.0,  objective: 35.21950968599074
      fs.econ.fcorrection_dp_shell: 2.9574259091643342
      fs.econ.fcorrection_dp_tube: 1.0087439453763214
      fs.econ.fcorrection_htc: 1.5268536248072353
    Bin number: 5.0,  objective: 30.91856463210748
      fs.econ.fcorrection_dp_shell: 1.4380544995158824
      fs.econ.fcorrection_dp_tube: 0.9908482353215727
      fs.econ.fcorrection_htc: 1.4567787137668406
    Bin number: 6.0,  objective: 25.175477857405607
      fs.econ.fcorrection_dp_shell: 1.6784930656455663
      fs.econ.fcorrection_dp_tube: 0.9708616867233275
      fs.econ.fcorrection_htc: 1.6105607872288856
    Bin number: 7.0,  objective: 37.12084834804889
      fs.econ.fcorrection_dp_shell: 1.6305437743245397
      fs.econ.fcorrection_dp_tube: 0.985273819266198
      fs.econ.fcorrection_htc: 1.409776694017239
    Bin number: 8.0,  objective: 30.156535947030978
      fs.econ.fcorrection_dp_shell: 1.4394598853980247
      fs.econ.fcorrection_dp_tube: 0.9951814482866638
      fs.econ.fcorrection_htc: 1.4268339124627714
    Bin number: 9.0,  objective: 19.20802998411669
      fs.econ.fcorrection_dp_shell: -4.030369420105185
      fs.econ.fcorrection_dp_tube: 1.0206238544140287
      fs.econ.fcorrection_htc: 1.4190343071600475
    Bin number: 10.0,  objective: 18.509771059164372
      fs.econ.fcorrection_dp_shell: 5.151170316917237
      fs.econ.fcorrection_dp_tube: 0.9760535814757234
      fs.econ.fcorrection_htc: 1.4599740953267906
    Bin number: 11.0,  objective: 14.443342904782563
      fs.econ.fcorrection_dp_shell: 0.18629348987212993
      fs.econ.fcorrection_dp_tube: 0.975897428007684
      fs.econ.fcorrection_htc: 1.5332937937876265
    Bin number: 12.0,  objective: 30.448947383732342
      fs.econ.fcorrection_dp_shell: 1.126538961266892
      fs.econ.fcorrection_dp_tube: 0.9871876690013094
      fs.econ.fcorrection_htc: 1.5620731361097653
    Bin number: 13.0,  objective: 29.772869119980953
      fs.econ.fcorrection_dp_shell: 0.022247628244675303
      fs.econ.fcorrection_dp_tube: 1.003799679518737
      fs.econ.fcorrection_htc: 1.4207904817280208
    Bin number: 14.0,  objective: 15.122137182360682
      fs.econ.fcorrection_dp_shell: 3.266422360203944
      fs.econ.fcorrection_dp_tube: 0.9904366373523035
      fs.econ.fcorrection_htc: 1.5063977260803123
    Bin number: 15.0,  objective: 24.573505418146837
      fs.econ.fcorrection_dp_shell: -1.3510851414408096
      fs.econ.fcorrection_dp_tube: 0.9749335407546094
      fs.econ.fcorrection_htc: 1.5359023609270999
    Bin number: 16.0,  objective: 38.8151170035052
      fs.econ.fcorrection_dp_shell: 1.7314079709652839
      fs.econ.fcorrection_dp_tube: 0.9936463239468311
      fs.econ.fcorrection_htc: 1.4424448588816015
    Bin number: 17.0,  objective: 44.077127826114065
      fs.econ.fcorrection_dp_shell: -0.28961077101303173
      fs.econ.fcorrection_dp_tube: 0.9792848546166794
      fs.econ.fcorrection_htc: 1.486249407473518
    Bin number: 18.0,  objective: 17.72456952067039
      fs.econ.fcorrection_dp_shell: 3.4876705148228955
      fs.econ.fcorrection_dp_tube: 1.032609242520667
      fs.econ.fcorrection_htc: 1.524213347039179
    Bin number: 19.0,  objective: 36.145032359093385
      fs.econ.fcorrection_dp_shell: 1.9094695066882954
      fs.econ.fcorrection_dp_tube: 0.9889213975877813
      fs.econ.fcorrection_htc: 1.39723224990267
    Bin number: 20.0,  objective: 18.596609920111813
      fs.econ.fcorrection_dp_shell: -1.7014772769788509
      fs.econ.fcorrection_dp_tube: 1.02145681107997
      fs.econ.fcorrection_htc: 1.5362982564124965
    Bin number: 21.0,  objective: 21.873052464241706
      fs.econ.fcorrection_dp_shell: -3.492852795619613
      fs.econ.fcorrection_dp_tube: 0.9686870688419849
      fs.econ.fcorrection_htc: 1.4516017885647503
    Bin number: 22.0,  objective: 17.530864638570304
      fs.econ.fcorrection_dp_shell: 0.7344647445786916
      fs.econ.fcorrection_dp_tube: 0.9928592961486812
      fs.econ.fcorrection_htc: 1.492194715617218
    Bin number: 23.0,  objective: 46.64682745383793
      fs.econ.fcorrection_dp_shell: -1.5661400011512785
      fs.econ.fcorrection_dp_tube: 0.9799353308244462
      fs.econ.fcorrection_htc: 1.5087758072503916
    Bin number: 24.0,  objective: 214.67011658875725
      fs.econ.fcorrection_dp_shell: 0.20185372095508283
      fs.econ.fcorrection_dp_tube: 1.0226643989739546
      fs.econ.fcorrection_htc: 1.438302889069996
    Bin number: 25.0,  objective: 58.698074551331516
      fs.econ.fcorrection_dp_shell: 2.9764729925057454
      fs.econ.fcorrection_dp_tube: 0.9521419368002099
      fs.econ.fcorrection_htc: 1.3955382542398527


.. code:: ipython3

    # Save results
    import json
    
    with open("econ_parmest_result.json", "w") as f:
        json.dump(parmest_results, f)

