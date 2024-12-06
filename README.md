# Efficiency Simulations

*Framework for efficiency simulations for HV-CMOS sensors in HEP applications*

- Goal: Test the readout mechanism of HV-CMOS MAPS in simulation to characterise their theoretical hit rate capabilities

## Versions

### v1.0
*First efficiency simulations conducted for MightyPix1 using small set of simulation data generated for the LHCb environment*
- Chip: MightyPix1
- Data: Small set of simulation data generated for the LHCb environment (500 events)

### v1.1

*Efficiency simulations conducted for MightyPix1 and test version of MightyPix2*
- Chips: MightyPix1 and MightyPix1 with an improved readout mechanism later implemented in LF-MightyPix and planned to be implemented in MightyPix2
- Data: Large set of random Poisson distributed simulation data generated for this purpose (500 000 events)

### v1.2
*Efficiency simulations conducted for AstroPix5*
- Chip: AstroPix5
- Data: Set of random Poisson distributed simulation data generated for this purpose

## Simulation instructions

### `1_root_data_extraction.py`

Generate the wanted input data from the **root files** containing physical simulation data. This data set is very limited in statistics!

### `2_run_simulation_data.py`

Main part of the simulation! Use the `config.yaml` file to choose what type of simulation you want to run and to set all the paramters. These are partly explained in the comments of the `config.yaml` file and partly in my thesis. When you are happy with the settings you chose in the `config.yaml` file, run the simulation with `python 2_run_simulation_data.py`.

Depending on whether you chose `trafo` or `comp` as your `sim` setting, the simulation will either _transform_ the raw data into the _input file_, meaning the format required for the MightyPix1 model, or _compare_ the already present input and output (what the MightyPix1 simulation gives you back when sending in data) files.

When choosing `scifi` as your data type, the input data will be generated from the physical simulation data contained in the root files. These contain only 500 events, so careful with the limited statistics. When choosing `random` as your data type, poisson distributed data will be generated for the `rate` your specified.

After you have created your _input file_ this can be sent through the MightyPix1 simulation, which is run using the Cadence tool xcelium. You will receive an output file from that simulation, which contains the data that MightyPix saw.

Finally the _input_ and _output file_ can be _compared_, again running the `2_run_simulation_data.py` script, but setting the `sim` type in the `config.yaml` file to `comp`. You can choose to generate a log file of the analysis results by setting `mute` to `False`. This file then gives you info on how many hits the MightyPix1 simulation detected, and what's up with the wrongly detected hits. If you set the `save` or `show` option of `plots` to `True`, you can see how long it took to read the hits, and how the missing hits are distributed.

### `3a_efficiency_plots.py`

This script is still a bit messy, with it you can plot the rate versus efficiency for your _randomly created data_.

### ``3b_sensor_plots.py``

This script is probably even messier, but with it you can (in principle) plot maps of all chips over the entire tracker region. Chips need to be numbered in the conventioned used in the root files (see top), which is why it can currently really only be used when using `scifi` input data.

**_Have fun!_**
