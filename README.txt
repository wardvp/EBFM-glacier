!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! EBFM QUICK START GUIDE !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Here is a quick guide on how to get started with EBFM in MATLAB. 

!!!!!!!!!!!!!!!!!!!!
!! BEFORE RUNNING !!
!!!!!!!!!!!!!!!!!!!!

Before running the model for one or several glaciers some preparations are needed. First of all, there are some modelling 
choices to make. Additionally, in order to run, the model requires information about the grid and meteorological
 conditions, and it needs to know which output variables to store and whether or not to read/write a restart file. The
 steps below describe how to get started with a new application of the model. 

1) Modelling choices:
The standard model physics and model parameter setup are described in Van Pelt et al. (2012, 2019, 2021) and references
therein. In INIT_parameters() choices can be made regarding the deep water percolation scheme (to use a 'bucket',
'normal', 'linear' or 'uniform' scheme; see Marchenko et al. 2017 for more details) and the snow/firn compaction scheme
(to use only firn physics 'firn-only' or to have a complementary seasonal snow model 'firn+snow'; see Kampenhout et al.
2017 for more details). Note that Van Pelt et al. (2019, 2021) used a 'normal' percolation scheme and a 'firn-only' 
compaction model.

In INIT_parameters() you will further find time, grid, model and input/output parameters that can be changed to your
needs. Explanations of the different parameters are given in the inline comments. Some notes: 
- A vertical grid should be specified by giving the number of vertical layers (grid.nl) and the maximum vertical spacing
  (grid.max_subZ). In case doubling of vertical layer thickness (by merging neighboring layers) at certain depths is 
  preferred, then grid.doubledepth should be set to 1 and it should be specified in grid.split around which layers this
  occurs. 
- The input/ouput parameter io.example_run is by default set to 1 to enable running an example application of the model
  for glaciers in Svalbard with arbitrary weather conditions. When setting up and running a new application of the model
  io.example_run should be set to 0. 

2) Creating a model grid:
EBFM performs calculations on a three-dimensional grid. The model therefore needs to know horizontal coordinates (input.x 
and input.y) and elevation (input.z) of all grid cells, which are provided by the user as input. Additionally, a mask 
(input.mask) indicating which cells are included in the simulation needs to be specified. The horizontal coordinates, 
elevation and mask are read in INIT_grid(). The horizontal coordinates should be in UTM coordinates (unit m) and the 
corresponding UTM zone should be specified in INIT_parameters() with grid.utmzone. The mask (input.mask) defines which 
grid cells are modelled (1) or not modelled (0). For example, you may want to set the mask for all non-glacier grid cells 
to 0. IMPORTANT: the model can only work with distributed grids that are exactly north-south oriented (i.e. not rotated), 
and have a regular grid spacing!

3) Specifying weather conditions:
In order to solve the energy-balance equation, the model needs to know current near-surface weather conditions. The 
required meteorological input fields include air temperature, precipitation, wind speed, relative humidity, cloud cover 
and air pressure. Data sources for meteorological conditions can be output from a regional climate model, a reanalysis 
dataset (e.g. ERA-40/ERA-INTERIM) or observational time-series from nearby weather stations. Regardless of the source, 
some post-processing of the meteorological input data is needed in order to generate a climate forcing on the model grid. 
Reading and processing of meteorological data is done in TIME_climate_forcing(). For example, if you want to use air 
temperature from a regional climate model in the model, you will need to regrid the raw data to the model grid. If the 
model grid has a finer resolution, interpolation is needed (e.g. using scatteredInterp), and it could be useful to 
additionally apply a temperature lapse rate to account for local topographic effects on temperature. Similarly, you may 
want to account for local topographic influences on precipitation and air pressure.   

The weather parameters temperature, relative humidity, precipitation, cloud cover, wind speed and air pressure are 
defined every time-step in the function TIME_climate_forcing_read_data() as the model variables input.T, input.RH, 
input.P, input.C, input.WS and input.Pres respectively. Note that these (and other) spatially distributed variables are 
(temporarily) stored as 1-D vectors rather than 2-D arrays when running the model. This significantly reduces numerical 
cost when performing calculations with the variables.

4) Saving model output to files:
At the end of every model time-step a selection of model variables is saved to output files. Output is usually saved to 
binary files (.bin) with one file per variable. Output files are saved to the directory specified in INIT_parameters() as 
io.outdir and with a frequency io.freq. Standard the whole list of output variables specified between lines 12-44 in 
TIME_write_to_file() is saved to output files, but the number of output files to be saved can be reduced by commenting 
out lines (using %) with variables to be excluded from saving. At the end of every model run, an additional 'run file' 
with time, grid, parameter and input/output settings for that run is saved to 'runinfo.mat'. 

5) Creating/reading a restart file:
It can be useful to save a restart file at the end of run, in case you want to use the final state of the model as the 
starting point for a next run. This is for example useful when performing one or several spinup runs. A restart file is 
generated by setting io.writebootfile (in INIT_parameters()) to 1 and a restart file with the name specified in 
io.bootfileout will be generated at the end of the run in the directory 'reboot'. A restart file is read when 
io.readbootfile is set to 1, in which case the restart file with the name specified in io.bootfilein is read.
 
!!!!!!!!!!!!!!!!!!!!!!!
!! RUNNING THE MODEL !!
!!!!!!!!!!!!!!!!!!!!!!!

Running the model is simply done by running the MAIN.m script in MATLAB. During the run, model timesteps will be 
displayed on the command line. In case you want to do an example run with the model without specifying any of the 
required user input described above, you can set io.example_run to 1 in INIT_parameters() and run the model. In that case 
the model will perform a simulation of all glaciers in Svalbard with arbitrary meteorological conditions (see 
TIME_climate_forcing_read_data for more details). Do not forget to set io.example_run to 0 when preparing your own model 
experiment.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VISUALIZING THE OUTPUT !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MATLAB provides an ideal platform to make figures of the output, e.g. to create timeseries plots, show distributed maps 
or to plot time-evolution of subsurface conditions. Example scripts for plotting are included in the directory 
'Visualization' and show how to make a time-series plot (Figures_time.m), a distributed map (Figures_distributed.m) and 
subsurface evolution plot (Figures_subsurface.m). For more detailed instructions, see the scripts and the inline comments.
 
!!!!!!!!!!!!!!!!
!! REFERENCES !!
!!!!!!!!!!!!!!!!

Marchenko, S., W.J.J. van Pelt, B. Claremar, H. Machguth, C.H. Reijmer, R. Pettersson and V.A. Pohjola (2017). 
Parameterizing deep water percolation improves subsurface temperature simulations by a multilayer firn model. Frontiers 
in Earth Science: Cryospheric Sciences, 5, 16. https://doi.org/10.3389/feart.2017.00016

Van Pelt, W.J.J., T.V. Schuler, V.A. Pohjola, and R. Pettersson (2021). Accelerating future mass loss of Svalbard 
glaciers from a multi-model ensemble. Journal of Glaciology, 67(263), 485-499. https://doi.org/10.1017/jog.2021.2

Van Pelt, W.J.J., V.A. Pohjola, R. Pettersson, S. Marchenko, J. Kohler, B. Luks, J.O. Hagen, T.V. Schuler, T. Dunse, B. 
Noël, and C.H. Reijmer (2019). A long-term dataset of climatic mass balance, snow conditions and runoff in Svalbard 
(1957–2018). The Cryosphere, 13, 2259-2280. https://doi.org/10.5194/tc-13-2259-2019

Van Pelt, W.J.J., J. Oerlemans, C.H. Reijmer, V.A. Pohjola, R. Pettersson and J.H. van Angelen (2012). Simulating melt, 
runoff and refreezing on Nordenskiöldbreen, Svalbard, using a coupled snow and energy balance model. The Cryosphere, 6, 
641-659. https://doi.org/10.5194/tc-6-641-2012


 
