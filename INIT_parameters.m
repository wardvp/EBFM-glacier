%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL PARAMETERS:
%%% User-defined run parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid, time, io, phys] = INIT_parameters()

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.ts = '01-Sep-2015 00:00';                                              % Date and time start run
time.te = '01-Sep-2016 00:00';                                              % Date and time end run
time.dt = 0.125;                                                            % Timestep (days)
time.tn = round((datenum(time.te)- ...                                      % Nr. of time-steps
    datenum(time.ts))/time.dt)+1;
time.dT_UTC = 1;                                                            % time difference relative to UTC (hours)          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grid parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid.utmzone = 33;                                                          % UTM zone
grid.max_subZ = 0.1;                                                        % Maximum first layer thickness (m)
grid.nl = 50;                                                               % Number of vertical layers
grid.doubledepth = 1;                                                       % Double vertical layer depth at layer grid.split (1=yes, 0=no)
grid.split = [15;25;35];                                                    % Vertical layer nr's at which layer depth doubles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys.percolation = 'normal';                                                % Water percolation scheme:
                                                                            %   - 'bucket': tipping-bucket method (all water added at the surface)
                                                                            %   - 'normal': normally distributed deep percolation
                                                                            %   - 'linear': linearly distributed deep percolation
                                                                            %   - 'uniform': uniformly distributed deep percolation
phys.snow_compaction = 'firn_only';                                         % Snow and firn compaction scheme:
                                                                            %   - 'firn_only': apply Ligtenberg et al. (2011) for all snow and firn layers
                                                                            %   - 'firn+snow': apply Ligtenberg et al. (2011) for firn layers and Kampenhout et al. (2017) for seasonal snow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input/output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
io.homedir = pwd;                                                           % Home directory
io.outdir = [io.homedir '/Output/Testrun/'];                                 % Output directory
io.rebootdir = [io.homedir '/Reboot/'];                                     % Restart file directory
io.example_run = 1;                                                         % Run example case (no user input required)
io.readbootfile = 0;                                                        % REBOOT: read initial conditions from file (1=yes, 0=no)
io.writebootfile = 1;                                                       % REBOOT: write file for rebooting (1=yes, 0=no)  
io.bootfilein = 'boot_init.mat';                                            % REBOOT: bootfile to be read  
io.bootfileout = 'boot_final.mat';                                          % REBOOT: bootfile to be written
io.freqout = 1;                                                             % OUTPUT: frequency of storing output (every n-th time-step)

end

