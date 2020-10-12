%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL PARAMETERS:
%%% User-defined run parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid, time, io, phys] = INIT_parameters()

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.ts = '1-Sep-2015 0:00';                                                % Date start run
time.te = '1-Sep-2016 0:00';                                                % Date end run
time.dt = 0.125;                                                            % Timestep (days)
time.tn = round((datenum(time.te)- ...                                      % Nr. of time-steps
    datenum(time.ts))/time.dt)+1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grid parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid.utmzone = 33;                                                          % UTM zone
grid.max_subZ = 0.1;                                                        % Maximum first layer thickness (m)
grid.nl = 50;                                                               % Number of vertical layers
grid.doubledepth = 1;                                                       % Double vertical layer depth at layer grid.split (1=yes, 0=no)
grid.split = [15;25;35];                                                    % Vertical layer nr's at which layer depth doubles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model physics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys.percolation = 2;                                                       % Deep percolation scheme (1 = bucket, 2 = normal dist., 3 = linear dist., 4 = uniform dist.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input/output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
io.homedir = 'G:\Projects\Project EBFM_GLACIER';                            % Home directory
io.outdir = [io.homedir '/Output/'];                                        % Output directory
io.example_run = 1;                                                         % Run example case (no user input required)
io.readbootfile = 0;                                                        % REBOOT: read initial conditions from file (1=yes, 0=no)
io.writebootfile = 0;                                                       % REBOOT: write file for rebooting (1=yes, 0=no)  
io.bootfilein = 'boot_init.mat';                                            % REBOOT: bootfile to be read  
io.bootfileout = 'boot_final.mat';                                          % REBOOT: bootfile to be written
io.out_surface = 1;                                                         % OUTPUT: write surface variables to files (1=yes, 0=no)
io.out_subsurface = 1;                                                      % OUTPUT: write subsurface variables to files (1=yes, 0=no)
io.freqout = 8;                                                             % OUTPUT: frequency of storing output (every n-th time-step)

end

