%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SURFACE ENERGY BALANCE - SNOW & FIRN MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model setup & initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[grid,time,io,phys]         = INIT_parameters();
[C]                         = INIT_constants();
[grid]                      = INIT_grid(grid,io);
[OUT,IN,OUTFILE]            = INIT_initial_conditions(C,grid,io);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time-loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:time.tn
    
    % Print time to screen
    [time] = TIME_print_time(t,time);

    % Read and prepare climate input
    [IN,OUT] = TIME_climate_forcing(C,grid,IN,t,time,OUT,io);

    % Surface energy balance model
    [OUT] = TIME_energy_balance_model(C,OUT,IN,time,grid);
    
    % Snow & firn model
    [OUT] = TIME_snow_model(C,OUT,IN,time.dt,grid,phys);
        
    % Mass balance
    [OUT] = TIME_mass_balance(OUT,IN,C);
    
    % Write output to files
    [io,OUTFILE] = TIME_write_to_file(OUTFILE,io,OUT,grid,t,time,C);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save restart-file at end run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FINAL_create_boot_file(OUT,io);
