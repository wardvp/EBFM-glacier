%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Energy-balance + snow model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%% Model setup & initialization
[grid,time,io,phys]         = func_init_params();
[C]                         = func_init_constants();
[grid]                      = func_init_grid(grid,io);
[OUT,IN,OUTFILE]            = func_init_arrays(C,grid,io);

%% Time-loop
for t=1:time.tn
    
    %% Print time to screen
    [time] = func_printtime(t,time);

    %% Read and prepare climate input
    [IN,OUT] = func_loadclimate(C,grid,IN,t,time,OUT,io);

    % Surface energy balance model
    [OUT] = func_energybalance(C,OUT,IN,t,time,grid);
    
    %% Snow/firn model
    [OUT] = func_snowmodel(C,OUT,IN,time.dt,grid,time,phys);
        
    %% Mass balance
    [OUT] = func_massbalance(OUT,IN,C);
    
    %% Write output to files
    [io,OUTFILE] = func_writetofile(OUTFILE,io,OUT,grid,t,time,C);
end

%% Save restart-file at end run
func_createbootfile(OUT,io);
