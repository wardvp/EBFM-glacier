function input = user_grid_input(io)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT: Provide grid information:
%%%         - input.x: 2-D array containing UTM easting coordinates (m)
%%%         - input.y: 2-D array containing UTM northing coordinates (m)
%%%         - input.z: 2-D array containing elevation (m)
%%%         - input.mask: 2-D array containing mask (0 = no glacier, 
%%%                       1 = glacier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if io.example_run
    load('G:\Projects\Project EBFM_GLACIER\Grid\dem_and_mask_fullgrid.mat');
    input = grid_svalbard;
else
    % SPECIFY USER INPUT HERE!
end

end