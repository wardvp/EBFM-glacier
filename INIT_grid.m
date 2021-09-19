%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GRID:
%%% Read surface heights and glacier mask from input file(s)
%%% Calculate derived grid parameters (e.g. slope, aspect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid] = INIT_grid(grid,io)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT: Provide grid information:
%%%         - input.x: 2-D array containing UTM easting coordinates (m)
%%%         - input.y: 2-D array containing UTM northing coordinates (m)
%%%         - input.z: 2-D array containing elevation (m)
%%%         - input.mask: 2-D array containing mask (0 = no glacier, 
%%%                       1 = glacier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if io.example_run
    load([io.homedir '\Grid\dem_and_mask.mat']);
    input = grid_svalbard;
else
    % SPECIFY USER INPUT HERE!
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid.x_2D = input.x;
grid.y_2D = input.y;
grid.z_2D = input.z;

grid.Lx = size(grid.x_2D,1);
grid.Ly = size(grid.y_2D,2);

grid.mask_2D = input.mask;

[~,FY] = gradient(grid.y_2D);
if FY(1)<0
    grid.x_2D = flipud(grid.x_2D);
    grid.y_2D = flipud(grid.y_2D);    
    grid.z_2D = flipud(grid.z_2D);
    grid.mask_2D = flipud(grid.mask_2D);
end
[FX,~] = gradient(grid.x_2D);
if FX(1)<0
    grid.x_2D = fliplr(grid.x_2D);
    grid.y_2D = fliplr(grid.y_2D);    
    grid.z_2D = fliplr(grid.z_2D);
    grid.mask_2D = fliplr(grid.mask_2D);
end

grid.gpsum = sum(grid.mask_2D(:)==1);
grid.mask = grid.mask_2D(grid.mask_2D(:)==1);

[grid.lat_2D,grid.lon_2D] = INIT_grid_utm2ll(grid.x_2D,grid.y_2D,grid.utmzone);

grid.x = grid.x_2D(grid.mask_2D(:)==1); 
grid.y = grid.y_2D(grid.mask_2D(:)==1);
grid.z = grid.z_2D(grid.mask_2D(:)==1);
grid.ind = find(grid.mask_2D==1);
[grid.xind, grid.yind] = find(grid.mask_2D==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grid slope and aspect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ASPECT, SLOPE, gradN, gradE] = gradientm(grid.lat_2D,grid.lon_2D,grid.z_2D);
grid.slope = tan(SLOPE*pi/180);
grid.slope_x = gradE;
grid.slope_y = gradN;
grid.aspect = ASPECT;

grid.slope = grid.slope(grid.mask_2D(:)==1);
grid.slope_x = grid.slope_x(grid.mask_2D(:)==1);
grid.slope_y = grid.slope_y(grid.mask_2D(:)==1);
grid.aspect = grid.aspect(grid.mask_2D(:)==1);
grid.lat = grid.lat_2D(grid.mask_2D(:)==1);
grid.lon = grid.lon_2D(grid.mask_2D(:)==1);

grid.slope_beta = atan(grid.slope);
grid.slope_gamma = atan(-grid.slope_x./grid.slope_y).* (grid.slope_y>=0)...
                    + (-pi + atan(-grid.slope_x./grid.slope_y)) .* ...
                    (grid.slope_y<0 & grid.slope_x>0) + (pi + ...
                    atan(-grid.slope_x./grid.slope_y)) .* ...
                    (grid.slope_y<0 & grid.slope_x<0);
grid.slope_gamma(grid.slope_x==0 & grid.slope_y<0) = pi;
grid.slope_gamma(grid.slope_x==0 & grid.slope_y==0) = 0;
grid.slope_gamma = -grid.slope_gamma; 

end

