%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INCOMING SHORTWAVE RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SWin,OUT] = TIME_flux_SWin(C,OUT,IN,grid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Direct, diffuse and total radiation after shading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.TOAdir = (0.2+0.65.*(1-IN.C)).*(1-OUT.shade).*OUT.TOA;                  % SOURCE: Oerlemans (1992)
OUT.TOAdiff = (0.8-0.65.*(1-IN.C)).*OUT.TOA;
OUT.TOAshade = OUT.TOAdir + OUT.TOAdiff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Atmospheric transmissivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmissivity after gaseous absorption / scattering
m = 35.*(IN.Pres./C.Pref).*(1224.*(OUT.TOAflat./OUT.I0).^2+1).^(-0.5);      % SOURCE: Meyers and Dale (1983)
t_rg = 1.021-0.084.*sqrt(m.*(949d0*(IN.Pres./1d3).*1d-5+0.051));            % SOURCE: Atwater and Brown Jr (1974) 

% Transmissivity after water vapor absorption
tempdew_kelvin = (1/273.15-(C.Rv/C.Ls).*log(IN.q.*IN.Pres ...               % SOURCE: McDonald (1960)
    ./(C.eps*C.VP0))).^(-1);		
tempdew_fahr = 32.0+1.8.*(tempdew_kelvin-273.15d0);
if abs(mean(grid.lat))<20                                              % SOURCE: Smith (1966)
    lambda = 2.91;
elseif abs(mean(grid.lat))<30
    lambda = 3.12;
elseif abs(mean(grid.lat))<40
    lambda = 3.00;
elseif abs(mean(grid.lat))<50
    lambda = 2.78;
elseif abs(mean(grid.lat))<60
    lambda = 2.79;
elseif abs(mean(grid.lat))<70
    lambda = 2.41;
elseif abs(mean(grid.lat))<80
    lambda = 2.03;
else
    lambda = 1.62;
end
u = exp(0.1133-log(lambda+1d0)+0.0393.*tempdew_fahr);
t_w = 1-0.077.*(u.*m).^0.3;                                                 

% Transmissivity after aerosol absorption
t_a = C.k_aer.^m;                                                           % SOURCE: Houghton (1954)

% Transmissivity after cloud absorption / scattering
t_cl = 1.0-0.128.*IN.C-0.346.*IN.C.^2;                                      % SOURCE: Van Pelt et al. (2012)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incoming solar radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWin = OUT.TOAshade .* t_rg .* t_w .* t_a .* t_cl;

end

