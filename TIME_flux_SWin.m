%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INCOMING SHORTWAVE RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SWin,OUT] = TIME_flux_SWin(C,OUT,IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Direct, diffuse and total radiation after shading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.TOAdir = (0.2+0.65.*(1-IN.C)).*(1-OUT.shade).*OUT.TOA;
OUT.TOAdiff = (0.8-0.65.*(1-IN.C)).*OUT.TOA;
OUT.TOAshade = OUT.TOAdir + OUT.TOAdiff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Atmospheric transmissivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmissivity after gaseous absorption / scattering
m = 35.*(IN.Pres./C.Pref).*(1224.*(OUT.TOAflat./OUT.I0).^2+1).^(-0.5);
t_rg = 1.021-0.084.*sqrt(m.*(949d0*(IN.Pres./1d3).*1d-5+0.051));

% Transmissivity after water vapor absorption
tempdew_kelvin = (1/273.15-(C.Rv/C.Ls).*log(IN.q.*IN.Pres ...
    ./(C.eps*C.VP0))).^(-1);		
tempdew_fahr = 32.0+1.8.*(tempdew_kelvin-273.15d0);						
u = exp(0.1133-log(C.lambda+1d0)+0.0393.*tempdew_fahr);
t_w = 1-0.077.*(u.*m).^0.3;

% Transmissivity after aerosol absorption
t_a = C.k_aer.^m; 

% Transmissivity after cloud absorption / scattering
t_cl = 1.0-0.128.*IN.C-0.346.*IN.C.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incoming solar radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWin = OUT.TOAshade .* t_rg .* t_w .* t_a .* t_cl;

end

