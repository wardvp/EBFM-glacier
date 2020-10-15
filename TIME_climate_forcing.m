%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METEOROLOGICAL FORCING:
%%% - Specify / read meteorological input for current time-step, including
%%%     air temperature, precipitation, relative humidity, cloud cover,
%%%     wind-speed and air pressure
%%% - Determine derived meteorological fields, including snowfall, 
%%%     rainfall, annual accumulation, specific humidity, vapor pressure,
%%%     air density, time since last snowfall and potential temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN,OUT] = TIME_climate_forcing(C,grid,IN,t,time,OUT,io)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY/READ METEO FORCING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IN = TIME_climate_forcing_read_data(IN,io,C,time,grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DERIVED METEOROLOGICAL FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snowfall & rainfall
IN.snow = IN.P .* (IN.T < C.rainsnowT-1);
IN.rain = IN.P .* (IN.T > C.rainsnowT+1);
IN.snow = IN.snow + IN.P .* (C.rainsnowT-IN.T+1)./2 .* ...
    (IN.T < C.rainsnowT+1 & IN.T > C.rainsnowT-1);
IN.rain = IN.rain + IN.P .* (1+IN.T-C.rainsnowT)./2 .* ...
    (IN.T < C.rainsnowT+1 & IN.T > C.rainsnowT-1);

% Annual snow accumulation
OUT.ys = (1.0-1.0/(C.yeardays/time.dt)).*OUT.ys + IN.P.*1d3;
logys = log(OUT.ys);
IN.yearsnow = repmat(OUT.ys,[1 grid.nl]);
IN.logyearsnow = repmat(logys,[1 grid.nl]);

% Vapor pressure & specific humidity
VPsat = C.VP0.*exp(C.Lv/C.Rv.*(1.0./273.15-1.0./IN.T)) .* (IN.T>=273.15)...
    + C.VP0.*exp(C.Ls/C.Rv.*(1.0./273.15-1.0./IN.T)) .* (IN.T<273.15);
IN.VP = IN.RH .* VPsat;
IN.q = IN.RH .* (VPsat .* C.eps ./ IN.Pres);
    
% Air density
IN.Dair = IN.Pres./C.Rd./IN.T;

% Time since last snow fall event
OUT.timelastsnow(IN.snow/(time.dt*24*3600)>C.Pthres) = time.TCUR;
if t==1
    OUT.timelastsnow(:) = time.TCUR; 
end

% Potential temperature & lapse rate
IN.Theta = IN.T.*(C.Pref./IN.Pres).^(C.Rd/C.Cp);
IN.f_Pottemp = fit(grid.z,IN.Theta,'poly1');
IN.Theta_lapse = max(IN.f_Pottemp.p1,0.0015);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STORE RELEVANT VARIABLES IN OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.climT = IN.T;
OUT.climP = IN.P;
OUT.climC = IN.C;
OUT.climRH = IN.RH;
OUT.climWS = IN.WS;
OUT.climPres = IN.Pres;
OUT.climsnow = IN.snow;
OUT.climrain = IN.rain;

end