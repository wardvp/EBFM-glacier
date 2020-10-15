%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT: Provide a climate forcing:
%%%     - IN.T: grid vector containing temperature at current time
%%%     - IN.P: grid vector containing precipitation at current time
%%%     - IN.C: grid vector containing cloud cover at current time
%%%     - IN.RH: grid vector containing relative humidity at current time
%%%     - IN.WS: grid vector containing wind speed at current time
%%%     - IN.Pres: grid vector containing pressure at current time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IN = TIME_climate_forcing_read_data(IN,io,C,time,grid)

if ~io.example_run
    % SPECIFY USER INPUT HERE!

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXAMPLE: Random weather conditions
    % - sinuisoidal seasonal temperature cycle (no diurnal cycle)
    % - one precipitation event per week
    % - alternating weeks with high and low cloud fraction and rel humidity
    % - random with speed between 0 and 10 m s-1
    % - fixed sea level pressure with elevation decline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Air temperature (K)
    yearfrac = day(time.TCUR_DT,'dayofyear')/C.yeardays;
    T_amplitude = 10;                                                       % Seasonal temperature amplitude (K)
    T_mean_sea_level = 268;                                                 % Mean sea level temperature (K)
    T_lapse_rate = -0.005;                                                  % Temperature lapse rate (K m-1)
    IN.T(:) = T_mean_sea_level + T_amplitude*sin(2*pi*yearfrac-0.65*pi);    % Apply sinusoidal seasonal temperature cycle
    IN.T(:) = IN.T(:) + T_lapse_rate*grid.z(:);                             % Apply elevation lapse rate
    
    % Precipitation (m w.e.)
    P_annual_sea_level = 0.5;                                               % Annual precipitation at sea level (m w.e.)
    P_z_gradient = 0.1;                                                     % Precipitation - elevation gradient (% m-1)
    if week(time.TCUR_DT)~=week(time.TPREV_DT)
        IN.P(:) = P_annual_sea_level/52 * (1 + P_z_gradient* ...
            grid.z(:)/100);
    else
        IN.P(:) = 0;
    end
    
    % Cloud cover (fraction)
    if mod(week(time.TCUR_DT),2)==0                                         % Alternating overcast and clear-sky weeks
        IN.C(:) = 1;
    else
        IN.C(:) = 0;
    end

    % Relative humidity (fraction)
    if mod(week(time.TCUR_DT),2)==0                                         % Alternating high- and low-humidity weeks
        IN.RH(:) = 0.8;
    else
        IN.RH(:) = 0.5;
    end
    
    % Wind speed (m s-1)
    max_WS = 10;                                                            % Maximum wind speed (m s-1)
    IN.WS(:) = max_WS*rand;                                                 % Random wind speed between 0-max_WS
    
    % Air pressure (Pa)
    Pres_sea_level = 1015d2;                                                % Sea level pressure (Pa)
    IN.Pres(:) = Pres_sea_level*exp(-1.244d-4*grid.z);                      % Fixed air pressure decreasing with elevation
end

end
