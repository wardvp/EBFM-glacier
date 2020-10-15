%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the mass balance and snow mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = TIME_mass_balance(OUT,IN,C)

% Climatic mass balance
OUT.cmb = IN.snow + IN.rain - OUT.runoff ...
    + OUT.moist_deposition + OUT.moist_condensation ...
    - OUT.moist_sublimation - OUT.moist_evaporation;
OUT.cmb_cumulative = OUT.cmb_cumulative + OUT.cmb;

% Snow mass
OUT.snowmass = max(OUT.snowmass + OUT.cmb,0);
OUT.snowmass(all(OUT.subD==C.Dice,2),1) = 0;

end


