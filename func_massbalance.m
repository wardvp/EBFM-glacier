%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the mass balance and snow mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = func_massbalance(OUT,IN,C)

%% Update the mass balance
OUT.smb = IN.snow + IN.rain - OUT.runoff ...
    + OUT.moist_deposition + OUT.moist_condensation ...
    - OUT.moist_sublimation - OUT.moist_evaporation;
OUT.mbal = OUT.mbal + OUT.smb;

OUT.mbal_snow = max(OUT.mbal_snow + OUT.smb,0);
OUT.mbal_snow(all(OUT.subD==C.Dice,2),1) = 0;

end


