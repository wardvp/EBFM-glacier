%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute energy balance and solve for surface temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = func_energybalance(C,OUT,IN,t,time,grid)

%% Compute fluxes SWin, SWout, LWin, Qrain and GHF components
OUT         = func_flux_insol(grid,time,OUT);
[SWin,OUT]  = func_flux_SWin(C,OUT,IN);
[SWout,OUT] = func_flux_SWout(C,time,OUT,SWin);
LWin        = func_flux_LWin(C,IN);

%% Iterative procedure to solve the energy balance for surface temperature:
% - start by setting temperature range
% - compute surface temperature dependent energy fluxes
% - sum all fluxes to obtain the energy budget
% - use bisection method to iteratively update surface temperature
% - stop when surface temperature changes less than predefined limit

tt1 = OUT.Tsurf-40d0;
tt2 = OUT.Tsurf+40d0;
dT = tt2-tt1;
for c=1:20																
	dT=0.5*dT;
	ttmid=(tt1+tt2)/2d0;
    
    
    LWout  = func_flux_LWout(C,tt1);
    LHF    = func_flux_LHF(C,tt1,IN,true(grid.gpsum,1));
    SHF    = func_flux_SHF(C,tt1,IN,true(grid.gpsum,1));
    GHF    = func_flux_GHF(tt1,OUT,true(grid.gpsum,1));

	ebal_tt1 = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    LWout  = func_flux_LWout(C,ttmid);
    LHF    = func_flux_LHF(C,ttmid,IN,true(grid.gpsum,1));
    SHF    = func_flux_SHF(C,ttmid,IN,true(grid.gpsum,1));
    GHF    = func_flux_GHF(ttmid,OUT,true(grid.gpsum,1));
    
	ebal_ttmid = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    cond_EB = ebal_ttmid.*ebal_tt1<0;
    tt2(cond_EB) = ttmid(cond_EB);
    tt1(~cond_EB) = ttmid(~cond_EB);
  
    if (max(dT)<C.dTacc), break; end
end

%% In case the computed surface temperature is higher than zero:
% - set surface temperature to melting point
% - excess energy is used for melting (Emelt)

ttmid(abs(ttmid-273.15)<C.dTacc) = ttmid(abs(ttmid-273.15)<C.dTacc)-C.dTacc;
ttmid = min(ttmid,273.15);

LWout   = func_flux_LWout(C,ttmid);
LHF     = func_flux_LHF(C,ttmid,IN,true(grid.gpsum,1));
SHF     = func_flux_SHF(C,ttmid,IN,true(grid.gpsum,1));
GHF     = func_flux_GHF(ttmid,OUT,true(grid.gpsum,1));

Emelt   = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
Emelt(ttmid<273.15) = 0.0;

OUT.melt = 24.*3600.*time.dt.*Emelt./C.Lm./1d3;

%% Moisture exchange & melt
OUT.moist_deposition = 24.*3600.*time.dt.*LHF./C.Ls./1d3 .* (ttmid<273.15 & LHF>0);
OUT.moist_condensation = 24.*3600.*time.dt.*LHF./C.Lv./1d3 .* (ttmid>=273.15 & LHF>0);
OUT.moist_sublimation = -24.*3600.*time.dt.*LHF./C.Ls./1d3 .* (ttmid<273.15 & LHF<0);
OUT.moist_evaporation = -24.*3600.*time.dt.*LHF./C.Lv./1d3 .* (ttmid>=273.15 & LHF<0);

%% Avoid evaporation of absent melt
max_evap = OUT.melt;
OUT.moist_evaporation = min(OUT.moist_evaporation,max_evap);

%% Store local variables in OUT 
OUT.Tsurf = ttmid;
OUT.Emelt = Emelt;
OUT.LHF = LHF;
OUT.SWin = SWin;
OUT.SWout = SWout;
OUT.LWin = LWin;
OUT.LWout = LWout;
OUT.SHF = SHF;
OUT.LHF = LHF;
OUT.GHF = GHF;

end