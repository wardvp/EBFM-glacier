%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SURFACE ENERGY BALANCE MODEL:
%%% - Calculates incoming & reflected solar radiation [SWin & SWout], 
%%%     incoming & outgoing longwave radiation[LWin & LWout], sensible heat
%%%     flux [SHF], latent heat flux [LHF] and the ground heat flux [GHF] 
%%% - Solves the surface energy balance to determine the surface 
%%%     temperature [Tsurf] and energy involved in melting [Emelt]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = TIME_energy_balance_model(C,OUT,IN,time,grid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOLVE THE SURFACE ENERGY BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute SWin, SWout, LWin and GHF (not surface temperature dependent)
OUT         = TIME_flux_insol(grid,time,OUT);
[SWin,OUT]  = TIME_flux_SWin(C,OUT,IN);
[SWout,OUT] = TIME_flux_SWout(C,time,OUT,SWin);
LWin        = TIME_flux_LWin(C,IN);

% Set initial temperature range
Tlow = OUT.Tsurf-40d0;
Thigh = OUT.Tsurf+40d0;
dT = Thigh-Tlow;
for c=1:20																
	dT=0.5*dT;
	Tmid=(Tlow+Thigh)/2d0;
    
    % Compute surface energy balance for Tsurf = Tlow
    LWout  = TIME_flux_LWout(C,Tlow);
    LHF    = TIME_flux_LHF(C,Tlow,IN,true(grid.gpsum,1));
    SHF    = TIME_flux_SHF(C,Tlow,IN,true(grid.gpsum,1));
    GHF    = TIME_flux_GHF(Tlow,OUT,true(grid.gpsum,1));

	ebal_Tlow = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    % Compute surface energy balance for Tsurf = (Thigh+Tlow)/2
    LWout  = TIME_flux_LWout(C,Tmid);
    LHF    = TIME_flux_LHF(C,Tmid,IN,true(grid.gpsum,1));
    SHF    = TIME_flux_SHF(C,Tmid,IN,true(grid.gpsum,1));
    GHF    = TIME_flux_GHF(Tmid,OUT,true(grid.gpsum,1));
    
	ebal_Tmid = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    % If Ebal(Tmid) and Ebal (Tlow) are of opposite sign: Thigh = Tmid, 
    %   else Tlow = Tmid 
    cond_EB = ebal_Tmid.*ebal_Tlow<0;
    Thigh(cond_EB) = Tmid(cond_EB);
    Tlow(~cond_EB) = Tmid(~cond_EB);
    
    % Stop when surface temperature changes less than predefined limit
    if (max(dT)<C.dTacc)
        break; 
    end
    
    if c==20
        error('Energy balance did not converge below limit C.dTacc');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SURFACE MELT:
%%% - In case the computed surface temperature is higher than zero: 
%%%     - set surface temperature to melting point
%%%     - excess energy is used for melting (Emelt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmid(abs(Tmid-C.T0)<C.dTacc) = Tmid(abs(Tmid-C.T0)<C.dTacc)-C.dTacc;
Tmid = min(Tmid,C.T0);

LWout   = TIME_flux_LWout(C,Tmid);
LHF     = TIME_flux_LHF(C,Tmid,IN,true(grid.gpsum,1));
SHF     = TIME_flux_SHF(C,Tmid,IN,true(grid.gpsum,1));
GHF     = TIME_flux_GHF(Tmid,OUT,true(grid.gpsum,1));

Emelt   = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
Emelt(Tmid<C.T0) = 0.0;

OUT.melt = C.dayseconds.*time.dt.*Emelt./C.Lm./1d3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOISTURE FLUXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.moist_deposition = C.dayseconds.*time.dt.*LHF./C.Ls./1d3 .* ...
    (Tmid<C.T0 & LHF>0);
OUT.moist_condensation = C.dayseconds.*time.dt.*LHF./C.Lv./1d3 .* ...
    (Tmid>=C.T0 & LHF>0);
OUT.moist_sublimation = -C.dayseconds.*time.dt.*LHF./C.Ls./1d3 .* ...
    (Tmid<C.T0 & LHF<0);
OUT.moist_evaporation = -C.dayseconds.*time.dt.*LHF./C.Lv./1d3 .* ...
    (Tmid>=C.T0 & LHF<0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AVOID EVAPORATION OF ABSENT MELT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_evap = OUT.melt;
OUT.moist_evaporation = min(OUT.moist_evaporation,max_evap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STORE RELEVANT VARIABLES IN OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.Tsurf = Tmid;
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