%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REFLECTED SHORTWAVE RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SWout,OUT] = TIME_flux_SWout(C,time,OUT,SWin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Albedo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ice_cond = OUT.subD(:,1)==C.Dice | OUT.snowmass==0;
snow_cond = OUT.subD(:,1)<C.Dice & OUT.snowmass>0;

% For a snow surface
OUT.tstar(OUT.Tsurf==C.T0 & snow_cond) = C.tstar_wet;
OUT.tstar(OUT.Tsurf<C.T0 & snow_cond) = C.tstar_dry + min(C.T0-...
    OUT.Tsurf(OUT.Tsurf<C.T0 & snow_cond),10)*C.tstar_K;
OUT.alb_snow(snow_cond & OUT.timelastsnow<time.TCUR) = OUT.alb_snow( ...
    snow_cond & OUT.timelastsnow<time.TCUR) - max(OUT.alb_snow(snow_cond...
    & OUT.timelastsnow<time.TCUR)-C.alb_firn,0.0)./OUT.tstar(snow_cond ...
    & OUT.timelastsnow<time.TCUR)*time.dt;
OUT.alb_snow(OUT.timelastsnow==time.TCUR) = C.alb_fresh;
OUT.alb_snow(ice_cond) = C.alb_ice;
OUT.alb(snow_cond,1) = OUT.alb_snow(snow_cond,1);

% For an ice surface
OUT.alb(ice_cond,1) = C.alb_ice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reflected SW radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWout = SWin .* OUT.alb;

end

