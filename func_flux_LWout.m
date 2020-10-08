%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute outgoing shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LWout] = func_flux_LWout(C,Tsurf)

% Blackbody emission of thermal radiation
LWout = C.boltz.*Tsurf.^4;

end