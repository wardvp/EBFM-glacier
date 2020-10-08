%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute sensible heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SHF] = func_flux_SHF(C,Tsurf,IN,cond)

% Sensible heat flux (bulk equations)
C_kat = max(C.k_turb .* (IN.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* IN.Theta_lapse .* C.Pr)),0);
C_turb = 0.5*(C.turb + C_kat);
SHF = IN.Dair(cond).*C.Cp.*C_turb.*(IN.T(cond)-Tsurf);

end