%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TURBULENT SENSIBLE HEAT FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SHF] = TIME_flux_SHF(C,Tsurf,IN,cond)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Turbulent exchange coefficient [Klok & Oerlemans, 2002]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_kat = max(C.k_turb .* (IN.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* ...        % SOURCE: Klok and Oerlemans (2002)
    IN.Theta_lapse .* C.Pr)),0);
C_turb = 0.5*(C.turb + C_kat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sensible heat flux (bulk equations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHF = IN.Dair(cond).*C.Cp.*C_turb.*(IN.T(cond)-Tsurf);                      % SOURCE: Klok and Oerlemans (2002)

end