%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute incoming longwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LWin] = func_flux_LWin(C,IN)

% Clear-sky emissivity
ecs = 0.23 + C.b.*(IN.VP./IN.T).^0.125;

% Sky emissivity
e = ecs.*(1.0-IN.C.^C.p) + C.ecl.*IN.C.^C.p;
		
% Incoming longwave radiation
LWin = e .* C.boltz .* IN.T.^4;	

end