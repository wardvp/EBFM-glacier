%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute subsurface heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GHF] = func_flux_GHF(Tsurf,OUT,cond)

GHF_k       = (0.138-1.01d-3.*OUT.subD+3.233d-6.*OUT.subD.^2);
GHF_C       = (GHF_k(:,1).*OUT.subZ(:,1)+GHF_k(:,2).*OUT.subZ(:,2)+0.5.*GHF_k(:,3).*OUT.subZ(:,3)) ...
                ./(OUT.subZ(:,1)+OUT.subZ(:,2)+0.5.*OUT.subZ(:,3)).^2;

GHF = GHF_C(cond).*(OUT.subT(cond,3)-Tsurf(:));

end