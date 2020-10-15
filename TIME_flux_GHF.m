%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBSURFACE HEAT FLUX:
%%%     equals k*dT/dz, calculated over the depth range from the surface to 
%%%     the mid-point of the 2nd subsurface layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GHF] = TIME_flux_GHF(Tsurf,OUT,cond)

GHF_k       = (0.138-1.01d-3.*OUT.subD+3.233d-6.*OUT.subD.^2);              % effective conductivity (SOURCE: Sturm et al. 1997)
GHF_C       = (GHF_k(:,1).*OUT.subZ(:,1)+0.5*GHF_k(:,2).*OUT.subZ(:,2)) ...
                ./(OUT.subZ(:,1)+0.5.*OUT.subZ(:,2)).^2;

GHF = GHF_C(cond).*(OUT.subT(cond,2)-Tsurf(:));

end