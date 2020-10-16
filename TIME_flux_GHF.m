%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBSURFACE HEAT FLUX:
%%%     equals k*dT/dz, calculated over the depth range from the surface to 
%%%     the mid-point of the 2nd subsurface layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GHF] = TIME_flux_GHF(C,Tsurf,OUT,cond)

% Effective conductivity                                                    % SOURCE: Sturm et al. (1997) [snow], Westermann et al. (2011) [soil]
GHF_k = (0.138-1.01d-3.*OUT.subD+3.233d-6.*OUT.subD.^2) ...                 
    .*(OUT.subSOIL==0) + C.soil_Kfrozen .* (OUT.subSOIL==1 & ...
    OUT.subT<C.T0) + C.soil_Kthawed .* (OUT.subSOIL==1 & OUT.subT>=C.T0); 

% Subsurface heat flux
GHF_C = (GHF_k(:,1).*OUT.subZ(:,1)+0.5*GHF_k(:,2).*OUT.subZ(:,2)) ...
    ./(OUT.subZ(:,1)+0.5.*OUT.subZ(:,2)).^2;
GHF = GHF_C(cond).*(OUT.subT(cond,2)-Tsurf(:));

end