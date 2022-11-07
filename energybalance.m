% Get energy balance residual
% 
% This function calculates the residual of the energy balance at the
% surface of a debris-covered glacier, given an estimated value of debris 
% surface temperature and debris temperatures from previous time step
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function psi = energybalance(T_s,T_sP,T_dP,nLayers,hLayer,S_in,L_in,T_a,...
    u,q_a,r,p_a,debalpha,debc,debepsilon,debk,debrho,debz_0,z_a,z_u,...
    timestep,constc_a,constc_w,constg,constk_vk,constL_v,constM_a,...
    constR,constrho_w,constsigma,constT_i,f_sv,shade,Z,Z_dash,A,A_dash)

% Calculate temperature profile given air temperature or previous
% surface temperature and previous temperature profile
tempProf = tempprof(T_s,T_sP,T_dP,nLayers,hLayer,constT_i,debk,timestep,...
    debrho,debc);

% Calculate heat fluxes
S = srflux(S_in,debalpha,shade,f_sv,Z,Z_dash,A,A_dash);
L = lrflux(L_in,T_s,debepsilon,constsigma,T_a,f_sv);
H = shfluxns(T_a,T_s,u,q_a,p_a,z_a,z_u,constM_a,constR,constg,...
    constc_a,constk_vk,debz_0);
LE = lhfluxns(T_a,T_s,u,q_a,p_a,z_a,z_u,constM_a,constR,constg,constL_v,...
    constk_vk,debz_0);
G = chflux(T_s,tempProf(1),hLayer,debk);
P = phflux(T_a,T_s,r,constrho_w,constc_w);

% Calculate psi, the debris surface energy flux, which must equal zero for
% each time step
psi = S+L+H+LE+G+P;

end
