% Calculate latent heat flux using Monin-Obukhov approach to get bulk
% transfer coefficient

% Notes
% - Script by Mike McCarthy (2019)

function LE = lhfluxns(T_a,T_s,u,q_a,p_a,z_a,z_u,constM_a,constR,constg,...
    constL_v,constk_vk,debz_0)

% Calculate the specific humidity at the debris surface
q_s = q_a*T_s/T_a;

% Calculate air density
rho_a = p_a*constM_a/(constR*T_a);

% Calculate bulk transfer coefficient
C_bt = constk_vk^2/(log(z_u/debz_0)*log(z_a/debz_0));

% Calculate latent heat flux
LE = rho_a*constL_v*u*(q_a-q_s)*C_bt;

end
