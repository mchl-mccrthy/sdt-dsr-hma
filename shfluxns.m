% Calculate sensible heat flux

% This function calculates the sensible heat flux at the surface of a
% debris-covered glacier, assuming neutral atmospheric stability
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function H = shfluxns(T_a,T_s,u,q_a,p_a,z_a,z_u,constM_a,constR,...
    constg,constc_a,constk_vk,debz_0)

% Calculate air density based on temperature
rho_a = p_a*constM_a/(constR*T_a);

% Calculate specific heat capacity of air
c_p = constc_a*(1+0.84*q_a);

% Calculate bulk transfer coefficient
C_bt = constk_vk^2/(log(z_u/debz_0)*log(z_a/debz_0));

% Calculate sensible heat flux
H = rho_a*c_p*u*(T_a-T_s)*C_bt;

end