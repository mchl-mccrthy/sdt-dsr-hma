% Calculate heat flux due to precipitation

% This function calculates the heat flux due to precipitation at the
% surface of a debris-covered glacier
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function P = phflux(T_a,T_s,r,constrho_w,constc_w)

P = constrho_w*r*constc_w*(T_a-T_s);

end