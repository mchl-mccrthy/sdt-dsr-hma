% Calculate conductive heat flux

% This function calculates the conductive heat flux through a layer of
% material 
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function G = chflux(T_s,T_x,hLayer,debk)

G = -debk*(T_s-T_x)/hLayer;

end