% Calculate shortwave radiation flux

% This function calculates the shortwave radiation flux at the surface of a
% debris-covered glacier. It uses equations from Arnold et al (2006), but 
% with solar zenith angle as input instead of solar elevation angle, e.g. 
% Oke (1987)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function S = srflux(S_in,debalpha,shade,f_sv,Z,Z_dash,A,A_dash)

% Daytime incoming solar radiation cannot be greater than I0*cosZ 
% (Oke, 1987; p340). This is an important correction. If it is not 
% performed, very large incoming shortwave radiation values can result
S_in = min(S_in,max(0,1368*cosd(Z)));

% Calculate diffuse fraction of input shortwave radiation
S_in_DIF = 0.15*S_in;

% Calculate direct fraction of input shortwave radiation 
S_in_DIR = S_in-S_in_DIF;

% Calculate direct shortwave radiation reaching a surface normal to the
% solar beam
S_in_bDIR = S_in_DIR/(cosd(Z));

% Calculate diffuse shortwave radiation at the grid cell that is reflected 
% from nearby terrain
S_in_gTER = 0.25*(1-f_sv)*S_in;

% Calculate diffuse shortwave radiation at the grid cell
S_in_gDIF = f_sv*S_in_DIF+S_in_gTER;

% Calculate direct shortwave radiation at the grid cell, accounting for
% shading and topography
if ~shade % shade = 0 is 'in shade', shade = 1 is 'in sun'
    S_in_gDIR = 0;
else
    S_in_gDIR = S_in_bDIR*(cosd(Z)*cosd(Z_dash)+sind(Z)*sind(Z_dash)*...
        cosd(A-A_dash));
    
    % Account for self-shading
    if S_in_gDIR < 0
        S_in_gDIR = 0;
    end
end

% Calculate shortwave radiation flux
S = (1-debalpha)*(S_in_gDIR+S_in_gDIF);

end