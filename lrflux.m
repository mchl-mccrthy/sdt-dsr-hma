% Calculate longwave radiation flux

% This function calculates the longwave radiation flux at the surface of a
% debris-covered glacier. Using approach of Arnold et al (2006) for 
% reflected longwave from terrain
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function L = lrflux(L_in,T_s,epsDeb,constsigma,T_a,f_sv)

% Calculate incoming longwave radiation from the sky
L_inSky = L_in*f_sv;

% Calculate longwave radiation from surrounding terrain. First term
% is emitted, second is reflected 
L_inTer = 0.96*constsigma*T_a^4*(1-f_sv)+...
    (1-0.96)*L_in*f_sv;

% Calculate outgoing longwave radiation from debris, where the first term
% is emitted, second is reflected
L_out = epsDeb*constsigma*T_s^4+...
    (1-epsDeb)*(L_inSky+L_inTer);

% Calculate longwave radiation flux
L = L_inSky+L_inTer-L_out; 

end