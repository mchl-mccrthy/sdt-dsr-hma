function [T_s,T_d,melt] = ebmodel(h,S_in,L_in,T_a,u,q_a,r,p_a,SD,...
    nLayers,debalpha,debc,debepsilon,debk,debrho,debz_0,z_a,z_u,...
    timestep,constc_a,constc_w,constg,constk_vk,constL_f,constL_v,...
    constM_a,constR,constrho_i,constrho_w,constsigma,constT_i,f_sv,...
    shade,Z,Z_dash,A,A_dash)

% This function solves the energy balance at the surface of a
% debris-covered glacier in order to calculate ice melt, debris surface 
% temperature and debris internal temperatures, broadly following the 
% approach of Reid and Brock (2010) (the DEB model).
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Get number of timesteps
nTimesteps = length(S_in);

% Get thickness of each debris sublayer
hLayer = h/nLayers;

% Preallocate space for debris surface and internal temperatures
T_s = nan(nTimesteps+1,1);
T_d = nan(nLayers+1,nTimesteps+1);

% Initial condition is linear temperature gradient between debris surface
% and ice surface
T_s(1) = T_a(1);
T_d(:,1) = linspace(T_a(1),constT_i,nLayers+1);

% Specify tolerances for Newton's method
rangeT_s = 0.5;
tolT_s = 0.01;
maxIter = 100;
maxStepT_s = 1;

% Loop through time steps solving energy balance for surface temperature. 
% Use air temperature or previous surface temperature as initial point for 
% solver. Aim is to find the surface temperature value for which the energy
% balance function f(Ts) = 0
for iTimestep = 2:nTimesteps+1
    
    % If snow depth is more than zero, surface temperature is zero
    if SD(iTimestep-1) > 0
        T_s(iTimestep) = constT_i;
        
    % Otherwise, solve for surface temperature
    else

        % Cannot pass indexed variables in function handle object, so...
        S_inP = S_in(iTimestep-1);
        T_aP = T_a(iTimestep-1);
        T_dP = T_d(2:end-1,iTimestep-1);
        T_sP = T_s(iTimestep-1);
        L_inP = L_in(iTimestep-1);
        uP = u(iTimestep-1);
        q_aP = q_a(iTimestep-1);
        rP = r(iTimestep-1);
        shadeP = shade(iTimestep-1);
        ZP = Z(iTimestep-1);
        AP = A(iTimestep-1);
        
        % Create a function handle for f(Ts)
        fun = @(Ts) energybalance(Ts,T_sP,...
            T_dP,nLayers,hLayer,S_inP,...
            L_inP,T_aP,uP,q_aP,rP,p_a,debalpha,debc,debepsilon,...
            debk,debrho,debz_0,z_a,z_u,timestep,constc_a,constc_w,...
            constg,constk_vk,constL_v,constM_a,constR,constrho_w,...
            constsigma,constT_i,f_sv,shadeP,ZP,Z_dash,AP,A_dash);
        
        % First guess of T_s for each timestep is T_a
        T_s0 = T_a(iTimestep-1);
        
        % Solve for surface temperature
        T_s(iTimestep) = newtonsmethod(fun,T_s0,tolT_s,maxIter,rangeT_s,...
            maxStepT_s);
    end
    
    % Calculate temperature profile, given surface temperature solution
    T_d(2:end-1,iTimestep) = tempprof(T_s(iTimestep),T_s(iTimestep-1),...
        T_d(2:end-1,iTimestep-1),nLayers,hLayer,constT_i,debk,timestep,...
        debrho,debc);
end

% Add surface temperature and ice temperature to temperature profiles, and
% get rid of initial condition
T_d(1,:) = T_s;
T_d(end,:) = constT_i;
T_s(1) = [];
T_d(:,1) = [];

% Calculate melt per timestep. Melt cannot be negative
Q_m = chflux(T_d(end,:),T_d(end-1,:),hLayer,debk);
melt = Q_m*timestep/(constrho_i*constL_f);
melt(melt <= 0) = 0;

end
