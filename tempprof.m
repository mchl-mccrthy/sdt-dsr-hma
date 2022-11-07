% Calculate debris temperature profile

% This function calculates the temperature profile through a debris layer
% when debris surface temperature and ice surface temperature are known,
% implementing equations A8 to A12 of Reid and Brock (2010)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function tempProf = tempprof(Ts_t1,Ts_t,Td,nLayers,h,constT_i,debk,timestep,...
    debrho,debc)

    % Calculate C
    C = debk*timestep/(2*debrho*debc*h^2);
    
    % Calculate a, b and c
    a = C;
    b = 2*C+1;
    c = C;
    
    % Calculate d
    d(1) = C*Ts_t1+C*Ts_t+(1-2*C)*Td(1)+C*Td(2);
    d(2:nLayers-2) = C*Td(1:nLayers-3)+(1-2*C)*Td(2:nLayers-2)+C*Td...
        (3:nLayers-1);
    d(nLayers-1) = 2*C*constT_i+C*Td(nLayers-2)+(1-2*C)*Td(nLayers-1);
    
    % Calculate A and S
    A = zeros(nLayers-1,1);
    S = zeros(nLayers-1,1);
    A(1) = b;
    S(1) = d(1);
    for iLayer = 2:nLayers-1
        A(iLayer) = b-a/A(iLayer-1)*c;
        S(iLayer) = d(iLayer)+a/A(iLayer-1)*S(iLayer-1);
    end
    
    % Calculate new temperature profile
    tempProf = zeros(nLayers-1,1);
    tempProf(nLayers-1) = S(nLayers-1)/A(nLayers-1);
    for iLayer = 1:nLayers-2
        tempProf(nLayers-1-iLayer) = 1/A(nLayers-1-iLayer)*...
            (S(nLayers-1-iLayer)+c*tempProf(nLayers-1-iLayer+1));
    end
return


