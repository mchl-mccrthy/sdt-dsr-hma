% Compute specific mass balance

% This function processes input data and runs a surface energy-balance 
% model for debris-covered glaciers to compute ice melt under debris as
% described by McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function smbTable = computesmb(dem,demR,demInfo,meteoData,msRowCol,...
    debProps,constants,mcVals,meteoParams,runPix)

% Unpack debris properties
debc = debProps.debc;
debrho = debProps.debrho;
debepsilon = debProps.debepsilon;

% Unpack constants
constg = constants.constg;
constk_vk = constants.constk_vk;
constsigma = constants.constsigma;
constR = constants.constR;
constM_a = constants.constM_a;
constp_0 = constants.constp_0;
constL_v = constants.constL_v;
constc_a = constants.constc_a;
constrho_i = constants.constrho_i;
constrho_w = constants.constrho_w;
constT_i = constants.constT_i;
constI_0 = constants.constI_0;
constpsi = constants.constpsi;
constT_0 = constants.constT_0;
constL_f = constants.constL_f;
constc_w = constants.constc_w;
constgamma = constants.constgamma;

% Unpack Monte Carlo ranges
debkSamps = mcVals.debkSamps;
debz_0Samps = mcVals.debz_0Samps;
debalphaSamps = mcVals.debalphaSamps;
gammaErrs = mcVals.gammaErrs;
dtSamps = mcVals.dtSamps;
S_inErrs = mcVals.S_inErrs;
L_inErrs = mcVals.L_inErrs;
T_aErrs = mcVals.T_aErrs;
uErrs = mcVals.uErrs;

% Unpack meteorological parameters
timestep = meteoParams.timestep;
z_a = meteoParams.z_a;
z_u = meteoParams.z_u;

% Unpack run rows and columns
runRow = real(runPix);
runCol = imag(runPix);

% Get meteorological data information
meteoDatetime = meteoData.dateTime;
hgt = meteoData.hgt;
nTimesteps = length(meteoDatetime);

% Get size of DEM pixels and UTM zone
utmZone = demInfo.Zone;
pixelSize = demR.CellExtentInWorldX;

% Get size of DEM, calculate slope and aspect, get DEM pixel coordinates
[nRows,nCols] = size(dem);
[slope,aspect] = slopeaspect(dem,pixelSize);
[xDem,yDem] = pixcenters(demR,nRows,nCols); 
[xDem,yDem] = meshgrid(xDem,yDem);

% Calculate atmospheric pressure at reanalysis elevation
p_aHgt = constp_0*((1-(constgamma*hgt/constT_0)).^(constg*constM_a/...
    (constR*constgamma)));

% Make arrays in which to store smbMod and xPix, yPix
nRunPix = length(runPix);
smbMod = nan(1,nRunPix);
xPix = smbMod;
yPix = smbMod;
zPix = smbMod;

% Run EB model, statistically downscaling meteorology and varying debris
% thickness and properties
for iRunPix = 1:nRunPix
    
    % Get run row and column
    row = runRow(iRunPix);
    col = runCol(iRunPix);
    
    % Find appropriate meteorological data pixel for DEM pixel
    meteoDataUniPix = table2array(meteoData.uniPix);
    msRowColiRunPix = table2array(msRowCol(iRunPix,:));
    [~,uniPix] = ismember(msRowColiRunPix,meteoDataUniPix,'rows');

    % Unpack meteorological data
    S_in = meteoData.S_in(uniPix,:);
    L_in = meteoData.L_in(uniPix,:);
    u = meteoData.u(uniPix,:);
    T_a = meteoData.T_a(uniPix,:);
    rh = meteoData.rh(uniPix,:);
    r = meteoData.r(uniPix,:);
    S = meteoData.S(uniPix,:);

    % Get debris properties, thickness etc from MC ranges
    debk = debkSamps(iRunPix);
    debalpha = debalphaSamps(iRunPix);
    debz_0 = debz_0Samps(iRunPix);
    dt = dtSamps(iRunPix);
    gammaMC = constgamma+gammaErrs(iRunPix);

    % Get x, y coordinates of pixel
    xPix(iRunPix) = xDem(row,col);
    yPix(iRunPix) = yDem(row,col);
    
    % Get elevation
    zPix(iRunPix) = dem(row,col);
    
    % Get slope and aspect of pixel
    Z_dash = slope(row,col);
    A_dash = aspect(row,col);

    % Get solar position in relation to DEM pixel
    [pixLat,pixLon] = utm2ll(xPix(iRunPix),yPix(iRunPix),utmZone);
    [A,E] = solarpositionfast(meteoDatetime,0,0,...
        pixLat,pixLon);
    Z = 90-E;

    % Downscale air temperature according to lapse rate
    T_a = T_a-gammaMC.*(zPix(iRunPix)-hgt(uniPix));

    % Calculate atmospheric pressure
    p_a = constp_0*((1-(gammaMC*zPix(iRunPix)/...
        constT_0))^(constg*constM_a/(constR*constgamma)));

    % Determine pixel shading
    scaleZ = 1/pixelSize; % Calculate elevation scale
    angXY = A-90; % Get solar height step, where 0 is East
    distStep = 3; % Choose step distance
    vec = complex(sind(angXY),cosd(angXY)); % Calc. dir'n vectors
    distStepVec = distStep*vec;
    tanSolEl = tand(E); % Calculate solar height step     
    shade = nan(nTimesteps,1);
    for iTimestep = 1:nTimesteps
        shade(iTimestep) = f_whether_shaded_pixel(dem,...
            complex(row,col),vec(iTimestep),...
            distStepVec(iTimestep),distStep,scaleZ,...
            tanSolEl(iTimestep),nRows,nCols);
    end
    
    % Get sky-view factor
    azInt = 12;
    f_sv = pixelskyviewfactor(dem,pixelSize,azInt,distStep,row,col);

    % Add error to meteorological variables 
    S_inMC = S_in+S_inErrs(iRunPix); S_inMC(S_inMC < 0) = 0;
    L_inMC = L_in+L_inErrs(iRunPix); L_inMC(L_inMC < 0) = 0;
    T_aMC = T_a+T_aErrs(iRunPix);
    uMC = u+uErrs(iRunPix); uMC(uMC < 0) = 0;

    % Convert precipitation rate from m/h to m/s
    r = r/timestep;
    
    % Convert wind speed to same height as air temperature
    uMC = uMC*(log(z_a)/debz_0)/(log(z_u)/debz_0);
    
    % Calculate saturation vapour pressure at the measurement height (using
    % Teten's equation, after Murray, 1966 - Journal of Applied Meteorology
    e_s_a = 610.78.*exp((17.27.*(T_aMC-273.15))./(T_aMC-35.86));

    % Calculate partial pressure of water at measurement height
    e_a = rh.*e_s_a./100;

    % Calculate specific humidity at measurement height, given partial 
    % pressure of water
    q_a = 0.622.*e_a./(p_a-(0.378.*e_a));

    % Run EB model to get modelled SMB
    nLayers = 10;
    [~,~,melt] = ebmodel(dt,S_inMC,L_inMC,T_aMC,uMC,...
        q_a,r,p_a,S,nLayers,debalpha,debc,debepsilon,...
        debk,debrho,debz_0,z_a,z_a,timestep,constc_a,...
        constc_w,constg,constk_vk,constL_f,constL_v,...
        constM_a,constR,constrho_i,constrho_w,...
        constsigma,constT_i,f_sv,shade,Z,Z_dash,A,A_dash);
    sumMelt = nansum(melt);
    smbMod(iRunPix) = -sumMelt;
end

% Write table of outputs (and useful inputs)
utmZones = ones(nRunPix,1)*utmZone;
dtSamps = dtSamps(:);
smbMod = smbMod(:);
xPix = xPix(:);
yPix = yPix(:);
zPix = zPix(:);
debkSamps = debkSamps(:);
debz_0Samps = debz_0Samps(:);
debalphaSamps = debalphaSamps(:);
gammaErrs = gammaErrs(:);
S_inErrs = S_inErrs(:);
L_inErrs = L_inErrs(:);
T_aErrs = T_aErrs(:);
uErrs = uErrs(:);
smbTable = table(dtSamps,smbMod,xPix,yPix,zPix,utmZones,debkSamps,...
    debz_0Samps,debalphaSamps,gammaErrs(:),S_inErrs(:),L_inErrs(:),...
    T_aErrs(:),uErrs(:));

end

