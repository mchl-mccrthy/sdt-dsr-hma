% Run energy-balance model

% This script runs a surface energy-balance model for debris-covered
% glaciers to compute ice melt under debris layers of varying thickness,
% and ultimately to calculate supraglacial debris thickness as described by
% McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Specify file and folder names
fnConstants = 'Input data/constants.mat';
fnDebProps = 'Input data/debris_properties.mat';
fnMcRanges = 'Input data/monte_carlo_ranges_Ta.mat';
fnMeteoParams = 'Input data/meteorological_parameters.mat';
foInput = 'Miles et al 2021';
foSmbMod = 'Modelled SMB';

% Specify date range over which to run EB model (middle year of 2000-2016)
simDatetime.start = datetime(2009,1,1,0,0,0);
simDatetime.end = datetime(2009,12,31,23,0,0);

% Get constants, debris properties, meteorological parameters and Monte
% Carlo simulation ranges
load(fnConstants);
load(fnDebProps);
load(fnMcRanges);
load(fnMeteoParams);

% Specify RGIID
RGIID = '15.03733';

% Generate filenames for datasets to be imported
fnDem = [foInput  '/' RGIID '_AW3D.tif'];
fnSurfType = [foInput  '/' RGIID '_debris.tif'];
fnSmbE = [foInput  '/' RGIID '_zSMBe.tif'];
fnSmb = [foInput  '/' RGIID '_zSMB.tif'];
    
% Get surface types and create debris mask
[surfType,~,~] = getgeotiff(fnSurfType);
debMask = surfType == 2;

% Get DEM
[dem,demR,demInfo] = getgeotiff(fnDem);

% Specify ELA or estimate as median elevation of glacier and create ELA
% mask
ela = 5315; % From Miles et al (2021)
if isnan(ela)
    ela = prctile(reshape(dem(surfType == 1 | debMask),...
        [],1),50);
end
elaMask = dem < ela;

% Get DEM row and column of all debris-covered pixels below ELA
[debRow,debCol] = find(debMask & elaMask);
debRowCol = complex(debRow,debCol);

% Get SMB
[smbDist,smbDistR,smbDistInfo] = getgeotiff(fnSmb);
[smbDistE,~,~] = getgeotiff(fnSmbE);

% Convert SMB and SMB error from w.e. to i.e.
smbDist = smbDist*constants.constrho_w/constants.constrho_i;
smbDistE = smbDistE*constants.constrho_w/constants.constrho_i;

% Get size of surface mass balance image
sizeSmbDist = size(smbDist);

% Get number of samples according to elevation range
demEla = dem;
demEla(demEla > ela) = NaN;
demEla(~debMask) = NaN;
zDebMin = min(demEla,[],'all','omitnan');
zDebMax = max(demEla,[],'all','omitnan');
zDebRange = zDebMax-zDebMin;
nSamps = max(round(zDebRange),100);

% Sample some number (nSamps) of pixels, uniformly at random, 
% from all pixels within the debris mask. These are the pixels 
% for which to run the EB model. In this example script we 
% simply load the run pixels rather than random sampling,
% because it is not possible to store all the potential meteo 
% data on GitHub
load('Input data/run_pixels.mat');

% Get relevant pixels of meteorological data for those pixels,
% and associated meteorological data
load('Input data/meteo_pixels.mat');
load('Input data/meteo_data.mat');

% Generate the Monte Carlo values for each glacier. In this example script 
% we load the samples we used for Khumbu Glacier
load('Input data/mc_samples.mat')

% Calculate SMB and write to table with debris thickness, 
% Monte Carlo and pixel lat-lon values
smbTable = computesmb(dem,demR,demInfo,meteoData,msRowCol,...
    debProps,constants,mcVals,meteoParams,runPix);

% Save table
writetable(smbTable,[foSmbMod '/smb_mod_' ...
    RGIID '.csv'])




