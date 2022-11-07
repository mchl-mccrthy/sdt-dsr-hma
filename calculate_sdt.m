% Calculate supraglacial debris thickness

% This script calculates supraglacial debris thickness from specific mass 
% balance data, as described by McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Specify file and folder names
foInput = 'Miles et al 2021';
foDT = 'Raw SDT';
foSmbMod = 'Modelled SMB';
fnMcRanges = 'Input data/monte_carlo_ranges_Ta.mat';
foOC = 'Ostrem curves';

% Specify RGIID
RGIID = '15.03733';

% Load Monte Carlo ranges
load(fnMcRanges);

% Generate filenames for datasets to be imported
fnDem = [foInput  '/' RGIID '_AW3D.tif'];
fnSurfType = [foInput  '/' RGIID '_debris.tif'];
fnSmbE = [foInput  '/' RGIID '_zSMBe.tif'];
fnSmb = [foInput  '/' RGIID '_zSMB.tif'];
fnDt = [foDT  '/' RGIID '_dt.tif'];
fnDtUp = [foDT  '/' RGIID '_dtup.tif'];
fnDtLow = [foDT  '/' RGIID '_dtlow.tif'];
fnSmbMod = [foSmbMod '/smb_mod_' RGIID '.csv'];
fnOC = [foOC '/oc_' RGIID '.csv'];

% Get surface types and create debris mask
[surfType,~,~] = getgeotiff(fnSurfType);
debMask = surfType == 2;
glacMask = surfType == 1 | surfType == 2;

% Get distributed input datasets
[dem,demR,demInfo] = getgeotiff(fnDem);
[smbDist,smbDistR,smbDistInfo] = getgeotiff(fnSmb);
[smbDistE,~,~] = getgeotiff(fnSmbE);

% Get size of surface mass balance image
sizeSmbDist = size(smbDist);

% Load modelled SMBs
smbTable = readtable(fnSmbMod);

% Get ELA or estimate as median elevation of glacier and create ELA
% mask
ela = 5315; % From Miles et al (2021)
if isnan(ela)
    ela = prctile(reshape(dem(surfType == 1 | debMask),...
        [],1),50);
end
elaMask = dem < ela;

% Get elevation range and number of samples
dem2 = dem;
dem2(~elaMask) = NaN;
dem2(~debMask) = NaN;
zDebMin = min(dem2,[],'all','omitnan');
zDebMax = max(dem2,[],'all','omitnan');
zDebRange = zDebMax-zDebMin;

% How many Ostrem curves to fit per glacier?
nOCs = max(round(zDebRange/100),1);

% Get elevation ranges for Ostrem curves
zRangeEdges = linspace(zDebMin,zDebMax,nOCs+1);

% Allocate space for debris thickness estimate
dtExt = nan(sizeSmbDist);
dtUp = dtExt;
dtLow = dtExt;

% Fit an Ostrem curve for each elevation range
startPt = [-6,0.02];
ocTable = table();
ocTable.zMin = zRangeEdges(1:end-1)';
ocTable.zMax = zRangeEdges(2:end)';
ocTable.c1 = nan(nOCs,1);
ocTable.c2 = nan(nOCs,1);
ocTable.r2 = nan(nOCs,1);
ocTable.nPts = nan(nOCs,1);
for iOC = 1:nOCs
        dt = smbTable.dtSamps(smbTable.zPix >= zRangeEdges(iOC) ...
            & smbTable.zPix < zRangeEdges(iOC+1));
        smb = smbTable.smbMod(smbTable.zPix >= zRangeEdges(iOC) ...
            & smbTable.zPix < zRangeEdges(iOC+1));
        ocTable.nPts(iOC) = length(dt);
    if length(dt) >= 30
        [dtFit,piLookupSMB,piLookupDT,ocTable.r2(iOC)] = ...
            fitostremcurve(smb,dt,startPt);
        ocTable.c1(iOC) = dtFit.c1;
        ocTable.c2(iOC) = dtFit.c2;
    end
    
    % Use curve to generate distributed DT
    if ocTable.r2(iOC) > 0.4

        % Create elevation mask
        elevMask = dem >= zRangeEdges(iOC) & dem < zRangeEdges(iOC+1);
        dtExt(elevMask & debMask) = dtFit.c1.*dtFit.c2./smbDist...
            (elevMask & debMask)-dtFit.c2;
        dtExt(dtExt < mcRanges.dtSamps(1)) = mcRanges.dtSamps(1);
        dtExt(dtExt > 5) = 5;

        % Get measured SMB uncertainty
        smbDistUp = smbDist+smbDistE;
        smbDistLow = smbDist-smbDistE;
        piIndUp = nearestpoint(smbDistUp(:),piLookupSMB(:,1));
        piIndUp = reshape(piIndUp,sizeSmbDist);
        dtUp(elevMask & debMask) = changem(piIndUp(elevMask & ...
            debMask),piLookupDT,1:length(piLookupDT));
        piIndLow = nearestpoint(smbDistLow(:),piLookupSMB(:,2));
        piIndLow = reshape(piIndLow,sizeSmbDist);
        dtLow(elevMask & debMask) = changem(piIndLow(elevMask & ...
            debMask),piLookupDT,1:length(piLookupDT));
        
        % If SMB is more than 0, make DT and uncertainties NaN
        dtExt(smbDist >= 0) = NaN;
        dtUp(smbDist >= 0) = NaN;
        dtLow(smbDist >= 0) = NaN;
        dtLow(smbDist < 0 & smbDistUp >= 0) = 0.01;
    end
end

% Write Ostrem curve table to file
writetable(ocTable,fnOC)

% Write post-processed debris thickness to geotiff
geotiffwrite([foDT '/' RGIID '_dt.tif'],...
    dtExt,smbDistR,'CoordRefSysCode',...
    smbDistInfo.GeoTIFFCodes.PCS)

% Write uncertainties to geotiff
geotiffwrite([foDT '/' RGIID '_dtlow.tif'],...
    dtLow,smbDistR,'CoordRefSysCode',...
    smbDistInfo.GeoTIFFCodes.PCS)
geotiffwrite([foDT '/' RGIID '_dtup.tif'],...
    dtUp,smbDistR,'CoordRefSysCode',...
    smbDistInfo.GeoTIFFCodes.PCS)
