% Calculate debris-supply rate and englacial debris content

% This script calculates the supply rate of debris to the surface of a 
% debris-covered glacier from its debris-supply slopes, as described by
% McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Construct filenames for input data
RGIID = '15.03733'; % RGIID of Khumbu Glacier
fnGlacMaskOrig = ['Miles et al 2021\' RGIID '_debris.tif'];
fnGlacMaskMod = ['Improved glacier masks\' RGIID '_Headwall.tif'];
fnDem = ['Miles et al 2021\' RGIID '_AW3D.tif'];
fnDt = ['Processed SDT\' RGIID '_dt.tif'];
fnDtUp = ['Processed SDT\' RGIID '_dtup.tif'];
fnDtLow = ['Processed SDT\' RGIID '_dtlow.tif'];
fnSmb = ['Miles et al 2021\' RGIID '_zSMB.tif'];
fnSmbErr = ['Miles et al 2021\' RGIID '_zSMBe.tif'];
fnV = ['Dehecq et al 2019\' RGIID '_v.tif'];
fnU = ['Dehecq et al 2019\' RGIID '_u.tif'];
fnDss = ['Debris-supply slopes/' RGIID '_dss.tif'];

% Load original glacier mask, modified glacier mask, debris supply area, 
% DEM, debris thickness etc
surfTypeMod = geotiffread(fnGlacMaskMod);
surfTypeOrig = geotiffread(fnGlacMaskOrig);
[dem,demR] = geotiffread(fnDem);
demInfo = geotiffinfo(fnDem);
dt = geotiffread(fnDt);
dtUp = geotiffread(fnDtUp);
dtLow = geotiffread(fnDtLow);
smb = geotiffread(fnSmb);
smbErr = geotiffread(fnSmbErr);
dssMask = geotiffread(fnDss);
u = geotiffread(fnU); % u component of velocity from Dehecq et al (2019)
v = geotiffread(fnV); % v compoment of velocity from Dehecq et al (2019)

% Create masks for headwall area computation 
surfTypeMod(isnan(surfTypeMod)) = 0;
iceMod = surfTypeMod == 1;
debMod = surfTypeMod == 2;

% Make glacier mask
glacMask = iceMod | debMod;
% glacMask = surfTypeOrig == 1 | surfTypeOrig == 2; % Can also use 
% unmodified Miles et al 2021 mask

% Get pixel centres and pixel size
[xg,yg] = pixcenters(demInfo,'makegrid');
pixelSize = mode(diff(xg(1,:)));

% Make debris mask
debMask = dt > 0;

% Load DEM and get slope and aspect
dem = double(dem);
[slope,aspect] = slopeaspect(dem,pixelSize);

% Compute debris flux with distance from terminus
d = terminus_distance(glacMask,dem,pixelSize);  
[~,x,Q_sd,xsArea,~,~,~,w_d] = deb_through_fluxes_distance...
    (debMask,d,u,v,dt,pixelSize,dem);
[~,~,Q_sdUp,~,~,~,~,~] = deb_through_fluxes_distance...
    (debMask,d,u,v,dtUp,pixelSize,dem); 
[~,~,Q_sdLow,~,~,~,~,~] = deb_through_fluxes_distance...
    (debMask,d,u,v,dtLow,pixelSize,dem);    

% Get maximum flux, distance from terminus, debris areas, and debris 
% emergence
Q_sdSmooth = smooth(Q_sd,0.1,'moving');
Q_sdSmoothUp = smooth(Q_sdUp,0.1,'moving');
Q_sdSmoothLow = smooth(Q_sdLow,0.1,'moving');
[Q_sdA_out,Q_sdA_outInd] = max(Q_sdSmooth,[],'omitnan');
[Q_sdA_outUp,~] = max(Q_sdSmoothUp,[],'omitnan');
[Q_sdA_outLow,~] = max(Q_sdSmoothLow,[],'omitnan');
Q_sdA_outErrUp = Q_sdA_outUp-Q_sdA_out;
Q_sdA_outErrLow = Q_sdA_out-Q_sdA_outLow;
dQ_sdA_out = x(Q_sdA_outInd);
Q_sdA_in = 0;
A_dA = sum(d > dQ_sdA_out & debMask,'all')*pixelSize^2;
A_dI = sum(d <= dQ_sdA_out & debMask,'all')*pixelSize^2;
A_d = sum(debMask,'all')*pixelSize^2;
q_edA = (Q_sdA_out-Q_sdA_in)/A_dA;

% Compute debris-supply-slope area
A_dsT = sum((dssMask*pixelSize^2)./cosd(slope),'all');

% Specify rock and debris densities
rho_d = 1842.3;
rho_r = 2700;

% Compute debris-supply rate (mm/yr) and englacial debris content (frac)
smb(smb < -10) = -10; % Remove unrealistic SMBs
smb = smb*999.7/915; % Convert m w.e. to i.e.
M_A = abs(nansum(smb(smb < 0 & d > dQ_sdA_out & debMask...
    ),'all'))/(A_dA/pixelSize^2);
M_Aerr = abs(nansum(smbErr(smb < 0 & d > dQ_sdA_out & debMask...
    ),'all'))/(A_dA/pixelSize^2);
c_edA = q_edA*rho_d/(M_A*rho_r+q_edA*rho_d);
c_ed = c_edA*915/850; % Density conversion (Huss, 2013)
M_I = abs(nansum(smb(smb < 0 & d <= dQ_sdA_out & debMask...
    ),'all'))/(A_dI/pixelSize^2);
q_edI = c_edA*M_I*rho_r/(rho_d-rho_d*c_edA);
q_edA_d = q_edA*A_dA+q_edI*A_dI;
q_dsT = 1000*rho_d*q_edA_d/(rho_r*A_dsT);

% Calculate debris-supply rate and englacial debris content uncertainties
[q_dsUp,c_edAUp,c_edUp] = getdsrerror(q_edA,...
    Q_sdA_outErrUp,Q_sdA_out,q_edI,...
    q_edA_d,A_dA,A_dI,rho_d,rho_r,M_A,...
    M_Aerr,q_dsT,c_edA,c_ed);
[q_dsLow,c_edALow,c_edLow] = getdsrerror(q_edA,...
    Q_sdA_outErrLow,Q_sdA_out,q_edI,...
    q_edA_d,A_dA,A_dI,rho_d,rho_r,M_A,...
    M_Aerr,q_dsT,c_edA,c_ed);


