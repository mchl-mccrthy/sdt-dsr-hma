% Fill gaps in supraglacial debris thickness data

% This function fills gaps in supraglacial debris thickness data as
% described by McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function sdt = sdtfilter2(sdt,sdt2,z,smb,smbE,wdw,sT)

% Note. sdt2 is actual supraglacial debris thickness. sdt can be e.g. the 
% uncertainties

% Filter SDT using SMB threshold
sdt(smb+smbE >= 0) = NaN;

% Make SDT 0 where there is ice
sdt(sT == 1) = 0;

% Create mask
mask = ~isnan(sdt2);

% Loop through elevations
zs = min(z(mask),[],'all'):1:max(z(mask),[],'all')-wdw+1;
nZs = length(zs);
[sdtAvg,zMid] = deal(nan(nZs,1));
for iZ = 1:nZs

    % Make elevation condition
    zCond = z >= zs(iZ) & z < zs(iZ)+wdw-1;

    % Get SDT in that elevation range
    sdtAvg(iZ) = nanmean(sdt2(zCond),'all');
    zMid(iZ) = (zs(iZ)+zs(iZ)+wdw-1)/2;
end

% Interpolate
if length(zMid) >= 2
    sdtAvg = interp1(zMid,sdtAvg,double(zs),'nearest','extrap');

    % Apply elevation-based filter
    for iZ = 1:nZs
        sdt(z == zs(iZ) & sdt2 > sdtAvg(iZ)*3 & sdt2 > 0.3) = NaN;
    end
end

% Make SDT NaN where there is ice
sdt(sT == 1) = NaN;

end