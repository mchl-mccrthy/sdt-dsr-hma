% Fill gaps in calculated supraglacial debris thickness data 

% This script fills gaps in calculated supraglacial debris thickness data,
% as described by McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Specify file and folder names
foInput = 'Miles et al 2021';
foDt = 'Raw SDT';
foDtProc = 'Processed SDT';

% Specify RGIID
RGIID = '15.03733';

% Load data
dem = geotiffread([foInput  '/' RGIID '_AW3D.tif']);
[dt,dtR] = geotiffread([foDt '/' RGIID '_dt.tif']);
[dtup,~] = geotiffread([foDt '/' RGIID '_dtup.tif']);
[dtlow,~] = geotiffread([foDt '/' RGIID '_dtlow.tif']);
dtInfo = geotiffinfo([foDt '/' RGIID '_dt.tif']);
surfType = geotiffread([foInput '/' RGIID '_debris.tif']);
segs = geotiffread([foInput  '/' RGIID '_zones.tif']);
smb = geotiffread([foInput  '/' RGIID '_zSMB.tif']);      
smbE = geotiffread([foInput  '/' RGIID '_zSMBe.tif']);           

% Get ELA or estimate as median elevation of glacier and create ELA
% mask
ela = 5315; % From Miles et al (2021)
if isnan(ela)
    ela = prctile(reshape(dem(surfType == 1 | surfType ==2),...
        [],1),50);
end

% Fill SDT gaps where necessary. Glacier 1814 has faulty DEM
if ~isnan(nanmean(dt,'all'))
    
    % Do debris thickness interpolation and average over elevation
    % bands
    [nRows,nCols] = size(dt);
    wdw = 100;
    dt = sdtfilter2(dt,dt,dem,smb,smbE,wdw,surfType);
    if ~isnan(nanmean(dt,'all')) % Check there is data to interpolate from
        dt(surfType == 1) = 0;
        dt = inpaint_nans(dt);
        dt(surfType ~= 2 | dem > ela) = NaN;
        dt(dt < 0.01) = 0.01;
        dt(dt > 5) = 5;    
        segNos = unique(segs(surfType == 2));
        nSegs = length(segNos); 
        for iSeg = 1:nSegs
            dt(segs == segNos(iSeg)) = mean(dt(segs == segNos(iSeg)),...
                'all');
        end

        dtup = sdtfilter2(dtup,dt,dem,smb,smbE,wdw,surfType);
        dtup(surfType == 1) = 0;
        dtup = inpaint_nans(dtup);
        dtup(surfType ~= 2 | dem > ela) = NaN;
        dtup(dtup < 0.01) = 0.01;
        dtup(dtup > 5) = 5;
        for iSeg = 1:nSegs
            dtup(segs == segNos(iSeg)) = mean(dtup(segs == ...
                segNos(iSeg)),'all');
        end

        dtlow = sdtfilter2(dtlow,dt,dem,smb,smbE,wdw,surfType);
        dtlow(surfType == 1) = 0;
        dtlow = inpaint_nans(dtlow);
        dtlow(surfType ~= 2 | dem > ela) = NaN;
        dtlow(dtlow < 0.01) = 0.01;
        dtlow(dtlow > 5) = 5;
        for iSeg = 1:nSegs
            dtlow(segs == segNos(iSeg)) = mean(dtlow(segs == ...
                segNos(iSeg)),'all');
        end

        % Model physics is poor for thinner debris as e.g. Evatt et
        % al (2015). Critical thickness in HMA is 0.036 m on average 
        % (Reznichenko et al. 2010), so where modelled SDT is < 0.05 m...
        dt(dt < 0.05) = 0.03;
        dtlow(dt < 0.05) = 0.01;
        dtup(dt < 0.05) = 0.05;

        % Save results to new folder
        geotiffwrite([foDtProc '/' RGIID '_dt.tif'],dt,dtR,...
            'CoordRefSysCode',dtInfo.GeoTIFFCodes.PCS)
        geotiffwrite([foDtProc '/' RGIID '_dtup.tif'],dtup,dtR,...
            'CoordRefSysCode',dtInfo.GeoTIFFCodes.PCS)
        geotiffwrite([foDtProc '/' RGIID '_dtlow.tif'],dtlow,dtR,...
            'CoordRefSysCode',dtInfo.GeoTIFFCodes.PCS)
    end
end
