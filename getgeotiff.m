% Get geotiff and associated information


function [GT,GTR,GTInfo] = getgeotiff(fnGT)

% Open DEM geotiff
[GT,GTR] = geotiffread(fnGT);
GT = double(GT);
GTInfo = geotiffinfo(fnGT);

end