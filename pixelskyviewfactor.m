% Sky-view factor

% This function calculates sky-view factor from an elevation model, 
% following Arnold et al (2006). The the bulk of the algorithm comes from 
% the function f_whether_shaded by Neil Arnold.
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function svf = pixelskyviewfactor(dem,pixelSize,azInt,distStep,row,col)

% Calculate some parameters
nAzs = 360/azInt; % Number of azimuths at which to calculate horizon angle
scaleZ = 1/pixelSize; % Get scale of elevation model in z direction

% Get dimensions of DEM
[dimX,dimY] = size(dem);

% Preallocate horizon angles
horAngs = zeros(1,nAzs);

% Get horizon angles at different azimuths
for iAz = 1:nAzs

    % Get azimuth for which to calculate horizon angle. Starts at
    % 90+azInt degrees because zero is east in MATLAB
    azAng = (iAz*azInt)-90;

    % Get elevation of pixel
    pixElev = scaleZ*dem(row,col);

    % Set initial position
    point = complex(col,row);
    dist = 0;
    done = false;

    % Calculate direction vectors
    vec = complex(sind(azAng),cosd(azAng));
    distStepVec = distStep*vec;

    % For each point on the transect
    while (~done)
        point = point+distStepVec;
        dist = dist+distStep;

        % Get DEM pixel to sample
        pixRowSamp = fix(imag(point));
        pixColSamp = fix(real(point));

        % Check if we've gone off the edge of the DEM
        done = ((pixColSamp <= 0   ) | ...
        +       (pixRowSamp <= 0   ) | ...
        +       (pixColSamp >= dimY) | ...
        +       (pixRowSamp >= dimX));

        % If we haven't gone off the edge
        if (~done)

            % Get elevation of sample pixel
            pixElevSamp = scaleZ*dem(pixRowSamp,pixColSamp);

            % Compute horizon angle of sample pixel
            tanHorAngSamp = (pixElevSamp-pixElev)/abs(dist);
            horAngSamp = atand(tanHorAngSamp);

            % If new horizon angle is greater than previous
            if horAngSamp > horAngs(iAz)
                horAngs(iAz) = horAngSamp;    
            end
        end
    end
end

% Get sky-view factor
svf = azInt/360*sum((cosd(horAngs)).^2);

end

        
        
        
        