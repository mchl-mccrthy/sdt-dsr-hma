% Calculate solar position 

% This function calculates the position of the sun following Walraven 
% (1978)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

% Notes
% - Set zone and dasvtm to zero for UTC
%
% Inputs
% - dateTime = date in MATLAB datetime format
% - zone = international time zone, where GMT is 0
% - dasvtm = daylight saving time, where 1 is yes, 0 is no
% - lat = latitude, -90 to 90, S is negative, N is positive (deg)
% - lon = longitude, -180 to 180, W is negative, E is positive (deg)
% 
% Outputs
% - solAz = solar azimuth (deg)
% - solEl = solar elevation (deg)
%
% % Test using example data in Walraven (1978)
% dateTime1 = datetime(1977,4,30,1,0,0);
% dateTime2 = datetime(1977,4,30,24,0,0);
% dateTime = dateTime1:hours(1):dateTime2;
% zone = 8;
% dasvtm = 1;
% lat = 38.538;
% lon = -121.758; % With minus sign because lon input here is east positive
% [solAz,solEl] = solarpositionfast(dateTime,zone,dasvtm,lat,lon);
% solAz = -solAz+180; % This because solAz output here is 0 to 360
% testData = readtable('walraven_data.csv');
% figure('Name','Solar azimuth (deg)')
% scatter(testData.solAz,solAz); hold on
% plot([-180,180],[-180,180])
% xlabel('Walraven (1978)'); ylabel('This function')
% figure('Name','Solar elevation (deg)')
% scatter(testData.solEl,solEl); hold on
% plot([-90,90],[-90,90])
% xlabel('Walraven (1978)'); ylabel('This function')

function [solAz,solEl] = solarpositionfast(dateTime,zone,dasvtm,lat,lon)

% Define constants
twoPi = 6.2831853;
rad = 0.017453293;

% Convert date to yr,doy,hr (where midnight = 24),min,sec, so doy is doy-1
% on days before March 1 in a leap year
[year,~,~,hour,minute,second] = datevec(dateTime);
tod = timeofday(dateTime);
idx1 = tod == 0;
dateTime(idx1) = dateTime(idx1)-days(1);
[year(idx1),~,~,hour(idx1),minute(idx1),second(idx1)] = datevec(dateTime...
    (idx1));
hour(idx1) = 24;
isleap = ~mod(year,4) & (mod(year,100) | ~mod(year,400));
doy = day(dateTime,'dayofyear');
doy(isleap) = doy(isleap)-1;

% Calculate time
delYr = year-1980;
leap = fix(delYr/4);
T = hour+(minute+second/60)/60;
time = delYr*365+leap+doy-1+T/24;
time(delYr == leap*4) = time(delYr == leap*4)-1;
time(delYr < 0 & delYr ~= leap*4) = time(delYr < 0 & delYr ~= leap*4)-1;

% Convert longitude to degrees west
lon = -lon;

% Calculate longitude of sun
theta = (360*time/365.25)*rad;
g = -0.031271-4.53963e-7*time+theta;
el = 4.900968+(3.67474e-7).*time+(0.033434-2.3e-9.*time).*sin(g)+...
    0.000349.*sin(2.*g)+theta;

% Calculate angle between plane of ecliptic and plane of celestial equator
eps = 0.409140-6.2149e-9*time;

% Calculate right ascenscion and declination of the sun
sel = sin(el);
a1 = sel.*cos(eps);
a2 = cos(el);
ra = atan2(a1,a2);
ra(ra < 0) = ra(ra < 0)+twoPi;
decl = asin(sel.*sin(eps));

% Calculate siderial and local siderial time
st = 1.759335+twoPi*(time/365.25-delYr)+3.649e-7*time;
st(st >= twoPi) = st(st >= twoPi)-twoPi;
s = st-lon.*rad+1.0027379.*(zone-dasvtm+T).*15.*rad;
s(s >= twoPi) = s(s >= twoPi)-twoPi;

% Calculate hour angle and local latitude
h = ra-s;
phi = lat.*rad;

% Calculate solar elevation and azimuth
solEl = asin(sin(phi).*sin(decl)+cos(phi).*cos(decl).*cos(h));
solAz = asin(cos(decl).*sin(h)./cos(solEl))/rad;

% Correct circular statistics
solAz(sin(solEl) < sin(decl)./sin(phi) & solAz < 0) = solAz(sin(solEl) <...
    sin(decl)./sin(phi) & solAz < 0)+360;
solAz(sin(solEl) < sin(decl)./sin(phi)) = 180-solAz(sin(solEl) < sin...
    (decl)./sin(phi));
solEl = solEl/rad;

% Change to 360 degrees
solAz = -solAz+180;

% Make sure solar azimuth is real number
solAz = real(solAz);

end

    