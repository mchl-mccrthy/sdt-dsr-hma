function [outdists]=terminus_distance(mask,DEM,dx)
% terminus_distance - Function to estimate the distance from any point
% within a mask to the glacier terminus (given by the lowest elevation)
% according to the geodesic distance along the glacier centerline
%
% Syntax:  [outdists]=terminus_distance(mask,DEM,dx)
%
% Inputs:
%    mask - binary raster mask of glacier extent
%    DEM - double-precision raster of elevation with same extent as mask
%    dx - pixel size in porjected units corresponding to veloctiy and DEM units (m)
%
% Outputs:
%    outdists - raster of geodesic distances from the terminus, determined
%    along the glacier centerline

% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% June 2020; Last revision: 01-July-2020distInt=1000; %interval of segmentation

%% set mask
mask2=double(mask);
mask2(mask2==0)=NaN;


%% smooth DEM for better contours
DEM1 = double(DEM).*mask2;

DEM2 = inpaint_nans(DEM1,0);

DEM3 = imgaussfilt(DEM2,2);
DEM3 = imgaussfilt(DEM3,2);

%% define terminus and skeleton
term=(DEM3==min(DEM3(mask)));
% 
% termDist=bwdist(term).*dmask2.*dx;
% termDist2=bwdistgeodesic(dmask,term)
% 
% edgeDist=bwdist(isnan(dmask)).*dmask2.*dx;
% imagesc(termDist);
% imagesc(edgeDist);

Skel=bwskel(mask);
Skel2=bwlabel(Skel);
if max(Skel2(:))>1
    Skel3=Skel;
    for iseg=1:max(Skel2(:))
        Da=bwdist(Skel2==iseg,'quasi-euclidean'); %distance from current seg
        Db=bwdist((Skel2>0) & (Skel2~=iseg),'quasi-euclidean'); %distance from other segs
        
        Dtot=round((Da+Db).*32)./32; %distance between current seg and other segs
        path = imregionalmin(Dtot); 
        path=bwmorph(path,'thin','inf'); %new path
        Skel3=Skel3 | path;
    end
    Skel=Skel3;
    clear Skel2 Skel3 D1 Db Dtot path
end

m=1:size(DEM,1);
n=1:size(DEM,2);
[n2,m2]=meshgrid(n,m);

%% setup centerline seeds
%extend skeleton to terminus pixel
D1=bwdist(term,'quasi-euclidean').*dx; %distance from the terminus
D2=bwdist(Skel,'quasi-euclidean').*dx; %distance from the skeleton
D=round((D1+D2).*32)./32;
paths = imregionalmin(D);
pathterm=bwmorph(paths,'thin','inf');
seeds=Skel|pathterm|term;

termDist=bwdistgeodesic(imdilate(seeds,strel('square',2)),term,'quasi-euclidean').*dx;
dists=double([termDist(seeds)]);
inds=find(seeds);
mi=m2(inds);ni=n2(inds);

% % nodedist=bwdist(Skel); %distance from skeleton
% 
% 
% [seeddist] = graydist(dmask,(seeds>0),'quasi-euclidean'); %cost-distance with heavy weight for mask
% 
% % zones = iseed.*uint32(mask);
% seeddist(isnan(seeddist))= 1000;
% zones = uint16(watershed(seeddist)).*uint16(mask);
% 
% [~, iseed] = bwdist((seeds>0).*mask);

% outdists=0.*mask;
F=scatteredInterpolant(ni(:),mi(:),dists(:),'linear','linear');
try
    outdists=reshape(F(n2(:),m2(:)),size(DEM)).*mask;
    outdists(isnan(mask2))=NaN;
catch ME
    outdists=NaN;
end
