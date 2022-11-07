function [tFL,Ds,Afluxes,ThxAs,vels,speeds,zs,wdth]=deb_through_fluxes_distance(mask,Tdist,U,V,DT,dx,dem)
% deb_through_fluxes - Function to calculate volumetric fluxes of debris with 
% respect to distance from terminus. 
%
% Syntax:  [ELA,FLout,cFLu,tFL,ELs,Afluxes]=through_fluxes(mask,DEM,U,V,DT,dx)
%
% Inputs:
%    mask - binary raster mask of glacier extent
%    dist - double-precision raster of elevation with same extent as mask
%    U - double-precision raster of x-component of surface velocity with same extent as mask
%    V - double-precision raster of y-component of surface velocity with same extent as mask
%    DT - double-precision raster of debris thickness with same extent as mask
%    dx - pixel size in porjected units corresponding to veloctiy and DEM units (m)
%
% Outputs:
%    tFL - raster of volumetric flux at each pixel
%    Ds - array of elevations at 25m intervals
%    Afluxes - array of volumetric debris fluxes through Ds
%    ThxAs - array of debris cross sections at Ds

% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% June 2020; Last revision: 01-July-2020

%% initialize distances
Tdist(Tdist==Inf)=NaN;
Tdist(Tdist==0)=NaN;
Ds=0:100:max(Tdist(:)); %set up segmentation

%% calculate total speed and flux
S=sqrt(U.^2+V.^2); %calculate speed
tFL =S.*DT; %total flux at any point, per meter cross section!! area accoutned for with trapz below
tFL(isnan(tFL))=0;
DT(isnan(DT))=0;

%% convert flux to perpendicular to 'distance' gradient - ie downglacier or not
[Fx,Fy] = gradient(Tdist);
Fx(mask==0)=NaN;
Fy(mask==0)=NaN;
fx = Fx./sqrt(Fx.^2+Fy.^2);
fy = Fy./sqrt(Fx.^2+Fy.^2);
cFLu = dot([U(:),-V(:)]',[-fx(:),-fy(:)]',1); %down-glacier velocity
cFLv = dot([U(:),-V(:)]',[fy(:),-fx(:)]',1); %cross-glacier velocity
cFLu=reshape(cFLu,size(U));
cFLv=reshape(cFLv,size(U));
cFLu(isnan(cFLu))=0;
    
%% contour distance to terminus as debris flux gates
Afluxes=nan(1,length(Ds));
ThxAs=Afluxes;
vels=Afluxes;
zs=Afluxes;
wdth=Afluxes;
for iD=1:length(Ds)
    M=contourc(Tdist,[Ds(iD),Ds(iD)]);
    if numel(M)>4
        [x,y]=C2xyz(M);
        fluxes=nan(1,length(x));
        DcxA=nan(1,length(x));
        vel=nan(1,length(x));
        speed=nan(1,length(x));
        segDist=nan(1,length(x));
        z=nan(1,length(x)); % ***
        for iL=1:length(x)
            curx=x{iL}';
            cury=y{iL}';
            if length(curx)>1
                segFL=interp2(cFLu.*DT,x{iL},y{iL});
                segS=interp2(S,x{iL},y{iL});
                segU=interp2(cFLu,x{iL},y{iL});
                segDT=interp2(DT,x{iL},y{iL});
                segZ=interp2(dem,x{iL},y{iL}); % ***
                delx = diff(x{iL});
                dely = diff(y{iL});
                dist = sqrt(delx.^2+dely.^2).*dx;
                cdist = [0,cumsum(dist)];
                fluxes(iL) = trapz(cdist,segFL);
                DcxA(iL) = trapz(cdist,segDT);
                speed(iL) = nanmean(segS); %output speed
                vel(iL) = nanmean(segU);
                segDist(iL) = max(cdist,[],'omitnan');
                z(iL) = nanmean(segZ); % ***
            else
                fluxes(iL)=NaN;
                DcxA(iL)=NaN;
                speed(iL)=NaN;
                vel(iL)=NaN;
                segDist(iL)=NaN;
                z(iL)=NaN; % ***
                cdist=NaN;
            end
        end
    else
        fluxes=NaN;
        DcxA=NaN;
        vel=NaN;
        speed=NaN;
        segDist=NaN;
        z=NaN; % ***
        cdist=NaN;
    end
    Afluxes(iD)=nansum(fluxes);
    ThxAs(iD)=nansum(DcxA);
    vels(iD)=nansum(vel.*segDist)./nansum(segDist);
    speeds(iD)=nansum(speed.*segDist)./nansum(segDist);
    zs(iD)=nanmean(z); % ***
    wdth(iD) = cdist(end);
end

end
