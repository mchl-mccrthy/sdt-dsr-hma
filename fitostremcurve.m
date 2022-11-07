% Fit Ostrem curve 

% This function fits an Ostrem curve to specific mass balance and
% supraglacial debris thickness data, as described by McCarthy et al (2022)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function [dtFit,piLookupSMB,piLookupDT,r2] = fitostremcurve(smbMod,...
    dtSamps,startPt)

% Delete NaNs... are there any?!
toDelete = isnan(smbMod(:)) | isnan(dtSamps(:));
smbMod(toDelete) = [];
dtSamps(toDelete) = [];

% Sort
[dtSamps,sortInd] = sort(dtSamps);
smbMod = smbMod(sortInd);

% Specify fit type and fit curve
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-12,0],...
    'Upper',[0,Inf]);
fitType = fittype('c1*(c2./(c2+x))','options',fo);
[dtFit,gofDt] = fit(dtSamps,smbMod,fitType,'StartPoint',startPt);
r2 = gofDt.rsquare;

% Get prediction intervals
piLookupDT = 0:0.001:5;
piLookupSMB = predint(dtFit,piLookupDT,0.68,'observation','off');

end