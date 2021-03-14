function cost = opt_dist_ACNmodel(x)

global im10k2M x2M y2M  Polarization Geometry smask lpeak q


 
photonenergy = 16.8; %in KeV

%detdist in um
detdist = x(1);% 

% pathlength = 8e3; % in um
% %pathlength = 0;%x(2); % in um
% jetwidth = 5e6; % in um; complete guess. should only affect absorption mask
% I2Mu = 211.8; % in cm^2/g for 9keV from CXRO
% I2weight = 126.9*2; % in g/mol
% SampleMu = I2Mu*I2weight*sampledensity*1e4; % should be linear attenuation. has converstion for um^-1 from cm^-1


[TwoTheta,Geometry,Polarization] = ...
    CorrectionsAngleDependant_noabs([x2M(:) y2M(:)],detdist,x(2), x(3));

%[TwoTheta,Geometry,Polarization] = ...
%    CorrectionsAngleDependant_xppl2816_noabs([px_x(:) px_y(:)],detdist,x(2),x(3));

QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3987;

mask = smask(:).*Geometry(:);
  
qbinedges=0.4:0.05:8;

[~,~,qbins] = histcounts(QQ.*smask(:),qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;



im2=im10k2M(:).*mask; 
for n=1:max(qbins)
  b0=im2(qbins==n);
  % b0a(n)=nanmean(b0(b0>0)); % mean intensity
  b0a(n)=trimmean(b0(b0>0),10);
end

%b0a(xq<0.35)=0;
 
range_fit=(xq>1.2  & xq<6) ;   

St=interp1(q(~isnan(lpeak)) ,lpeak(~isnan(lpeak)),xq,'spline');
   
%range_fit=(xq>1.45 & xq<2.8)  ;
cost= nansum( abs( b0a(range_fit)-(St((range_fit)).*x(end-1)+x(end) ) ) )./numel(range_fit);

%cost= sum( abs( b0a(range_fit)./b0a(range_fit)-St0(range_fit)./sum(St0,(range_fit))   ) )./numel(range_fit);
