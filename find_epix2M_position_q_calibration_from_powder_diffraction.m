clear all

%% read h5:
fpath='./h5/'
%xcslv9618_Run0010
runnum=9;
file=[fpath 'xcslv9618_Run' num2str(runnum,'%0.4i') '.h5'];
im10k2M=h5read(file,'/Sums/epix10k2M_calib');
x2M=h5read(file,'/UserDataCfg/epix10k2M/x');
y2M=h5read(file,'/UserDataCfg/epix10k2M/y');
z2M=h5read(file,'/UserDataCfg/epix10k2M/z');
photon_energy=h5read(file,'/ebeam/photon_energy');
diodeGon=h5read(file,'/diodeGon/channels');

ped=h5read(file,'/UserDataCfg/epix10k2M/ped');
mask0=h5read(file,'/UserDataCfg/epix10k2M/mask');
x2M=double(x2M);
y2M=double(y2M);
z2M=double(z2M);

px=x2M;
py=y2M;
pz=z2M;

%% quick masking
smask=ones(size(im10k2M));
% central lines of modules \ tiles
smask(192:193,:,:)=NaN;
smask(:,176:177,:)=NaN;
smask(end-4:end,:,:)=NaN;
smask(:,end-4:end,:)=NaN;
smask(1:5,:,:)=NaN;
smask(:,1:5,:)=NaN;
% apply mask:
im=(double(im10k2M)).*double(mask0).*smask;

% some manual masks for shadows in spesific modules
im(25:50,177:end,5)=0;
im(80:110,177:end,5)=0;
im(160:190,177:end,5)=0;
im(300:end,1:150,11)=0;
im(im<1)=NaN;

% prep image to threshold
imc=zeros(size(im));

% local threshold is established looking at at blurred avg intensity in
% each of the banks. The idea is that the powder diffraction is more
% intenstese and localized and so a blurred wider neighboroud will capture
% the background.
% the blurring is done via convolution with a filter such as a gaussian.


%filt=ones(19); % size for moving window avg
filt=fspecial('gaussian',29,5); % size for moving window avg
filt=filt./sum(filt(:)); %

% subtract local avg to level the intensity bkg
for k=1:size(im,3)
    for n=1:8
        clear banksL banksR
        % def banks (left and right, because it's easier)
        banksL(:,:,n)= im(n*48-47:n*48 ,1:176,k);
        banksR(:,:,n)= im(n*48-47:n*48 ,177:end,k);

        % Left side
        d2=nanconv(banksL(:,:,n),filt, 'nanout','edge' );
        banksL(:,:,n)=banksL(:,:,n)-d2;
        imc(n*48-47:n*48 ,1:176,k)=banksL(:,:,n);

        % Right side
        d2=nanconv(banksR(:,:,n),filt, 'nanout','edge' );
        banksR(:,:,n)=banksR(:,:,n)-d2;
        imc(n*48-47:n*48 ,177:end,k)=banksR(:,:,n);
    end
end

% threshold globaly
imc(imc<1)=NaN;
imct=imc(:)>nanmean(imc(:))+0.1*std(imc(~isnan(imc(:)))) ;
imct=reshape(imct,size(imc));

% plotting:
%------------------------------------------------
binratio=100;
nim0=get2Ddata(px(:),py(:),im(:) ,binratio);  % masked img
nim=get2Ddata(px(:),py(:),imc(:) ,binratio);  % local-bkg subt img
nimBW=get2Ddata(px(:),py(:),imct(:)  ,binratio); % thresholded img

imx=get2Ddata(px(:),py(:),px(:),binratio);
imy=get2Ddata(px(:),py(:),py(:),binratio);
imz=get2Ddata(px(:),py(:),pz(:),binratio);

% mask negetive values:
nim0(nim0<0)=NaN; nim(nim<0)=NaN; imx(imx==0)=NaN; imy(imy==0)=NaN;

% x-y vectors for plotting
nx=nanmean(imx,2); ny=nanmean(imy,1);

figure(1)
subplot(1,3,1);
imagesc(ny,nx, nim0 )
inten_limit=10*mad(nim0(:));
caxis([0 inten_limit]);
title('im')
axis square

subplot(1,3,2);  imagesc(ny,nx, nim )
caxis([0 inten_limit])
title('imc')
axis square

subplot(1,3,3);  imagesc(ny,nx, nimBW )
axis square
title('imct')
caxis([0 1 ])
colormap([1 1 1; jet(256) ])
%------------------------------------------------

%% find circles in thesholded img via ransac to get initial estimate of center:
clear c  M
% set ransac parametres
sample_size = 3;
max_distance = 1000;
pointmask=0.*imct(:);
figure
% we scan for circles that contain more inliers, then take these points
% out to scan for circles that contain less inliers...
mnpic=1000; % min_num_points_in_circle vector

tic % this can take a bit time... (60 sec)
for k=1:numel(mnpic)
    min_num_points_in_circle = mnpic(k);
    hold on

    M = find(imct(:).*~pointmask(:));
    c{k} = FindCirclesRansac(py(M),px(M),sample_size,max_distance,min_num_points_in_circle);

    %remove inlier points from found circles range for the next round of  min_num_points_in_circle
    if ~isempty(c{k})
        for m=1:size(c,1)
            pointmask= pointmask | abs(sqrt ( (px(:)-c{k}(m,2)).^2+(py(:)-c{k}(m,1)).^2)-c{k}(m,3))<2*max_distance;
        end
    end
end
toc

% plot what we found:
M = find(imct(:));
scatter(py(M),px(M),1,[0.5 0.5 0.5]); hold on

cc=vertcat(c{:});
[val, id]=sort(cc(:,3));

for n=1:numel(id)
    viscircles(cc(id(n),1:2),cc(id(n),3),'LineWidth',0.5,'Color',[1 0 0] );hold on
    text(cc(id(n),1)-cc(id(n),3), -1200+cc(id(n),2),num2str(n),'Color',[1 0 0],'FontSize',7);
end


% write found x,y on screen
yl=get(gca,'ylim');
xl=get(gca,'xlim');

[xc, std_xc]= fractrimmean(cc(:,2),37); % fractrimmean allows to calc the std as well
[yc, std_yc]= fractrimmean(cc(:,1),37);

% note that the x and y are in the counter-intutive order ,
% this is to conform with the LCLS on screen python plots
% in the beam time vs the matlab way to flatten arrays.
% The text on the plot is correct.

found_circles_str=...
    {['x = ' num2str(xc,'%10.1f\n'), ' \pm '  num2str(std_xc,'%10.1f\n')  ,'\mum' ];
    ['y = ' num2str(yc,'%10.1f\n'), ' \pm '  num2str(std_yc,'%10.1f\n') ,'\mum'] };

dim = [.65  .6  .3 .3];
annotation('textbox',dim,'String',found_circles_str,'FitBoxToText','on','FontSize',15);

xlim([min(nx) max(nx)])
ylim([min(ny) max(ny)])

ylim([xc-1.1*max(cc(:,3)) xc+1.1*max(cc(:,3))])
xlim([yc-1.1*max(cc(:,3)) yc+1.1*max(cc(:,3))])

title('run 10 - 2M detector')
uitable(gcf,'Data',cc,'Position',[1250,110,220,850],'ColumnName',{'Yc','Xc','R'});

axis square

%% read calculated powder diffraction signal
d=dlmread('LaB6_mp-2680_computed_18kev.txt'); % read ref file
photonenergy = 18; %in KeV

clear PDref q_vs_e
uq=0.001:0.001:15; % q vector
q_vs_e   =  (4*pi / (12.3984193/photonenergy)) .*  sind(d(:,1)./2); % q's of
PDref =interp1( q_vs_e ,d(:,2),uq,'spline');

% just pick out the peak locations and amplitudes:
[PDpks,PDlocsq] = findpeaks(PDref,uq,'MinPeakProminence',0.015);
PDpks(1)=[];
PDlocsq(1)=[];
plot(uq,PDref);hold on;
plot(PDlocsq,PDpks,'x');

xlim([0  15])
ylim([0 max(PDref) ])
Ang =char(197); % symbol for angstrom
xlabel(['Q (' Ang '^{-1})'])
title('LaB6 calculated diffraction')

%% find detector ceter and distance, hence Q and the detector corrections
% using bounded minimization of the radial image vs a powder reference:
photonenergy=trimmean(photon_energy-224.0935,50)*1e-3;
clc
% increase the tolerances to check convergance
options = optimset('MaxFunEvals',10,'MaxIter',10,'TolX',1e-1,'TolFun',1e-1,'Display' , 'iter');
%intial values for the optimization are:
% [det distance in um, detector center (x, y) in um, scaling,  broadeing factors for
% the reference:
init_vals=[80.1e3 xc yc 1e4 0.025];
LB = [ 79.5e3     xc-100  yc-250    0 0.01 ];
UB = [ 80.5e3     xc+100  yc+250  5e4 0.05 ];
pxyz =[px(:) , py(:),  (nanmean(pz(:))-pz(:))];
dq=0.01;
% prep parameters to pass to the costfun
optparams.pxyz=pxyz;
optparams.im=im;
optparams.uq=uq;
optparams.PDref=PDref;
optparams.photonenergy=photonenergy;
optparams.dq=dq;

tic
[x,fval,exitflag] = fminsearchbnd(@costfun ,init_vals,LB,UB,options,optparams);
toc

%% plot found fit:
clear imcorrected xq qbinedges refi b0a b0 im2 p  QQ x1

p = ([(pxyz(:,1)-x(2)) pxyz(:,2)-x(3)  pxyz(:,3)+x(1) ]);
Lp2=(p(:,1).^2+p(:,2).^2+p(:,3).^2); % squared scattering vector

% Two theta angle to Q
TwoTheta=acosd(p(:,3)./realsqrt(Lp2));
QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3984193;

% Polarisation correction
Pol=0; % Polarisation 0:1, 0=vertical polarisation
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(Lp2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));

% Geometry correction ( Eqs. 6 in Boesecke & Diat)
Geometry=1./( (p(:,3).^2./Lp2).*(p(:,3)./realsqrt(Lp2))  );

% make correction mask:
mask = ~isnan(im(:)) .*Geometry(:).*Polarization;

% bin detector QQ array in the radial q (called xq):
qbinedges=0.9:dq:8.4;
[~,~,qbins] = histcounts(QQ ,qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;
% detector intensity in xq bins prep:
b0a=zeros(1,numel(xq));

%apply correction mask
im2=im(:).*mask.*~isnan(im(:));

% add detector intensities to xq bins and average
for n=1:max(qbins)
    b0=im2(qbins==n);
    if all(isnan(b0)) % check for completly missing intensity at q bin
        b0a(n)=NaN;
    else
        [b0a(n), b0a_std(n)]= fractrimmean(b0(b0>0),37); %take the trimmed mean
    end
end

if sum(isnan(b0a))
    nr=~isnan(b0a);
    b0a(find(~nr))=interp1(find(nr),b0a(nr),find(~nr),'nearest','extrap');
end

%remove baseline:
imcorrected= b0a(:)-arPLS_baseline(b0a,1e5,1e-6);
xqi=xq;

% broaden peaks of reference powder to the resolution of the measurment
sigmaK=x(5);
K=exp(-(uq-mean(uq)).^2/(2*sigmaK^2));
K=K./sum(K(:));
PDref_broad=conv(PDref,K,'same')*x(4);
refi=interp1(uq' ,PDref_broad,xqi,'makima');

% plot
figure('Position',[100 100 800 900])

subplot(2,1,1)
plot(xq,b0a,xq,arPLS_baseline(b0a,1e5,1e-6),'LineWidth',2)
legend('Raw integrated signal','baseline (arPLS)','FontSize',18);
set(gca,'FontSize',18)

xlim([0  10])
subplot(2,1,2)
plot(xqi,imcorrected, 'LineWidth',2) ;hold on;
plot(xqi,refi, 'r:','LineWidth',2)
set(gca,'FontSize',18)
ylim([0 2e5])
xlim([0  8.7])


xlabel('Q [A^{-1}]')
h=legend('Measured signal (baseline subtracted)','Reference','FontSize',18);
dim=[h.Position(1),h.Position(2)-h.Position(4),h.Position(3),h.Position(4)];

words={['x = ' num2str(x(2),'%10.1f\n'), ' \pm '  num2str(std_xc,'%10.1f\n' ) ,'\mum' ] , ...
    ['y = ' num2str(x(3),'%10.1f\n'), ' \pm '  num2str(std_yc,'%10.1f\n') ,'\mum'] ,  ...
    ['Z = ' num2str(x(1),'%10.1f\n') ' \mum' ]};
annotation('textbox',dim,'String',words,'FitBoxToText','on','FontSize',14);
set(gca,'FontSize',18)


%----------------------------------------------------------------------
% End of code
%----------------------------------------------------------------------


%% Aux functions used in the code:
function cost=costfun(x,optparams)


dq=optparams.dq;
pxyz=optparams.pxyz;
im=optparams.im;
uq=optparams.uq;
PDref=optparams.PDref;
photonenergy=optparams.photonenergy;


p = ([(pxyz(:,1)-x(2)) pxyz(:,2)-x(3)  pxyz(:,3)+x(1)]);
Lp2=(p(:,1).^2+p(:,2).^2+p(:,3).^2); % squared scattering vector


% Two theta angle to Q
TwoTheta=acosd(p(:,3)./sqrt(Lp2));
QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3984193;

% Polarisation correction
Pol=0; % Polarisation 0:1, 0=vertical polarisation
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(Lp2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));

% Geometric correction ( Eqs. 6 in Boesecke & Diat)
Geometry=1./( (p(:,3).^2./Lp2).*(p(:,3)./realsqrt(Lp2))  );

mask = ~isnan(im(:)) .*Geometry(:).*Polarization;
qbinedges=0.9:dq:8.4;
[~,~,qbins] = histcounts(QQ ,qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;
b0a=zeros(1,numel(xq));

im2=im(:).*mask.*~isnan(im(:));

for n=1:max(qbins)
    b0=im2(qbins==n);
    if all(isnan(b0))
        b0a(n)=NaN;
    else
        %b0a(n)=trimmean(b0(b0>0),37);
        [b0a(n), b0a_std(n)]= fractrimmean(b0(b0>0),37); %take the trimmed mean
    end
end

if sum(isnan(b0a))
    nr=~isnan(b0a);
    b0a(find(~nr))=interp1(find(nr),b0a(nr),find(~nr),'nearest','extrap');
end

%remove baseline:
imcorrected= b0a(:)-arPLS_baseline(b0a,1e5,1e-6);
xqi=xq;

%broaden ref fit
sigmaK=x(5);
K=exp(-(uq-mean(uq)).^2/(2*sigmaK^2));
K=K./sum(K(:));
PDref_broad=conv(PDref,K,'same')*x(4);
refi=interp1(uq' ,PDref_broad,xqi,'spline');
qrange=xqi>1.45 & xqi<8.4;
weights =imcorrected(:)./b0a_std(:);
cost= 0.5*nansum( weights(qrange).*abs(refi(qrange)' - imcorrected(qrange) ).^2 )./nansum(qrange);

end

function circles = FindCirclesRansac(x,y,sample_size,max_distance,min_num_points_in_circle)
% parameters (may need some modification depends on data)
% sample_size = 3;
% max_distance = 1e3;
% min_num_points_in_circle = 30;

circles = {};
counter = 0;
while numel(x) > 1000 * sample_size && counter < 10
    try
        [circle, inlierIdx] = ransac([x, y], @fit_circle, @distFcn, ...
            sample_size, max_distance, 'MaxNumTrials',10000);
    catch
        break
    end
    % refit using only inliers
    circle = fit_circle([x(inlierIdx) y(inlierIdx)]);
    distf = distFcn(circle, [x y]);
    founfit = distf < max_distance;
    if sum(founfit) > min_num_points_in_circle
        % this model fits enough points
        circles{end+1} = circle;
        x(founfit) = [];
        y(founfit) = [];
    else
        counter = counter + 1;
    end
end
circles = vertcat(circles{:});
end

function crc = fit_circle(xy)  % xy is n-by-2 matrix of point coordinates
% fit in least squares sens
x = xy(:, 1);
y = xy(:, 2);
X = [-2*x, -2*y, ones(size(x))];
Y = -x.^2 - y.^2;
crc = (X\Y).';  % least squares solution
r2 = -crc(3) +crc(1).^2 + crc(2).^2;
if r2 <= 0
    crc(3) = 0;
else
    crc(3) = sqrt(r2);
end
% output crc is a 3 vector (cx, cy, r)
end

function distf = distFcn(crc, xy)
% how good a fit circle for points
x = xy(:, 1) - crc(1);
y = xy(:, 2) - crc(2);
distf = abs(sqrt(x.^2 + y.^2) - crc(3));
end

function d = get2Ddata(x, y, vec,N)
%This function loads an input of the format [x y int] and
%generates a matrix of int spanned by x and y.

%Resetting zeros & Dividing by N
x = ceil( (x - (min(x(:))))/N)+1;
y = ceil( (y - (min(y(:))))/N)+1;

d = zeros(  max(x(:))-min(x(:))+1 , max(y(:))-min(y(:))+1);

for n=1:numel(vec)
    d(x(n),y(n))=vec(n);
end
end


function baseline = arPLS_baseline(signal, smoothness_param, min_diff)
% implements Baseline correction using adaptive iteratively reweighted
% penalized leastsquares (10.1039/B922045C)

% defaults
if nargin<2
    smoothness_param = 1e3;
    min_diff = 1e-6;
end
signal = signal(:);

order = 2; % Difference filter order
N = 100; % Maximum # of iterations

L = numel(signal);
difference_matrix = diff(speye(L), order);
minimization_matrix = (smoothness_param*difference_matrix')*difference_matrix;
penalty_vector = ones(L,1);

for count = 1:N
    penalty_matrix = spdiags(penalty_vector, 0, L, L);
    % Cholesky decomposition
    C = chol(penalty_matrix + minimization_matrix);
    baseline = C \ (C'\(penalty_vector.*signal));
    d = signal - baseline;
    % make d-, and get penalty_vector^t with mean and std
    dn = d(d<0);
    penalty_vector_temp = 1./(1+exp(2*(d-(2*std(dn)- mean(dn)))/std(dn)));
    % check exit condition and backup
    if norm(penalty_vector-penalty_vector_temp)/norm(penalty_vector) < min_diff
        %         count
        break;
    end
    penalty_vector = penalty_vector_temp;
end
end