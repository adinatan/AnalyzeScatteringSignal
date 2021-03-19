clear all
close all
clc

load('./data/xcsx39718_Run0089.mat')
 
x2M=double(x2M)*1e3;
y2M=double(y2M)*1e3;

smask=ones(size(im10k2M));
smask(192:193,:,:)=NaN;
smask(:,176:177,:)=NaN;

smask(end-4:end,:,:)=NaN;
smask(:,end-4:end,:)=NaN;

smask(1:5,:,:)=NaN;
smask(:,1:5,:)=NaN;
%% read and analyze data
im=im10k2M(:).*smask(:); 
im=im./nanmax(im(:));
px=x2M(:);
py=y2M(:);

imn=im./max(im);
norm_i_edges=linspace(0.1,1,5e2);
yhist=histcounts(im,norm_i_edges);
xhist=norm_i_edges(1:end-1)+mean(diff(norm_i_edges))/2;

[pks, locs, pwid, prom]=findpeaks(smooth(yhist,10),xhist,'MinPeakProminence',6,'Annotate','extents');

%% find circles according to the intensity theshold:
w=locs(end)-1*max(pwid):min(pwid):locs(end)+1*max(pwid);
k=0;
clear circles c
 
c=[];
for w=locs(end)-1*max(pwid):min(pwid):locs(end)+1*max(pwid)
    k=k+1;
    M{k} = im>w & im<w+max(pwid)   ; % select point within intensity range 
    circles = FindCirclesRansac(py(M{k}),px(M{k}));
    c=[c;circles]; % stack found circles
end
 

%% plot the analysis above:
figure('Position',[50 50 700 700])
subplot(2,2,1)
plot(xhist,yhist,locs(end),pks(end),'s');
xlabel('norm. I')
ylabel('counts')
legend('normalized intensity of image','assesed intensity of ring')
axis square

subplot(2,2,2)
for k=1:numel(M)
    scatter(py(M{k}),px(M{k}),1); hold on
end

cmp=jet(size(c,1));
for n=1:size(c,1)
    viscircles(c(n,1:2),c(n,3),'LineWidth',0.5,'Color',cmp(n,:) );
end
axis square
hold off


subplot(2,2,3)
binratio=100;
nim=get2Ddata(px(:),py(:),im(:),binratio);
imx=get2Ddata(px(:),py(:),px(:),binratio);
imy=get2Ddata(px(:),py(:),py(:),binratio);
nim(nim==0)=NaN; imx(imx==0)=NaN; imy(imy==0)=NaN;
nx=nanmean(imx,2); ny=nanmean(imy,1);
imagesc(ny,nx,nim); hold on;colormap(flipud(bone)) ; hold on
for n=1:size(c,1)
    viscircles(c(n,1:2),c(n,3),'LineWidth',0.5,'Color','r' );
end
axis square


% write found x,y on screen
yl=get(gca,'ylim');
xl=get(gca,'xlim');
if any(std(c(:,1:2))./mean(c(:,1:2))<1e-3)
    xc= mean(c(:,2));
    yc= mean(c(:,1));
else
    xc= trimmean(c(:,2),33);
    yc= trimmean(c(:,1),33);
end

found_circles_str=...
    {['x = ' num2str(xc,'%10.1f\n'), ' \pm '  num2str(std(c(:,2)),'%10.1f\n' ) ,'\mum' ];
    ['y = ' num2str(yc,'%10.1f\n'), ' \pm '  num2str(std(c(:,1)),'%10.1f\n') ,'\mum'] };

text(xl(2)*1.5,mean(yl),found_circles_str,'FontSize',12);

% note that the x and y are in the counter-intutive order ,
% this is to conform with the LCLS on screen python plots
% in the beam time vs the matlab way to flatten arrays.
% The text on the plot is correct.


%% find Q (as well as P,G) with bounded minimization of the image vs a reference:

options = optimset('MaxFunEvals',300,'MaxIter',300,'TolX',1e-3,'TolFun',1e-3,'Display' , 'iter');
%options = optimset( 'Display' , 'iter');
tic 
% init_vals=[init_det_dist scale_fac1 scale_fac2];
init_vals=[39e3 300  0];

LB = [ 38e3 250 -10];
UB = [ 40e3 350  10];

load('./RefDataSets/ref_acn.mat');
p = ([(py(:)-yc) px(:)-xc]);
photonenergy = 16.8; %in KeV

tic
[x,fval,exitflag] = fminsearchbnd(@costfun ,init_vals,LB,UB,options,p,im,ref,photonenergy);
toc

%% plot fit and params:
close all
 
%x(2)=; x(3)=;
p3= x(1)*ones(size(p,1),1);
Lp2=(p(:,1).^2+p(:,2).^2+p3.^2); % squared scattering vector
 
% Polarisation correction
Pol=0; % Polarisation 0:1, 0=vertical polarisation
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(p(:,1).^2+p(:,2).^2+p3.^2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));
 
% Geometric correction ( Eqs. 6 in Boesecke & Diat)
Geometry=1./( (p3.^2./Lp2).*(p3./sqrt(Lp2))  );
 
% Two theta angle to Q
TwoTheta=acosd(p3./sqrt(Lp2));
QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3987;

mask = ~isnan(im(:)).*Geometry(:);
qbinedges=0.4:0.05:8;
[~,~,qbins] = histcounts(QQ.*~isnan(im(:)),qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;

im2=im(:).*mask; 
for n=1:max(qbins)
  b0=im2(qbins==n);
  b0a(n)=nanmean(b0(b0>0)); % mean intensity
 % b0a(n)=trimmean(b0(b0>0),10);
end 
imcorrected=b0a.*x(end-1)+x(end);

refi=interp1(ref(:,1) ,ref(:,2),xq,'spline');
%%
figure
plot(xq,imcorrected,xq,refi,'LineWidth',1) 
xlim([0 7.7])
xlabel('Q [A^{-1}]')
h=legend('Measured signal','Reference','FontSize',14);
dim=[h.Position(1),h.Position(2)-h.Position(4),h.Position(3),h.Position(4)];
annotation('textbox',dim,'String',['Z = ' num2str(x(1),'%10.1f\n') ' \mum'],'FitBoxToText','off','FontSize',14);


%----------------------------------------------------------------------
%% Aux functions used in the code:


function circles = FindCirclesRansac(x,y)
% parameters (may need some modification depends on data)
sample_size = 3;
max_distance = 1e6;
min_num_points_in_circle = 10;

circles = {};
counter = 0;
while numel(x) > 100 * sample_size && counter < 10
    try
        [circle, inlierIdx] = ransac([x, y], @fit_circle, @distFcn, ...
            sample_size, max_distance);
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

%%%%% fminsearchbnd part:
function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
  x0u(k:n) = [];
end

% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end

% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end

% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);

% undo the variable transformations into the original space
x = xtransform(xu,params);

% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);

% Use a nested function as the OutputFcn wrapper
  function stop = outfun_wrapper(x,varargin);
    % we need to transform x first
    xtrans = xtransform(x,params);
    
    % then call the user supplied OutputFcn
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    
  end

end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function

% transform
xtrans = xtransform(x,params);

% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});

end % sub function intrafun end

% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end

%% cost 


function cost=costfun(x,p,im,ref,photonenergy)
 

p3= x(1)*ones(size(p,1),1); % x(1) -detdist in um

Lp2=(p(:,1).^2+p(:,2).^2+p3.^2); % squared scattering vector
 
% Polarisation correction
Pol=0; % Polarisation 0:1, 0=vertical polarisation
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(p(:,1).^2+p(:,2).^2+p3.^2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));
 
% Geometric correction ( Eqs. 6 in Boesecke & Diat)
Geometry=1./( (p3.^2./Lp2).*(p3./sqrt(Lp2))  );
 
% Two theta angle to Q
TwoTheta=acosd(p3./sqrt(Lp2));
QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3987;

mask = ~isnan(im(:)).*Geometry(:);
  
qbinedges=0.4:0.05:8;

[~,~,qbins] = histcounts(QQ.*~isnan(im(:)),qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;

im2=im(:).*mask; 
for n=1:max(qbins)
  b0=im2(qbins==n);
  b0a(n)=nanmean(b0(b0>0)); % mean intensity
 % b0a(n)=trimmean(b0(b0>0),10);
end 

refi=interp1(ref(:,1) ,ref(:,2),xq,'spline');

range_fit=(xq>1.4  & xq<3) ;   

   
%range_fit=(xq>1.45 & xq<2.8)  ;
cost= nansum( abs( (b0a(range_fit).*x(end-1)+x(end))-refi(range_fit) ) ) ./sum(range_fit);


end



