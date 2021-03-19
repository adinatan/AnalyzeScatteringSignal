load('./data/xcsx39718_Run0165.mat')
 
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
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
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

text(xlim(2)*1.5,mean(ylim),found_circles_str,'FontSize',12);

% note that the x and y are in the counter-intutive order ,
% this is to conform with the LCLS on screen python plots
% in the beam time vs the matlab way to flatten arrays.
% The text on the plot is correct.





%----------------------------------------------------------------------
% Aux functions used in the code:


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


