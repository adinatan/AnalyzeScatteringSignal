
clear all
close all
clc

global im10k2M x2M y2M  Polarization Geometry smask lpeak q



load('xcsx39718_Run0089.mat')
load('ref_acn.mat')

x2M=double(x2M)*1e3;
y2M=double(y2M)*1e3;


smask=ones(size(im10k2M));
smask(192:193,:,:)=NaN;%smask(end-4:end,:,:)=NaN;
smask(:,176:177,:)=NaN;%smask(:,end-4:end,:)=NaN; 
%%
close all

 T=im10k2M;
 for n=1:size(T,3)
     %T(:,:,n)=  (medfilt2(T(:,:,n).*smask(:,:,n),[3 3]).*smask(:,:,n));
      Te(:,:,n) =  entropyfilt(T(:,:,n),[0 1 0; 1 1 1; 0 1 0]);%ones(3));
     
     Ts(:,:,n) = stdfilt(T(:,:,n),[0 1 0; 1 1 1; 0 1 0]);
     Tm=nanmean(reshape(T(:,:,n),1,[]));
     %Tm(:,:,n) = conv2(T(:,:,n).*smask(:,:,n),ones(3)./9,'same');
     %stdmask(:,:,n) = T(:,:,n)<(Tm(:,:,n)+3.*Ts(:,:,n)) & ...
     %             T(:,:,n)>(Tm(:,:,n)-3.*Ts(:,:,n)) ;
              
      Toutlier(:,:,n) = T(:,:,n) >= (Tm + 3*Ts(:,:,n)) | T(:,:,n) <= (Tm - 3*Ts(:,:,n)); 
     
      statmask(:,:,n)= (  Toutlier(:,:,n) &  ~Te(:,:,n) );

 end
 
 
 
 %%
 for n=1:size(T,3)
    
  
 end
 
 % T=T(:);
%     
% clear H
% xr=linspace(min(x2M(:)),max(x2M(:)),667);
% yr=linspace(min(y2M(:)),max(y2M(:)),667);
% 
%   
%  [im,I,J]=hist3d(x2M(:).*smask(:), y2M(:).*smask(:),T(:).*smask(:),xr,yr);
%  imagesc(im)
%  
 
% find center:
%%
close all
binratio=100;
im=AN_get2Ddata2(x2M(:),y2M(:),im10k2M(:).*smask(:).*statmask(:),binratio);
im(im==0)=NaN;
imx=AN_get2Ddata2(x2M(:),y2M(:),x2M(:),binratio);
imy=AN_get2Ddata2(x2M(:),y2M(:),y2M(:),binratio);

imx(imx==0)=NaN;
imy(imy==0)=NaN;

norm_i_edges=linspace(0.1,1,5e2);
yhist=histcounts(im(:)./max(im(:)),norm_i_edges);
xhist=norm_i_edges(1:end-1)+mean(diff(norm_i_edges))/2;

subplot(1,2,1)
imagesc(im);colormap(flipud(bone))

subplot(1,2,2)
plot(xhist,yhist);
xlabel('norm. I')
ylabel('counts')

% looking at the hist intensities we can test thresholds around 0.4-0.6

%%


%%
imn=im./max(im(:));

k=0;
clear circles c  
figure; 
nx=nanmean(imx,2);
ny=nanmean(imy,1);
imagesc(imn); hold on;colormap(flipud(bone))
  
for w=0.6:0.025:0.8
    k=k+1;
    clear M Mx My
    M=imn<w+0.025 & imn>w   ;
    [My, Mx]=find(M);
    
    
    circles = find_circles_ransac(Mx,My);
    c(k,:)=circles;
    
    %c2(:,k)=  AGSM_circle_fitting(M)';
    
    if 1
       figure
imagesc(imn); hold on;colormap(flipud(bone))

        %cmp=colormap(flag(size(circles,1)));
        for n=1:size(circles,1)
            viscircles(circles(n,1:2),circles(n,3),'LineWidth',0.5,'Color','r' );
            text( [circles(n,1)-circles(n,3),circles(n,1)-circles(n,3)],...
                [circles(n,2),circles(n,2)],num2str(n),'Color','b','FontSize',10);
            
            %     viscircles([AGSMcircles(1) AGSMcircles(2)],AGSMcircles(3),'LineWidth',0.5,'Color','g' );
            
        end
        
    end
    
end
disp( '.............................................................')
disp(['x = ' num2str(mean(c(:,1))), '   std = '  num2str(std(c(:,1)))])
disp(['y = ' num2str(mean(c(:,2))), '   std = '  num2str(std(c(:,2)))])
disp(['r = ' num2str(mean(c(:,3))), '   std = '  num2str(std(c(:,3)))])
 

%% play in pixel units:
x0=round(mean(c(:,1)));
y0=round(mean(c(:,2)));
rx=min(x0,size(im,1)-x0)-1;
ry=min(y0,size(im,2)-y0)-1;

imc=im(x0-rx:x0+rx,y0-ry:y0+ry);
qira=cart2tripolar(im,1:4); %cart2tripolar is a newer version of tripol
qI=nanmean(nanmean(qira,3)' )';
figure;
plot(qI);
 


%% translate units:

% idxgood=~(isnan(imx)); 
% [Ax,Ay]  = meshgrid(1:size(im,2),1:size(im,1));  
% Cx = griddata( Ax(idxgood),Ay(idxgood),imx(idxgood),  mean(c(:,1)),  mean(c(:,2)) ) ;
% Cy = griddata( Ax(idxgood),Ay(idxgood),imy(idxgood),  mean(c(:,1)),  mean(c(:,2)) ) ;
%  
% disp(['x = ' num2str(Cx)])
% disp(['y = ' num2str(Cy)]) 

Cx=interp1(1:numel(nx),nx,mean(c(:,1)),'spline');
Cy=interp1(1:numel(ny),ny,mean(c(:,2)),'spline');


%% find z:

clear init_vals
options = optimset('MaxFunEvals',2000,'MaxIter',2000,'TolX',1e-5,'TolFun',1e-5,'Display' , 'iter');
%options = optimset( 'Display' , 'iter');
tic 
% init_vals=[init_det_dist init_x init_y init_scale];
init_vals=[39e3 Cx Cy 0.7e5 0];

LB = [ 38e3 Cx-100 Cy-100  0.6e5 -5e5];
UB = [ 40e3 Cx+100 Cy+100 0.8e5 5e5 ];

[x,fval,exitflag] = fminsearchbnd(@opt_dist_ACNmodel ,init_vals,LB,UB,options);
toc

%% check what we got:
%x=[39.1e3 Cx  Cy  0.706e5 -1e5];
detdist = x(1);% in mm
 
photonenergy = 16.8; %in KeV
clear b0 b0a
 
[TwoTheta,Geometry,Polarization] = ...
    CorrectionsAngleDependant_noabs([x2M(:) y2M(:)],detdist,x(2),x(3));

QQ = 4*pi*sind(TwoTheta/2).*photonenergy./12.3987;

QQ2=AN_get2Ddata2(x2M(:),y2M(:),QQ(:).*smask(:),binratio);
QQ2(QQ2==0)=NaN;

% Polarization2=AN_get2Ddata2(x2M(:),y2M(:),Polarization(:).*smask(:),binratio);
% Polarization2(QQ2==0)=NaN;
% figure
% imagesc(Polarization2)
% 
% Geometry2=AN_get2Ddata2(x2M(:),y2M(:),Geometry(:).*smask(:),binratio);
% Geometry2(QQ2==0)=NaN;
% figure
% imagesc(Geometry2)
%mask = smask(:).*Geometry(:).*Polarization(:);

mask = smask(:).*Geometry(:);
 
qbinedges=0.4:0.05:8;


[~,~,qbins] = histcounts(QQ.*smask(:),qbinedges);
xq=qbinedges(1:end-1)+mean(diff(qbinedges))/2;



im2=im10k2M(:).*mask; 
for n=1:max(qbins)
  b0=im2(qbins==n);
   b0a(n)=nanmean(b0(b0>0)); % mean intensity
 % b0a(n)=trimmean(b0(b0>0),1);
end

St=interp1(q(~isnan(lpeak)) ,lpeak(~isnan(lpeak)),xq,'spline');
%
figure
plot(xq,St.*x(end-1)+x(end),'LineWidth',1); hold on
plot(xq,b0a-4e5,'LineWidth',1);
xlim([0 7.72])
legend('Ref','Measured Run 79')
xlabel('Q [A^{-1}]');
set(gca,'FontSize',14)

%end
return
%% 
Q2d=AN_get2Ddata2(x2M(:)+x(2),y2M(:)+x(3),QQ(:).*smask(:),binratio);
Q2d(im==0)=NaN;


q0=Q2d;
q0n=q0;
q0n(q0n==0)=NaN;
q0n=inpaint_nans(q0n,1);
[x0 y0]=find(q0n==min(q0n(:))); % assume min Q is center
%crop around center:
rx=min(x0,size(q0n,1)-x0)-1;
ry=min(y0,size(q0n,2)-y0)-1;



qc=q0n(x0-rx:x0+rx,y0-ry:y0+ry);
qira=cart2tripolar(qc,1:4); %cart2tripolar is a newer version of tripol
qira0=cart2tripolar(q0n,1:4); %cart2tripolar is a newer version of tripol

clear p qr qcr
qr=nanmean(nanmean(qira,3)' )';
r1=1;
r2=20;
p=polyfit( [ r1:numel(qr)-r2] , [ qr(r1:end-r2)]'  ,3);
qcr=polyval(p,1:numel(qr));
plot(qcr);hold on;plot(qr); title(num2str(qcr(1)))


%
% 
% x0=round(mean(c(:,1)));
% y0=round(mean(c(:,2)));
% rx=min(x0,size(im,1)-x0)-1;
% ry=min(y0,size(im,2)-y0)-1;

%im=AN_get2Ddata2(x2M(:)-x(2),y2M(:)-x(3),im10k2M(:).*smask(:).*Polarization(:).*Geometry(:),binratio);
im=AN_get2Ddata2(x2M(:)-x(2),y2M(:)-x(3),im10k2M(:).*smask(:),binratio);

im(im==0)=NaN;

imc=im(x0-rx:x0+rx,y0-ry:y0+ry);
ira=cart2tripolar(im,1:4); %cart2tripolar is a newer version of tripol
ir1=nanmean(nanmean(ira,3)' )';
%%
figure;
plot(qcr,ir1(1:numel(qcr))); hold on
plot(xq,St.*x(end)); hold on
xlim([0.5 7])


%%

x0=ceil(size(im,1)/2); y0=ceil(size(im,2)/2); % Assuming image is centered
L = min([x0,y0]);
RR=(0:L);
PPR=(floor(0.5*pi*(RR+1))-1); % the # of pixels per radius for cartesian quadrant
AngleInc = (0.5*pi./PPR'); % the angle increment per radius
AngleInc(1)=0; % avoid inf at origin
% create image quadrants

%set radius of square matrix around center
L = min([x0,y0,(size(im) - [x0,y0])]);

%create  image quadrants
if mod(size(im,1),2)==1 %case for even image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0:y0+L-1);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0:x0+L-1   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0:x0+L-1   ,y0:y0+L-1);
else %case for odd image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0+1:y0+L);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0+1:x0+L   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0+1:x0+L   ,y0+1:y0+L);
end
