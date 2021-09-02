clear all
close all
clc

% creating a data vector in the range [0,2*pi]

% % generate random vec up to beta 4 + noise
% 
N=64; % # of bins in angle from 0 to 2*pi
rv=rand(1,3);
w=randi(1e3,1,N);
vars=1./(sqrt(w));
vec=1000*(rv(1)-0.5)+... % baseline
    (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
    (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
    vars.*randn(1,N); % noise term

vec2=1000*(rv(1)-0.5)+... % baseline
    (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
    (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
    0*randn(1,N); % noise term

errorbar(vec,vars,'o'); hold on; plot(vec2)

%%
%vec(45:50)=NaN;
%plot(vec)
% say we want to fit up to beta6 :
bParams=[2:2:6];



% helper function to replace the use the leg in the original code:
paren = @(x, varargin) x(varargin{:});


PPR=numel(vec)-1;
AngleInc = 2*pi./PPR'; % angle increment per radius
betas = zeros(numel(bParams)+1,1);


npr=PPR;
qp=0:npr;

% for very small radii we reduce the maximal beta value
% allowed to avoid divergences at the origin

%if npr/2 <= max(bParams)
%    bParams=2:2:npr/2;
%end

y = vec;%ira(r,qp+1); % assign row of all pixels in radius rp
% this may have NaN values, one way to treat such a case is to
% interpolate \ extrapolate NaN values using splines or AR methods.
% We will ignore this solution for now and keep the NaNs as is.

B = zeros(1,numel(bParams)+1);
% this is beta0, including if there are NaN values
notnan=~isnan(y);
B(1)=nansum(y)/sum(notnan)*(npr+1);

% one fit coefficient for each B param
fitCoefs = ones(numel(bParams) + 1, sum(notnan));
for ii=(1:numel(bParams))
    % assign relevant Legendre polynomials to fitCoef
    %   fitCoefs(ii+1,:) = leg(bParams(ii), cos(AngleInc*qp(notnan)));
    fitCoefs(ii+1,:) = paren( legendre(bParams(ii), cos(AngleInc*qp(notnan))), 1,1:sum(notnan));
    B(ii+1) = y(notnan)*fitCoefs(ii+1,:)'; % fill B vector
end

A = fitCoefs * fitCoefs'; % matrix A for least square fitting
%A(1,1)=npr+1;
lastwarn(''); % reset last warning

Ain=A\eye(size(A)); % instead of inv(A)

b2=  A \ ( y(notnan)*fitCoefs')' ; % alternative articulation
 
[bw,sew_b,msew] = lscov(fitCoefs',y(notnan)',w(notnan)); % weighted LS
 
d = @(ww) w(notnan).*(yfun(ww)-y(notnan));
b3 = lsqnonlin(@(w) sqrt(2*(sqrt(1+abs(d(w)))-1)) ,ones(numel(bParams)+1,1)); %  the new cost function


[~, msgid] = lastwarn % capture warning in case Ain is close to singular
if strcmp(msgid,'MATLAB:nearlySingularMatrix');
    % switch to Moore-Penrose pseudoinverse in case of nearly singular
    Ain =  pinv(A);
end

Beta = zeros(1,numel(bParams)+1);
Beta(1)=B*Ain(:,1);
betas(1)=(Beta(1));   % Beta0 is just the Intensity Factor


for ii=1:numel(bParams)
    Beta(ii+1)=B*Ain(:,ii+1)/Beta(1); % generate beta matrix
    %  Beta(ii+1)=B*Ain(:,ii+1); % generate beta matrix
    betas(ii+1) = Beta(ii+1); % copy for output betas
end
bb=betas; 
alphaMat = (AngleInc * (0:PPR));
rrpCosAlphaMat = ones(1,PPR+1).*cos(alphaMat);
clear beta_contr

for ii=1:numel(bParams) %  add each beta contribution
    %  beta_contr(:,ii)= betas(ii+1)*leg(bParams(ii), rrpCosAlphaMat);
    beta_contr(:,ii)= betas(ii+1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:numel(y));
    beta_contr2(:,ii)= b2(ii+1)./b2(1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:numel(y));
    beta_contr3(:,ii)= b3(ii+1)./b3(1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:numel(y));

    beta_contrw(:,ii)= bw(ii+1)./bw(1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:numel(y));

end

factMat = betas(1).*(ones(1, PPR+1) + sum(beta_contr,2));
factMat2 = b2(1).*(ones(1, PPR+1) + sum(beta_contr2,2));

factMat3 = b3(1).*(ones(1, PPR+1) + sum(beta_contr3,2));

factMatw = bw(1).*(ones(1, PPR+1) + sum(beta_contrw,2));

% note that the /sqrt(r) is omitted because the function is suppose
% to work such as LD2(beta2cart(betas))=betas;
reco_vec(1:PPR+1)=factMat(1:PPR+1,1);% /sqrt(r);
reco_vec2(1:PPR+1)=factMat2(1:PPR+1,1);% /sqrt(r);
reco_vec3(1:PPR+1)=factMat3(1:PPR+1,1);% /sqrt(r);

reco_vecw(1:PPR+1)=factMatw(1:PPR+1,1);% /sqrt(r);

close all

subplot(1,2,1)
errorbar(vec,vars,'o','Color',[0.5 0.5 0.5]); hold on; plot(vec2,'k')
plot( reco_vec2,'b'); 
plot( reco_vec3,'g'); 
plot( reco_vecw,'r'); 
 
legend('vec','ground truth','LS','WLS (soft L1)','WLS');

subplot(1,2,2)
%errorbar(vec,vars,'o','Color',[0.5 0.5 0.5]); hold on; plot(vec2,'k')
plot( vec2-reco_vec2,'b'); hold on
plot( vec2-reco_vec3,'g'); 
plot( vec2-reco_vecw,'r'); 


%%
% for n=1:numel(betas)
%     if n>1
%       %  S{n}=['\beta_' num2str(2*n-2) ' = ' num2str(betas(n)*betas(1))];
%           S{n}=['\beta_' num2str(2*n-2) ' = ' num2str(b2(n))  ];
%        S2{n}=['\alpha_' num2str(2*n-2) ' = ' num2str(bw(n))  ];
%     else
%         %S{n}=['\beta_' num2str(2*n-2) ' = ' num2str( betas(n))];
%           S{n}=['\beta_' num2str(2*n-2) ' = ' num2str( b2(n))];
%         S2{n}=['\alpha_' num2str(2*n-2) ' = ' num2str( bw(n))]
%    
%     end
% end
%  [S ; S2];
%  
 
%% 
%clc



 
%text(1.9*pi,betas(1),S)
%text(0.9*pi,betas(1),S2)

