function [betas varargout]=LDSD(vec, bParams)
warning('off','all')

% Apply Legendre decomposition for a single data vector spanning [0,2*pi] 
% Inputs:
%  vec      - a vector that contain information that spans 0 to 2*pi
%  bParams  - a vector with the (even) beta parameter (legenadre order)
%             to be used in the fitting procedure (excluding 0th-order
%             which is always included), e.g. [2 4] for beta_2 & beta_4.
%
% Outputs:
%  betas      - a vector containing each beta_n such that positive
%               orders are normalized by beta_0 (intensity of vec),
%               for example: betas(1) is beta_0, betas(2) is beta_2/beta_0,
%               betas(3) is beta4/beta_0, ...
%
% reco_vec    - if a second output is asked,  the reconstructed vector 
%               using the betas is given 
%   Ver 1 (2019-11-14)
%   Adi Natan (natan@stanford.edu)
%
%
% % example:
% 
% % generate random vec up to beta 4 + noise
% N=64; % # of bins in angle from 0 to 2*pi
% vec=1000*(rand(1)-0.5)+...
%     (rand(1)-0.5)*cos(linspace(0,2*pi,N)).^2+...
%     (rand(1)-0.5)*cos(linspace(0,2*pi,N)).^4+...
%     0.1*(rand(1)-0.5)*rand(1,N);
% 
% %apply Legendre decomposition up to a higher order and see that indeed
% %higher order dont contribute
% [betas, reco_vec]=LDSD(vec, 2:2:6);
% 
% plot(linspace(0,2*pi,N),vec,'x');hold on; plot(linspace(0,2*pi,N),reco_vec); legend('vec','reco vec');
% for n=1:numel(betas)
%     if n>1
%         S{n}=['\beta_' num2str(2*n-2) ' = ' num2str(betas(n)*betas(1))]
%     else
%         S{n}=['\beta_' num2str(2*n-2) ' = ' num2str( betas(n))]
%     end
% end
% text(1.9*pi,betas(1),S)


%% defaults

if (nargin < 2);  bParams=[2 4]; end

% Check that the beta Parameters are in the scope of the code
if any(mod(bParams,2)) || any(bParams<=0) || any(bParams>42)
    error('Only even positive beta parameters of <=42 orders supported! Beyond that order there are floating point accuracy errors and vpa treatment is needed');
end

PPR=numel(vec)-1;
AngleInc = 2*pi./PPR'; % angle increment per radius
betas = zeros(numel(bParams)+1,1);


npr=PPR;
qp=0:npr;

% for very small radii we reduce the maximal beta value
% allowed to avoid divergences at the origin
if npr/2 <= max(bParams)
    bParams=2:2:npr/2;
end

y = vec;%ira(r,qp+1); % assign row of all pixels in radius rp
% this may have NaN values, one way to treat such a case is to
% interpolate \ extrapolate NaN values using splines or AR methods.
% We will ignore this solution for now and keep the NaNs as is.

B = zeros(1,numel(bParams)+1);
% this is beta0, including if there are NaN values
B(1)=nansum(y)/sum(~isnan(y))*(npr+1);

% one fit coefficient for each B param
fitCoefs = ones(numel(bParams) + 1, sum(~isnan(y)));
for ii=(1:numel(bParams))
    % assign relevant Legendre polynomials to fitCoef
    fitCoefs(ii+1,:) = leg(bParams(ii), cos(AngleInc*qp(~isnan(y))));
    B(ii+1) = y(~isnan(y))*fitCoefs(ii+1,:)'; % fill B vector
end

A = fitCoefs * fitCoefs'; % matrix A for least square fitting
A(1,1)=npr+1;
lastwarn(''); % reset last warning

Ain=A\eye(size(A)); % instead of inv(A)

[~, msgid] = lastwarn; % capture warning in case Ain is close to singular
if strcmp(msgid,'MATLAB:nearlySingularMatrix')
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



%% reconstruct signal
if nargout>1
alphaMat = (AngleInc * (0:PPR));
        rrpCosAlphaMat = ones(1,PPR+1).*cos(alphaMat);
        clear beta_contr
        
        for ii=1:numel(bParams) %  add each beta contribution
            beta_contr(:,ii)= betas(ii+1)*leg(bParams(ii), rrpCosAlphaMat);
        end
        
        factMat = betas(1).*(ones(1, PPR+1) + sum(beta_contr,2));
        % note that the /sqrt(r) is omitted because the function is suppose
        % to work such as LD2(beta2cart(betas))=betas;
        reco_vec(1:PPR+1)=factMat(1:PPR+1,1);% /sqrt(r);
        varargout{1}=reco_vec;
end
    
function p=leg(m,x)
%  This function returns Legendre polynomial P_m(x) where m is the degree
%  of polynomial and X is the variable.
  
        P = legendre(m,x);
        p = squeeze(P(1,:,:));  
end
