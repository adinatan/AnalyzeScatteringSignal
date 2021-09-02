function out = yfun(ww)   
N=64;
bParams=[2:2:6];
PPR=N-1;
AngleInc = 2*pi./PPR'; % angle increment per radius

alphaMat = (AngleInc * (0:PPR));
rrpCosAlphaMat = ones(1,PPR+1).*cos(alphaMat);
paren = @(x, varargin) x(varargin{:});
 

for ii=1:numel(bParams) %  add each beta contribution
      beta_contrw(:,ii)= ww(ii+1)./ww(1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:N);
end

factMatw = ww(1).*(ones(1, PPR+1) + sum(beta_contrw,2));
out(1:PPR+1)=factMatw(1:PPR+1,1);

