function out=ANhankel_transform(sig,vec,pOrder)
% perfom Hanekl transform
% sig - the signal to be transfromed
% vec -  the coordinate vector of sig
% pOrder -  the transform order
%
% The function hankel_matrix is based on https://doi.org/10.1364/JOSAA.21.000053
% The function bessel_zeros is based on: doi.org/10.1016/0021-9991(79)90134-7 

h = hankel_matrix(pOrder, vec(end),numel(vec), 1e-13);
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
% the inverse transform is:
% IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax; 
 
out = HT(sig);


function s_HT = hankel_matrix(order, rmax, Nr, eps_roots)
%HANKEL_MATRIX: Generates data to use for Hankel Transforms
%
%	s_HT = hankel_matrix(order, rmax, Nr, eps_roots)
%
%	s_HT	=	Structure containing data to use for the pQDHT
%	order	=	Transform order
%	rmax	=	Radial extent of transform
%	Nr		=	Number of sample points
%	eps_roots	=	Error in estimation of roots of Bessel function (optional)
%
%	s_HT:
%		order, rmax, Nr =	As above
%		J_roots		=	Roots of the pth order Bessel fn.
%		J_roots_N1	=	(N+1)th root
%		r			=	Radial co-ordinate vector
%		v			=	frequency co-ordinate vector
%		kr			=	Radial wave number co-ordinate vector
%		vmax		=	Limiting frequency
%					=	roots_N1 / (2*pi*rmax)
%		S			=	rmax * 2*pi*vmax (S product)
%		T			=	Transform matrix
%		J			=	Scaling vector
%					=	J_(order+1){roots}
%
%	
if (~exist('eps_roots', 'var')||isempty(eps_roots))
	s_HT.eps_roots = 1e-6;
else
	s_HT.eps_roots = eps_roots;
end
s_HT.order = order;
s_HT.rmax = rmax;
s_HT.Nr = Nr;
%	Calculate N+1 roots:
J_roots = bessel_zeros(1, s_HT.order, s_HT.Nr+1, s_HT.eps_roots);
s_HT.J_roots = J_roots(1:end-1);
s_HT.J_roots_N1 = J_roots(end);
%	Calculate coordinate vectors
s_HT.r = s_HT.J_roots * s_HT.rmax / s_HT.J_roots_N1;
s_HT.v = s_HT.J_roots / (2*pi * s_HT.rmax);
s_HT.kr = 2*pi * s_HT.v;
s_HT.vmax = s_HT.J_roots_N1 / (2*pi * s_HT.rmax);
s_HT.S = s_HT.J_roots_N1;
%	Calculate hankel matrix and vectors
%	I use (p=order) and (p1=order+1)
Jp = besselj(s_HT.order, (s_HT.J_roots) * (s_HT.J_roots.') / s_HT.S);
Jp1 = abs(besselj(s_HT.order+1, s_HT.J_roots));
s_HT.T = 2*Jp./(Jp1 * (Jp1.') * s_HT.S);
s_HT.J = Jp1;
end

function z = bessel_zeros(d, a, n, e)
%BESSEL_ZEROS: Finds the first n zeros of a bessel function
%
%	z = bessel_zeros(d, a, n, e)
%
%	z	=	zeros of the bessel function
%	d	=	Bessel function type:
%			1:	Ja
%			2:	Ya
%			3:	Ja'
%			4:	Ya'
%	a	=	Bessel order (a>=0)
%	n	=	Number of zeros to find
%	e	=	Relative error in root
%
%	This function uses the routine described in:
%		"An Algorithm with ALGOL 60 Program for the Computation of the
%		zeros of the Ordinary Bessel Functions and those of their
%		Derivatives".
%		N. M. Temme
%		Journal of Computational Physics, 32, 270-279 (1979)
z = zeros(n, 1);
function FI = fi(y)
	c1 = 1.570796;
	if (~y)
		FI = 0;
	elseif (y>1e5)
		FI = c1;
	else
		if (y<1)
			p = (3*y)^(1/3);
			pp = p^2;
			p = p*(1 + pp*(pp*(27 - 2*pp) - 210)/1575);
		else
			p = 1/(y + c1);
			pp = p^2;
			p = c1 - p*(1 + pp*(2310 + pp*(3003 + pp*(4818 + pp*(8591 + pp*16328))))/3465);
		end
		pp = (y+p)^2;
		r = (p - atan(p+y))/pp;
		FI = p - (1+pp)*r*(1 + r/(p+y));
	end
end
function Jr = bessr(d, a, x)
	switch (d)
		case (1)
			Jr = besselj(a, x)./besselj(a+1, x);
		case (2)
			Jr = bessely(a, x)./bessely(a+1, x);
		case (3)
			Jr = a./x - besselj(a+1, x)./besselj(a, x);
		otherwise
			Jr = a./x - bessely(a+1, x)./bessely(a, x);
	end
end
aa = a^2;
mu = 4*aa;
mu2 = mu^2;
mu3 = mu^3;
mu4 = mu^4;
if (d<3)
	p = 7*mu - 31;
	p0 = mu - 1;
	if ((1+p)==p)
		p1 = 0;
		q1 = 0;
	else
		p1 = 4*(253*mu2 - 3722*mu+17869)*p0/(15*p);
		q1 = 1.6*(83*mu2 - 982*mu + 3779)/p;
	end
else
	p = 7*mu2 + 82*mu - 9;
	p0 = mu + 3;
	if ((p+1)==1)
		p1 = 0;
		q1 = 0;
	else
		p1 = (4048*mu4 + 131264*mu3 - 221984*mu2 - 417600*mu + 1012176)/(60*p);
		q1 = 1.6*(83*mu3 + 2075*mu2 - 3039*mu + 3537)/p;
	end
end
if (d==1)|(d==4)
	t = .25;
else
	t = .75;
end
tt = 4*t;
if (d<3)
	pp1 = 5/48;
	qq1 = -5/36;
else
	pp1 = -7/48;
	qq1 = 35/288;
end
y = .375*pi;
if (a>=3)
	bb = a^(-2/3);
else
	bb = 1;
end
a1 = 3*a - 8;
% psi = (.5*a + .25)*pi;
for s=1:n
	if ((a==0)&(s==1)&(d==3))
		x = 0;
		j = 0;
	else
		if (s>=a1)
			b = (s + .5*a - t)*pi;
			c = .015625/(b^2);
			x = b - .125*(p0 - p1*c)/(b*(1 - q1*c));
		else
			if (s==1)
				switch (d)
					case (1)
						x = -2.33811;
					case (2)
						x = -1.17371;
					case (3)
						x = -1.01879;
					otherwise
						x = -2.29444;
				end
			else
				x = y*(4*s - tt);
				v = x^(-2);
				x = -x^(2/3) * (1 + v*(pp1 + qq1*v));
			end
			u = x*bb;
			v = fi(2/3 * (-u)^1.5);
			w = 1/cos(v);
			xx = 1 - w^2;
			c = sqrt(u/xx);
			if (d<3)
				x = w*(a + c*(-5/u - c*(6 - 10/xx))/(48*a*u));
			else
				x = w*(a + c*(7/u + c*(18 - 14/xx))/(48*a*u));
			end
		end
		j = 0;
		
while ((j==0)|((j<5)&(abs(w/x)>e)))
			xx = x^2;
			x4 = x^4;
			a2 = aa - xx;
			r0 = bessr(d, a, x);
			j = j+1;
			if (d<3)
				u = r0;
				w = 6*x*(2*a + 1);
				p = (1 - 4*a2)/w;
				q = (4*(xx-mu) - 2 - 12*a)/w;
			else
				u = -xx*r0/a2;
				v = 2*x*a2/(3*(aa+xx));
				w = 64*a2^3;
				q = 2*v*(1 + mu2 + 32*mu*xx + 48*x4)/w;
				p = v*(1 + (40*mu*xx + 48*x4 - mu2)/w);
			end
			w = u*(1 + p*r0)/(1 + q*r0);
			x = x+w;
		end
		z(s) = x;
	end
end
end

