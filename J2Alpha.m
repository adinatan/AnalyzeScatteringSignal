% let psi be a solution of a TDSE, psi_v(J,R), a func of vibration, rotation, int dist.
% for simplicity we only include J,R (seperating vib from rot)
clear all
%representing psi(R,alpha) we need to project J-states to Legendre:
rgrid=linspace(3,28,2^8)'; % in atomic units
jmax=39; % just high enough # for strong fields

PSI=zeros(jmax+1,numel(rgrid));
PSI(2,:)=exp(-(rgrid-10).^2);
PSI=PSI./sum(PSI(:));

n_alpha=200; % sampling etc
alpha_grid=linspace(0,pi,n_alpha);
cos_alpha_grid=cos(alpha_grid);
% the assumption here is for a parallel cos^2 transition
% if we want to emulate a perpendicular transition then change
% cos_alpha_grid=sin(alpha_grid);

Y=zeros(jmax+1,n_alpha);
Yo=Y;
for n=0:jmax+1-1
    L=legendre(n,cos_alpha_grid);
    Y(n+1,:)=L(1,:);
    Yo(n+1,:)=L(1,:).*sqrt(sin(alpha_grid));
end

rescaling_factor=diag(Yo*Yo.');
Y=Y./sqrt(repmat(rescaling_factor,[1 n_alpha]));
YT=Y.';

%consider for representation:
%alpha_grid=linspace(-pi/2,pi/2,n_alpha);
%Y=circshift(Y,[0 n_alpha/2]);
%YT=circshift(YT,[0 n_alpha/2]);

rho=abs(YT*PSI).^2;%+abs(YT*PSI_e).^2; % for ground and excited state case

rho=rho./sum(rho(:));

% plot to check
subplot(2,1,1); imagesc(rgrid,0:jmax,abs(PSI).^2); title('|\psi(J,R)|^2'); ylabel('J state');
subplot(2,1,2); imagesc(rgrid,alpha_grid,rho); title('|\psi(\alpha,R)|^2');xlabel('R [au]'); ylabel('\alpha [rad]')
 