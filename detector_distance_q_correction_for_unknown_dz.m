clear all
close all
% assume a round detector, where the otigin of the axes is at the center of
% the detector. We express the radial position of the detector by:
detector_r_grid=linspace(eps,82.5,830); % in mm similar to ePix10k 2M
EkeV = 9; % photon energy in keV

% q as function of detector distance z:
a=(4*pi / (12.3984193/EkeV));
qrz= @(z) a/sqrt(2).*sqrt(1-z./sqrt(detector_r_grid.^2+z.^2))+eps;

% assume we have a "sample" at some estimated distance z0 from the detector:
z0=60; % in mm
% so we define q to be:
q = qrz(z0);
Name = {'C','I','I','H','H'};
xyz  = [    0.142149937     -0.113392611     -0.010383833   %C
    1.804868889     -1.135903532      0.001393572   %I1
    -1.821501675     -1.141064312     -0.001561729   %I2
    0.191382262      0.584389635      0.898176095   %H1
    -0.052475117      0.636542526     -0.844064941]; %H2

% The simulated radial scattering of this sample, at these parameters is:
sig_at_z0 = SDS({EkeV z0 eps},Name,xyz);

% the theoretical \ simulated scattering curve at z0:
S0_theo=sig_at_z0.SN;

% However, while we thought we measured our sample at z0,
% we were wrong by some unknown amount nonzero distance shift dz:
dz=sign(rand(1)-0.5).*3*rand(1)+5;  

% the "measured" scattering signal will be (added some noise to it):
sig_at_z0dz = SDS({EkeV z0-dz eps},Name,xyz);
S0_measured = sig_at_z0dz.SN+0.005*max(sig_at_z0dz.SN)*randn(numel(sig_at_z0dz.SN),1);

% the correction term from z0 to z0-dz is:
ct= @(z0,dz) a*(sind( acotd((z0-dz)./(z0*tand(2*asind(q./a))  ))/2));

% assume some q-range to do the fitting to find dz:
qrang=q>0.8 & q<4;
errloss=@(x) sum( abs( S0_theo(qrang)' -(x(1)+x(2).*interp1(ct(z0*1e-3,x(3)*1e-3), S0_measured  ,q(qrang),'makima' ) )))./sum(qrang);
% assume some genreal signal fitting parameters in addition to dz:
%           "DC"     scaling    dz
init_vals= [ 0          1       .1  ];
LB =       [-1e-1     0.1      -10  ];
UB =       [1e-1       10       10  ];
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolX',1e-13,'TolFun',1e-13);
[x_fit,fval,exitflag] = fminsearchbnd(errloss,init_vals,LB,UB,options);

% applying the correction term (ct) to the scattering signal based on the fit:
S0_corrected=x_fit(1)+x_fit(2).*interp1(ct(z0*1e-3,x_fit(3)*1e-3),S0_measured,q);

% plot it all:
plot(q(qrang),S0_theo(qrang),'r','LineWidth',2); hold on
plot(q(qrang),S0_measured(qrang),'b','LineWidth',1);
plot(q(qrang),S0_corrected(qrang),':k','LineWidth',2);
AA = char(197); % symobl for Angstrom for plotting
xlabel(['q [' AA '^{-1}]']); ylabel('I [e^2]')
legend('S_0(q) theoretical',['S_0(q) measured with unknown error dz=' num2str(dz)],['S_0(q) corrected with dz fit =' num2str(x_fit(3))])
