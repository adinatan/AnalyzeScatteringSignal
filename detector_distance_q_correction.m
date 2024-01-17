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
S0_z0=sig_at_z0.SN;

plot(q,S0_z0,'LineWidth',2); hold on
AA = char(197); % symobl for Angstrom for plotting
xlabel(['q [' AA '^{-1}]']); ylabel('I [e^2]')

% However, if we estimated the distance from the detecytor to be different
% for example by z0-dz, with:
dz=5; % in mm
% the scattering signal will be:
sig_at_z0dz = SDS({EkeV z0-dz eps},Name,xyz);
S0_z0dz=sig_at_z0dz.SN;

% but with the incorrect estimated q:
plot(q,S0_z0dz,'LineWidth',2); 

% the correction term from z0 to z0-dz is:
ct= @(z0,dz) a*(sind( acotd((z0-dz)./(z0*tand(2*asind(q./a))  ))/2)); 

% applying to the incorrectly assume distance scattering signal:
S0_corrected=interp1(ct(z0*1e-3,dz*1e-3),S0_z0dz,q);
plot(q,S0_corrected,':','LineWidth',2); 

legend('S_0(q) at z_0','S_0(q) at z_0-dz','S_0(q) at z_0-dz corrected')