function [Q,Geometry,Polarization]=DetectorCorrections(photonenergy,xyz0,pxyz)
% a function that generates the Q, scattering angle, and detector corrections:
%
% 1. Polarization correction based on Hura et al https://doi.org/10.1063/1.1319614 
% 2. Geometry correction for a plane detector based on Boesecke & Diat https://doi.org/10.1107/S0021889897001647  (Eq 6.)
%
%Coordinate convention, cartesian lab frame: incoming beam propagates in
%positive z-direction: (0 0 1), y-axis positive upwards. The detector is at z=0, perpendicular to the incoming beam.

% Inputs:
% photonenergy - a scalar in keV
% xyz0 -  a 1x3 vector for the estiamted x,y,z shifts to center the detector and find its distance to the sample
% pxyz -  an Nx3 array with coordinates for all detector pixel positions in the lab frame (N pixels, xyz postion per pixel)

% Outputs:
% Q - an Nx3 array of radial Q magnitudes per detector pixels 
% Polarization - an Nx3 array of polarization scaling per detector pixels
% Geometry - an Nx3 array of geomerty scaling per detector pixels

% apply shifts to xyz
p = ([(pxyz(:,1)-xyz0(1)) pxyz(:,2)-xyz0(2)  pxyz(:,3)+xyz0(3)]);
Lp2=(p(:,1).^2+p(:,2).^2+p(:,3).^2); % squared scattering vector

% Two theta angle to Q
TwoTheta=acosd(p(:,3)./sqrt(Lp2));
Q = 4*pi*sind(TwoTheta/2).*photonenergy./12.3984193;

% Pol. correction
Pol=0; % 0 for vertical, 1 for horizontal
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(Lp2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));
 
% Geo. correction ( Eqs. 6 in Boesecke & Diat)
Geometry=1./( (p(:,3).^2./Lp2).*(p(:,3)./realsqrt(Lp2))  );
 
