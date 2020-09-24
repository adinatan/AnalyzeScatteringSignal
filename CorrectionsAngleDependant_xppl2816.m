function [TT,Geometry,Polarization,Absorption]=CorrectionsAngleDependant_xppl2816(xy,L0,bcx,bcy,JetThickness,JetWidth,SampleMu)
% Polarisation based on Hura 
% J chem Phys 2000 A high-quality x-ray scattering experiment on liquid
% water at ambient conditions.
% Tim Brandt van Driel, August 2011
%
% Geometry correction for a plane detector based on Boesecke & Diat, J.
% Appl. Cryst 1997, p.867, eqs 6.
% Martin Meedom Nielsen, August 2011
%
% Absorption correction calculated for a tilted liquid sheet.
% Martin Meedom Nielsen, Aug 2011

% Coordinate convention, cartesian lab frame: incoming beam propagates in
% positive z-direction: (0 0 1), y-axis positive upwards. 
%
% Detector frame: coinciding with the labframe, for a detectorplan
% perpendicular to the incoming beam, and shifted positive along the y-axis
% corresponding to the sample-detector distance
%
% xy: a Nx2 matrix with coordinates for all pixel positions in the lab frame
% dx,dy pixel size in x and z direction (in the detector frame) 
% bcx,bcy: pixel coordinates of the direct beam (in the detctor frame)
% L0: smallest distance between detector and sample
%
% Modified from original for xppl2816. 6/2016 Robert Hartsock
%% Construct variables used for the calculations
px=xy(:,1);
py=xy(:,2);
% List of pixel coordinates in lab frame(real distances):
p=[(px(:)-bcx) py(:)-bcy L0*ones(size(px(:)))];
Lp2=(p(:,1).^2+p(:,2).^2+p(:,3).^2); % squared scattering vector
Lp=sqrt(Lp2); % length of scatering vector

%% Polarisation
Pol=0; % Polarisation 0:1, 0=vertical polarisation
sint=realsqrt(p(:,1).^2+p(:,2).^2)./realsqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
sinp=p(:,2)./realsqrt(p(:,1).^2+p(:,2).^2);
cosp=p(:,1)./realsqrt(p(:,1).^2+p(:,2).^2);
Polarization=1./(Pol.*(1-(sinp.*sint).^2)+(1-Pol).*(1-(cosp.*sint).^2));
%% Geometric correction
% Eqs. 6 in Boesecke & Diat
dOmegaNormalized=(p(:,3).^2./Lp2).*(p(:,3)./sqrt(Lp2));

Geometry=1./dOmegaNormalized;
%% Absorption correction through the liquid sheet
% Calculate the transmission through a liquid sheet of thickness JetThickness and
% width JetWidth. The solution is discretized in npoints points and avaraged.
% 
% NB! The thickness of the jet is assumed to be very small compared to the
% distance from the jet to the detector (i.e. JetThickness << abs(L0)), allowing
% the direction vectors to the pixels to be treated as independent of the
% sample point through the jet
npoints=1000;
% Some help parameters:
JetThicknessh=JetThickness/2;
JetWidthh=JetWidth/2;
% Coordinates p are in the Lab frame 
point=zeros(npoints,3);
% create npoints samples along the beam inside the jet. 
point(:,3)=linspace(-JetThicknessh,JetThicknessh,npoints);
%
d=zeros(size(p));
nor=sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
d(:,1)=p(:,1)./nor;
d(:,2)=p(:,2)./nor;
d(:,3)=p(:,3)./nor;
%
% rotate direction vectors to the jet system - a rotation about the y-axis
% in the lab frame:
b=0; % forgive the extraneous lines of code. this is a patch work correction. RWH 6/2016
Rb=[cos(b) 0.0 -sin(b); 0.0 1.0 0.0; sin(b) 0.0 cos(b)];
c=Rb*d';
%Rotate the line of starting points to the jet system in the same way:
pp=Rb*point';
% In the Jet coordinate system, the jet is bounded in the x-z plane by the
% lines: x=-JetWidth/s, x=+JetWidth/2, z=-JetThickness/2, z=+JetThickness/2. Note that negative
% values of c(3,:) corresponds to rays to the detector exiting from the
% 'backside' of the jet
t=zeros(1,length(d));
for n=1:npoints
    l1=(JetWidthh-pp(1,n))./c(1,:);
    idx=l1<0; % find the negative solution for x
    l1(idx)=(-JetWidthh-pp(1,n))./c(1,idx); % they cut the other face of the jet
    l2=(JetThicknessh-pp(3,n))./c(3,:); 
    idx=l2<0; % find the negative solution for z
    l2(idx)=(-JetThicknessh-pp(3,n))./c(3,idx); % they cut the other face of the jet
    % the miniSampleMum distance travelled until cutting at either x or z
    % determines the distance spend in the jet, when adding the distance
    % travelled along the beam in the jet until the point of scattering
    l=min(l1,l2)+(point(n,3)+JetThicknessh);
    t=t+exp(-SampleMu*l);
end
Totaltransmission=t'/npoints; %return the average value
Absorption=1./Totaltransmission;
%% Two theta angle
TT=acos(p(:,3)./Lp)*180/pi;
%% Plot all corrections seperately and combined.
% figure
% Plotfig(xyz(:,1), xyz(:,2), Polarization,777)
% title('Polarisation Correction')
% colorbar
% figure
% Plotfig(xyz(:,1), xyz(:,2), Geometry,778)
% title('Geometric Correction')
% colorbar
% figure
% Plotfig(xyz(:,1), xyz(:,2), Absorption,779)
% title('Absorption Correction')
% colorbar

