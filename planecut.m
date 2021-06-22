function  varargout  = planecut(X,Y,Z,X1,Y1,Z1,translation_vec,rotation_vec,cutTolerance)
% extracts an arbitray cut from a plane (X1,Y1,Z1) and a Z(X,Y) manifold
% by first rotating around the origin by rotation_vec=(alpha,beta,gamma) [deg]
% and translating to some point given by translation_vec (xt,yt,zt);
% natan@stanford.edu Ver 1 , 11 April 2020
if nargin<9
    cutTolerance = 5e-2;
end

    
XYZ=[X1(:),Y1(:),Z1(:)];

Rx=R3d(rotation_vec(1),[1 0 0]);
Ry=R3d(rotation_vec(2),[0 1 0]);
Rz=R3d(rotation_vec(3),[0 0 1]);

C= XYZ*Ry*Rx*Rz ; % Rotation transform

xd=reshape(C(:,1),size(X1));
yd=reshape(C(:,2),size(Y1));
zd=reshape(C(:,3),size(Z1));

% translation transform
subX = xd + translation_vec(1);
subY = yd + translation_vec(2);
subZ = zd + translation_vec(3);


cut=0;
while max(sum(cut))<2 % cut should be less than 2 pixels wide
    cutTolerance = cutTolerance+1e-2;
    Zi=interp2(X,Y,Z,subX,subY,'cubic',0);
    cut=abs(subZ-Zi)<cutTolerance;
    if cutTolerance>1e-1 && sum(cut(:))<1 % check if there is a cut
        disp('no cut found')
        break
    end
    
end
varargout{1}= cut;
varargout{2}= Zi;
varargout{3}= subX;
varargout{4}= subY;
varargout{5}= subZ;

end

function R=R3d(x,u)
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%
% R=R3d(x,u)
%
% x: The counter-clockwise rotation about the axis in degrees.
% u: A 3-vector specifying the axis direction. Must be non-zero

R=eye(3);
u=u(:)/norm(u);

for ii=1:3
    v=R(:,ii);  
    %Rodrigues' formula
    R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
end

end
