function [ira , varargout]=cart2tripolarQuad(Q)
% Cartesian to polar pixel transform - with compact representation
% The code starts with a Cartesian image of square pixel size dx*dy=1, and
% outputs a polar image of polar pixels that has a similar amount of
% information in polar coordinates dr*dtheta~=1.  The signal in each polar
% pixel is determined via its fractional overlap with the four surrounding
% Cartesian pixels. The result is a triangular polar representation,
% because the number of polar pixels per radius per angle increase with
% radius, i.e. the larger the radius the higher the angular resolution.
% The code was originally used in the context of analyzing images with
% symmetry around a quadrant so that functionality was kept.
% The code support NaN valued pixels for masking.
%
% Inputs:
%  Q       - a quadrant of an image, assumed to be centered and corrected for deformations, tilt, etc.
%             the image may contain NaN, and the function will ignore them.
%
% Optional Outputs:
%
% ira      - Intensity per Radius per Angle - a 2d triangular array in
%            polar coordinates of each im quadrant (pi/2 angle spread),
%            where the quadrants are represented in the 3 dimension (as
%            pages). Area outside the triangular array is set to NaN.
% ira_tot  -  Same as ira but of the entire image that is included (up to
%            0 to 2*pi)
% spol     - taking ira_tot and stretching it to a standard 2d square array
%             using linear interpolation.
%
%


Q=double(Q);

% set geometry parameters
include_corners=0;
if include_corners
    % make image square if it is not padding NaN
    ims=size(Q);
    if ims(1)>ims(2)
        Q(:,end+1:ims(1),:)=NaN;
    elseif ims(1)<ims(2)
        Q(end+1:ims(2),:,:)=NaN;
    end
    % pad NaNs to include corners in the transfromation
    Q=padarray(Q, ceil(size(Q)*(sqrt(2)-1)/2),NaN);
end

L = min( size(Q) );

RR=(0:L);
PPR=(floor(0.5*pi*(RR+1))-1); % the # of pixels per radius for cartesian quadrant
AngleInc = (0.5*pi./PPR'); % the angle increment per radius
AngleInc(1)=0; % avoid inf at origin
%% create image quadrants

%set radius of square matrix around center

 
%%
ira = NaN(L-2,PPR(L)); % initialize the  matrix
     a4=Q;
    ira(1,1) = a4(1,1);      % origin pixel remains the same
    
    for r=2:L-2
        npr=PPR(r); % determine # polar pix in radius
        angincr=AngleInc(r);    % the angular increment per radius
        qp=0:npr;
        xp=r*sin(angincr*qp)+1;  % polar x-coordinate
        yp=r*cos(angincr*qp)+1;  % polar y-coordinate
        
        % the signal in a polar pixel at (r,qp) is determined by its
        %  fractional overlap with the four surrounding Cartesian pixels
        %  define scale fractional weight of cartesian pixels in polar pixels
        xc=round(xp);
        yc=round(yp);
        xd=1-abs(xc-xp);
        yd=1-abs(yc-yp);
        
        xpxc=(-1).^(xp(2:npr)<xc(2:npr));
        ypyc=(-1).^(yp(2:npr)<yc(2:npr));
        
        a4r=zeros(4,npr+1);
        c=a4r;
        
        a4r=[    (a4(xc+(yc-1)*L));
            0    (a4(xc(2:npr)+(yc(2:npr)-1+ypyc)*L))      0;
            0    (a4(xc(2:npr)+xpxc+yc(2:npr)*L))          0;
            0    (a4(xc(2:npr)+xpxc+L*(yc(2:npr)+ypyc)))   0];      

        c=[xd.*yd;
            0    xd(2:npr).*(1-yd(2:npr))      0;
            0    (1-xd(2:npr)).*yd(2:npr)      0;
            0    (1-xd(2:npr)).*(1-yd(2:npr))  0];
        
        % only when all 4 cartesian pixels are NaN the polar pixel will be NaN
        % otherwise, only the not NaN values will be considered according
        % to their fractional contribution
        
        c=bsxfun(@rdivide,c.*(~isnan(a4r)),sum(c.*(~isnan(a4r))));
        
        % gather intensity per radius per angle (ira)
        ira(r,1:npr+1,1) = nansum(c.*a4r);
        ira(r,all(isnan(c.*a4r)),1)=NaN;
        %PESR(r,1:npr+1) = r - 1 + ira(r,1:npr+1);
    end
 

%% additional outputs
if max(nargout,1)>1
  
        spol=t2s(ira);
        varargout{1} = spol;
    end
    
end



function squmat=t2s(ira)
% transfer from triangular polar matrix to square polar matrix
PPR=(floor(0.5*pi*((0:size(ira,1))+1))-1); % calc the  # of pixels per radius in quadrant
%qParams=find(squeeze(~all(all(isnan(ira)))));
totPPR = numel(1)*(PPR+1); % total # of polar pixels per radius considering qParams

ira_tot=  ira ; % get the ira_tot

% create a matrix that "squares" the triangular polar matrix
squmat=zeros(size(ira_tot)); % preallocate
% the first pixel is just the intensity at the origin
squmat(1,:)=ones(1,size(ira_tot,2)).*ira_tot(1,1);
for k=2:size(ira_tot,1)
    try % in case data is all zero, then interp1 wont work
        
        
        %first map pixels in given range to new range
        smap=interp1(1:totPPR(k) , single(~isnan(ira_tot(k,1:totPPR(k)))) , linspace(1,totPPR(k),size(ira_tot,2)) ,'nearest');
        % interp only on non NaN pixels
        nr=~isnan(ira_tot(k,1:totPPR(k)));
        squmat(k,:)=interp1(find(nr),ira_tot(k,nr),linspace(1,totPPR(k),size(ira_tot,2)),'nearest');
        squmat(k,~smap)=NaN;
        
    catch
    end
end
end