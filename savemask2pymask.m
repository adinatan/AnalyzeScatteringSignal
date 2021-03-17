function[pymask]=savemask2pymask(ScanNumber,mask)

%% save pymask, 185 lines  size(a)= 5920         388
%% reshape mask to 32x185x388, ones are masked so pymask=~mask;
%a=textread('X:\1603_IH_Anisotropy\Test\pixel_mask\183-end.data');
%size(a)
% pymask=mask;% no longer oposite in python
% %Images=rdCSPADdata('X:\1603_IH_Anisotropy\HDF5\xpp00316-r0255.h5',1,1:1);
% %Images=permute(Images,[3 2 1]);
% %size(Images)
% pymask=reshape(pymask,[388,185,32]);
% pymask=reshape(pymask,[388,185*32]);
% pymask=permute(pymask,[2 1]);
%pymask=reshape(pymask,[185*32,388]);


% back transform
% size(a)=        5920         388
% a=permute(a,[2 1]);
% size(a)=        388        5920
% a=reshape(a,[388,185,32]);
% size(a)=        388   185    32
% mask=logical(a(:));


%%
%size(pymask);
%sum(pymask(:));
%size(a)
%sum(a(:))
%%
%a=textread('X:\1603_IH_Anisotropy\Test\pixel_mask\183-end.data');
%size(a)

[pymask] =mask2pymask(mask);
%% 20e+00
%dlmwrite(['X:\1603_IH_Anisotropy\Test\pixel_mask\' num2str(ScanNumber) '-end.data'],pymask,'delimiter',' ','precision','%0.18e') % 18 float
dlmwrite([num2str(ScanNumber) '-end.data'],pymask,'delimiter',' ','precision','%0.18e') % 18 float