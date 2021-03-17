function [pymask] =mask2pymask(mask)
%%
pymask=mask;
pymask=reshape(pymask,[388,185,32]);
pymask=reshape(pymask,[388,185*32]);
pymask=logical(permute(pymask,[2 1]));