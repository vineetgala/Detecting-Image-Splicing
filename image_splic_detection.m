I = imread('canong3_nikond70_sub_03.tif');
% input image
aacc = local_noise_var(I);
% local noise variance of I
aacc = single(aacc);
% converting to single for K-means
m = mean(mean(aacc));
% mean of local noise variance
[L,C] = imsegkmeans(aacc,2);
% Image segmentation on local noise variance using K-means and K=2
B = labeloverlay(aacc,L);
% B = segmented image base on local noise variance
imshow(B);




