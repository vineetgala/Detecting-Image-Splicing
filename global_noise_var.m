function gnv = global_noise_var(I)
    %takes a *.tif (rgb) image as input and returns the calculated global
    %noise variance
    %I = imread('sample.tif');
    %I = im2double(I);
    I = rgb2gray(I);
    T = dctmtx(8);
    dct = @(block_struct) T * block_struct.data * T';
    dct_full = blockproc(I,[8 8],dct);
    %Has size same as input and each 8*8 block is 2D DCT transform of 8*8
    %input image
    
    [height, wid] = size(I);
    
    dct_patches = zeros(8,8,height*wid/64);
    %3d matrix with each 2d plane being a 8*8 2D DCT patch (64 channels)
    
    for y_coord = 1:height/8
        for x_coord = 1:wid/8
            dct_patches(:,:,x_coord + ((wid/8)*(y_coord-1))) = dct_full(y_coord*8-7:y_coord*8, x_coord*8-7:x_coord*8);
        end
    end
    kurt_mtx = kurtosis(dct_patches,1,3);
    var_mtx = var(dct_patches,0,3);
    avg_root_kurt = (sum(kurt_mtx.^0.5,'all')-kurt_mtx(1,1)^0.5)/63;
    avg_var_inv_sq = (sum(1./(var_mtx.^2),'all')-1/(var_mtx(1,1)^2))/63;
    avg_root_kurt_var_inv = (sum((kurt_mtx.^0.5)./var_mtx,'all')-(kurt_mtx(1,1)^0.5)/var_mtx(1,1))/63;
    avg_var_inv = (sum(1./(var_mtx),'all')-1/(var_mtx(1,1)))/63;
    root_kurt_conc = ((avg_root_kurt*avg_var_inv_sq) - (avg_root_kurt_var_inv*avg_var_inv))/(avg_var_inv_sq - avg_var_inv^2);
    gnv = 1/avg_var_inv - avg_root_kurt/(avg_var_inv*root_kurt_conc);
end