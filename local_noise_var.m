function aacc = local_noise_var(I)
%     I = imread("sample_s3.png");
    I = rgb2gray(I);
    I = double(I);
    %I(200:300,200:300) = I(200:300,200:300) + 10*randn(101,101);
    [width, height] = size(I);
    width = width -7;
    height = height - 7;
    dct_pixelwise = zeros(width, height, 64);
    T = dctmtx(8);
    for i = 1:width
        for j = 1:height
            dct_pixelwise(i,j,:) = reshape(T*I(i:(i+7),j:(j+7))*T',[64,1]);
        end
    end

    orders = zeros(width,height,64,4);
    orders(:,:,:,1) = dct_pixelwise;
    orders(:,:,:,2) = dct_pixelwise.^2;
    orders(:,:,:,3) = dct_pixelwise.^3;
    orders(:,:,:,4) = dct_pixelwise.^4;

    presum = zeros(width,height,64,4);
    for i = 1:width
        for j = 1:height
            if i > 1
                presum(i,j,:,:) = presum(i,j,:,:) + presum(i-1,j,:,:);
            end
            if j > 1
                presum(i,j,:,:) = presum(i,j,:,:) + presum(i,j-1,:,:);
            end
            if i >1
                if j > 1
                    presum(i,j,:,:) = presum(i,j,:,:) - presum(i-1,j-1,:,:);
                end
            end
            presum(i,j,:,:) = presum(i,j,:,:) + orders(i,j,:,:);
        end
    end
    window_size = 32;
    pic = zeros(width-window_size+1,height-window_size+1);
    for i = 1:(width-window_size+1)
        for j = 1:(height-window_size+1)
            mu = presum(i+31,j+31,:,:);
            if i>1
                mu = mu- presum(i-1,j+31,:,:);
            end
            if j>1
                mu = mu - presum(i+31,j-1,:,:);
            end
            if i>1
                if j>1
                    mu = mu + presum(i-1,j-1,:,:);

                end
            end
            mu = reshape(mu,[64,4]);
            mu = mu./(window_size*window_size);
            mu = mu(2:64,:);
            sig2 = mu(:,2) - mu(:,1).^2;
            kur = (mu(:,4).^2 -4*mu(:,3).*mu(:,1)+6*mu(:,2).*(mu(:,1).^2) -3*(mu(:,1).^4))./(mu(:,2).^2 - 2*mu(:,2).*(mu(:,1).^2) + (mu(:,1).^4)) - 3;

            sqrt_kur = sqrt(max(kur,0));
            inv_sig2 = 1./sig2;

            c_kur = (mean(sqrt_kur,"all")*mean(inv_sig2.^2,"all") - mean(sqrt_kur.*inv_sig2,"all")*mean(inv_sig2,"all"))/(mean(inv_sig2.^2,"all") - mean(inv_sig2,"all").^2);
            noise_var = (1-mean(sqrt_kur,"all")/c_kur)/(mean(inv_sig2,"all"));
            pic(i,j) = sqrt(max(noise_var,0));
        end
    end
    aacc = pic;
end

