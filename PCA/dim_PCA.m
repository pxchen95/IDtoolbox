function [d] = dim_PCA(X,T)
% Performs global PCA on the data, explaining (1-T)*100% of data's variance
% INPUTS:
%        X = DxN matrix, each col = 1 data point
%        T = threshold for PCA, must be in [0,1]
% OUTPUT:
%        d = intrinsic dimension estimate (\in {0,1,2,...,D})

if T < 0 || T > 1
    disp('ERROR: T must be in [0,1]')
else
    X = X'; % transpose so each row = 1 data pt

    [~,~,latent] = pca(X);       % compute eigenvalues (=var. explained) of cov. matrix 
    latent = latent/sum(latent); % normalize so sums to 1
    [D,~] = size(latent);        % extrinsic dimension

    for j = 0:D                      % loop over possible dimensions
        if sum(latent(1:j)) >= (1-T) % find # e.v.s needed to account for (1-T) of variance
            d = j;
            break
        end
    end
end

end