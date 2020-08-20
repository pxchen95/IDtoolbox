function [T] = local_pca_2k(X,k,n)
% Use neighborhoods of size 2k - perform PCA on every other k neighbor, 
% fit to (i.e. calculate error w/) remaining k neighbors
% INPUTS: 
%        X = DxN manifold data, each column = 1 data pt
%        k = # of neighbors per neighborhood (2k <= total sample size)
%        n = # of neighborhoods to consider (default = N)
% OUTPUT: T = average projection error per dimension

[D,N] = size(X); % [embedding dimension, # of samples]

if nargin < 3
    n = N;
end

if 2*k > N
    disp('ERROR: k must be <= (sample size)/2');
else
    I = datasample(1:N,n,'Replace',false); % randomly pick samples

    I_nn = knnsearch(X',X','K',2*k); % indices of 2k nearest neighbors (incl. yourself)
    
    mu_global = mean(X,2); % mean of entire dataset
    mu_global = repmat(mu_global,1,k);

    j = 1;            % counter for number of re-samples (# of center points for PCA)
    t = zeros(n,D+1); % projection error
    for i = I % loop over resampled samples
        Y = X(:,I_nn(i,1:2:(2*k-1))); % k neighbors to be used for PCA
        Z = X(:,I_nn(i,2:2:(2*k)));   % k neighbors to be fitted
        
        % perform PCA  
        if k == 1                                 % covariance matrix of neighborhood
            C = cov(repmat(Y',2,1));
        else
            C = cov(Y');                              
        end
        [V, e] = eig(C);                          % V = cols = eigenvectors of C, D = diagonals = eigenvalues of C in ascending order
        [~, I_lambda] = sort(diag(e), 'descend'); % reorder eigenvalues in descending order
        V = V(:,I_lambda);                        % reorder corresponding eigenvectors
        mu = mean(Y,2);                           % neighborhood mean
        mu = repmat(mu,1,k);
        
        % calculate projection error in fitting neighbors
        Xd = zeros(D,k); % projection onto 1st d dimensions of basis
        for d = 1:D+1 % test every possible ID (d = 0:D)
            if d == 1
                t(j,d) = sum(vecnorm(Z-mu,2,1).^2)/sum(vecnorm(Z-mu_global,2,1).^2); % projection error
            else
                Xd = Xd + V(:,d-1)*dot(Z-mu,repmat(V(:,d-1),1,k),1);                      % recentered projection
                t(j,d) = sum(vecnorm(Z - Xd-mu,2,1).^2)/sum(vecnorm(Z-mu_global,2,1).^2); % projection error
            end     
        end
        j = j+1;
    end

    T = mean(t,1); % average projection error per dimension
    
    % subplot(1,2,1)
    % plot(0:D,T); xlabel('d'); ylabel('Average Projection Error')
    % subplot(1,2,2)
    % plot(1:D,T(2:end)); xlabel('d'); ylabel('Average Projection Error')
    % title('Zoom In (ignoring 1st data point)')
    % xlim([1 D])
end
end