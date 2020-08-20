function [T] = local_pca(X,k,n)
% Perform PCA on local neighborhoods 
% INPUTS: 
%        X = DxN manifold data, each column = 1 data pt
%        k = # of neighbors per neighborhood (incl. center pt)
%        n = # of neighborhoods to consider (default = N)
% OUTPUT: T = average projection error per dimension

[D,N] = size(X); % [embedding dimension, # of samples]

if nargin < 3
    n = N;
end

I = datasample(1:N,n,'Replace',false); % randomly pick samples

I_nn = knnsearch(X',X','K',k); % indices of nearest neighbors (incl yourself)
mu_global = mean(X,2);         % mean of entire dataset
mu_global = repmat(mu_global,1,k);
    
j = 1;            % counter for number of re-samples
t = zeros(n,D+1); % projection error
for i = I % loop over resampled samples
    
    Y = X(:,I_nn(i,:)); % k-nearest neighbors of X(:,i)
    
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
    
    Xd = zeros(D,k); % projection onto 1st d dimensions of basis
    for d = 1:D+1 % test every possible ID (d = 0:D)
        if d == 1
            t(j,d) = sum(vecnorm(Y-mu,2,1).^2)/sum(vecnorm(Y-mu_global,2,1).^2); % projection error
        else
            Xd = Xd + V(:,d-1)*dot(Y-mu,repmat(V(:,d-1),1,k),1);                      % recentered projection
            t(j,d) = sum(vecnorm(Y - Xd-mu,2,1).^2)/sum(vecnorm(Y-mu_global,2,1).^2); % projection error
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