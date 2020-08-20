function [d] = twoNN(X)
% Two Nearest Neighbor algorithm as described in:
%    Facco, E., d’Errico, M., Rodriguez, A., & Laio, A. (2017). Estimating 
%    the intrinsic dimension of datasets by a minimal neighborhood
%    information. Scientific reports, 7(1), 12140.
% INPUT:  X = DxN matrix, each col = 1 data pt
% OUTPUT: d = ID estimator

[~,N] = size(X); % sample size

% use this instead of pdist2 if it crashes due to N too large
% D = zeros(N,N); % pairwise euclidean distances
% for i = 1:N
%     for j = 1:N
%         D(i,j) = norm(X(:,i)-X(:,j));
%     end
% end
% D = sort(D,2); % sort in ascending order

D = pdist2(X',X','euclidean','Smallest',3); % pairwise distances

r1 = D(2,:); % distance to closest neighbor
r2 = D(3,:); % distance to 2nd closest neighbor

mu = r2./r1;   % ratio between two shortest distances
mu = sort(mu); % sort in ascending order

F_emp = (1:N)/N; % empirical cdf

I_out = ceil(.9*N);                        % only consider 1st 90% of points
mu(I_out:end) = []; F_emp(I_out:end) = []; % discard outliers

x = log(mu);       % x-coordinate for linear regression
y = -log(1-F_emp); % y-coordinate for linear regression
d = x'\y';         % ID estimate = slope of linear regression passing through origin

% generate same plots as in Figure 1 of the paper
% plot(x,y,'.'); hold on; plot([0,max(x)],[0,d_hat*max(x)]);

end

