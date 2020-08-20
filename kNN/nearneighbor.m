function [d_hat] = nearneighbor(X,k,tol,maxiter)
% Near neighbor algorithm as described in:
%   Pettis, K. W., Bailey, T. A., Jain, A. K., & Dubes, R. C. (1979). An 
%   intrinsic dimensionality estimator from near-neighbor information. IEEE 
%   Transactions on pattern analysis and machine intelligence, (1), 25-37.
% INPUTS:
%        X       = DxN matrix, each col = 1 data point
%        k       = neighborhood size
%        tol     = tolerance for min diff b/t iterates (sugg. 0.01)
%        maxiter = maximum number of iterations allowed (sugg. 4)
% OUPUTS:
%        d_hat   = ID estimate (may not be an integer)

X = X';          % transpose so each row = 1 data pt
[N,~] = size(X); % sample size

% find nearest neighbors
[~, r] = knnsearch(X,X,'K',k+1,'Distance','euclidean'); % Euclidean distance to k nearest neighbors
r(:,1) = [];                                            % don't include yourself as a neighbor

% remove outliers
m_max = sum(r(:,k))/N;                                     % measure of avg distance
s_max = sum((r(:,k)-m_max).^2)/(N-1); s_max = sqrt(s_max); % measure of avg spread of distances
out = (r(:,k) <= m_max + s_max);                           % identify outliers
r = r(out,:); N = sum(out);                                % remove outliers
r_bar = sum(r,1)/N;                                        % r(k) = avg distance to kth neighbor

% iterative method to calculate estimate d_hat
G = @(k,d) k.^(1./d).*gamma(k)./gamma(k+1./d);
x = log(1:k); y = log(r_bar);                                  % compute est using slope of least squares regression
d0 = (k*sum(x.^2) - (sum(x))^2)/(k*sum(x.*y) - sum(x)*sum(y)); % initialize estimate
i = 0;                                                         % iteration counter
eps = tol+1;                                                   % diff b/t estimate iterates
while eps > tol && i <= maxiter
    y = log(r_bar) + G(1:k,d0);                                    % new y-coord for least squares regression
    d1 = (k*sum(x.^2) - (sum(x))^2)/(k*sum(x.*y) - sum(x)*sum(y)); % current estimate
    eps = abs(d1 - d0);                                            % diff b/t last two estimates
    d0 = d1;                                                       % reset previous estimate
    i = i + 1;                                                     % iterate number
end

d_hat = d1; %round(d1); % final ID estimate

end

