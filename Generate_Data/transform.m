function [X] = transform(X,D)
% Random linear transform to embed d-dimensional data into D-dimensions
% INPUTS:
%        X = D0 x N matrix, N points in D0-dimensions (each col = 1 pt)
%        D = desired embedding dimension
% OUTPUTS:
%        X = D x N matrix, transformed data (each col = 1 pt) 

[D0,~] = size(X); % original dimension of the data

if D0 ~= D
    T = rand(D,D0); % random lin. transform
    T = orth(T);    % orthogonalize so distance-preserving
    X = T*X;        % embed in D-space
end

end