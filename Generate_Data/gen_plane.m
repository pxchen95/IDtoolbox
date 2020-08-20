function [X] = gen_plane(d,D,N,pad_zero)
% Generates N points uniformly sampled from a d-dimensional hyperplane 
% embedded in D-dimensional space
% INPUTS:  d = intrinsic dimension
%          D = embedding dimension
%          N = # of sample points
%          pad_zero = true  - embed by padding with zeros, 
%                     false - embed via random linear transform
% OUTPUTS: X = D-by-N matrix, each column = 1 data point

if pad_zero == true || d == D % embed in higher dim by padding with zeros
    X = zeros(D,N); X(1:d,:) = rand(d,N); % randomly sample from unit D-hypercube
    
else % embed in higher-dim via random linear transform
    T = rand(D,d); % lin. transform to randomly generate a plane
    T = orth(T);   % orthogonalize so distance-preserving
    X = rand(d,N); % randomly sample from unit D-hypercube
    X = T*X;       % embed in D-space
end

end