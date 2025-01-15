function [X] = Pro2TraceNorm_hsi(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S, V, D] = svd(Z,'econ');
V = max(diag(V) - tau', 0);
n = sum(V > 0);
X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
% [m, n] = size(Z);
% if 2*m < n
%     AAT = Z*Z';
%     [S, Sigma, D] = svd(AAT,'econ');
%     Sigma = diag(Sigma);
%     V = sqrt(Sigma);
%     tol = max(size(Z)) * eps(max(V));
%     n = sum(V > max(tol, tau));
%     mid = max(V(1:n)-tau, 0) ./ V(1:n) ;
%     X = S(:, 1:n) * diag(mid) * S(:, 1:n)' * Z;
%     return;
% end
% if m > 2*n
%     [X, n, Sigma] = Pro2TraceNorm(Z', tau);
%     X = X';
%     return;
% end
% [S,V,D] = svd(Z);
% Sigma = diag(V);
% n = sum(diag(V) > tau);
% X = S(:, 1:n) * max(V(1:n,1:n)-tau, 0) * D(:, 1:n)';




