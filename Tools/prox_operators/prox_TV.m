function sol = prox_TV(b, lambda, param)
if nargin<3, param=struct; end

if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end

% Initializations
[r, s] = gradient_op(b*0);
pold = r; qold = s;
told = 1; prev_obj = 0;

% Main iterations
if param.verbose > 1
    fprintf('  Proximal TV operator:\n');
end
for iter = 1:param.max_iter
    
    % Current solution
    sol = b - lambda*div_op(r, s);
    
    % Objective function value
    obj = .5*norm(b(:)-sol(:), 2)^2 + lambda * TV_norm(sol);
    rel_obj = abs(obj-prev_obj)/obj;
    prev_obj = obj;
    
    % Stopping criterion
    if param.verbose>1
        fprintf('   Iter %i, obj = %e, rel_obj = %e\n', ...
            iter, obj, rel_obj);
    end
    if rel_obj < param.rel_obj
        crit_TV = 'TOL_EPS'; break;
    end
    
    % Udpate divergence vectors and project
    [dx, dy] = gradient_op(sol);
    r = r - 1/(8*lambda) * dx; s = s - 1/(8*lambda) * dy;
    weights = max(1, sqrt(abs(r).^2+abs(s).^2));
    p = r./weights; q = s./weights;
    
    % FISTA update
    t = (1+sqrt(4*told^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    s = q + (told-1)/t * (q - qold); qold = q;
    told = t;
    
end

% Log after the minimization
if ~exist('crit_TV', 'var'), crit_TV = 'MAX_IT'; end
if param.verbose >= 1
    fprintf(['  Prox_TV: obj = %e, rel_obj = %e,' ...
        ' %s, iter = %i\n'], obj, rel_obj, crit_TV, iter);
end

end