function fApop = compute_fApop(Apaf0, Cytc_star, lambda, fitParams)
% Compute fApop for given Apaf0, Cytc_star, and lambda.
% If fitParams is not given, the function will call fit_fApop_params(lambda) (costly).
%
% Usage:
%   P = fit_fApop_params(1);        % run once to get fitParams (expensive)
%   f = compute_fApop(100, 100, 1, P);
%
% Or let compute_fApop compute fits automatically:
%   f = compute_fApop(100,100,1);

if nargin < 4 || isempty(fitParams)
    fprintf('No fitParams supplied. Running fit_fApop_params(%g) - this may take a while...\n', lambda);
    fitParams = fit_fApop_params(lambda);
end

% get nondimensional c0
c0 = Cytc_star ./ Apaf0;

% extract fitted params
b1 = fitParams.beta1;
b2 = fitParams.beta2;
x7_1 = fitParams.x7_1;
x7_inf = fitParams.x7_inf;

% piecewise g(c0)
if c0 <= 1
    g = ( exp(b1 * c0) - 1 ) / ( exp(b1) - 1 ) * x7_1;
else
    g = ( x7_1 - x7_inf ) * exp( b2 * (c0 - 1) ) + x7_inf;
end

% a0 = [Apaf]/[Apaf]_0, in our nondimensionalization we used a0 = 1 (Apaf/Apaf0),
% but if Apaf (current) differs from Apaf0 you'd include it. Here assume a0 = 1
a0 = 1;
fApop = a0 * g;

end
