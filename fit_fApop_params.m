function fitParams = fit_fApop_params(lambda)
% Fit beta1, beta2, and extract x7_ss(1) and x7_ss(inf) for given lambda.
% Returns struct with fields: beta1, beta2, x7_1, x7_inf, cgrid, x7_ss
%
% WARNING: this integrates many ODEs until steady state; it can take some time.
% Typical usage:
%   P = fit_fApop_params(1);

if nargin < 1
    lambda = 1;%given in Harrington model
end

% c0 grid (nondimensional: c0 = Cytc* / Apaf0). Use log-spaced grid covering <<1 to >>1
cgrid = logspace(-3, 3, 60);  % 0.001 ... 1000

x7_ss = zeros(size(cgrid));

% integration options
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

for k = 1:length(cgrid)
    c0 = cgrid(k);
    % initial conditions: a(0)=1, c(0)=c0, x1..x7 = 0
    Y0 = zeros(9,1);
    Y0(1) = 1;    % a(0)
    Y0(2) = c0;   % c(0)
    % integrate until steady state in nondimensional tau
    tspan = [0 1e4];   % long time; ODE is stiff in range -> ode15s can be used externally
    try
        [~, Y] = ode15s(@(t,y) apaf_ode(t,y,lambda), tspan, Y0, opts);
    catch
        % fallback to ode45 if ode15s missing or fails
        [~, Y] = ode45(@(t,y) apaf_ode(t,y,lambda), tspan, Y0, opts);
    end
    Yend = Y(end,:)';
    x7_ss(k) = Yend(9); % x7
end

% interpolate x7 at c0 = 1
x7_1 = interp1(cgrid, x7_ss, 1, 'pchip');

% estimate x7_inf as the value at the largest cgrid (should converge)
x7_inf = x7_ss(end);

% Fit beta1 for c0 <= 1 using g1(c0) = ((exp(beta1*c0)-1)/(exp(beta1)-1)) * x7_1
idx_low = find(cgrid <= 1);
c_low = cgrid(idx_low);
x_low = x7_ss(idx_low);

% objective for beta1
obj1 = @(b) sum( ( x_low - ( (exp(b*c_low)-1)./(exp(b)-1) ) * x7_1 ).^2 );

% initial guess
b0 = 1.0;
b1 = fminsearch(obj1, b0, optimset('TolX',1e-8,'MaxFunEvals',2000));

% Fit beta2 for c0 > 1 using g2(c0) = (x7_1 - x7_inf)*exp(beta2*(c0-1)) + x7_inf
idx_high = find(cgrid > 1);
c_high = cgrid(idx_high);
x_high = x7_ss(idx_high);

obj2 = @(b) sum( ( x_high - ( (x7_1 - x7_inf).*exp(b*(c_high-1)) + x7_inf ) ).^2 );
b2 = fminsearch(obj2, -0.01, optimset('TolX',1e-8,'MaxFunEvals',2000)); % decay constant typically negative

% assemble results
fitParams.beta1 = b1;
fitParams.beta2 = b2;
fitParams.x7_1 = x7_1;
fitParams.x7_inf = x7_inf;
fitParams.cgrid = cgrid;
fitParams.x7_ss = x7_ss;
fitParams.lambda = lambda;

% optional save
save(sprintf('fApop_fit_lambda_%g.mat',lambda),'fitParams');

end
