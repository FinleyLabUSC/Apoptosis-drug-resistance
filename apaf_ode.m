function dYdt = apaf_ode(t, Y, lambda)
% Nondimensional Apaf oligomerization ODEs (n = 7).
% Y = [a; c; x1; x2; x3; x4; x5; x6; x7]
% lambda is the ratio k2/k1 used by the paper.

a = Y(1);
c = Y(2);
x = Y(3:9); % x1..x7

dYdt = zeros(9,1);

% da/dtau = -a*c
dYdt(1) = - a * c;

% dc/dtau = -a*c
dYdt(2) = - a * c;

% dx1/dtau = a*c - lambda * x1 * (2*x1 + x2 + x3 + x4 + x5 + x6)
dYdt(3) = a * c - lambda * x(1) * ( 2*x(1) + x(2) + x(3) + x(4) + x(5) + x(6) );

% for i = 2..7:
% dx_i/dtau = lambda * ( sum_{j=1}^{floor(i/2)} x_j * x_{i-j} - x_i * sum_{j=1}^{7-i} (1+delta_{j,i-j}) x_j )
% Implementation carefully follows combinatorics. delta_{ij} acts only in second sum when j == i-j.

for i = 2:7
    % index in x: i -> x(i)
    % compute first sum: sum_{j=1}^{floor(i/2)} x_j * x_{i-j}
    idx = i; % i from 2..7
    m = floor(idx/2);
    s1 = 0;
    for j = 1:m
        s1 = s1 + x(j) * x(idx - j);
    end

    % compute second sum: sum_{j=1}^{7-idx} (1 + delta_{j, idx-j}) * x_j
    s2 = 0;
    maxj = 7 - idx;
    for j = 1:maxj
        if j == (idx - j)
            mult = 1 + 1; % delta = 1
        else
            mult = 1 + 0; % delta = 0
        end
        s2 = s2 + mult * x(j);
    end

    dYdt(2 + idx) = lambda * ( s1 - x(idx) * s2 );
end

end

