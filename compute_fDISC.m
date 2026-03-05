function [rss, lss, fDISC] = compute_fDISC(FasL0, FasR0, KDisc)
    % Computes rss from quartic, then lss, then fDISC = c2ss + c3ss.

    % --- Compute coefficients a, b, y ---
    a = 3*FasL0 - FasR0 + 3*KDisc;
    b = 3*KDisc*(2*FasL0 - FasR0 + KDisc);
    y = KDisc^2 * (3*FasL0 - 3*FasR0 + KDisc);

    % --- Solve quartic for rss ---
    % Polynomial: rss^4 + a*rss^3 + b*rss^2 + y*rss - r0*(KDisc^3) = 0
    polyCoeff = [1, a, b, y, -FasR0*(KDisc^3)];

    roots_r = roots(polyCoeff);

    % Keep only real, positive roots (physical requirement)
    real_roots = roots_r(abs(imag(roots_r)) < 1e-9); % tolerance for numerical noise
    real_roots = real_roots(real_roots >= 0);

    if isempty(real_roots)
        error('No non-negative real solution found for rss.');
    end

    % If multiple valid roots exist, pick the smallest positive one (usually the physical solution)
    rss = min(real_roots);

    % --- Compute lss ---
    ratio = rss / KDisc;
    lss = FasL0 / (1 + 3*ratio + 3*ratio^2 + 3*ratio^3);

    % --- Compute fDISC = c2ss + c3ss ---
    c2ss = 3*lss*(ratio^2);
    c3ss = 3*lss*(ratio^3);

    fDISC = c2ss + c3ss;
    fprintf('fDISC = %g\n', fDISC);
end