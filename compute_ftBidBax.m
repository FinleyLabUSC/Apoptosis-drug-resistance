function [rss, lss, ftBidBax] = compute_ftBidBax(tBid, Bax0, KtBidBax)
% Computes rss, lss, and ftBidBax from the cubic steady-state equations.

    K = KtBidBax;

    % --- Compute cubic coefficients ---
    a = 2*tBid - Bax0 + 2*K;
    b = K*(2*tBid - 2*Bax0 + K);

    % Cubic equation:
    % rss^3 + a*rss^2 + b*rss - r0*(K^2) = 0
    polyCoeff = [1, a, b, -Bax0*(K^2)];

    % Solve cubic
    roots_r = roots(polyCoeff);

    % Keep real, non-negative roots
    real_roots = roots_r(abs(imag(roots_r)) < 1e-9);
    real_roots = real_roots(real_roots >= 0);

    if isempty(real_roots)
        error('No valid real positive root for rss in compute_ftBidBax.');
    end

    % Choose smallest positive real root (biophysical branch)
    rss = min(real_roots);

    % --- Compute lss ---
    ratio = rss / K;
    lss = tBid / (1 + ratio^2);

    % --- Compute c2ss = ftBidBax ---
    ftBidBax = lss * ratio^2;
    fprintf('ftBidBax = %g\n', ftBidBax);

end
