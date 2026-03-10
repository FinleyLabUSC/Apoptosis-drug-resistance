function sol = run_model_with_scaled_params() %Supplemental Figure 2
% Build params (use your existing vector values)
params = [
    1e-4      %  1 k1
    1e-4      %  2 k2 
    5.8e-4     %  3 k3 4
    3e-3      %  4 k4 
    0.035      %  5 k4r
    3e-3       %  6 k5
    5e-3       %  7 k6
    0.035      %  8 k6r
    5e-4       %  9 k7
    1e-3      % 10 k8
    1e-3      % 11 k9
    7e-3       % 12 k10
    2.21e-3    % 13 k10r
    2e-4       % 14 k11
    2e-4       % 15 k12
    5e-5       % 16 k13
    1.06e-4    % 17 k14
    1e-3      % 18 k14r
    1.3e-3     % 19 sC8 216.67  *-- adjusted
    1.3e-3         % 20 sC3 35  *-- adjusted
    1.111e-3      % 21 sXIAP 66.67 *-- adjusted
    1.111e-3   % 22 sBAR  66.67
    4.168e-4   % 23 sBid  25
    1e-3      % 24 sCytCm 100
    0.0167     % 25 sSMACm 0.0167 100
    1.3e-3     % 26 sC9 1.3e-3 20
    8.807e-3   % 27 dDISC
    6.5e-5     % 28 dC8
    9.667e-5   % 29 dC8a
    6.5e-5     % 30 dC3
    9.667e-5   % 31 dC3a
    1.933e-4   % 32 dIAP
    2.883e-4   % 33 dC3aIAP
    1.667e-5   % 34 dBAR
    1.933e-4   % 35 dC8aBAR
    1.667e-5   % 36 dBid
    1.667e-5   % 37 dtBid
    0.0264     % 38 dtBid_Bax
    1e-5      % 39 dCytcm
    1e-5      % 40 dCytc
    1.667e-5   % 41 dSMACm
    1.667e-5   % 42 dSMAC
    1.933e-4   % 43 dSMACIAP
    1.487e-5   % 44 dApop
    6.5e-5     % 45 dC9
    9.667e-5   % 46 dC9a
    2.883e-4   % 47 dC9aXIAP
    2         % 48 FasL0  (also used as tBid in ftBidBax)
    10         % 49 FasR0
    1.032      % 50 KDisc
    83.33      % 51 Bax0
    100        % 52 KtBidBax
    100        % 53 Apaf0
    1          % 54 lambda
    0          % 55 placeholder for fDISC (filled below)
    0          % 56 placeholder for ftBidBax
    0          % 57 placeholder for fApop (baseline)    
    0;%9e-3         %Tocopheryloxybutyrate -must turn on with = 1e-3(on) =0(off)
    0;%9         %Narciclasine = 5(on) =0(off)
    1;%9           %Celecoxib  = 5(on) =0(off)
    1;%10          %beta WITH Toco (for bell shape) = 10(on) =1(off)
];


% --- Precompute fDISC (quartic solve) ---
FasL0 = params(48);
FasR0 = params(49);
KDisc = params(50);
[rss_fDISC, lss_fDISC, fDISC] = compute_fDISC(FasL0, FasR0, KDisc);

params(55) = fDISC;

% --- Precompute ftBidBax (cubic solve) ---
tBid = params(48);       % we use same slot for tBid (adjust if you want separate param)
Bax0  = params(51);
KtBidBax = params(52);
[rss_tBid, lss_tBid, ftBidBax] = compute_ftBidBax(tBid, Bax0, KtBidBax);

params(56) = ftBidBax;

% --- Load or fit Apaf params (P) needed by compute_fApop ---
lambda = params(54);
fitFile = sprintf('fApop_fit_lambda_%g.mat', lambda);
if exist(fitFile, 'file')
    S = load(fitFile,'fitParams');
    P = S.fitParams;
    fprintf('Loaded fApop fit file: %s\n', fitFile);
else
    fprintf('Fit file %s not found — running fit_fApop_params(%g). This may take some time...\n', fitFile, lambda);
    P = fit_fApop_params(lambda);  % expensive: integrates many ODEs
    %saved inside fit_fApop_params already
end

% compute a baseline fApop using initial Cytc* = (use initial Cytc from model)
Cytc_star_guess = 100; 
Apaf0 = params(53);
fApop0 = compute_fApop(Apaf0, Cytc_star_guess, lambda, P);
params(57) = fApop0;


% Initial conditions y0 (21 species) 
y0 = zeros(21,1);
y0(2,1) = 216.67; %C8 
y0(3,1) = 0;
y0(4,1) = 35; %C3
y0(5,1) = 0;
y0(6,1) = 66.67; %XIAP
y0(8,1) = 66.67; %BAR
y0(10,1) = 25; %Bid 
y0(13,1) = 10; %CytochromeC-Mitochondria *assumption to avoid negatives
y0(14,1) = 100; %CytochromeC
y0(16,1) = 100; %SMAC
y0(19,1) = 20; %Caspase9
y0(20,1) = 0;



% time span
tspan = [0 172800]; %seconds timescale consistent with Harrington's params

% Solve ODE and pass params and P to ODE function
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
sol = ode15s(@(t,y) coreFile_DISC(t,y,params, P), tspan, y0, opts);


% === Plot all species ===
t = sol.x;
Y = sol.y;

t_end = sol.x(end);
t = linspace(0, t_end, 2000);
Y = deval(sol, t);

disp(sol.x(end)/86400)   % days


speciesNames = {
    'DISC'          % 1
    'C8'            % 2
    'C8a'           % 3
    'C3'            % 4
    'C3a'           % 5
    'XIAP'          % 6
    'C3a:XIAP'      % 7
    'BAR'           % 8
    'C8a:BAR'       % 9
    'Bid'           % 10
    'tBid'          % 11
    'tBid:Bax'      % 12
    'Cytc_mito'     % 13
    'Cytc'          % 14
    'SMAC_mito'     % 15
    'SMAC'          % 16
    'SMAC:XIAP'     % 17
    'Apoptosome'    % 18
    'C9'            % 19
    'C9a'           % 20
    'C9a:XIAP'      % 21
};

figure;
set(gcf,'Position',[100 100 1400 900])
tiledlayout(5,5,'Padding','compact','TileSpacing','compact');

% Highlighted species indices
highlightSpecies = [3 4 5 20];
numHighlight = length(highlightSpecies);

% Base baby blue color
baseColor = [0.54, 0.81, 0.94];

% Shades for species
babyBlueShades = linspace(0.7, 1.0, numHighlight)'; 
colors = babyBlueShades .* baseColor; 
colors(colors>1) = 1;

% Create tiled layout
tiledlayout(numHighlight,1,'TileSpacing','compact','Padding','compact');

for idx = 1:numHighlight
    speciesIdx = highlightSpecies(idx);
    nexttile
    hold on

    % Condition 1: reduced parameters
    params1 = params;             
    params1(1:54) = params1(1:54)*0.02; 
    sol1 = ode15s(@(t,y) coreFile_DISC(t,y,params1, P), tspan, y0, opts);
    t1 = sol1.x; Y1 = sol1.y;
    plot(t1/60, Y1(speciesIdx,:), 'LineWidth', 3, 'Color', colors(idx,:), 'LineStyle','--') % dashed

    % Condition 2: baseline parameters
    params2 = params; % baseline
    sol2 = ode15s(@(t,y) coreFile_DISC(t,y,params2, P), tspan, y0, opts);
    t2 = sol2.x; Y2 = sol2.y;
    plot(t2/60, Y2(speciesIdx,:), 'LineWidth', 3, 'Color', colors(idx,:), 'LineStyle','-') % solid

    % Formatting
    title(speciesNames{speciesIdx}, 'Interpreter','none','FontWeight','bold')
    xlabel('Time (minutes)','FontWeight','bold')
    ylabel('Amount (nM)','FontWeight','bold')
    set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off
    if idx == 1
        legend({'0.02x','baseline'},'Location','best')
    else
        legend off
    end

    hold off
end

sgtitle('Caspase Species Dynamics under Two Parameter Conditions','FontSize',14,'FontWeight','bold')

end

