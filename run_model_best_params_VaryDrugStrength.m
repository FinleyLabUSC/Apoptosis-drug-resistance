%% Figures 7 and 8 - vary drug strength

%% ===============================
% Base parameter vector
% ===============================
params = [
    1e-4       %  1 k1   
    1e-4       %  2 k2 
    5.8e-4     %  3 k3
    3e-3       %  4 k4
    0.035      %  5 k4r
    3e-3       %  6 k5
    5e-3       %  7 k6
    0.035      %  8 k6r
    5e-4       %  9 k7
    1e-3       % 10 k8
    1e-3       % 11 k9
    7e-3       % 12 k10
    2.21e-3    % 13 k10r
    2e-4       % 14 k11
    2e-4       % 15 k12 
    5e-5       % 16 k13
    1.06e-4    % 17 k14
    1e-3       % 18 k14r
    1.3e-3     % 19 sC8   
    1.3e-3     % 20 sC3   
    1.111e-3   % 21 sXIAP 
    1.111e-3   % 22 sBAR  
    4.168e-4   % 23 sBid  
    1e-3       % 24 sCytCm 
    0.0167     % 25 sSMACm 
    1.3e-3     % 26 sC9 
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
    1e-5       % 39 dCytcm
    1e-5       % 40 dCytc
    1.667e-5   % 41 dSMACm
    1.667e-5   % 42 dSMAC
    1.933e-4   % 43 dSMACIAP
    1.487e-5   % 44 dApop
    6.5e-5     % 45 dC9
    9.667e-5   % 46 dC9a
    2.883e-4   % 47 dC9aXIAP
    0.0385572  % 48 fDISC (derived from function)
    0.0159727  % 49 ftBidBax (derived from function)
    0.0601891  % 50 fApop (derived from function)
    0          % 51 Tocopheryloxybutyrate = 9e-3(on) = 0(off)
    0          % 52 Narciclasine = 9(on) = 0(off)
    1          % 53 Celecoxib  = 9(on); = 1(off)
    1          % 54 beta WITH Toco (for bell shape) = 10(on); = 1(off)
];


params(1:19) = params(1:19)*0.02; 
params(20:26) = params(20:26)*0.2; %rescaling only synthesis reactions
params(27:47) = params(27:47)*0.02;    


%% ===============================
% Initial conditions
% ===============================
y0 = zeros(21,1);
y0(2)  = 216.67; %C8
y0(3)  = 10;     %C8a
y0(4)  = 35;     %C3    
y0(5)  = 10;     %C3a
y0(6)  = 66.67;  %XIAP
y0(8)  = 66.67;  %BAR
y0(10) = 25;     %Bid    
y0(13) = 10;     %Cytochrome C_mitochondria    
y0(14) = 100;    %Cytochrome C   
y0(16) = 100;    %SMAC   
y0(19) = 20;     %C9    
y0(20) = 10;     %C9a

%% ===============================
% Load best-fit PSO parameters
% ===============================
load('gp_best.mat');
load('gf_best.mat');
load('gp_best_All.mat');
load('gf_best_All.mat');

bestSet = gp_best;

paramIdx = [2 3 15 16 20 26 31];
icIdx    = [2 4 19];

params(paramIdx) = bestSet(1:length(paramIdx));
y0(icIdx)        = bestSet(length(paramIdx)+1:end);

%% ===============================
% Time span
tspan = [0 172800]; % seconds
opts  = odeset('RelTol',1e-6,'AbsTol',1e-9);

drugNames = {'Tocopheryloxybutyrate','Narciclasine','Celecoxib'};

drugGrads = {1e-3:1e-3:9e-3, 1:1:9, 1:1:9};

figure
tiledlayout(3,1,'TileSpacing','compact')

for d = 1:3

    nexttile
    hold on   % keep all curves

    strengths = drugGrads{d};
    n = numel(strengths);

    % --- color maps ---
    if d == 1 % Toco (pink)
        cmap = [linspace(1.00,0.60,n)', linspace(0.80,0.00,n)', linspace(0.90,0.60,n)'];
    elseif d == 2 % Narc (blue)
        cmap = [linspace(0.70,0.00,n)', linspace(0.85,0.45,n)', linspace(1.00,0.60,n)'];
    else % Cele (purple)
        cmap = [linspace(0.80,0.40,n)', linspace(0.70,0.20,n)', linspace(1.00,0.50,n)'];
    end

    maxC3a_vals = zeros(n,1);

    for s = 1:n

        % --- clone params for this iteration ---
        params_local = params;

        % --- reset other drugs ---
        params_local(51) = 0;
        params_local(52) = 0;
        params_local(53) = 1; 
        params_local(54) = 1; 

        % --- activate only the current drug ---
        if d == 1
            params_local(51) = strengths(s);
            params_local(54) = 10;
        elseif d == 2
            params_local(52) = strengths(s);
        elseif d == 3
            params_local(53) = strengths(s);
        end

        % --- solve ODE ---
        sol = ode15s(@(t,y) coreFile_DISC(t,y,params_local), tspan, y0, opts);
        t = sol.x;
        Y = sol.y;

        % --- plot C3a ---
        plot(t/60, Y(5,:)*1e4, ...
            'LineWidth', 2, ...
            'Color', cmap(s,:), ...
            'DisplayName', sprintf('%.3g', strengths(s)));

        % --- debug output ---
        maxC3a_vals(s) = max(Y(5,:)*1e4);

    end

 % ===== COMPUTE SUMMARY (AFTER LOOP — FIXED) =====
    lowVal  = min(maxC3a_vals);
    highVal = max(maxC3a_vals);

    absChange  = highVal - lowVal;
    foldChange = highVal / lowVal;

    % ===== ADD TO PLOT =====
    text(0.02, 0.85, ...
        sprintf('\\DeltaC3a = %.2f\nFold = %.2fx', absChange, foldChange), ...
        'Units','normalized', ...
        'FontSize', 11, ...
        'FontWeight','bold');

    % ===== formatting =====
    title([drugNames{d} ' effect on C3a'])
    xlabel('Time (minutes)')
    ylabel('Normalized C3a')

    legend('show','Location','bestoutside')
    set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

    hold off
end

%% Same figure but C9a this time

drugNames = {'Tocopheryloxybutyrate','Narciclasine','Celecoxib'};

drugGrads = {1e-3:1e-3:9e-3, 1:1:9, 1:1:9};

figure
tiledlayout(3,1,'TileSpacing','compact')

for d = 1:3

    nexttile
    hold on   % keep all curves

    strengths = drugGrads{d};
    n = numel(strengths);

    % --- color maps ---
    if d == 1 % Toco (pink)
        cmap = [linspace(1.00,0.60,n)', linspace(0.80,0.00,n)', linspace(0.90,0.60,n)'];
    elseif d == 2 % Narc (blue)
        cmap = [linspace(0.70,0.00,n)', linspace(0.85,0.45,n)', linspace(1.00,0.60,n)'];
    else % Cele (purple)
        cmap = [linspace(0.80,0.40,n)', linspace(0.70,0.20,n)', linspace(1.00,0.50,n)'];
    end

    maxC9a_vals = zeros(n,1);

    for s = 1:n

        % --- clone params for this iteration ---
        params_local = params;

        % --- reset other drugs ---
        params_local(51) = 0;
        params_local(52) = 0;
        params_local(53) = 1;  % Cele OFF by default

        % --- activate only the current drug ---
        if d == 1
            params_local(51) = strengths(s);
        elseif d == 2
            params_local(52) = strengths(s);
        elseif d == 3
            params_local(53) = strengths(s);
        end

        % --- solve ODE ---
        sol = ode15s(@(t,y) coreFile_DISC(t,y,params_local), tspan, y0, opts);
        t = sol.x;
        Y = sol.y;

        % --- plot C3a ---
        plot(t/60, Y(20,:)*1e4, ...
            'LineWidth', 2, ...
            'Color', cmap(s,:), ...
            'DisplayName', sprintf('%.3g', strengths(s)));

        % --- debug output ---
        maxC9a_vals(s) = max(Y(20,:)*1e4);

    end

 % ===== COMPUTE SUMMARY (AFTER LOOP — FIXED) =====
    lowVal  = min(maxC9a_vals);
    highVal = max(maxC9a_vals);

    absChange  = highVal - lowVal;
    foldChange = highVal / lowVal;

    % ===== ADD TO PLOT =====
    text(0.02, 0.85, ...
        sprintf('\\DeltaC9a = %.2f\nFold = %.2fx', absChange, foldChange), ...
        'Units','normalized', ...
        'FontSize', 11, ...
        'FontWeight','bold');

    % ===== formatting =====
    title([drugNames{d} ' effect on C9a'])
    xlabel('Time (minutes)')
    ylabel('Normalized C9a (×10^4 scaled)')

    legend('show','Location','bestoutside')
    set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

    hold off
end

%% ==========================================
% Normalized C3a (relative to highest dose) *
% ==========================================

figure
tiledlayout(3,1,'TileSpacing','compact')

for d = 1:3

    nexttile
    hold on

    strengths = drugGrads{d};
    n = numel(strengths);

    % --- colormaps (same as before) ---
    if d == 1
        cmap = [linspace(1.00,0.60,n)', linspace(0.80,0.00,n)', linspace(0.90,0.60,n)'];
    elseif d == 2
        cmap = [linspace(0.70,0.00,n)', linspace(0.85,0.45,n)', linspace(1.00,0.60,n)'];
    else
        cmap = [linspace(0.80,0.40,n)', linspace(0.70,0.20,n)', linspace(1.00,0.50,n)'];
    end

    C3a_traj = cell(n,1);
    t_traj   = cell(n,1);
    maxVals  = zeros(n,1);

    % --- simulate all strengths first ---
    for s = 1:n

        params_local = params;
        params_local(51) = 0;
        params_local(52) = 0;
        params_local(53) = 1;
        params_local(54) = 1;

        if d == 1
            params_local(51) = strengths(s);
            params_local(54) = 10;
        elseif d == 2
            params_local(52) = strengths(s);
        else
            params_local(53) = strengths(s);
        end

        sol = ode15s(@(t,y) coreFile_DISC(t,y,params_local), tspan, y0, opts);
        t = sol.x;
        Y = sol.y;

        C3a = Y(5,:) * 1e4;

        C3a_traj{s} = C3a;
        t_traj{s}   = t;
        maxVals(s)  = max(C3a);
    end

    % --- normalization factor = highest dose ---
    normFactor = maxVals(end);

    % --- plot normalized ---
    for s = 1:n
        plot(t_traj{s}/60, ...
             C3a_traj{s} / normFactor, ...
             'LineWidth',2, ...
             'Color',cmap(s,:), ...
             'DisplayName',sprintf('%.3g', strengths(s)));
    end

    title([drugNames{d}])
    xlabel('Time (minutes)')
    ylabel({'C3a (normalized', 'to max value)'})
    legend('show','Location','bestoutside')
    set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

    hold off
end

%% ==========================================
% Normalized C9a (relative to highest dose)
% ==========================================

figure
tiledlayout(3,1,'TileSpacing','compact')

for d = 1:3

    nexttile
    hold on

    strengths = drugGrads{d};
    n = numel(strengths);

    % --- colormaps (same as before) ---
    if d == 1
        cmap = [linspace(1.00,0.60,n)', linspace(0.80,0.00,n)', linspace(0.90,0.60,n)'];
    elseif d == 2
        cmap = [linspace(0.70,0.00,n)', linspace(0.85,0.45,n)', linspace(1.00,0.60,n)'];
    else
        cmap = [linspace(0.80,0.40,n)', linspace(0.70,0.20,n)', linspace(1.00,0.50,n)'];
    end

    C9a_traj = cell(n,1);
    t_traj   = cell(n,1);
    maxVals  = zeros(n,1);

    % --- simulate all strengths first ---
    for s = 1:n

        params_local = params;
        params_local(51) = 0;
        params_local(52) = 0;
        params_local(53) = 1;
        params_local(54) = 1;

        if d == 1
            params_local(51) = strengths(s);
            params_local(54) = 10;
        elseif d == 2
            params_local(52) = strengths(s);
        else
            params_local(53) = strengths(s);
        end

        sol = ode15s(@(t,y) coreFile_DISC(t,y,params_local), tspan, y0, opts);
        t = sol.x;
        Y = sol.y;

        C9a = Y(20,:) * 1e4;

        C9a_traj{s} = C9a;
        t_traj{s}   = t;
        maxVals(s)  = max(C9a);
    end

    % --- normalization factor = highest dose ---
    normFactor = maxVals(end);

    % --- plot normalized ---
    for s = 1:n
        plot(t_traj{s}/60, ...
             C9a_traj{s} / normFactor, ...
             'LineWidth',2, ...
             'Color',cmap(s,:), ...
             'DisplayName',sprintf('%.3g', strengths(s)));
    end

    title([drugNames{d}])
    xlabel('Time (minutes)')
    ylabel({'C9a (normalized', 'to max value)'})
    legend('show','Location','bestoutside')
    set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

    hold off
end


%% Figure 10

% Max C3a Range Across Drug Treatments (Normalized to No Drug Case) ---

% Drug strengths
drugGrads = {1e-3:1e-3:9e-3, 1:1:9, 1:1:9};

% Case labels
caseNames = {'Untreated','Toco','Narc','Cele', ...
             'Toco+Narc','Toco+Cele','Narc+Cele', ...
             'All Three'};
nCases = numel(caseNames);

% Colors for plotting
colors = [
    1.0 1.0 1.0;     % No Drugs (white)
    1.0 0.5 0.7;     % Toco (pink)
    0.2 0.5 1.0;     % Narc (blue)
    0.6 0.4 1.0;     % Cele (purple)
    0.8 0.8 0.8;     % Toco+Narc
    0.5 0.5 0.5;     % Toco+Cele
    0.3 0.3 0.3;     % Narc+Cele
    0.1 0.1 0.1      % All Three
];

minVals = zeros(nCases,1);
maxVals = zeros(nCases,1);

% --- Precompute "No Drug" baseline ---
params_base = params;
params_base(51) = 0; % Toco off
params_base(52) = 0; % Narc off
params_base(53) = 1; % Cele baseline

sol_base = ode15s(@(t,y) coreFile_DISC(t,y,params_base), tspan, y0, opts);
Y_base = sol_base.y;
maxC3a_base = max(Y_base(5,:)*1e4);  % store baseline max C3a

figure
hold on

for c = 1:nCases
    maxC3a_all = [];

    for s = 1:9
        params_local = params;

        % --- assign strengths safely ---
        toco = drugGrads{1}(min(s,end));
        narc = drugGrads{2}(min(s,end));
        cele = drugGrads{3}(min(s,end));

        % --- activate drugs based on case ---
        switch c
            case 1 % No Drug
                params_local(51) = 0; params_local(52) = 0; params_local(53) = 1;
            case 2 % Toco
                params_local(51) = toco;
            case 3 % Narc
                params_local(52) = narc;
            case 4 % Cele
                params_local(53) = cele;
            case 5 % Toco + Narc
                params_local(51) = toco; params_local(52) = narc;
            case 6 % Toco + Cele
                params_local(51) = toco; params_local(53) = cele;
            case 7 % Narc + Cele
                params_local(52) = narc; params_local(53) = cele;
            case 8 % All Three
                params_local(51) = toco; params_local(52) = narc; params_local(53) = cele;
        end

        % --- simulate ---
        sol = ode15s(@(t,y) coreFile_DISC(t,y,params_local), tspan, y0, opts);
        Y = sol.y;

        % --- store normalized max C3a ---
        maxC3a_all(end+1) = max(Y(5,:)*1e4) / maxC3a_base;  % normalized
    end

    % --- store range for plotting ---
    minVals(c) = min(maxC3a_all);
    maxVals(c) = max(maxC3a_all);

    % --- vertical range line ---
    plot([c c], [minVals(c) maxVals(c)], 'Color', colors(c,:), 'LineWidth', 5)

    % --- bottom dot ---
    plot(c, minVals(c), 'o', 'MarkerSize',8,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','k')
    % --- top dot ---
    plot(c, maxVals(c), 'o', 'MarkerSize',8,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','k')

    % --- summary print ---
    fprintf('\n=== %s ===\n', caseNames{c});
    fprintf('Min max C3a (norm): %.4f\n', minVals(c));
    fprintf('Max max C3a (norm): %.4f\n', maxVals(c));
    fprintf('Range: %.4f\n', maxVals(c) - minVals(c));
end

% --- formatting ---
xlim([0 nCases+1])
xticks(1:nCases)
xticklabels(caseNames)
xtickangle(30)
ylabel('Max C3a (normalized to untreated)')
yscale("log");
title('Normalized Max C3a Across Drug Treatments')
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off