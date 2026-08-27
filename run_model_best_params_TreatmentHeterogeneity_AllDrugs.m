%% Figures 9 and 10

rng(0); % to get reproducible random number sequence

% ===============================
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

%% Simulate results for Figure 9

% Time span
tspan = [0 172800]; %seconds timescale consistent with Harrington's params

% Solve ODE and pass params and P to ODE function
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
sol = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0, opts);

% Proteins
y0_names = {'C8','C3','XIAP','BAR','Bid','CytoC','SMAC','C9'};
y0_idx   = [2,4,6,8,10,14,16,19];
y0_base  = [y0(2), y0(4), y0(6), y0(8), y0(10), y0(14), y0(16), y0(19)];
nProteins = numel(y0_idx);
nSamples = 500;

% Initialize storage
C3a_max_noDrug = nan(nSamples, nProteins);
C3a_max_combo  = nan(nSamples, nProteins);
allProteins_init = nan(nSamples, nProteins);

% Drug indices
toco_idx = 51;
narc_idx = 52;
cele_idx = 53;
tocoBeta_idx = 54;

% Doses
drug_doses = [9e-3, 1, 1, 10];   % Toco, Narc, Cele, Toco(beta)
noDrug_doses = [0, 0, 1, 1];     % Toco, Narc, Cele, Toco(beta)

% Simulate for each protein
% -----------------------------
for p = 1:nProteins
    i = y0_idx(p);
    base = y0_base(p);

    % Log-uniform samples ±2 orders of magnitude
    y0_samples = base * 10.^(-2 + 4*rand(1,nSamples));
    allProteins_init(:,p) = y0_samples;

    for s = 1:nSamples
        y0_mod = y0;
        y0_mod(i) = y0_samples(s);

        % --- No drug ---
        params([toco_idx narc_idx cele_idx tocoBeta_idx]) = noDrug_doses;
        sol = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0_mod, opts);
        Y_noDrug = sol.y; %Y(Y<0) = 0;
        C3a_max_noDrug(s,p) = max(Y_noDrug(5,:));

        % --- Combination therapy ---
        params([toco_idx narc_idx cele_idx tocoBeta_idx]) = drug_doses;
        sol = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0_mod, opts);
        Y_withDrugCombo = sol.y; %Y(Y<0) = 0;
        C3a_max_combo(s,p) = max(Y_withDrugCombo(5,:));
    end


end

% Reset drugs
params([toco_idx narc_idx cele_idx tocoBeta_idx]) = noDrug_doses;

% Normalize to no-drug baseline
% -----------------------------
C3a_max_noDrug_norm = C3a_max_noDrug./ mean(C3a_max_noDrug,1);   % fold-change = 1 for baseline
C3a_max_combo_norm  = C3a_max_combo./ mean(C3a_max_noDrug,1);

%% Figure 9: Plot paired boxplots

% Prepare boxplot data
% -----------------------------
allData = [];
group   = [];

for p = 1:nProteins
    allData = [allData; C3a_max_noDrug_norm(:,p); C3a_max_combo_norm(:,p)];
    group   = [group; repmat((p-1)*2+1, nSamples,1); repmat((p-1)*2+2, nSamples,1)];
end


% -----------------------------
figure('Color','w','Position',[200 200 1400 600])
h = boxplot(allData, group, 'Colors','k','Symbol','');

% Colors for boxes: alternating light/dark grey
lightGrey = [0.8 0.8 0.8];
darkGrey  = [0.3 0.3 0.3];
colorsBox = repmat([lightGrey; darkGrey], nProteins,1);

boxObjs = findobj(gca,'Tag','Box');
xCenters = arrayfun(@(b) mean(get(b,'XData')), boxObjs);
[~, sortIdx] = sort(xCenters, 'ascend');
boxObjs = boxObjs(sortIdx);

% Apply colors
for j = 1:length(boxObjs)
    patchX = get(boxObjs(j),'XData');
    patchY = get(boxObjs(j),'YData');
    patch(patchX, patchY, colorsBox(j,:), 'FaceAlpha',1, 'EdgeColor','k');
end

% X-axis labels
xticks_pos = mean(reshape(1:nProteins*2,2,nProteins));
set(gca,'XTick', xticks_pos, 'XTickLabel', y0_names)

ylabel('Maximum C3a fold-change (compared to untreated)')
title('Effect of Heterogenous Initial Protein Amounts on C3a Response to Tocopheryloxybutyrate')
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

% Add baseline reference line at fold-change = 1
% -----------------------------
hold on
yline(1, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);

% Add jittered points only (no connecting lines)
% -----------------------------
for p = 1:nProteins
    x1 = (p-1)*2+1; x2 = (p-1)*2+2;
    jitterX1 = x1 + 0.1*(rand(nSamples,1)-0.5);
    jitterX2 = x2 + 0.1*(rand(nSamples,1)-0.5);
    scatter(jitterX1, C3a_max_noDrug_norm(:,p), 20, 'k', 'filled', 'MarkerFaceAlpha',0.5)
    scatter(jitterX2, C3a_max_combo_norm(:,p), 20, 'k', 'filled', 'MarkerFaceAlpha',0.5)
end


% Set y-axis limits 
% -----------------------------
ymax = (max(max([C3a_max_noDrug, C3a_max_combo])))+1;
ymin = (min(min([C3a_max_noDrug, C3a_max_combo])));
ylim([0 ymax]);
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off
set(gca, 'YScale', 'log')
hold off;

% KS tests
% -----------------------------
fprintf('Kolmogorov-Smirnov tests (No drug vs Combo therapy):\n');
for p = 1:nProteins
    [h_ks, p_ks] = kstest2(C3a_max_noDrug(:,p), C3a_max_combo(:,p));
    fprintf('%s: p = %.3g, h = %d\n', y0_names{p}, p_ks, h_ks);
end

%% Simulate results for Figure 10

% Heterogeneous population
%==============================
nCells = 2000; %nSamples;

% Species indices
idx_C8  = 2;   % C8
idx_C3  = 4;   % C3
idx_C3a = 5;   % C3a
idx_IAP = 6;   % XIAP
idx_BAR = 8;   % BAR
idx_C9  = 20;  % C9

% Drug parameter indices
toco_idx = 51;
narc_idx = 52;
cele_idx = 53;
tocoBeta_idx = 54;

% Identify non-zero initial conditions
nonzero_idx = find(y0 > 0);
y0_base = y0(nonzero_idx);

% Storage
C8_init     = nan(nCells,1);
C3_init     = nan(nCells,1);
C3a_max     = nan(nCells,1);
C3a_max_rel = nan(nCells,1);
IAP_init    = nan(nCells,1);
BAR_init    = nan(nCells,1);
C9_init     = nan(nCells,1);
all_init    = nan(nCells,length(y0));

simsToPlot = 50;
darkGrey  = [0.3 0.3 0.3];
Grey = [0.75 0.75 0.75];

% Run simulations with NO DRUG
% -----------------------------
figure;

for c = 1:nCells
    % Create heterogeneous initial state
    y0_mod = y0;
    scaleFactors = 10.^(-2 + 4*rand(numel(nonzero_idx),1));
    y0_mod(nonzero_idx) = y0_base .* scaleFactors;

    % Store initial values
    all_init(c,:) = y0_mod;
    C8_init(c)  = y0_mod(idx_C8);
    C3_init(c)  = y0_mod(idx_C3);
    IAP_init(c) = y0_mod(idx_IAP);
    BAR_init(c) = y0_mod(idx_BAR);
    C9_init(c)  = y0_mod(idx_C9);

    sol = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0_mod, opts);
    Y_base = sol.y; Y_base(Y_base<0) = 0;

    % Max C3a
    C3a_max_base(c) = max(Y_base(idx_C3a,:));


    %plot time courses for the first 100 cells
    subplot(2,1,1)
    if c < simsToPlot
        plot(sol.x/60,Y_base(5,:),'LineWidth',0.5,'Color',darkGrey)
        hold on;
    end

    % Format axes, labels, title
    set(gca, ...
        'Box','off', ...
        'LineWidth',2, ...        
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'FontSize',14, ...
        'Layer','top')
        grid off
    xlim([0 3000])
    ylim([0 1000])
    xlabel('Time (minutes)','FontSize',14)
    ylabel('C3a (nM)','FontSize',14)
    title('Untreated')
    
    % Axes formatting
    grid on
    set(gca, ...
        'Box','off', ...
        'LineWidth',2, ...        
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'FontSize',14, ...
        'Layer','top')
        grid off
end
hold off


% Run simulations drug ON
% -----------------------------

params([toco_idx narc_idx cele_idx tocoBeta_idx]) = drug_doses;

for c = 1:nCells
    % specify initial values (uses the same values as without drug)
    y0_mod(idx_C8) = C8_init(c);
    y0_mod(idx_C3) = C3_init(c);
    y0_mod(idx_IAP) = IAP_init(c);
    y0_mod(idx_BAR) = BAR_init(c);
    y0_mod(idx_C9)  = C9_init(c);
    % y0_mod(idx_C3a) = C3a_init(c);

    sol = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0_mod, opts);
    Y_withDrugCombo = sol.y; Y_withDrugCombo(Y_withDrugCombo<0) = 0;

    % Max C3a
    C3a_max_withDrug(c) = max(Y_withDrugCombo(idx_C3a,:));

    %plot time courses for the first 100 cells
    subplot(2,1,2)
    if c < simsToPlot
        plot(sol.x/60,Y_withDrugCombo(5,:),'LineWidth',0.5,'Color',Grey)
        hold on;
    end

    % Format axes, labels, title
    set(gca, ...
        'Box','off', ...
        'LineWidth',2, ...        
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'FontSize',14, ...
        'Layer','top')
        grid off
    xlim([0 3000])
    ylim([0 1000])
    xlabel('Time (minutes)','FontSize',14)
    ylabel('C3a (nM)','FontSize',14)
    title('Treated with all three drugs')

    % Axes formatting
    grid on
    set(gca, ...
        'Box','off', ...
        'LineWidth',2, ...        
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'FontSize',14, ...
        'Layer','top')
        grid off
end
hold off;

%% Figure S6

% Make boxplot
C3a_max_hetero = [C3a_max_base; C3a_max_withDrug];

figure;
h = boxplot(C3a_max_hetero', [1;2], 'Colors','k','Symbol','');


% Colors for boxes: alternating light/dark grey
lightGrey = [0.8 0.8 0.8];
darkGrey  = [0.3 0.3 0.3];
colorsBox = repmat([lightGrey; darkGrey], nProteins,1);

boxObjs = findobj(gca,'Tag','Box');
xCenters = arrayfun(@(b) mean(get(b,'XData')), boxObjs);
[~, sortIdx] = sort(xCenters, 'ascend');
boxObjs = boxObjs(sortIdx);

% Apply colors
for j = 1:length(boxObjs)
    patchX = get(boxObjs(j),'XData');
    patchY = get(boxObjs(j),'YData');
    patch(patchX, patchY, colorsBox(j,:), 'FaceAlpha',1, 'EdgeColor','k','FaceAlpha',0.5);
end

% X-axis labels
set(gca, 'XTickLabel', {'Untreated','Treated with all three drugs'})

ylabel('Maximum C3a (nM)')
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off

hold on;

% Plot individual values
jitterX1 = 1 + 0.22*(rand(nCells,1)-0.5);
jitterX2 = 2 + 0.22*(rand(nCells,1)-0.5);
scatter(jitterX1', C3a_max_hetero(1,:), 20, 'k', 'filled', 'MarkerFaceAlpha',0.5)
scatter(jitterX2', C3a_max_hetero(2,:), 20, 'k', 'filled', 'MarkerFaceAlpha',0.5)


medLines = findobj(gca,'Tag','Median');
set(medLines, 'LineWidth', 3, 'Color', 'k');
uistack(medLines,'top');


% KS test
% -----------------------------
fprintf('Kolmogorov-Smirnov tests (No drug vs Combo therapy):\n');
[h_ks, p_ks] = kstest2(C3a_max_hetero(1,:), C3a_max_hetero(2,:))


%% Figure 10, panel A

% Scatter plot (initial C3 vs Max C3a)
figure;
scatter(C3_init, C3a_max_withDrug, 20, 'k', 'filled', 'MarkerFaceAlpha',0.45)
set(gca,'XScale','log')

ymin = min(C3a_max_withDrug);
ymax = max(C3a_max_withDrug);

ylim([ymin*0.95 ymax*1.05])
xlim([min(C3_init)*0.9 max(C3_init)*1.1])

xlabel('Initial C3 (nM)')
ylabel('Max C3a (nM)')

set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',16, ...
    'Layer','top')
    grid off

hold on
yline(200, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'Label','Max C3a threshold');


%% Prepare colormap for scatter plot

% Custom green gradient
nColors = 256;

% Define three anchor colors
lightYG   = [0.9 1 0.6];   % yellowish green
forestG   = [0 0.4 0];     % forest green
darkGreen = [0 0.1 0];     % near black-green

% Split colormap into two gradients
n1 = round(nColors/2);
n2 = nColors - n1;

% Gradient 1: yellow-green → forest green
map1 = [linspace(lightYG(1), forestG(1), n1)', ...
        linspace(lightYG(2), forestG(2), n1)', ...
        linspace(lightYG(3), forestG(3), n1)'];

% Gradient 2: forest green → dark green
map2 = [linspace(forestG(1), darkGreen(1), n2)', ...
        linspace(forestG(2), darkGreen(2), n2)', ...
        linspace(forestG(3), darkGreen(3), n2)'];

% Combine
greenMap = [map1; map2];
colormap(greenMap)


%% Figure 10, panel B

% Scatter plot for initial conditions vs C3a max with treatment
% C3 vs XIAP and BAR vs XIAP

allProteins_init_hetero = [C8_init C3_init BAR_init IAP_init C9_init];
y0_names = {'C8','C3','BAR','XIAP','C9'};

min_init = min(allProteins_init_hetero)
max_init = max(allProteins_init_hetero)

% make the size of the dot correlated to max C3a
dotSize = C3a_max_withDrug/5; 

figure;
subplot(2,1,1)
colormap(greenMap)
scatter(allProteins_init_hetero(:,2),allProteins_init_hetero(:,4), dotSize, C3a_max_withDrug, 'filled', 'MarkerFaceAlpha',0.6)
set(gca,'XScale','log','YScale','log')
xlim([min_init(2)*.5 max_init(2)]*1.5)
ylim([min_init(4)*.5 max_init(4)]*1.5)

colorbar
caxis([min(C3a_max_withDrug(:)) max(C3a_max_withDrug(:))])
cb = colorbar;
cb.Label.String = 'Max C3a (nM)';
cb.FontSize = 12;

% Format axes, labels, title
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...       
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',16, ...
    'Layer','top')
    grid off
set(gca,'XScale','log','YScale','log')
xlabel('Initial C3 (nM)','FontSize',14)
ylabel('Initial XIAP (nM)','FontSize',14)

% Axes formatting
grid on
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',16, ...
    'Layer','top')
    grid off


subplot(2,1,2)
colormap(greenMap)
scatter(allProteins_init_hetero(:,3),allProteins_init_hetero(:,4), dotSize, C3a_max_withDrug, 'filled', 'MarkerFaceAlpha',0.6)
set(gca,'XScale','log','YScale','log')
xlim([min_init(3)*.5 max_init(3)]*1.5)
ylim([min_init(4)*.5 max_init(4)]*1.5)

colorbar
cb = colorbar;
cb.Label.String = 'Max C3a (nM)';
cb.FontSize = 12;

% Format axes, labels, title
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',16, ...
    'Layer','top')
    grid off
set(gca,'XScale','log','YScale','log')
xlabel('Initial BAR (nM)','FontSize',14)
ylabel('Initial XIAP (nM)','FontSize',14)

% Axes formatting
grid on
set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',16, ...
    'Layer','top')
    grid off

    
%% Additional data - pairwise scatter plots for initial conditions

% Scatter plot for initial conditions vs C3a max with treatment - for all pairs of proteins
figure; 
counter = 0;
for i = 1:5
    for j = 1:5
        counter = counter + 1;
        if i>j

            subplot(5,5,counter)
            colormap(greenMap)
    
            scatter(allProteins_init_hetero(:,i),allProteins_init_hetero(:,j), C3a_max_withDrug/10, C3a_max_withDrug, 'filled', 'MarkerFaceAlpha',0.6)
            
            set(gca,'XScale','log','YScale','log')
            xlim([min_init(i)*.5 max_init(i)]*1.5)
            ylim([min_init(j)*.5 max_init(j)]*1.5)
            xlabel(y0_names(i),'FontSize',8)
            ylabel(y0_names(j),'FontSize',8)
            set(gca, ...
                'Box','on', ...
                'LineWidth',1.5, ...        
                'TickDir','out', ...
                'TickLength',[0.015 0.015], ...
                'FontSize',10, ...
                'Layer','top')
                grid off
        end

    end
end


%% Figure 10, panel C

% Select cells with  C3a above threshold
above_threshold_idx = [];
below_threshold_idx = [];
threshold1 = 200;

for i = 1:nCells
    if C3a_max_withDrug(i) > threshold1 
        above_threshold_idx = [above_threshold_idx; i];
    else
        below_threshold_idx = [below_threshold_idx; i];
    end
end


% ---------------------%BAR
proteinOfInterest = 3; %BAR

data1 = allProteins_init_hetero(below_threshold_idx,proteinOfInterest);
data2 = allProteins_init_hetero(above_threshold_idx,proteinOfInterest);
combinedData = [data1; data2];

median_BAR_below_threshold = median(data1);
median_BAR_above_threshold = median(data2);

grouping = [repmat({'Max C3a BELOW threshold'}, length(data1), 1); ...
            repmat({'Max C3a ABOVE threshold'}, length(data2), 1)];

figure;
h = boxplot(combinedData, grouping, 'OutlierSize',0.1, 'Colors','k', 'Symbol','');

boxObjs = findobj(gca,'Tag','Box');
set(boxObjs, 'LineWidth', 2, 'Color', 'k');

medLines = findobj(gca,'Tag','Median');
set(medLines, 'LineWidth', 3, 'Color', 'k');
uistack(medLines,'top');
ylabel(append('Initial ', y0_names(proteinOfInterest), ' (nM)'), 'FontSize', 14)


% Colors for boxes: alternating light/dark grey
lightGrey = [0.8 0.8 0.8];
darkGrey  = [0.3 0.3 0.3];
colorsBox = repmat([lightGrey; darkGrey], nProteins,1);

boxObjs = findobj(gca,'Tag','Box');
xCenters = arrayfun(@(b) mean(get(b,'XData')), boxObjs);
[~, sortIdx] = sort(xCenters, 'ascend');
boxObjs = boxObjs(sortIdx);

% Apply colors
for j = 1:length(boxObjs)
    patchX = get(boxObjs(j),'XData');
    patchY = get(boxObjs(j),'YData');
    patch(patchX, patchY, colorsBox(j,:), 'FaceAlpha',1, 'EdgeColor','k','FaceAlpha',0.5);
end


set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off
set(gca,'YScale','log')

hold on

% Add jittered points only (no connecting lines)
% -----------------------------
for p = 1:2
    x1 = (p-1)*2+1; 
    x2 = (p-1)*2+2;
    jitterX1 = x1 + 0.1*(rand(length(data1),1)-0.5);
    jitterX2 = x2 + 0.1*(rand(length(data2),1)-0.5);
    scatter(jitterX1, allProteins_init_hetero(below_threshold_idx,proteinOfInterest), 20, lightGrey,'filled','MarkerFaceAlpha',0.5)
    scatter(jitterX2, allProteins_init_hetero(above_threshold_idx,proteinOfInterest), 20, darkGrey,'filled','MarkerFaceAlpha',0.5)
end


% ---------------------%XIAP
proteinOfInterest = 4; %XIAP

data1 = allProteins_init_hetero(below_threshold_idx,proteinOfInterest);
data2 = allProteins_init_hetero(above_threshold_idx,proteinOfInterest);
combinedData = [data1; data2];

median_XIAP_below_threshold = median(data1)
median_XIAP_above_threshold = median(data2)

grouping = [repmat({'Max C3a BELOW threshold'}, length(data1), 1); ...
            repmat({'Max C3a ABOVE threshold'}, length(data2), 1)];

figure;
h = boxplot(combinedData, grouping, 'OutlierSize',0.1, 'Colors','k', 'Symbol','');

boxObjs = findobj(gca,'Tag','Box');
set(boxObjs, 'LineWidth', 2, 'Color', 'k');

medLines = findobj(gca,'Tag','Median');
set(medLines, 'LineWidth', 3, 'Color', 'k');
uistack(medLines,'top');
ylabel(append('Initial ', y0_names(proteinOfInterest), ' (nM)'), 'FontSize', 14)

% Colors for boxes: alternating light/dark grey
lightGrey = [0.8 0.8 0.8];
darkGrey  = [0.3 0.3 0.3];
colorsBox = repmat([lightGrey; darkGrey], nProteins,1);

boxObjs = findobj(gca,'Tag','Box');
xCenters = arrayfun(@(b) mean(get(b,'XData')), boxObjs);
[~, sortIdx] = sort(xCenters, 'ascend');
boxObjs = boxObjs(sortIdx);

% Apply colors
for j = 1:length(boxObjs)
    patchX = get(boxObjs(j),'XData');
    patchY = get(boxObjs(j),'YData');
    patch(patchX, patchY, colorsBox(j,:), 'FaceAlpha',1, 'EdgeColor','k','FaceAlpha',0.5);
end


set(gca, ...
    'Box','off', ...
    'LineWidth',2, ...        
    'TickDir','out', ...
    'TickLength',[0.02 0.02], ...
    'FontSize',14, ...
    'Layer','top')
    grid off
set(gca,'YScale','log')

hold on

% Add jittered points only (no connecting lines)
% -----------------------------
for p = 1:2
    x1 = (p-1)*2+1; 
    x2 = (p-1)*2+2;
    jitterX1 = x1 + 0.1*(rand(length(data1),1)-0.5);
    jitterX2 = x2 + 0.1*(rand(length(data2),1)-0.5);
    scatter(jitterX1, allProteins_init_hetero(below_threshold_idx,proteinOfInterest), 20, lightGrey,'filled','MarkerFaceAlpha',0.5)
    scatter(jitterX2, allProteins_init_hetero(above_threshold_idx,proteinOfInterest), 20, darkGrey,'filled','MarkerFaceAlpha',0.5)
end