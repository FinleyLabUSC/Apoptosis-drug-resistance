%% Baseline model simulations (no drug treatment)

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


% ===============================
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

% ===============================
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

% ===============================
% time span
tspan = [0 172800]; %seconds timescale consistent with Harrington's params


% Solve ODEs
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t, Y] = ode15s(@(t,y) coreFile_DISC(t,y,params), tspan, y0, opts);

% === Plot all species ===

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
    'Cytc-mito'     % 13
    'Cytc'          % 14
    'SMAC-mito'     % 15
    'SMAC'          % 16
    'SMAC:XIAP'     % 17
    'Apoptosome'    % 18
    'C9'            % 19
    'C9a'           % 20
    'C9a:XIAP'      % 21
};

figure;
speciesToPlot = length(y0);

% Create tiled layout
tiledlayout(6,4,'TileSpacing','compact','Padding','compact');

for idx = 1:length(y0)
    speciesIdx = y0(idx);
    nexttile
    hold on
    plot(t/60,Y(:,idx),'LineWidth', 3)

    % Formatting
    title(speciesNames{idx}, 'Interpreter','none','FontWeight','bold')
    xlabel('Time (minutes)','FontWeight','bold')
    ylabel('Amount (nM)','FontWeight','bold')
    set(gca, ...
    'Box','off', ...
    'LineWidth',1.5, ...        % <-- makes axes thick
    'TickDir','out', ...
    'TickLength',[0.01 0.01], ...
    'FontSize',10, ...
    'Layer','top')
    grid off

end

sgtitle('Baseline model (no treatment)','FontSize',14,'FontWeight','bold')

