%% Figure 6 - PSO Fits

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

param_errors = [all_gf_best; all_gp_best]';
params_sortedByError = sortrows(param_errors,1);

pos = params_sortedByError(1:60, 2:11)'; %top 60 parameter sets that do not have a param that hits the bounds 

std_mult = 1; % number of standard deviations


%--- specify the experimental Caspase data (relative to last time point available)
Toco_Caspase3 = [0	0.769
1440	3.850
2880	1]; 

Toco_Caspase9 = [0	0.111
1440	0.678
2880	1]; 

Narc_proCaspase3 = [0	1.968776
960	    1.1249
1440	1];

Narc_Caspase8 = [0	0.3615
960     0.7040
1440	0.7230
2880	1];

Cele_Caspase3 = [0	0.3409200
720     0.3517476
1440	0.8987100
2880	1];

%Specify Experimental x and y for all points 

% --- Toco ---
x_Toco_C3   = Toco_Caspase3(:,1);   % convert to minutes if needed
y_Toco_C3   = Toco_Caspase3(:,2);

x_Toco_C9   = Toco_Caspase9(:,1);
y_Toco_C9   = Toco_Caspase9(:,2);

% --- Narc ---
x_Narc_C3   = Narc_proCaspase3(:,1);
y_Narc_C3   = Narc_proCaspase3(:,2);

x_Narc_C8   = Narc_Caspase8(:,1);
y_Narc_C8   = Narc_Caspase8(:,2);

% --- Cele ---
x_Cele_C3   = Cele_Caspase3(:,1);
y_Cele_C3   = Cele_Caspase3(:,2);%% Set details
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
Time =[0 43200 57600 86400 172800];

tstep = 1;
tsim = [0:tstep:172800];

%Tocopheryloxybutyrate
figure

% --- Activate drug ---
params(51) = 9e-3;
params(54) = 10;

numRuns = size(pos,2);

Toco_pred_C3_full = zeros(numRuns, length(tsim));
Toco_pred_C9_full = zeros(numRuns, length(tsim));

for i = 1:numRuns
    
% Assign parameters 
params(2)  = pos(1,i); %k2
params(3)  = pos(2,i); %k3
params(15) = pos(3,i); %k12
params(16) = pos(4,i); %k13
params(20) = pos(5,i); %sC3
params(26) = pos(6,i );%sC9
params(31) = pos(7,i); %dC3a
y0(2)  = pos(8,i); %iC8
y0(4)  = pos(9,i); %iC3
y0(19) = pos(10,i); %iC3

    % Simulate
    [tsim, pred] = ode15s(@(t,y) coreFile_DISC(t,y,params), tsim, y0, options);

    % Normalize
    C3 = pred(:,5)/pred(end,5);
    C9 = pred(:,20)/pred(end,20);


    Toco_pred_C3_full(i,:) = C3';
    Toco_pred_C9_full(i,:) = C9';
end

mean_TocoC3_full = mean(Toco_pred_C3_full,1);
std_TocoC3_full  = std(Toco_pred_C3_full,0,1);

mean_TocoC9_full = mean(Toco_pred_C9_full,1);
std_TocoC9_full  = std(Toco_pred_C9_full,0,1);

% Normalize data
Toco_normData = Toco_Caspase3(:,2)/Toco_Caspase3(end,2);
Toco_normData2 = Toco_Caspase9(:,2)/Toco_Caspase9(end,2);

% --- Plot C3a ---
subplot(1,2,1); hold on

x = tsim / 60; % convert to minutes

fill([x; flipud(x)], ...
     [mean_TocoC3_full'+std_mult*std_TocoC3_full'; flipud(mean_TocoC3_full'-std_mult*std_TocoC3_full')],...
      [1 0.6 0.8],'EdgeColor','none','FaceAlpha',0.4); 

plot(x,mean_TocoC3_full,'-','Color',[1 0 0.6],'LineWidth',2);
plot(x_Toco_C3, y_Toco_C3, 's','MarkerEdgeColor','k','LineWidth', 4,'MarkerSize',10);
ylim([0, max(mean_TocoC3_full' + std_mult*std_TocoC3_full')*1.1])
title('Tocopheryloxybutyrate'); ylabel('Normalized C3a'); xlabel('Time (min)');


subplot(1,2,2); hold on % C9a

x = tsim / 60;

fill([x; flipud(x)], ...
     [mean_TocoC9_full'+std_mult*std_TocoC9_full'; flipud(mean_TocoC9_full'-std_mult*std_TocoC9_full')],...
      [1 0.6 0.8],'EdgeColor','none','FaceAlpha',0.4); 


plot(x,mean_TocoC9_full,'-','Color',[1 0 0.6],'LineWidth',2);
plot(x_Toco_C9, y_Toco_C9, 's','MarkerEdgeColor','k','LineWidth', 4,'MarkerSize',10);
ylim([0, max(mean_TocoC9_full'+std_mult*std_TocoC9_full')*1.1])
title('Tocopheryloxybutyrate'); ylabel('Normalized C9a'); xlabel('Time (min)');

%START added block for additional cell lines 
% --- NEW CELL LINE DATA (C3a only) ---

% Base timepoints (minutes)
t_new = [0 1440 2880];

% Small offsets to avoid overlap (in minutes)
offset = 40;  % tweak if needed for spacing

t_PC3        = t_new - offset;   % shift left
t_LNCaP_TS   = t_new;            % center
t_LNCaP_TOB  = t_new + offset;   % shift right

% Data
PC3_alphaTOB    = [0.205, 0.857, 1];
LNCaP_alphaTS   = [0.233, 0.395, 1];
LNCaP_alphaTOB  = [0.270, 0.378, 1];

% Plot
plot(t_PC3, PC3_alphaTOB, 's', ...
    'MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor','k', 'MarkerSize',10);

plot(t_LNCaP_TS, LNCaP_alphaTS, '^', ...
    'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize',10);

plot(t_LNCaP_TOB, LNCaP_alphaTOB, '^', ...
    'MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor','k', 'MarkerSize',10);
%END added block for additional cell lines


subplot(1,2,2); hold on % C9a

x = tsim / 60;

fill([x; flipud(x)], ...
     [mean_TocoC9_full'+std_mult*std_TocoC9_full'; flipud(mean_TocoC9_full'-std_mult*std_TocoC9_full')],...
      [1 0.6 0.8],'EdgeColor','none','FaceAlpha',0.4); 


plot(x,mean_TocoC9_full,'-','Color',[1 0 0.6],'LineWidth',2);
plot(x_Toco_C9, y_Toco_C9, 's','MarkerEdgeColor','k','LineWidth', 4,'MarkerSize',10);
ylim([0, max(mean_TocoC9_full'+std_mult*std_TocoC9_full')*1.1])
title('Tocopheryloxybutyrate'); ylabel('C9a'); xlabel('Time (min)');

%START
% --- NEW CELL LINE DATA (C9a only) ---

% Base timepoints (minutes)
t_new = [0 1440 2880];

% Small offsets to avoid overlap (in minutes)
offset = 40;  

t_PC3        = t_new - offset;   % shift left
t_LNCaP_TS   = t_new;            % center
t_LNCaP_TOB  = t_new + offset;   % shift right

% Data
C9a_PC3_alphaTOB    = [0.185, 1.019, 1];
C9a_LNCaP_alphaTS   = [0.1, 0.6, 1];
C9a_LNCaP_alphaTOB  = [0.154, 0.415, 1];

% Plot
plot(t_PC3, C9a_PC3_alphaTOB, 's', ...
    'MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor','k', 'MarkerSize',10);

plot(t_LNCaP_TS, C9a_LNCaP_alphaTS, '^', ...
    'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize',10);

plot(t_LNCaP_TOB, C9a_LNCaP_alphaTOB, '^', ...
    'MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor','k', 'MarkerSize',10);
%END

legend_handles = [
    plot(nan,nan,'s','MarkerEdgeColor','k','LineWidth', 4), ...
    plot(nan,nan,'s','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','k'), ...
    plot(nan,nan,'^','MarkerFaceColor','k','MarkerEdgeColor','k'), ...
    plot(nan,nan,'^','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','k')
];

legend_labels = {
    'PC3-\alphaTS (Experimental Data)', ...
    'PC3-\alphaTOB', ...
    'LNCaP-\alphaTS', ...
    'LNCaP-\alphaTOB'
};

legend(legend_handles, legend_labels, ...
    'Position',[0.35 0.02 0.3 0.05], ... % bottom center
    'Orientation','horizontal', ...
    'Box','off');

% Turn off drug
params(51) = 0;
params(54) = 1;

%% Narciclasine
params(52)=9; % activate

numRuns = size(pos,2);
numTime1 = size(Narc_proCaspase3,1);
numTime2 = size(Narc_Caspase8,1);


Narc_pred_C3_full = zeros(numRuns,length(tsim));
Narc_pred_C8_full = zeros(numRuns,length(tsim));

for i=1:numRuns

% Assign parameters 
params(2)  = pos(1,i); %k2
params(3)  = pos(2,i); %k3
params(15) = pos(3,i); %k12
params(16) = pos(4,i); %k13
params(20) = pos(5,i); %sC3
params(26) = pos(6,i );%sC9
params(31) = pos(7,i); %dC3a
y0(2)  = pos(8,i); %iC8
y0(4)  = pos(9,i); %iC3
y0(19) = pos(10,i); %iC3

    [tsim,pred]=ode15s(@(t,y) coreFile_DISC(t,y,params), tsim, y0, options);

    C3 = pred(:,4)/pred(86400,4);
    C8 = pred(:,3)/pred(end,3);

    Narc_pred_C3_full(i,:) = C3';
    Narc_pred_C8_full(i,:) = C8';
end


mean_NarcC3_full = mean(Narc_pred_C3_full,1);
std_NarcC3_full  = std(Narc_pred_C3_full,0,1);

mean_NarcC8_full = mean(Narc_pred_C8_full,1);
std_NarcC8_full  = std(Narc_pred_C8_full,0,1);


% Normalize data
Narc_normData = Narc_proCaspase3(:,2)/Narc_proCaspase3(end,2);
Narc_normData2 = Narc_Caspase8(:,2)/Narc_Caspase8(end,2);


params(52)=0; %turn off drug


%% Celecoxib
params(53)=9;

numRuns=size(pos,2);
numTime=size(Cele_Caspase3,1);

Cele_pred_full=zeros(numRuns,length(tsim));

for i=1:numRuns

% Assign parameters 
params(2)  = pos(1,i); %k2
params(3)  = pos(2,i); %k3
params(15) = pos(3,i); %k12
params(16) = pos(4,i); %k13
params(20) = pos(5,i); %sC3
params(26) = pos(6,i );%sC9
params(31) = pos(7,i); %dC3a
y0(2)  = pos(8,i); %iC8
y0(4)  = pos(9,i); %iC3
y0(19) = pos(10,i); %iC3

    [tsim,pred]=ode15s(@(t,y) coreFile_DISC(t,y,params), tsim, y0, options);

    C3a = pred(:,5)/pred(end,5);

    Cele_pred_full(i,:) = C3a';
end

mean_Cele_full = mean(Cele_pred_full,1);
std_Cele_full  = std(Cele_pred_full,0,1);

Cele_normData = Cele_Caspase3(:,2)/Cele_Caspase3(end,2);

params(53)=1; %turn off drug


%% --- PSO predictions with experimental points ---
% Figure 6

% Tocopherylbutyrate 

figure;
subplot(3,2,1); hold on
% Full dynamics: mean ± std
fill([tsim(:)/60; flipud(tsim(:)/60)], ...
     [(mean_TocoC3_full'+std_mult*std_TocoC3_full'); flipud(mean_TocoC3_full'-std_mult*std_TocoC3_full')], ...
     [1 0.6 0.8], 'EdgeColor','none','FaceAlpha',0.4);
plot(tsim(:)/60, mean_TocoC3_full, '-', 'Color',[1 0 0.6],'LineWidth',2);
% Experimental points
plot(Toco_Caspase3(:,1), Toco_normData, 's','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');
title('Tocopherylbutyrate'); ylabel('C3a'); xlabel('Time (min)');
ylim([0, max(mean_TocoC3_full' + std_mult*std_TocoC3_full')*1.1])

subplot(3,2,2); hold on
fill([tsim(:)/60; flipud(tsim(:)/60)], ...
     [(mean_TocoC9_full'+std_mult*std_TocoC9_full'); flipud(mean_TocoC9_full'-std_mult*std_TocoC9_full')], ...
     [1 0.6 0.8], 'EdgeColor','none','FaceAlpha',0.4);
plot(tsim(:)/60, mean_TocoC9_full, '-', 'Color',[1 0 0.6],'LineWidth',2);
% Experimental points
plot(Toco_Caspase9(:,1), Toco_normData2, 's','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');
title('Tocopherylbutyrate'); ylabel('Normalized C9a'); xlabel('Time (min)');
ylim([0, max(mean_TocoC9_full'+std_mult*std_TocoC9_full')*1.1])

% Narciclasine 
subplot(3,2,3); hold on
fill([tsim(:)/60; flipud(tsim(:)/60)], ...
     [(mean_NarcC3_full'+std_mult*std_NarcC3_full'); flipud(mean_NarcC3_full'-std_mult*std_NarcC3_full')], ...
     [0.6 0.7 1], 'EdgeColor','none','FaceAlpha',0.4);
plot(tsim(:)/60, mean_NarcC3_full, '-', 'Color',[0.25 0.41 0.88],'LineWidth',2);
% Experimental points
plot(Narc_proCaspase3(:,1), Narc_normData, 's','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');
title('Narciclasine'); ylabel('Normalized C3'); xlabel('Time (min)');
ylim([0, max(mean_NarcC3_full'+std_mult*std_NarcC3_full')*1.1])

subplot(3,2,4); hold on
fill([tsim(:)/60; flipud(tsim(:)/60)], ...
     [(mean_NarcC8_full'+std_mult*std_NarcC8_full'); flipud(mean_NarcC8_full'-std_mult*std_NarcC8_full')], ...
     [0.6 0.7 1], 'EdgeColor','none','FaceAlpha',0.4);
plot(tsim(:)/60, mean_NarcC8_full, '-', 'Color',[0.25 0.41 0.88],'LineWidth',2);
plot(Narc_Caspase8(:,1), Narc_normData2, 's','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');
title('Narciclasine'); ylabel('Normalized C8a'); xlabel('Time (min)');
ylim([0, max(mean_NarcC8_full'+std_mult*std_NarcC8_full')*1.1])

% Celecoxib 
subplot(3,2,5); hold on
fill([tsim(:)/60; flipud(tsim(:)/60)], ...
     [(mean_Cele_full'+std_mult*std_Cele_full'); flipud(mean_Cele_full'-std_mult*std_Cele_full')], ...
     [0.8 0.6 0.9], 'EdgeColor','none','FaceAlpha',0.4);
plot(tsim(:)/60, mean_Cele_full, '-', 'Color',[0.58 0.27 0.82],'LineWidth',2);
% Experimental points
plot(Cele_Caspase3(:,1), Cele_normData, 's','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');
title('Celecoxib'); ylabel('Normalized C3a'); xlabel('Time (min)');
ylim([0, max(mean_Cele_full'+std_mult*std_Cele_full')*1.1])

% No second readout for Celecoxib
subplot(3,2,6); axis off

sgtitle('PSO Predictions: Mean ± Std with Experimental Points');

