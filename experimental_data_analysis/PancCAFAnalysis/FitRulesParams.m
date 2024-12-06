clearvars
addpath("./src") % two helper functions

%% load the file 
file = "./data/motility_features_30.csv";
X = readtable(file);

%% get the relevant info
speed = X.avg_speed;
cell_type = string(X.Cell_type);
ecm_conc = string(X.ECM_conc);

unique_cell_types = unique(cell_type);
unique_ecm_conc = unique(ecm_conc);

n_cell_types = length(unique_cell_types);
n_ecm_concs = length(unique_ecm_conc);

%% get mean and std of avg_speed for each condition
mean_speeds = zeros(n_cell_types, n_ecm_concs);
std_speeds = zeros(n_cell_types, n_ecm_concs);
for i = 1:n_cell_types
    for j = 1:n_ecm_concs
        index = cell_type==unique_cell_types(i) & ecm_conc==unique_ecm_conc(j);
        s = speed(index);
        mean_speeds(i,j) = mean(s);
        std_speeds(i,j) = std(s);
    end
end

%% reproduce Yoseph's bar plots
figure;
for i = 1:n_cell_types
    nexttile
    bar(mean_speeds(i,:))
    title(regexprep(unique_cell_types(i),"_"," "))
end


%%
signals = 2:7; % concentrations of collagen
optimal_pars = zeros(7,n_cell_types);
for i = 1:n_cell_types
    behaviors = mean_speeds(i,:);
    errors = std_speeds(i,:);
    fn = @(x) sum(((hillResponse(signals, x) - behaviors)./errors).^2, "omitnan");

    base_val = 0.2;
    min_val = 0;
    max_val = 3;
    decrease_ec50 = 4;
    decrease_pow = 4;
    increase_ec50 = 1;
    increase_pow = 5;

    x0 = [base_val; min_val; max_val; decrease_ec50; decrease_pow; increase_ec50; increase_pow];
    lb = zeros(length(x0),1);
    lb(4) = 1; % decrease_ec50 >= 1
    lb(5) = 1; % decrease_pow >= 1
    ub = [1; 0.4; 4; 4; 6; 10; 10];


    % Ax <= b
    A = zeros(0,7);
    b = zeros(0,1);

    A(1,[1,2]) = [-1,1]; % min value <= base value
    b(1) = 0;

    A(2,[1,3]) = [1,-1]; % base value <= max value
    b(2) = 0;

    optimal_pars(:,i) = fmincon(fn, x0, A, b, [], [], lb, ub);
end

%% normalize YLims function
% note: if using MATLAB Version < R2024a, either
%    i) move this local function to the bottom of the script OR
%   ii) put it in a separate file named normalizeYLims.m in this folder
function normalizeYLims(in)

% take a figure in and set all the axes to have the same y limits
% if input is set of axes, set them to all have the same y limits


switch class(in)
    case 'matlab.ui.Figure'
        h = findall(in,'type','axes');
        yls = arrayify(h,'YLim');
        YL = [min(yls(:,1)),max(yls(:,2))];
        set(h,'YLim',YL)

    case 'matlab.graphics.axis.Axes'
        yls = arrayify(in(:),'YLim');
        YL = [min(yls(:,1)),max(yls(:,2))];
        set(in,'YLim',YL)
end
end

%% check fits
f=figure;
cols = lines(2);
for i = 1:n_cell_types
    nexttile; hold on
    hR = hillResponse(signals, optimal_pars(:,i));
    errorbar(signals, mean_speeds(i,:), std_speeds(i,:), "Color", cols(2,:), "LineWidth", 2, "Marker","^","MarkerFaceColor",cols(2,:),"MarkerSize",8)
    plot(signals, hR, "Color", cols(1,:), "Marker","^","MarkerFaceColor",cols(1,:),"LineWidth",2,"MarkerSize",8)
    legend({"Data","Fit"})
    xlabel("ECM Density")
    ylabel("Migration Speed (micron/min)")
    set(gca, "FontSize", 16)
    title(regexprep(unique_cell_types(i),"_"," "))
end
normalizeYLims(f)

%% make table for export
T = array2table(optimal_pars);
for i = 1:size(T,2)
    T.Properties.VariableNames(i) = unique_cell_types(i);
end

T.Properties.RowNames = ["base_val (um/min)", "min_val (um/min)", "max_val (um/min)", "decrease_ec50 (mg/mL)", "decrease_pow", "increase_ec50 (mg/mL)", "increase_pow"];

writetable(T, "FittedRulesParams.csv", "WriteRowNames",true)
