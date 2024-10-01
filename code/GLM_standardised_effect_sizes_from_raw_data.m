%% load data

clear;

data = readtable("data.csv");

%% set variables

n = size(data, 1); % sample size
X = ones(n, 1); % design matrix, intercept-only model
y = data.Magnitude; % outcomes
w = data.Weight; % outcome weights

%% fit weighted regression model

mdl = fitlm(X, y, Weights=w, Intercept=false); % turn off intercept to stop matlab automatically adding intercept to X

%% 

t = mdl.Coefficients.tStat; % t-statistic

df = mdl.DFE; % degrees of freedom for error


%% calculate standardised effect sizes

pearsons_r = sign(t) * sqrt(t^2 / (length(y) - rank(X) + t^2)); % t-to-r from brainder (https://brainder.org/2015/03/04/all-glm-formulas/)

cohens_f2 = pearsons_r^2 / (1 - pearsons_r^2); % r-to-f from G*Power 3.1 paper (https://doi.org/10.3758/BRM.41.4.1149)

cohens_d = (2*pearsons_r) / sqrt(1 - pearsons_r^2); % r-to-d from textbook 'Introduction to meta-analysis' Ch7 pg45 eq7.5

j = 1 - (3 / ((4*df) - 1)); % Hedges's g correction factor from textbook 'Introduction to meta-analysis' Ch4 pg27 eq4.22

hedges_g = j * cohens_d; % d-to-g from textbook 'Introduction to meta-analysis' Ch4 pg27 eq4.23
