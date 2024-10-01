%% make data
clear;

% set parameters for distributions from which data will be drawn
% these parameters are set to those used for Petal trial sample size
% calculation
mu_1 = 1.07;
sigma_1 = 0.66;
mu_2 = 0.642;
sigma_2 = 0.66;

n_1 = 51;
n_2 = 51;

%%

% O19 = mu_1;
% O20 = sigma_1;
% O21 = n1;
% 
% Q19 = mu_2;
% Q20 = sigma_2;
% Q21 = n2;

tstat_numerator = mu_1 - mu_2;
tstat_denominator = sqrt(((((n_1-1)*sigma_1^2)+((n_2-1)*sigma_2^2))/(n_1+n_2-2))*((1/n_1+1/n_2)));
tstat = tstat_numerator / tstat_denominator;

cohensD = (mu_1-mu_2) / (sqrt((((n_1-1)*sigma_1^2)+((n_2-1)*sigma_2^2))/(n_1+n_2-2)));
cohensD_var = ((n_1+n_2)/(n_1*n_2)) + (cohensD^2)/(2*(n_1+n_2));
cohensD_se = sqrt(cohensD_var);

j = (1-(3/(4*(n_2+n_1-2)-1)));
hedgesG2 = j * cohensD;
hedgesG2_var = (j^2) * cohensD_var;
hedgesG2_se = sqrt(hedgesG2_var);

a = ((n_1 + n_2)^2) / (n_1 * n_2);
r_from_d = cohensD / sqrt(cohensD^2 + a);
r2_from_d = r_from_d^2;
cohensf2_from_cohensD = r2_from_d / (1 - r2_from_d);


%%

% set sample sizes to be large (e.g. 10000) so the result will closely
% match the parameters above and thus the sample size calculation should
% match the Petal trial (n=102 for alpha=0.5, power=0.9, two-tailed,
% allocation ration = 1
% set sample sizes to be small (e.g. 100) to get variable results
n_1_sim = 10000;
n_2_sim = 10000;

% draw dataset 1 and get mean and standard deviation
d_1 = normrnd(mu_1, sigma_1, n_1_sim, 1);
mean_1 = mean(d_1);
std_1 = std(d_1);

% draw dataset 2 and get mean and standard deviation
d_2 = normrnd(mu_2, sigma_2, n_2_sim, 1);
mean_2 = mean(d_2);
std_2 = std(d_2);

%% calculate Cohen's D

% formula for Cohen's D is equation 1 from here: https://doi.org/10.3389/fpsyg.2013.00863 
numerator = (mean_1 - mean_2);
denominator = sqrt((((n_1_sim - 1) * std_1^2) + ((n_2_sim - 1) * std_2^2)) / (n_1_sim + n_2_sim - 2));

Cohens_D = numerator / denominator;

%% calculate Cohen's f2

y = [d_1; d_2];
x1 = ones(size(y));
x2 = [ones(size(d_1)); ones(size(d_2))*(-1)];
X = [x1, x2];
mdl = fitlm(X, y, Intercept=false);
t = mdl.Coefficients.tStat(2);

% formula for partial correlation coefficient from t-statistic is from
% here: https://brainder.org/2015/03/04/all-glm-formulas/
r = sign(t) * sqrt(t^2 / (length(y) - rank(X) + t^2));
r2 = r^2;

% formula for Cohen's f2 is from here: https://doi.org/10.3758/BRM.41.4.1149
Cohens_f2 = r2 / (1 - r2);

%% Convert between Cohen's D and Cohen's f2

Cohens_f2_from_Cohens_D = Cohens_D^2 / 4;

Cohens_D_from_Cohens_f2 = 2 * sqrt(Cohens_f2);

Cohens_D_from_r = (2 * r) / sqrt(1 - r^2);
