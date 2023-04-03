%% prep data

% real data
T = readtable('data.csv'); % data should have variables as columns with the first row being names of variables

variable_names = T.Properties.VariableNames;
data = table2array(T);

% % fake data
% desired_correlation_between_x_and_y = 0.5;
% mu = 50;
% sigma = 5;
% M = mu + sigma*randn(1000,2);
% R = [1 desired_correlation_between_x_and_y; desired_correlation_between_x_and_y 1];
% L = chol(R);
% M = M*L;
% x = M(:,1);
% y = M(:,2);
% data = cat(2,x,y);
% data = zscore(data);


%%

alpha = 0.05;


%% get correlation matrix eigenvalues (lambda)

M = size(data, 2);
corr_mat = corr(data);


% convert the correlation matrix (corr_mat) to positive-definite (corr_mat_pd)
[~,flag] = chol(corr_mat,'lower');

if flag == 0
    corr_mat_pd = corr_mat;
else
    [V,D] = eig(corr_mat);  % Calculate the eigendecomposition of your matrix (A = V*D*V') where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
    d= diag(D);             % Get the eigenvalues in a vector "d"
    d(d <= 1e-7) = 1e-7;    % Set any eigenvalues that are lower than threshold "TH" ("TH" here being equal to 1e-7) to a fixed non-zero "small" value (here assumed equal to 1e-7)
    D_c = diag(d);          % Built the "corrected" diagonal matrix "D_c"
    corr_mat_pd = V*D_c*V'; % Recalculate your matrix "A" in its PD variant "A_PD"
end

lambda = eig(corr_mat_pd);


%% get effective number of tests (M_e)

% Nyholt method
% https://doi.org/10.1086/383251
% don't use
% M_e = 1 + (M - 1) * (1 - var(lambda) / M);

% Li method
% https://doi.org/10.1007/s00439-011-1118-2
M_e = M - sum( (lambda > 1) .* (lambda - 1) );


%% get adjusted significance level (alpha_prime)

alpha_prime_bonferroni = alpha / M;

alpha_prime_li = alpha / M_e;
