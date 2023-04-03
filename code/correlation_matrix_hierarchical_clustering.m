%% load data

T = readtable('data.csv'); % data should have variables as columns with the first row being names of variables

variable_names = T.Properties.VariableNames;
data = table2array(T);
num_variables = size(data,2);

%% perform cluster analysis on correlation matrix

correlationMatrix = corr(data);

correlationMatrix_maxValue = max(max(abs(tril(correlationMatrix, -1))));

correlationMatrix_distancesMatrix = 1 - correlationMatrix;

correlationMatrix_distancesVector = squareform(correlationMatrix_distancesMatrix);

hierarchical_cluster_tree = linkage(correlationMatrix_distancesVector, 'ward');


%% plot results

% general plot variables

dendrogram_height=0.8;

correlationMatrixPlot_height = 0.7;

correlationMatrixPlot_width = 0.5 / (num_variables + 1);

overall_fig_size_x = 7;

overall_fig_size_y = 7;


figure('position',[10 10 1000 500]);


% dendrogram

subplot('position',[0 dendrogram_height 1 1-dendrogram_height]);

colour_threshold = 'default';
[handles_to_lines, ~, variable_order] = dendrogram(hierarchical_cluster_tree, 0, 'colorthreshold', colour_threshold);

xticklabels = variable_names(variable_order);
set(gca, 'ytick', [], 'xticklabel', xticklabels, 'TickLength', [0,0]);
set(handles_to_lines, 'LineWidth', 3);


% correlation matrix

subplot('position',[correlationMatrixPlot_width 0 1-2*correlationMatrixPlot_width correlationMatrixPlot_height-0.01]);

correlationMatrix_clusteredForPlotting = correlationMatrix;
correlationMatrix_clusteredForPlotting = correlationMatrix_clusteredForPlotting(variable_order, variable_order);
correlationMatrix_clusteredForPlotting(eye(length(correlationMatrix_clusteredForPlotting))>0) = Inf; % set diagonal to infinity for nicer plot

colormap('jet');
grotcolourmapHandle = colormap;
grotcolourmapHandle(end,:) = [.8 .8 .8];
colormap(grotcolourmapHandle);

imagesc(correlationMatrix_clusteredForPlotting,[(correlationMatrix_maxValue + 0.1)*(-1), (correlationMatrix_maxValue + 0.1)]);

axis off;
daspect('auto');

set(gcf,...
    'Units', 'Inches', ...
    'Position', [0, 0, overall_fig_size_x, overall_fig_size_y], ...
    'PaperPositionMode', 'auto');

saveas(gcf, 'correlation_clustering.jpeg');
