%% load data

T = readtable('DatasetQ4.csv');

variable_names = T.Properties.VariableNames;
data = table2array(T);
num_variables = size(data,2);

%% perform cluster analysis on correlation matrix

correlationMatrix = corr(data);

correlationMatrix_maxValue = max(max(abs(tril(correlationMatrix, -1))));

correlationMatrix_distancesMatrix = sqrt(2) * sqrt(1 - correlationMatrix);

hierarchical_cluster_tree = linkage(correlationMatrix_distancesMatrix, 'ward');


%% plot results

% general plot variables

dendrogram_height=0.8;

correlationMatrixPlot_height = 0.7;

correlationMatrixPlot_width = 0.5 / (num_variables + 1);

overall_fig_size_x = 7;

overall_fig_size_y = 7.5;


figure('position',[10 10 1000 500]);


% dendrogram

% subplot('position',[0 dendrogram_height 1 1-dendrogram_height]);
subplot('position',[0 dendrogram_height 0.91 1-dendrogram_height]);

colour_threshold = 'default';
[handles_to_lines, ~, variable_order] = dendrogram(hierarchical_cluster_tree, 0, 'colorthreshold', colour_threshold);
set(handles_to_lines, 'LineWidth', 3);
ax = gca;
ax.YTick = [];
ax.YColor = 'none';
ax.XRuler.Axle.Visible = 'off';
ax.XTickLabel = variable_names(variable_order);


% correlation matrix

subplot('position',[correlationMatrixPlot_width 0 1-2*correlationMatrixPlot_width correlationMatrixPlot_height-0.01]);

correlationMatrix_clusteredForPlotting = correlationMatrix;
correlationMatrix_clusteredForPlotting = correlationMatrix_clusteredForPlotting(variable_order, variable_order);
correlationMatrix_clusteredForPlotting(eye(length(correlationMatrix_clusteredForPlotting))>0) = Inf; % set diagonal to infinity for nicer plot

colormap('jet');
colourmapHandle = colormap;
colourmapHandle(end,:) = [.8 .8 .8];
colormap(colourmapHandle);

heatmap_colorRange = [(correlationMatrix_maxValue + 0.1)*(-1), (correlationMatrix_maxValue + 0.1)];
% heatmap_colorRange = [-1, 1];

imagesc(correlationMatrix_clusteredForPlotting, heatmap_colorRange);

c = colorbar('eastoutside');
c = colorbar;
c.Label.String = 'Pearson correlation';

axis off;
daspect('auto');
set(gcf,...
    'Units', 'Inches', ...
    'Position', [0, 0, overall_fig_size_x, overall_fig_size_y], ...
    'PaperPositionMode', 'auto');

leftShift_to_centre_corValues = 0.15;
for cormat_row_position = 1: num_variables
    for cormat_column_position = 1 : num_variables
        text( (cormat_column_position - leftShift_to_centre_corValues), cormat_row_position, num2str(correlationMatrix_clusteredForPlotting(cormat_row_position, cormat_column_position), '%.2f'), 'Color', [0.8 0.8 0.8]);
    end
end

saveas(gcf, 'correlation_clustering.jpeg');
