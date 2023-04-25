% set variables

group1_X = rand(50,1);
group1_y = rand(50,1);

group2_X = rand(48,1);
group2_y = rand(48,1);

num_bootsamples = 100;

% do bootstrap regression

group1_sample_size = length(group1_y);
group1_index_orig = (1 : group1_sample_size)';
[~, group1_index_boot] = bootstrp(num_bootsamples, [], group1_index_orig);
group1_index_boot_complete = cat(2, group1_index_orig, group1_index_boot);


group2_sample_size = length(group2_y);
group2_index_orig = (1 : group2_sample_size)';
[~, group2_index_boot] = bootstrp(num_bootsamples, [], group2_index_orig);
group2_index_boot_complete = cat(2, group2_index_orig, group2_index_boot);

boot_ES = zeros(num_bootsamples+1, 1);

warning('off', 'stats:regress:RankDefDesignMat')
for boot_counter = 1 : length(boot_ES)

    tmp_group1_idx = group1_index_boot_complete(:,boot_counter);
    tmp_group1_y = group1_y(tmp_group1_idx, :);
    tmp_group1_X = group1_X(tmp_group1_idx, :);

    tmp_group2_idx = group2_index_boot_complete(:,boot_counter);
    tmp_group2_y = group2_y(tmp_group2_idx, :);
    tmp_group2_X = group2_X(tmp_group2_idx, :);

    tmp_X = cat(1, tmp_group1_X, tmp_group2_X);
    tmp_y = cat(1, tmp_group1_y, tmp_group2_y);
    
    tmp_col_ones = ones(length(tmp_y), 1);
    tmp_X = cat(2, tmp_col_ones, tmp_X);

    tmp_boot_ES_all = regress(tmp_y, tmp_X);
    boot_ES(boot_counter) = tmp_boot_ES_all(2);

end
warning('on', 'stats:regress:RankDefDesignMat')

bootstrap_CI90 = prctile(boot_ES, [5 95])
bootstrap_CI95 = prctile(boot_ES, [2.5 97.5])
bootstrap_CI99 = prctile(boot_ES, [0.5 99.5])
