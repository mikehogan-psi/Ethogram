master_directory = 'Z:\Mike\Data\Psilocybin Fear Conditioning\Cohort 4_06_05_25 (SC PAG Implanted Animals)';

session = 'Extinction';

% Select only mouse data folders
mouse_files = dir(fullfile(master_directory, 'Mouse*'));

mice_to_analyse = [2 3 4 5 6 7];

num_mice = numel(mice_to_analyse);

combined_data_folders = cell(num_mice, 1);

for mouse = mice_to_analyse
        mouse_name = mouse_files(mouse).name;
        mouse_path = mouse_files(mouse).folder;
        combined_data_folders{mouse} = fullfile(mouse_path, mouse_name, ...
            session, 'Combined Data'); 
end

X_data = cell(num_mice, 1);
cell_ids_all = cell(num_mice, 1);

for mouse = mice_to_analyse
    current_data = dir(fullfile(combined_data_folders{mouse}, '*mouse*'));
    current_data_file = fullfile(current_data.folder, current_data.name);
    current_variables = load(current_data_file, 'X', 'cell_ids');
    X_data{mouse} = current_variables.X;
    cell_ids_all{mouse} = current_variables.cell_ids;
end

X_data_cat = [];
cell_ids_cat = [];

for mouse = mice_to_analyse
    X_data_reshaped = reshape(X_data{mouse}, 10, 2640, []);
    X_data_cat = cat(3, X_data_cat, X_data_reshaped);
    cell_ids_reshaped = reshape(cell_ids_all{mouse}, 2640, []);
    cell_ids_cat = cat(2, cell_ids_cat, cell_ids_reshaped);
end

cell_ids_cat = unique(cell_ids_cat(1, :));

load("C:\Users\G71044MH\OneDrive - The University of Manchester\Documents\GitHub\Ethogram\Neural Data Analysis\PoissonNN\results\results.mat")

%%
valid_cells = pseudoR2 > 0.05;
n_valid = sum(valid_cells);

good_weights = weights(valid_cells, :);
good_biases = biases(valid_cells);
good_cells = cell_ids_cat(valid_cells);
good_cell_idxs = find(valid_cells);

% Last 4 covariates are binary, no need to normalise
isBinary = [false, false, false, false, false, false, true, true, true, true];

norm_X = zeros(size(X_data_cat));

for neuron = 1:size(X_data_cat, 3)
    for covariate = 1:size(X_data_cat, 1)
        if isBinary(covariate)
            current_row = X_data_cat(covariate, :, neuron);         
        else
            current_row = zscore(X_data_cat(covariate, :, neuron));
            current_row(~isfinite(current_row)) = 0;
        end
        norm_X(covariate, :, neuron) = current_row;
    end
end
%%
pseudoR2 = pseudoR2';
weights = weights';
biases = biases';
cell_ids = cell_ids';
X = X';
y = y';
y_pred = double(y_pred');

valid_cells = pseudoR2 > 0.05;
n_valid = sum(valid_cells);

good_weights = weights(:, valid_cells);
good_biases = biases(valid_cells);
good_cells = cell_ids(valid_cells);

% Last 4 covariates are binary, no need to normalise
isBinary = [false, false, false, false, false, false, true, true, true, true];

norm_X = zeros(size(X));

for neuron = 1:size(X, 2)
    for covariate = 1:size(X, 3)
        if isBinary(covariate)
            current_row = X(:, neuron, covariate);         
        else
            current_row = zscore(X(:, neuron, covariate));
            current_row(~isfinite(current_row)) = 0;
        end
        norm_X(:, neuron, covariate) = current_row;
    end
end
    
valid_idx = find(valid_cells);

glm_output = nan(n_valid, size(X, 3));

for neuron = 1:n_valid
    current_neuron = valid_idx(neuron);
    current_y_pred = y_pred(:, current_neuron);
    current_X = squeeze(norm_X(:, current_neuron, :));
    current_X = double(current_X);
    b = glmfit(current_X, current_y_pred, 'poisson', 'link', 'log');
    glm_output(neuron, :) = b(2:end)';
end


%%

X = reshape(X, 2640, 202, 10);

% Last 4 covariates are binary, no need to normalise
isBinary = [false, false, false, false, false, false, true, true, true, true];

% Z-score covariates (except binary ones)
norm_X = zeros(size(X));

for neuron = 1:size(X, 2)
    for covariate = 1:size(X, 3)
        if isBinary(covariate)
            current_row = X(:, neuron, covariate);         
        else
            current_row = zscore(X(:, neuron, covariate));
            current_row(~isfinite(current_row)) = 0;
        end
        norm_X(:, neuron, covariate) = current_row;
    end
end

valid_idx = find(valid_cells);

glm_output = nan(n_valid, size(X, 3));

for neuron = 1:n_valid
    current_neuron = valid_idx(neuron);
    current_y_pred = y_pred(:, current_neuron);
    current_X = squeeze(norm_X(:, current_neuron, :));
    current_X = double(current_X);
    b = glmfit(current_X, current_y_pred, 'poisson', 'link', 'log');
    glm_output(neuron, :) = b(2:end)';
end

%%

glm_output_z = zscore(glm_output, 0, 2);  % z-score across covariates per neuron

imagesc(glm_output_z)
colorbar
xlabel('Covariates'); ylabel('Neurons');
title('GLM beta weights for each neuron');
%%
[coeff, score, ~, ~, explained] = pca(glm_output_z);
scatter(score(:,1), score(:,2), 40, pseudoR2(valid_cells), 'filled')
xlabel(['PC1 (' num2str(round(explained(1))) '%)'])
ylabel(['PC2 (' num2str(round(explained(2))) '%)'])
title('Neurons projected onto first two tuning PCs')
colorbar; ylabel(colorbar,'pseudoR^2');

%%
k = 3;  % try different k’s
[idx, C] = kmeans(score(:,1:3), k, 'Replicates',10);

figure;
gscatter(score(:,1), score(:,2), idx);
xlabel('PC1'); ylabel('PC2');
title('Neuron clusters based on behavioural tuning');

%%
nCov = size(glm_output,2);
mean_betas = zeros(k, nCov);

for c = 1:k
    mean_betas(c,:) = mean(glm_output_z(idx==c,:), 1);
end

figure;
imagesc(mean_betas)
colorbar
yticks(1:k)
ylabel('Cluster')
xticks(1:nCov)
xticklabels({'Grooming','Rearing','Darting','Freezing','Velocity','TimeBin','TrialNum', 'TrialID', 'LoomOn','FlashOn'})
xtickangle(45)
title('Average tuning pattern per cluster (z-scored betas)')
