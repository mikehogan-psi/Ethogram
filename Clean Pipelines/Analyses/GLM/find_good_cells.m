% FIND_GOOD_CELLS
% 
% Calculates empirical p value for each cell based on shuffled and real
% pseudo-R2 values and then applies 5% false discovery rate correction to
% generate a mask for good cells to keep (that perform better than shuffled
% model)
%
% INPUT: pseudoR2- n_cells x 1 double of real pR2 values
%        pseudoR2_shuff_all- n_cells x 1 cell array of shuffled pR2 values
%
% OUTPUT: keep_mask- logical mask with indicies of cells to keep

function[keep_mask] = find_good_cells(pseudoR2, pseudoR2_shuff_all)

n_cells = numel(pseudoR2);
p_emp = nan(n_cells, 1);

for neuron = 1:n_cells
    current_n = pseudoR2_shuff_all{neuron};
    current_n = current_n(isfinite(current_n));   % drop NaN/Inf
    
    if isempty(current_n) || ~isfinite(pseudoR2(neuron))
        continue
    end
    p_emp(neuron) = (1 + sum(current_n >= pseudoR2(neuron)))/(1 + numel(current_n));
end

q = mafdr(p_emp, 'BHFDR', true);
keep_mask = (q <= 0.05) & isfinite(p_emp) & isfinite(pseudoR2);

end
