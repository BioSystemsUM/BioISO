parent_parent_dir = dir ('fast_gap_fill_models');

times = zeros(2,5);
models_names = cell(2, 5);
names = cell(2,5);
rxns = zeros(2,5);
mets = zeros(2,5);
failled = cell(2,5);

for i = 1:size(parent_parent_dir)
    
    if ~strcmp(parent_parent_dir(i).name,'.') && ~strcmp(parent_parent_dir(i).name,'..')
        
        parent_parent_dir(i).name
        concat_dir = strcat('fast_gap_fill_models/', parent_parent_dir(i).name);
        parent_dir = dir (concat_dir);
        
        for j = 1:size(parent_dir)
            
            if ~strcmp(parent_dir(j).name,'.') && ~strcmp(parent_dir(j).name,'..')
                
                parent_dir(j).name
                model_dir = strcat(concat_dir, '/', parent_dir(j).name);
                
                model = readCbModel(model_dir, 'fileType', 'SBML');
                try
                    tic;
                    dead_end_metabolites = gapFind(model, 'true');
                    blocked_reactions = findBlockedReaction(model);
                    final_t = toc;
                    fail = 'no';
                catch
                    dead_end_metabolites = cell(1,1);
                    blocked_reactions = cell(1,1);
                    final_t = 0;
                    fail = 'yes';
                end
                times(i, j) = final_t;
                failled{i, j} = fail;
                models_names{i, j} = parent_parent_dir(i).name;
                names{i, j} = parent_dir(j).name;
                [mets_n_rows, mets_n_cols] = size(dead_end_metabolites);
                [rxns_n_rows, rxns_n_cols] = size(blocked_reactions);
                mets(i, j) = mets_n_rows;
                rxns(i, j) = rxns_n_cols;
                
            end
        end
    end
end









