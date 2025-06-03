%% Pipeline of individual networks mapping
% 3R-BRAIN Project
% Written by Xue-Ru Fan 2025-06-03 @ Beijing Normal University

%% Step 12 Agreement(%) of Consensus maps

% n=1;t=1;i=1;
for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));
        
        output_dir = fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation');
        consensus_file = fullfile(output_dir, 'Consensus.mat');

        load(consensus_file);

        matFiles = dir(fullfile(output_dir, 'sub-*.mat'));
        
        lh_count = zeros(size(lh_consensus));
        rh_count = zeros(size(rh_consensus));
        
        for i = 1:nsub

            sub_data = load(fullfile(output_dir, matFiles(i).name));
            
            lh_match = (sub_data.lh_labels == lh_consensus);
            lh_count = lh_count + double(lh_match);
            
            rh_match = (sub_data.rh_labels == rh_consensus);
            rh_count = rh_count + double(rh_match);
        end
        
        lh_percentage = (lh_count / nsub) * 100;
        rh_percentage = (rh_count / nsub) * 100;
        
        result_file = fullfile(output_dir, 'Agreement.mat');
        save(result_file, 'lh_percentage', 'rh_percentage');
        
    end
end

disp("Step 12 DONE!")

%% Step 13 Make agreement map visualizations

for n = 1:xnet
    n_str = num2str(networks(n));
    for t = 1:xthe
        t_str = num2str(threshold(t));

        vis_dir = fullfile(data_dir, 'visual', ['t-' t_str '_n-' n_str]);
        load(fullfile(data_dir, 'output', ['t-' t_str '_n-' n_str], 'ind_parcellation', 'Agreement.mat'));
    
        temp = cifti_read(fullfile(temp_dir, "DU15NET_AgreeProb_fsLR_32k.dscalar.nii"));
        temp.metadata = [];
        temp.cdata = [lh_percentage; rh_percentage];

        cifti_write(temp, fullfile(vis_dir, 'Agreement.dscalar.nii'));

    end
end

disp("Step 13 DONE!")