% Calculate tSNR
% Xueru Fan Jun 3 2025 @ BNU

clc; clear; close;
data_dir = 'D:\TravelClub\data';
resu_dir = 'D:\TravelClub\tsnr';
sites = {'A', 'B', 'C', 'D', 'E', 'G'};
sub_ids = [2, 4, 10];

temp = cifti_read('C:\TravelClub\templet\DU15NET_AgreeProb_fsLR_32k.dscalar.nii');
temp.metadata = [];

%% Step 1 计算每个OC-ME的tSNR并保存成dscalar文件

for s = 1:length(sites)
    site = sites{s};
    site_dir = fullfile(data_dir, ['3RB2_', site]);
    
    for p = 1:length(sub_ids)
        sub_id = sub_ids(p);
        sub_str = sprintf('sub-%03d', sub_id);
        sub_dir = fullfile(site_dir, sub_str);
        
        if ~exist(sub_dir, 'dir')
            fprintf('站点 %s 中被试 %s 目录不存在: %s\n', site, sub_str, sub_dir);
            continue;  % 跳过这个被试
        end

        outputFolder = fullfile(resu_dir, ['3RB2_', site], sub_str);  

        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder); 
        end

        run_folders = dir(fullfile(sub_dir, '*PA_run*'));
        
        for r = 1:length(run_folders)
            run_folder = run_folders(r).name;
            run_dir = fullfile(sub_dir, run_folder);

            boldts = ft_read_cifti(fullfile(run_dir, [run_folder, '_Atlas_s0.dtseries.nii']));
            mask = ismember(boldts.brainstructure, [1, 2]);
            data = boldts.dtseries(mask, :);
            
%             boldts = cifti_read(fullfile(run_dir, [run_folder, '_Atlas_s0.dtseries.nii']));
%             data = boldts.cdata(1:64984, :);

            mean_ts = mean(data, 2);
            std_ts = std(data, 0, 2);
            tSNR = mean_ts ./ std_ts;

            data_to_write = temp;
            data_to_write.cdata = tSNR;
            cifti_write(data_to_write, fullfile(outputFolder, [run_folder, '_tSNR.dscalar.nii']));

        end
    end
end

disp('Step 1 DONE!');

%% Step 2 计算每个被试跨run和站点的OC-ME的平均tSNR

ave_dir = fullfile(resu_dir, 'ave_subwise');
if ~exist(ave_dir, 'dir')
    mkdir(ave_dir);
end

for p = 1:length(sub_ids)
    sub_id = sub_ids(p);
    sub_str = sprintf('sub-%03d', sub_id);
    
    all_tSNR = [];

    for s = 1:length(sites)
        site = sites{s};
        site_dir = fullfile(resu_dir, ['3RB2_', site], sub_str);
        
        if ~exist(site_dir, 'dir')
            continue;  % 跳过不存在的目录
        end
        
        tSNR_files = dir(fullfile(site_dir, '*_tSNR.dscalar.nii'));
        
        for f = 1:length(tSNR_files)
            file_path = fullfile(site_dir, tSNR_files(f).name);
            tSNR_data = cifti_read(file_path);
            
            if isempty(all_tSNR)
                all_tSNR = zeros(size(tSNR_data.cdata, 1), length(tSNR_files)*length(sites));
                file_count = 0;
            end
            
            file_count = file_count + 1;
            all_tSNR(:, file_count) = tSNR_data.cdata;
        end
    end
    
    mean_tSNR = mean(all_tSNR, 2);
    
    output_data = temp;
    output_data.cdata = mean_tSNR;
    cifti_write(output_data, fullfile(ave_dir, [sub_str, '_mean_tSNR.dscalar.nii']));

end

disp('Step 2 DONE!');