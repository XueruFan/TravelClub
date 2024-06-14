function TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh, target, output_file1, output_file2,...
    varargin_file1, varargin_file2, outlier_text, split_data)

s_mesh = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', ...
    'fs_LR_32k_downsample_900', 'fslr_downsample_900mesh_parcellation.dlabel.nii'));

% read in both left and right text files.
fid = fopen(varargin_file1, 'r');
i = 0;
while(1);
tmp = fscanf(fid, '%s\n', 1);
    if(isempty(tmp))
        break
    else
        i = i + 1;
        varargin1{i} = tmp;
    end
end
fclose(fid);

% Compute profile for left hemi
for i = 1:length(varargin1)
    if(exist('outlierin', 'var'))
        % {0,1} vector, where uncensored time points are 0
        outliers = dlmread(outlierin{i}); 
    end
    
    input = varargin1{i};
    
    % read input file
    [input_series, t_series, input_size] = read_fmri(input);
    t_series = t_series';
    if(exist('outliers', 'var'))
        t_series = t_series(outliers==1, :);
    end
    
    lh_seed_series = t_series;
    seed_ind = find(input_series.brainstructure==1 | input_series.brainstructure==2);
    seed_ind = seed_ind(s_mesh.dlabel~=0);
    s_series = lh_seed_series(:, seed_ind);
    
    % normalize series (note that series are now of dimensions: T x N) and compute correlation
    if(split_data~=0 && length(varargin1)==1)
        for split = 1:split_data
            rang_start = floor(size(s_series,1)/split_data) * (split - 1) + 1;
            rang_end = floor(size(s_series,1)/split_data) * split;
            if(split == split_data)
                rang_end = size(s_series,1);
            end
            curr_s_series = s_series(rang_start:rang_end,:);
            curr_t_series = t_series(rang_start:rang_end,:);
            curr_corr_mat1 = CBIG_corr(curr_s_series, curr_t_series);

            disp(['Fake run ' num2str(split) ', isnan: ' num2str(sum(sum(isnan(curr_corr_mat1)))) ' out of '...
            num2str(numel(curr_corr_mat1))]);
            tmp_corr = curr_corr_mat1;
            tmp_corr(isnan(curr_corr_mat1)) = 0;
            corr_mat1(:,:,split) = tmp_corr;
        end
        if(i == 1)
            output = corr_mat1;
        else
            output = output + corr_mat1;
        end
    else        
        corr_mat1 = CBIG_corr(s_series, t_series);
        corr_mat1(isnan(corr_mat1)) = 0;
        if(i == 1)
            output = corr_mat1;
        else
            output = output + corr_mat1;
        end
    end
    clear outliers
end

output = output / length(varargin1);
corr_mat1 = output;
corr_mat1(corr_mat1 > 0) = 1;
corr_mat1(corr_mat1 < 0) = 0;
write_fmri(output_file1, input_series, corr_mat1', [input_size(1:end-1) size(corr_mat1, 1)]);


function [fmri, vol, vol_size] = read_fmri(fmri_name)
fmri = ft_read_cifti(fmri_name);
vol = fmri.dtseries;
vol = vol(fmri.brainstructure==1|fmri.brainstructure==2, :);
vol_size = size(vol);
fmri.dtseries = [];


function write_fmri(fmri_name, fmri, vol, vol_size)
profile_mat = single(reshape(vol, vol_size));
save(fmri_name, 'profile_mat', '-v7.3');


