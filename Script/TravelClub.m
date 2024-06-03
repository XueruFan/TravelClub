%% Step 1
project_dir = 'E:\PhDproject\HCP\HCP_test4MSHBM'
TravelClub_MSHBM_Prep(project_dir)
%% Step 2
% see Kong2019_MSHBM README.md for input suggestions
TravelClub_MSHBM_SetInput

%% Step 3
% Generating profiles and initialization parameters
% make sure you have added the CBIG folder into your path

for sub = 1:str2num(nsub)
 for sess = 1:str2num(maxsess)
    CBIG_MSHBM_generate_profiles(seed_mesh, targ_mesh, project_dir, ...
        num2str(sub), num2str(sess), '0');
 end
end

%% Step 4
% To obtain the group averaged profiles
CBIG_MSHBM_avg_profiles(seed_mesh, targ_mesh, project_dir, nsub, maxsess);

%% Step 5
% To run Yeo2011 clustering algorithm for generate or own group prior
% 这里用的就是自己的数据
num_clusters = 15;
num_initialization = 2; % use 1000 for real analysis
CBIG_MSHBM_generate_ini_params(seed_mesh, targ_mesh, num_clusters, ...
    num_initialization, project_dir)

% 这一步Du用HCP的40个人的数据，从2~21个网络都run了一遍，各做了一个组水平的
% 先验，最终是结合后面基于ROI的连接图，确定了15个网络
%% Step 6



