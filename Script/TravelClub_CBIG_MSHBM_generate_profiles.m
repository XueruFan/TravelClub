function TravelClub_CBIG_MSHBM_generate_profiles(seed_mesh,targ_mesh,out_dir,sub,sess,split_flag)
out_profile_dir = fullfile(out_dir,'profiles',['sub' sub],['sess' sess]);
if(~exist(out_profile_dir))
    mkdir(out_profile_dir);
end
profile_file = fullfile(out_profile_dir,['sub' sub '_sess' sess '_' targ_mesh '_roi' seed_mesh '.surf2surf_profile.mat']);
fMRI_list = fullfile(out_dir,'data_list','fMRI_list',['sub' sub '_sess' sess '.txt']);
censor_list = fullfile(out_dir,'data_list','censor_list',['sub' sub '_sess' sess '.txt']);
if(~exist(censor_list))
    censor_list = 'NONE';
end
TravelClub_CBIG_ComputeCorrelationProfile(seed_mesh,targ_mesh, profile_file, 'NONE', fMRI_list, 'NONE', censor_list, split_flag);
end

