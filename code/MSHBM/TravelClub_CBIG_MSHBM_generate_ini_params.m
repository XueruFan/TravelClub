function process_network_threshold(targ_mesh, network, output_file, avg_profile_file, niter)

    temp_file = [tempname '.mat'];
    
    try
        CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(...
            targ_mesh, '', network, temp_file, ...
            avg_profile_file, 'NONE', 0, niter, 0, 1000, 1);
        
        if exist(temp_file, 'file')
            data = load(temp_file);
            clustered.mtc = data.mtc;
            clustered.lowerbound = data.lowerbound;
            clustered.lambda = data.lambda;
            
            save(output_file, 'lh_labels', 'rh_labels', 'clustered', '-v7.3');
            delete(temp_file);
        end
    catch ME
        if exist(temp_file, 'file')
            delete(temp_file); 
        end
        warning('Error: %s\n错误: %s', output_file, ME.message);
    end
end