function Params = TravelClub_CBIG_MSHBM_estimate_group_priors(project_dir,mesh,num_sub,num_sess,site,num_clusters,varargin)
pnames = {'max_iter' 'conv_th' 'save_all'};
dflts =  {'50' '1e-5' '0'};
[max_iter, conv_th, save_all] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%% setting parameters
setting_params.num_sub = str2double(num_sub); 
setting_params.num_session = str2double(num_sess);
setting_params.num_clusters = str2double(num_clusters);
setting_params.epsilon = 1e-4;
setting_params.conv_th = str2double(conv_th);
setting_params.max_iter = str2double(max_iter);
setting_params.save_all = str2double(save_all);
setting_params.mesh = mesh;

%% read in data
data = fetch_data(project_dir, setting_params.num_session, site, setting_params.num_sub, setting_params.mesh);
setting_params.dim = size(data.series,2) - 1;
if(setting_params.dim < 1200) 
    setting_params.ini_concentration = 500; 
elseif(setting_params.dim >= 1200 && setting_params.dim < 1800 )
    setting_params.ini_concentration = 650;                                          
else
    error(['Data dimension is higher than 1800,' ...
          'please set the minimal concentration parameter value ini_concentration,' ...
          'where besseli((D/2-1),ini_concentration) > 0 and relatively small\n']);
end    

%% paramter initialization
% load group parcellation and parameters
group = load(fullfile(project_dir, 'group', 'group.mat'));
setting_params.g_mu = transpose(group.clustered.mtc);
Params.mu = setting_params.g_mu;
Params.epsil = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
Params.s_psi = repmat(setting_params.g_mu, 1, 1, setting_params.num_sub);
Params.sigma = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
Params.s_t_nu = repmat(setting_params.g_mu, 1, 1, setting_params.num_session, setting_params.num_sub);
Params.kappa = setting_params.ini_concentration*ones(1, setting_params.num_clusters);

log_vmf = permute(Params.s_t_nu, [1,2,4,3]);
log_vmf = mtimesx(data.series, log_vmf);%NxLxSxT
log_vmf = bsxfun(@times, permute(log_vmf,[2,1,3,4]), transpose(Params.kappa));%LxNxSxT
log_vmf = bsxfun(@plus,Cdln(transpose(Params.kappa), setting_params.dim), log_vmf);%LxNxSxT
log_vmf = CBIG_nansum(log_vmf, 4);
log_vmf = permute(log_vmf, [2,1,3]);
s_lambda = bsxfun(@minus, log_vmf, max(log_vmf,[],2));
mask = repmat((sum(s_lambda,2)==0), 1, setting_params.num_clusters, 1);
s_lambda = exp(s_lambda);
Params.s_lambda = bsxfun(@times, s_lambda, 1./sum(s_lambda,2));
Params.s_lambda(mask) = 0;

Params.theta = mean(Params.s_lambda, 3);
theta_num = sum(Params.s_lambda ~= 0, 3);
Params.theta(Params.theta ~= 0) = Params.theta(Params.theta ~= 0)./theta_num(Params.theta ~= 0);

%% EM
stop_inter = 0;
cost_inter = 0;
Params.iter_inter = 0;
while(stop_inter == 0)
    Params.iter_inter = Params.iter_inter + 1;
    %% Intra subject variability
    cost = 0;
    iter_intra_em = 0;
    stop_intra_em = 0;
    Params.sigma = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
    Params.s_psi = repmat(setting_params.g_mu, 1, 1, setting_params.num_sub);
    while(stop_intra_em == 0)
        iter_intra_em = iter_intra_em + 1;

        Params.kappa = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
        Params.s_t_nu = repmat(setting_params.g_mu, 1, 1, setting_params.num_session, setting_params.num_sub);
        Params = vmf_clustering_subject_session(Params, setting_params, data);

        Params = intra_subject_var(Params, setting_params);

        update_cost = bsxfun(@times, Params.s_psi, permute(Params.s_t_nu, [1,2,4,3]));
        update_cost = sum(bsxfun(@times, Params.sigma, update_cost), 1);
        update_cost = bsxfun(@plus, Cdln(Params.sigma, setting_params.dim), update_cost);
        update_cost = CBIG_nansum(sum(sum(update_cost, 2), 3), 4) ...
                    + sum(sum(bsxfun(@plus, ...
                      sum(bsxfun(@times, bsxfun(@times, Params.mu, Params.s_psi),Params.epsil), 1), ...
                      Cdln(Params.epsil, setting_params.dim)), 2), 3);
        update_cost = update_cost + sum(Params.cost_em);
        if(abs(abs(update_cost - cost)./cost) <= setting_params.epsilon)
            stop_intra_em = 1;
            Params.cost_intra = update_cost;
        end
        if(iter_intra_em >= 50)
            stop_intra_em = 1;
            Params.cost_intra = update_cost;
        end
        cost = update_cost;
    end

    %% Inter subject variability
    Params = inter_subject_var(Params, setting_params);

    update_cost_inter = Params.cost_intra;    
    Params.Record(Params.iter_inter) = update_cost_inter;
    if(abs(abs(update_cost_inter - cost_inter)./cost_inter) <= setting_params.conv_th || ...
            Params.iter_inter >= setting_params.max_iter)
        stop_inter = 1;
        Params.cost_inter = update_cost_inter;

        % set s_lambda, s_psi, s_t_nu to be empty to save space and time
        if(save_all==0)
            Params.s_lambda = [];
            Params.s_psi = [];
            Params.s_t_nu = [];
        end
        if(~exist(fullfile(project_dir, 'priors')))
            mkdir(fullfile(project_dir, 'priors'));
        end
        save(fullfile(project_dir, 'priors', 'Params_Final.mat'), 'Params');
    end
    cost_inter = update_cost_inter;
    if(~exist(fullfile(project_dir, 'priors')))
        mkdir(fullfile(project_dir, 'priors'));
    end
    save(fullfile(project_dir, 'priors', ['Params_iteration',num2str(Params.iter_inter),'.mat']), 'Params');
end

end


%% sub-functions

function Params=inter_subject_var(Params,setting_params)

% Inter-subject variability level

% update mu
mu_update = sum(Params.s_psi, 3);
mu_update = bsxfun(@times, mu_update, 1./sqrt(sum((mu_update).^2)));

% update epsil
epsil_update = bsxfun(@times, Params.s_psi, mu_update);
epsil_update = sum(sum(epsil_update, 1), 3);
epsil_update = epsil_update./setting_params.num_sub;
for i = 1:setting_params.num_clusters
    epsil_update(i) = invAd(setting_params.dim, epsil_update(i));
    if(epsil_update(i) < setting_params.ini_concentration)
        epsil_update(i) = setting_params.ini_concentration;
        fprintf('[WARNING] epsil of %d is less than the minimal value %d \n',i, setting_params.ini_concentration);
    end
    if(isinf(epsil_update(i)))
        epsil_update(i) = Params.epsil(i);
        fprintf('[WARNING] epsil of %d is Inf\n',i);
    end
end
Params.mu = mu_update;
Params.epsil = epsil_update;
end

function Params = intra_subject_var(Params,setting_params)

% Intra-subject variability level

flag_psi = zeros(setting_params.num_sub, 1);
stop_intra = 0;
iter_intra = 0;
while(stop_intra == 0)
    iter_intra = iter_intra + 1;
    fprintf('It is inter interation %d intra iteration %d..update s_psi and sigma..\n',Params.iter_inter,iter_intra);
    % update s_psi
    s_psi_update = CBIG_nansum(bsxfun(@times,Params.s_t_nu,repmat(Params.sigma,size(Params.s_t_nu,1), ...
                   1,size(Params.s_t_nu,3),size(Params.s_t_nu,4))),3);
    s_psi_update = reshape(s_psi_update,...
                   size(s_psi_update,1),size(s_psi_update,2),size(s_psi_update,3)*size(s_psi_update,4));
    s_psi_update = bsxfun(@plus,s_psi_update,bsxfun(@times,Params.epsil,Params.mu));
    s_psi_update = bsxfun(@times,s_psi_update,1./sqrt(sum((s_psi_update).^2)));
   
    for s = 1:setting_params.num_sub
        checkpsi = diag(s_psi_update(:,:,s)'*Params.s_psi(:,:,s));
        checkpsi_flag = (sum(1-checkpsi < setting_params.epsilon) < setting_params.num_clusters);
        if(checkpsi_flag < 1)
            flag_psi(s,1) = 1;
        end
    end
    Params.s_psi = s_psi_update;
    
    % update sigma
    sigma_update = bsxfun(@times,Params.s_psi,permute(Params.s_t_nu,[1,2,4,3]));
    sigma_update = CBIG_nanmean(mean(sum(sigma_update,1),3),4);%1xLxSxT=>1xL
    for i = 1:setting_params.num_clusters
        sigma_update(i) = invAd(setting_params.dim,sigma_update(i));
    end
    
    if((sum(flag_psi) == setting_params.num_sub) && (mean(abs(Params.sigma-sigma_update)./Params.sigma) ...
       < setting_params.epsilon))
        stop_intra = 1;
    end
    Params.sigma = sigma_update;
end
end

function Params = vmf_clustering_subject_session(Params,setting_params,data)

% Inter-region level

stop_em = 0;
iter_em = 0;
cost = zeros(1, setting_params.num_sub);

while(stop_em == 0)
    iter_em = iter_em + 1;
    fprintf('It is EM iteration.. %d..\n',iter_em);
    %% Mstep
    flag_nu = zeros(setting_params.num_sub,setting_params.num_session);
    stop_m = 0;
    iter_m = 0;
    fprintf('M-step..\n');
    while(stop_m == 0)
        iter_m = iter_m + 1;
        
        % update kappa
        s_lambda = Params.s_lambda;%NxLxS
        kappa_update = mtimesx(data.series,permute(Params.s_t_nu,[1,2,4,3]));
        kappa_update = bsxfun(@times,s_lambda,kappa_update);
        kappa_update = sum(sum(CBIG_nanmean(sum(kappa_update,1),4),3));
        kappa_update = kappa_update./sum(sum(sum(s_lambda,1),3));
        kappa_update = invAd(setting_params.dim,kappa_update);
        kappa_update = repmat(kappa_update,1,setting_params.num_clusters);

        if (sum(kappa_update == Inf) ~= 0)
            fprintf('[WARNING] kappa is Inf !\n')
            kappa_update(kappa_update == Inf) = Params.kappa(kappa_update == Inf);
        end
        if (sum(kappa_update < setting_params.ini_concentration) ~= 0)
            kappa_update(kappa_update < setting_params.ini_concentration) = setting_params.ini_concentration;
            fprintf('[WARNING] kappa is less than the minimal value %d \n', setting_params.ini_concentration);
        end
       
        for s = 1:setting_params.num_sub
            for t = 1:setting_params.num_session
                
                % update s_t_nu
                checknu = [];
                X = data.series(:,:,s,t);
                s_lambda = Params.s_lambda(:,:,s);
                lambda_X = bsxfun(@times,kappa_update,X'*s_lambda) + bsxfun(@times,Params.sigma,Params.s_psi(:,:,s));
                s_t_nu_update = bsxfun(@times,lambda_X,1./sqrt(sum((lambda_X).^2)));
                checknu = diag(s_t_nu_update'*Params.s_t_nu(:,:,t,s));
                checknu(isnan(checknu)) = 1;
                checknu_flag = (sum(1-checknu < setting_params.epsilon) < setting_params.num_clusters);
                Params.s_t_nu(:,:,t,s) = s_t_nu_update;

                if(checknu_flag < 1)
                    flag_nu(s,t) = 1;
                end
                
            end
        end
        if((sum(sum(flag_nu)) == setting_params.num_sub*setting_params.num_session) ...
           && (mean(abs(Params.kappa-kappa_update)./Params.kappa) < setting_params.epsilon))
            stop_m=1;
        end
        Params.kappa = kappa_update;
    end
    %% Estep
    fprintf('Estep..\n');
    
    % estimate s_lambda
    log_vmf = permute(Params.s_t_nu,[1,2,4,3]);
    log_vmf = mtimesx(data.series,log_vmf);%NxLxSxT
    log_vmf = bsxfun(@times,permute(log_vmf,[2,1,3,4]),transpose(Params.kappa));%LxNxSxT
    log_vmf(:,sum(log_vmf==0,1)==0) = bsxfun(@plus,Cdln(transpose(Params.kappa),setting_params.dim), ...
                                      log_vmf(:,sum(log_vmf==0,1)==0));%LxNxSxT
    log_vmf =   CBIG_nansum(log_vmf,4);%NxLxS
    idx = sum(log_vmf == 0,1) ~= 0;
    
    log_vmf = bsxfun(@plus,permute(log_vmf,[2,1,3]),log(Params.theta));
    s_lambda = bsxfun(@minus,log_vmf,max(log_vmf,[],2));
    s_lambda = exp(s_lambda);
    Params.s_lambda = bsxfun(@times,s_lambda,1./sum(s_lambda,2));
    Params.s_lambda = permute(Params.s_lambda,[2,1,3]);
    Params.s_lambda(:,idx) = 0;
    Params.s_lambda = permute(Params.s_lambda,[2,1,3]);
    
    % estimate theta
    Params.theta = mean(Params.s_lambda,3);

    %% em stop criteria
    for s = 1:setting_params.num_sub    
        for t = 1:setting_params.num_session
            X = data.series(:,:,s,t);
            setting_params.num_verts = size(X,1);
            log_vmf = vmf_probability(X,Params.s_t_nu(:,:,t,s), Params.kappa);
            if(t == 1)
                log_lambda_prop = log_vmf;
            else
                log_lambda_prop = CBIG_nansum(cat(3,log_lambda_prop,log_vmf),3);
            end
        end
        theta_cost = Params.theta;
        log_theta_cost = log(theta_cost);
        log_theta_cost(isinf(log_theta_cost)) = log(eps.^20);

        s_lambda_cost = Params.s_lambda(:,:,s);
        log_s_lambda_cost = log(s_lambda_cost);
        log_s_lambda_cost(isinf(log_s_lambda_cost)) = log(eps.^20);

        update_cost(:,s) = sum(sum(s_lambda_cost.*log_lambda_prop)) ...
                           + sum(sum(s_lambda_cost.*log_theta_cost)) ...
                           - sum(sum(s_lambda_cost.*log_s_lambda_cost));    
    end
    sub_set=find((abs(abs(update_cost-cost)./cost) > setting_params.epsilon) == 0);
    if(length(sub_set) == setting_params.num_sub)
        stop_em = 1;
        Params.cost_em = cost;
    end
    if(iter_em > 100)
        stop_em = 1;
        Params.cost_em = update_cost;
        warning('vem can not converge');
    end
    cost = update_cost;
end
fprintf('EM..Done\n');
end


function log_vmf = vmf_probability(X,nu,kap)
dim = size(X,2) - 1;
log_vmf = bsxfun(@plus,Cdln(kap,dim),bsxfun(@times,kap,X*nu));
end

function out = Ad(in,D)
out = besseli(D/2,in) ./ besseli(D/2-1,in);
end

function out = Cdln(k,d,k0)
k = double(k);
sizek = size(k);
k = k(:);

out = (d/2-1).*log(k)-log(besseli((d/2-1)*ones(size(k)),k));
if(d<1200)
    k0 = 500;
elseif(d>=1200 && d<1800)
    k0 = 650;
else
    error('dimension is too high, need  to specify k0');
end
fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
nGrids = 1000;

maskof = find(k>k0);
nkof = length(maskof);

if nkof > 0

    kof = k(maskof);

    ofintv = (kof - k0)/nGrids;
    tempcnt = (1:nGrids) - 0.5;
    ks = k0 + repmat(tempcnt,nkof,1).*repmat(ofintv,1,nGrids);
    adsum = sum( 1./((0.5*(d-1)./ks) + sqrt(1+(0.5*(d-1)./ks).^2)) ,2);

    out(maskof) =  fk0 - ofintv .* adsum;

end

out = single(reshape(out,sizek));
end

function [outu] = invAd(D,rbar)

rbar = double(rbar);

outu = (D-1).*rbar./(1-rbar.^2) + D/(D-1).*rbar;

[i] = besseli(D/2-1,outu);


if ((i == Inf)||(isnan(i)) || (i==0))
    out = outu - D/(D-1)*rbar/2;
    exitflag = Inf;
else
    [outNew, fval exitflag]  = fzero(@(argum) Ad(argum,D)-rbar,outu);
    if exitflag == 1
        out = outNew;
    else
        out = outu - D/(D-1)*rbar/2;
    end
end
end


function data = fetch_data(project_dir,num_session,site,num_sub,mesh)    

% read in input functional connectivity profiles
load(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib', ...
    'fs_LR_32k_medial_mask.mat'));
for t = 1:num_session
    data_profile = fullfile(project_dir,'profile_list','training_set',['sess' site{t} '.txt']);
    profile_name = table2cell(readtable(data_profile,'Delimiter',' ','ReadVariableNames',false));
    if(t==1)
        if(strcmp(profile_name{1,1},'NONE'))
            error('The first session of the first subject can not be empty');
        end
    end
    for i = 1:num_sub
        avg_file = profile_name{i,1};
        if(strcmp(avg_file,'NONE'))
            data.series(:,:,i,t) = zeros(size(data.series,1),size(data.series,2))*NaN;
        else

            [~, series, ~] = CBIG_MSHBM_read_fmri(avg_file);       
            series(~medial_mask, :) = 0;

            series = bsxfun(@minus,series,mean(series, 2));
            series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), ...
                                         sqrt(sum(series(all(series,2)~=0,:).^2,2)));
            data.series(:,:,i,t) = series;
        end
    end
end 
end
