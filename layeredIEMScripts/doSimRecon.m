function mean_recon = doSimRecon(parcombs,thisdat,nconds,nvox,vrf_w_std,...
    sim_response_std,trn_condList,test_list,trn_all,tst_all,trnX,basis_set,...
    xx,yy)

% set MATLAB to catch singular matrix warnings as errors
warning('error', 'MATLAB:nearlySingularMatrix');
warning('error', 'MATLAB:singularMatrix');

all_recons = []; all_rpos = []; all_rcond = []; all_parcomb = [];
mean_recon = [];

for p = 1:size(parcombs,1)
    vrf_dat = thisdat;
    
    nochange = find(~parcombs(p,:));
    for ni = 1:length(nochange)
        vrf_dat(:,1:2,nochange(ni)) = repmat(vrf_dat(:,3,nochange(ni)),1,2,1);
    end

    %% generate vRF weights
    for cc = 1:nconds
        for vv = 1:nvox
            % baseline + 2dcos(x,y,sz)*amp
            vrfw(cc,vv,:) = vrf_dat(vv,cc,5) + make2dcos(xx,yy,...
                vrf_dat(vv,cc,1),vrf_dat(vv,cc,2),...
                vrf_dat(vv,cc,3),7)*vrf_dat(vv,cc,4);
        end
    end

    % add noise to the weights
    vrfwn = vrfw + randn(size(vrfw))*vrf_w_std;

    %% generate simulated data from weights

    sim_response = nan(size(trn_condList,1),nvox);
    sim_test = nan(size(test_list,1),nvox);

    % response = sum(channel weight .* channel response to the stimulus)
    % X_all is n_trials x n_channels (with matching channels...)

    for cc = 1:nconds
        tridx = trn_condList==cc;
        sim_response(tridx,:) = trn_all(tridx,:) * reshape(vrfwn(cc,:,:),...
            nvox,[])';
        testidx = test_list(:,2)==cc;
        sim_test(testidx,:) = tst_all(testidx,:) * reshape(vrfwn(cc,:,:),...
            nvox,[])';
    end

    sim_response_n = sim_response + randn(size(sim_response))*sim_response_std;
    sim_test_n = sim_test + randn(size(sim_test))*sim_response_std;

    %% compute channel responses

    % make design matrix
    trn = sim_response_n;
    w = trnX\trn;

    for cc = 1:nconds-1
        % apply the EM weight matrix to the test data (L or R)
        tstset = test_list(:,2)==cc;
        tst = sim_test_n(tstset,:);
        try
            x = inv(w*w')*w*tst';
        catch ME
            switch ME.identifier
                case 'MATLAB:nearlySingularMatrix'
                    return;
                case 'MATLAB:singularMatrix'
                    return;
            end
        end
        chan_resp(tstset,:) = x';
        clear x;
    end

    %% weighted sum of basis set = representation

    % basis_set is numChan x res^2
    % chan_resp is nTrials x nChan

    recon_vec = chan_resp*basis_set;

    all_rpos = cat(1,all_rpos,test_list(:,1));
    all_rcond = cat(1,all_rcond,test_list(:,2));
    all_parcomb = cat(1,all_parcomb,repmat(p,size(test_list,1),1));
    all_recons = cat(1,all_recons,recon_vec);

    clear recon_vec w vrfw vrfwn chan_resp x
    clear sim_response sim_response_n sim_test sim_test_n
    clear vrf_dat
end
clear subdat stim tststim trial_list test_list X_all tst_all

%% avg these and calculate error measurements

nstim = length(unique(all_rpos));
clim = [min(all_recons(:)) max(all_recons(:))];

for p = 1:size(parcombs,1)
    for cc = 1:nconds-1
        for ss = 1:nstim
            thisidx = all_rpos==ss & all_rcond==cc & all_parcomb == p;
            mean_recon(p,cc,ss,:) = nanmean(all_recons(thisidx,:));
            clear thisidx;
        end
    end
end