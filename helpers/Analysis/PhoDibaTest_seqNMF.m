% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded
% Stage 3 of the processing pipeline

addpath(genpath('../../helpers'));
addpath(genpath('../../libraries/buzcode/'));

clear temp

if ~exist('data_config','var')
    Config;
end


if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading data from %s...\n', data_config.output.intermediate_file_paths{2});
    load(data_config.output.intermediate_file_paths{2}, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end

if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'results_array', 'general_results');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('PhoDibaTest_seqNMF ready to process!\n');

plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;


current_binning_index = 1;
active_binning_resolution = processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = timesteps_array{current_binning_index};
temp.curr_processed = active_processing.processed_array{current_binning_index};
active_results = results_array{current_binning_index};

fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);


% Binned Spike rates per time:
[PhoDibaTest_seqNMF_temp.activeMatrix] = fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts);  % 35351x126 double

PhoDibaTest_seqNMF_config.training_subset_start_stop_seconds = [11512.6074066973, 12000];
PhoDibaTest_seqNMF_config.training_subset_start_stop_bins = PhoDibaTest_seqNMF_config.training_subset_start_stop_seconds ./ active_binning_resolution;
PhoDibaTest_seqNMF_config.training_subset_start_stop_bins = [floor(PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(1)), ...
    ceil(PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(2))];

if processing_config.showOnlyAlwaysStableCells
    isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix(:, isAlwaysStable);
end

PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix'; % should be units x time
% Extract only the portion specified as the training portion
PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix(:, PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(1):PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(2));


%% Procedure for choosing K (Number of factors)
PhoDibaTest_seqNMF_config.L = 100; % L: Length (timebins) of each factor exemplar, units of seconds

Ws = {};
Hs = {};
numfits = 3; %number of fits to compare
for k = 1:10
    display(sprintf('running seqNMF with K = %i',k))
    for ii = 1:numfits
        [Ws{ii,k},Hs{ii,k}] = seqNMF(PhoDibaTest_seqNMF_temp.activeMatrix,'K',k, 'L', PhoDibaTest_seqNMF_config.L,'lambda', 0,'maxiter',100,'showplot',0); 
        % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    end
    inds = nchoosek(1:numfits,2);
    parfor i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    end
    
end
%% Plot Diss and choose K with the minimum average diss.
figure,
plot(1:10,Diss,'ko'), hold on
h1 = plot(1:10,median(Diss,1),'k-','linewidth',2);
h2 = plot([3,3],[0,0.5],'r--');
legend([h1 h2], {'median Diss','true K'})
xlabel('K')
ylabel('Diss')






%% Main Procedure:
PhoDibaTest_seqNMF_config.K = 5; % K: Number of factors
PhoDibaTest_seqNMF_config.L = 50; % L: Length (timebins) of each factor exemplar, units of seconds
PhoDibaTest_seqNMF_config.lambda = .001; % lambda: Regularization parameter
shg; clf
display('Running seqNMF on simulated data (2 simulated sequences + noise)')
[active_results.seqNMF.W, active_results.seqNMF.H] = seqNMF(PhoDibaTest_seqNMF_temp.activeMatrix,'K',PhoDibaTest_seqNMF_config.K, 'L', PhoDibaTest_seqNMF_config.L,'lambda', PhoDibaTest_seqNMF_config.lambda);

%% Look at factors
figure; SimpleWHPlot(active_results.seqNMF.W, active_results.seqNMF.H); title('SeqNMF reconstruction')
figure; SimpleWHPlot(active_results.seqNMF.W, active_results.seqNMF.H, PhoDibaTest_seqNMF_temp.activeMatrix); title('SeqNMF factors, with raw data')




% [W,H,cost,loadings,power] = seqNMF(active_binned_spike_data_matrix,'K',K,'L',L,'lambda',0.01);




fprintf('PhoDibaTest_seqNMF complete!\n');
