%
% Code written by Stefan Arnau, July 2020
% Email: arnau@ifado.de
%

% Better...
clear all;

% Set paths
PATH_FIELDTRIP       = '/home/plkn/fieldtrip-master/';
PATH_EEGLAB          = '/home/plkn/eeglab2021.1/';
PATH_IN              = '/mnt/data_fast/schroeger_2fac/0_resplocked_data/';
PATH_SOMETHING_DONE  = '/mnt/data_fast/schroeger_2fac/1_something_has_been_done_resplocked/';
PATH_PLOT            = '/mnt/data_fast/schroeger_2fac/2_plots_resplocked/';
PATH_VEUSZ           = '/mnt/data_fast/schroeger_2fac/3_veusz_resplocked/';

% Part switch
dostuff = {'thing1'};

% THING 1: prepare the data
if ismember('thing1', dostuff)

    % Initialize eeglab for topoplots
    addpath(PATH_EEGLAB);
    eeglab;
        
    % Read data
    performance_data = readtable([PATH_IN, 'vs_schroeger_behavior.csv']);
    n_subjects = size(performance_data, 1);
    erps_resplocked = [];
    for f = 1 : n_subjects
        load([PATH_IN, performance_data.VPCode{f}, '_resplocked.mat']);
        erps_resplocked(f, 1, 1, :, :) = squeeze(erps(1, :, :)); % std short
        erps_resplocked(f, 1, 2, :, :) = squeeze(erps(2, :, :)); % std long
        erps_resplocked(f, 2, 1, :, :) = squeeze(erps(3, :, :)); % dev short
        erps_resplocked(f, 2, 2, :, :) = squeeze(erps(4, :, :)); % dev long
    end

    aa=bb;

    subj=1;
    fprintf('%s', performance_data.VPCode{subj})
    figure
    plot(squeeze(erps_resplocked(subj, 1, 1, 13, :)));

    % New electrode order
    new_order = [1, 29, 30,...          % FP1 FPz Fp2 
                 3, 4, 31, 27, 28,...   % F7 F3 Fz F4 F8
                 2, 5, 26,...           % FC3 FCz FC4 
                 7, 8, 32, 23, 24,...   % T7 C3 Cz C4 T8
                 6, 9, 25,...           % CP3 CPz CP4
                 11, 12, 13, 19, 20,... % P7 P3 Pz P4 P8
                 14, 22, 18,...         % PO3 POz PO4
                 15, 16, 17];           % O1 Oz O2

    % Link the mastoids
    for x1 = 1 : size(erps_resplocked, 1)
        for x2 = 1 : size(erps_resplocked, 2)
            for x3 = 1 : size(erps_resplocked, 3)
                ref_val = squeeze(mean(erps_resplocked(x1, x2, x3, [10, 21], :), 4)); % Re-reference to linked mastoids
                for x4 = 1 : size(erps_resplocked, 4)
                    erps_resplocked(x1, x2, x3, x4, :) = squeeze(erps_resplocked(x1, x2, x3, x4, :)) - ref_val;
                end
            end
        end
    end
    
    % Remove mastoid channels
    erps_resplocked(:, :, :, [10, 21], :) = [];

    % Adjust new order
    new_order(new_order > 21) = new_order(new_order > 21) - 1;
    new_order(new_order > 10) = new_order(new_order > 10) - 1;

    % Re-order channels
    erps_resplocked = erps_resplocked(:, :, :, new_order, :);

    % Save erps_resplocked
    save([PATH_SOMETHING_DONE 'erps_resplocked'], 'erps_resplocked');

    % Save erp time
    save([PATH_SOMETHING_DONE 'erp_times'], 'erp_times');

end % End part1

% THING 2: Test conditions
if ismember('thing2', dostuff)

    % Load
    load([PATH_SOMETHING_DONE 'erps_resplocked']); % erps
    load([PATH_SOMETHING_DONE 'erp_times']); % erp_times
    load([PATH_SOMETHING_DONE 'neighbours']); % neighbours
    load([PATH_SOMETHING_DONE 'chanlabs']); % chanlabs
    load([PATH_SOMETHING_DONE 'chanlocs']); % chanlocs

    erps = erps_resplocked;

    % Init fieldtrip
    addpath(PATH_FIELDTRIP);
    ft_defaults;

    % Init eeglab
    addpath(PATH_EEGLAB);
    eeglab;

    % A template for GA structs
    cfg=[];
    cfg.keepindividual = 'yes';
    d = [];
    d.dimord = 'chan_time';
    d.label = chanlabs;
    prune_idx = erp_times >= -700 & erp_times <= 400;
    d.time = erp_times(prune_idx);

    % erps dimord: subject x oddball x length x channel x time

	% Average erps for conditions at FCz and Pz
	idx_fcz = find(strcmp({chanlocs.labels}, 'FCz'));
	idx_pz = find(strcmp({chanlocs.labels}, 'Pz'));
	erps_fcz = [];
	erps_pz = [];
	counter = 0;
	for odd = 1 : size(erps, 2)
        for tone = 1 : size(erps, 3)
			counter = counter + 1;
			erps_fcz(counter, :) = squeeze(mean(squeeze(erps(:, odd, tone, idx_fcz, prune_idx)), 1));
			erps_pz(counter, :) = squeeze(mean(squeeze(erps(:, odd, tone, idx_pz, prune_idx)), 1));
		end
	end
	dlmwrite([PATH_VEUSZ, 'erps_fcz.csv'], erps_fcz);
	dlmwrite([PATH_VEUSZ, 'erps_pz.csv'], erps_pz);
	dlmwrite([PATH_VEUSZ, 'erp_times.csv'], erp_times(prune_idx));

    % Tone duration difference wave
    duration_fcz_long = squeeze(mean(squeeze(erps(:, :, 2, idx_fcz, prune_idx)), 2));
    duration_fcz_short = squeeze(mean(squeeze(erps(:, :, 1, idx_fcz, prune_idx)), 2));
    duration_fcz_diff = duration_fcz_long - duration_fcz_short;
    diff_erp_fcz_duration = [mean(duration_fcz_long)', std(duration_fcz_long)', mean(duration_fcz_short)', std(duration_fcz_short)', mean(duration_fcz_diff)', std(duration_fcz_diff)'];
    dlmwrite([PATH_VEUSZ, 'diff_erp_fcz_duration.csv'], diff_erp_fcz_duration, 'delimiter', '\t');

    duration_pz_long = squeeze(mean(squeeze(erps(:, :, 2, idx_pz, prune_idx)), 2));
    duration_pz_short = squeeze(mean(squeeze(erps(:, :, 1, idx_pz, prune_idx)), 2));
    duration_pz_diff = duration_pz_long - duration_pz_short;
    diff_erp_pz_duration = [mean(duration_pz_long)', std(duration_pz_long)', mean(duration_pz_short)', std(duration_pz_short)', mean(duration_pz_diff)', std(duration_pz_diff)'];
    dlmwrite([PATH_VEUSZ, 'diff_erp_pz_duration.csv'], diff_erp_pz_duration, 'delimiter', '\t');

    % Deviant standard difference wave
    oddball_fcz_deviant = squeeze(mean(squeeze(erps(:, 2, :, idx_fcz, prune_idx)), 2));
    oddball_fcz_standard = squeeze(mean(squeeze(erps(:, 1, :, idx_fcz, prune_idx)), 2));
    oddball_fcz_diff = oddball_fcz_deviant - oddball_fcz_standard;
    diff_erp_fcz_oddball = [mean(oddball_fcz_deviant)', std(oddball_fcz_deviant)', mean(oddball_fcz_standard)', std(oddball_fcz_standard)', mean(oddball_fcz_diff)', std(oddball_fcz_diff)'];
    dlmwrite([PATH_VEUSZ, 'diff_erp_fcz_oddball.csv'], diff_erp_fcz_oddball, 'delimiter', '\t');

    oddball_pz_deviant = squeeze(mean(squeeze(erps(:, 2, :, idx_pz, prune_idx)), 2));
    oddball_pz_standard = squeeze(mean(squeeze(erps(:, 1, :, idx_pz, prune_idx)), 2));
    oddball_pz_diff = oddball_pz_deviant - oddball_pz_standard;
    diff_erp_pz_oddball = [mean(oddball_pz_deviant)', std(oddball_pz_deviant)', mean(oddball_pz_standard)', std(oddball_pz_standard)', mean(oddball_pz_diff)', std(oddball_pz_diff)'];
    dlmwrite([PATH_VEUSZ, 'diff_erp_pz_oddball.csv'], diff_erp_pz_oddball, 'delimiter', '\t');
    
    % Interaction fcz
    fcz_standard_long = squeeze(squeeze(erps(:, 1, 2, idx_fcz, prune_idx)));
    fcz_standard_short = squeeze(squeeze(erps(:, 1, 1, idx_fcz, prune_idx)));
    fcz_deviant_long = squeeze(squeeze(erps(:, 2, 2, idx_fcz, prune_idx)));
    fcz_deviant_short = squeeze(squeeze(erps(:, 2, 1, idx_fcz, prune_idx)));

    fcz_duration_in_dev = fcz_deviant_long - fcz_deviant_short;
    fcz_duration_in_std = fcz_standard_long - fcz_standard_short;
    fcz_oddball_in_long = fcz_deviant_long - fcz_standard_long;
    fcz_oddball_in_short = fcz_deviant_short - fcz_standard_short;

    fcz_duration_diff = fcz_duration_in_dev - fcz_duration_in_std;
    fcz_oddball_diff = fcz_oddball_in_long - fcz_oddball_in_short;

    interaction_fcz = [mean(fcz_standard_short)', std(fcz_standard_short)',...
                       mean(fcz_standard_long)', std(fcz_standard_long)',...
                       mean(fcz_deviant_short)', std(fcz_deviant_short)',...
                       mean(fcz_deviant_long)', std(fcz_deviant_long)',...
                       mean(fcz_duration_in_dev)', std(fcz_duration_in_dev)',...
                       mean(fcz_duration_in_std)', std(fcz_duration_in_std)',...
                       mean(fcz_oddball_in_long)', std(fcz_oddball_in_long)',...
                       mean(fcz_oddball_in_short)', std(fcz_oddball_in_short)',...
                       mean(fcz_duration_diff)', std(fcz_duration_diff)',...
                       mean(fcz_oddball_diff)', std(fcz_oddball_diff)',...
                       ];

    dlmwrite([PATH_VEUSZ, 'interaction_fcz.csv'], interaction_fcz, 'delimiter', '\t');

    % Interaction pz
    pz_standard_long = squeeze(squeeze(erps(:, 1, 2, idx_pz, prune_idx)));
    pz_standard_short = squeeze(squeeze(erps(:, 1, 1, idx_pz, prune_idx)));
    pz_deviant_long = squeeze(squeeze(erps(:, 2, 2, idx_pz, prune_idx)));
    pz_deviant_short = squeeze(squeeze(erps(:, 2, 1, idx_pz, prune_idx)));

    pz_duration_in_dev = pz_deviant_long - pz_deviant_short;
    pz_duration_in_std = pz_standard_long - pz_standard_short;
    pz_oddball_in_long = pz_deviant_long - pz_standard_long;
    pz_oddball_in_short = pz_deviant_short - pz_standard_short;

    pz_duration_diff = pz_duration_in_dev - pz_duration_in_std;
    pz_oddball_diff = pz_oddball_in_long - pz_oddball_in_short;

    interaction_pz = [mean(pz_standard_short)', std(pz_standard_short)',...
                       mean(pz_standard_long)', std(pz_standard_long)',...
                       mean(pz_deviant_short)', std(pz_deviant_short)',...
                       mean(pz_deviant_long)', std(pz_deviant_long)',...
                       mean(pz_duration_in_dev)', std(pz_duration_in_dev)',...
                       mean(pz_duration_in_std)', std(pz_duration_in_std)',...
                       mean(pz_oddball_in_long)', std(pz_oddball_in_long)',...
                       mean(pz_oddball_in_short)', std(pz_oddball_in_short)',...
                       mean(pz_duration_diff)', std(pz_duration_diff)',...
                       mean(pz_oddball_diff)', std(pz_oddball_diff)',...
                       ];

    dlmwrite([PATH_VEUSZ, 'interaction_pz.csv'], interaction_pz, 'delimiter', '\t');

    % Test long versus short trials
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(mean(squeeze(erps(s, :, 1, :, prune_idx)), 1));
        D{s} = d;
    end 
    GA_short = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(mean(squeeze(erps(s, :, 2, :, prune_idx)), 1));
        D{s} = d;
    end 
    GA_long = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('long-vs-short', GA_long, GA_short, 'long', 'short', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test standard versus deviant trials
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(mean(squeeze(erps(s, 1, :, :, prune_idx)), 1));
        D{s} = d;
    end 
    GA_std = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(mean(squeeze(erps(s, 2, :, :, prune_idx)), 1));
        D{s} = d;
    end 
    GA_dev = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('dev-vs-std', GA_dev, GA_std, 'dev', 'std', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test deviant versus standard short trials only
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 1, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_std_short = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_dev_short = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('dev-vs-std-in-short', GA_dev_short, GA_std_short, 'deviant', 'standard', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test deviant versus standard long trials only
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 1, 2, :, prune_idx));
        D{s} = d;
    end 
    GA_std_long = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 2, :, prune_idx));
        D{s} = d;
    end 
    GA_dev_long = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('dev-vs-std-in-long', GA_dev_long, GA_std_long, 'deviant', 'standard', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test short versus long in std only
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 1, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_short_std = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 1, 2, :, prune_idx));
        D{s} = d;
    end 
    GA_long_std = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('long-vs-short-in-std', GA_long_std, GA_short_std, 'long', 'short', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test short versus long in dev only
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_short_dev = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 2, :, prune_idx));
        D{s} = d;
    end 
    GA_long_dev = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('long-vs-short-in-dev', GA_long_dev, GA_short_dev, 'long', 'short', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test interaction based on oddball differences
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 1, :, prune_idx)) - squeeze(erps(s, 1, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_diff_short = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 2, :, prune_idx)) - squeeze(erps(s, 1, 2, :, prune_idx));
        D{s} = d;
    end 
    GA_diff_long = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('interaction-oddball-diff', GA_diff_long, GA_diff_short, 'diff-in-long', 'diff-in-short', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

    % Test interaction based on length differences
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 1, 2, :, prune_idx)) - squeeze(erps(s, 1, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_diff_std = ft_timelockgrandaverage(cfg, D{1, :});
    D = {};
    for s = 1 : size(erps, 1)
        d.avg = squeeze(erps(s, 2, 2, :, prune_idx)) - squeeze(erps(s, 2, 1, :, prune_idx));
        D{s} = d;
    end 
    GA_diff_dev = ft_timelockgrandaverage(cfg, D{1, :});
    a_cluster_test('interaction-length-diff', GA_diff_dev, GA_diff_std, 'diff-in-dev', 'diff-in-std', chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ);

end % End part2















% Helper functions
function[outerp] = znorm_erp(inerp, tvec)
    m = mean(inerp(tvec < 0));
    s = std(inerp(tvec < 0));
    outerp = (inerp - m) / s;
end

function[outerp] = vecnorm_erp(inerp)
    outerp = inerp / sqrt(sum(inerp .^ 2));
end

% Function that plots correlation clusters
function[] =  plot_correlation_clusters(stat, thresh, title_string, clim, PATH_VEUSZ, PATH_PLOT, chanlabs, chanlocs)

    mkdir([PATH_PLOT, title_string])
    PATH_OUTPUT = [PATH_PLOT, title_string, '/'];

    try
        sig_pos = find([stat.posclusters.prob] <= thresh);
    catch
        sig_pos = [];
    end
    try
        sig_neg = find([stat.negclusters.prob] <= thresh);
    catch
        sig_neg = [];
    end

    if sig_pos
        for cl = 1 : length(sig_pos)
            idx = stat.posclusterslabelmat == sig_pos(cl);
            pval = round(stat.posclusters(sig_pos(cl)).prob, 3);
            dlmwrite([PATH_VEUSZ, title_string, '_poscluster_', num2str(sig_pos(cl)), '.csv'], idx);
            chans_sig = find(sum(idx, 2));
            times_sig = find(sum(idx, 1));
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            figure('Visible', 'off'); clf;
            subplot(2, 1, 1)
            pd = stat.rho;
            contourf(stat.time, [1 : numel(chanlabs)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(chanlabs)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;
            subplot(2, 1, 2)
            pd = stat.rho;
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['poscluster #' num2str(sig_pos(cl)) ' - p=' num2str(pval)])
            saveas(gcf, [PATH_OUTPUT, title_string, '_poscluster_', num2str(sig_pos(cl)), '.png']);
        end
    end
    if sig_neg
        for cl = 1 : length(sig_neg)
            idx = stat.negclusterslabelmat == sig_neg(cl);
            pval = round(stat.negclusters(sig_neg(cl)).prob, 3);
            dlmwrite([PATH_VEUSZ, title_string, '_negcluster_', num2str(sig_neg(cl)), '.csv'], idx);
            chans_sig = find(sum(idx, 2));
            times_sig = find(sum(idx, 1));
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            figure('Visible', 'off'); clf;
            subplot(2, 1, 1)
            pd = stat.rho;
            contourf(stat.time, [1 : numel(chanlabs)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(chanlabs)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;        
            subplot(2, 1, 2)
            pd = stat.rho;
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['negcluster #' num2str(sig_neg(cl)) ' - p=' num2str(pval)])
            saveas(gcf, [PATH_OUTPUT, title_string, '_negcluster_', num2str(sig_neg(cl)), '.png']);
        end
    end
end

% Function that performs test and creates cluster-plots
function[] = a_cluster_test(titlestring, cond1, cond2, cond1string, cond2string, chanlocs, neighbours, PATH_PLOT, PATH_VEUSZ)

    % Create output directory
    mkdir([PATH_PLOT, titlestring])
    PATH_OUTPUT = [PATH_PLOT, titlestring, '/'];

    % Testparams
    testalpha  = 0.025;
    voxelalpha  = 0.01;
    nperm = 1000;

    % Set config
    cfg = [];
    cfg.tail             = 0;
    cfg.statistic        = 'depsamplesT';
    cfg.alpha            = testalpha;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = voxelalpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = nperm;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;
    cfg.uvar             = 2;
    cfg.design           = [ones(1, size(cond1.individual, 1)), 2 * ones(1, size(cond1.individual, 1)); 1 : size(cond1.individual, 1), 1 : size(cond1.individual, 1)];
    
    % The test
    [stat] = ft_timelockstatistics(cfg, cond1, cond2);  

    % Calculate effect sizes
    apes = [];
    n_subjects = size(cond1.individual, 1);
    for ch = 1 : numel(chanlocs)
        petasq = (squeeze(stat.stat(ch, :)) .^ 2) ./ ((squeeze(stat.stat(ch, :)) .^ 2) + (n_subjects - 1));
        apes(ch, :) = petasq - (1 - petasq) .* (1 / (n_subjects - 1));
    end

    % Plot effect sizes
    cmap = 'jet';
    clim = [-0.6, 0.6];
    figure('Visible', 'off'); clf;
    contourf(stat.time, [1 : numel(cond1.label)], apes, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', [clim])
    colorbar;
    title(['effect sizes: ', titlestring])
    saveas(gcf, [PATH_OUTPUT, 'effect_sizes_', titlestring, '.png']);

    % Plot conditions
    pd1 = squeeze(mean(cond1.individual, 1));
    pd2 = squeeze(mean(cond2.individual, 1));
    cmap = 'jet';
    lim = max(abs([pd1(:); pd2(:)]));
    clim = [-lim, lim];
    figure('Visible', 'off'); clf;
    subplot(1, 2, 1)
    contourf(stat.time, [1 : numel(cond1.label)], pd1, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', [clim])
    title(cond1string)
    subplot(1, 2, 2)
    contourf(stat.time, [1 : numel(cond1.label)], pd2, 40, 'linecolor','none')
    colormap(cmap)
    set(gca, 'clim', [clim])
    title(cond2string)
    saveas(gcf, [PATH_OUTPUT, 'ersp_', titlestring, '.png']);

    % Save effect sizes
    dlmwrite([PATH_VEUSZ, titlestring, '_effect_sizes.csv'], apes);

    % Save averages
    dlmwrite([PATH_VEUSZ, titlestring, '_' cond1string '_average.csv'], squeeze(mean(cond1.individual, 1)));
    dlmwrite([PATH_VEUSZ, titlestring, '_' cond2string '_average.csv'], squeeze(mean(cond2.individual, 1)));
    dlmwrite([PATH_VEUSZ, titlestring, '_difference.csv'], squeeze(mean(cond1.individual, 1)) - squeeze(mean(cond2.individual, 1)));

    % Set threshold to 0.025
    sig_pos = find([stat.posclusters.prob] <= testalpha);
    sig_neg = find([stat.negclusters.prob] <= testalpha);

    % Plot clusters
    if sig_pos
        for cl = 1 : length(sig_pos)

            % Indices of the cluster
            idx = stat.posclusterslabelmat == sig_pos(cl);
            pval = round(stat.posclusters(sig_pos(cl)).prob, 3);

            % Save contour of cluster
            dlmwrite([PATH_VEUSZ, titlestring, '_contour_poscluster_', num2str(sig_pos(cl)), '.csv'], idx);

            % Identify significant channels and time points
            chans_sig = find(sum(idx, 2));
            times_sig = find(sum(idx, 1));

            % Plot a topo of effect sizes
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            clim = [-0.5, 0.5];
            figure('Visible', 'off'); clf;
            pd = mean(apes(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            saveas(gcf, [PATH_OUTPUT, [titlestring, '_effectsize_topo_pos_'], num2str(sig_pos(cl)), '.png']);

            % Plot
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            clim = [-5, 5];
            figure('Visible', 'off'); clf;

            subplot(2, 2, 1)
            pd = squeeze(mean(cond1.individual, 1));
            contourf(stat.time, [1 : numel(cond1.label)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(cond1.label)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;
            title(cond1string)

            subplot(2, 2, 2)
            pd = squeeze(mean(cond2.individual, 1));
            contourf(stat.time, [1 : numel(cond1.label)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(cond1.label)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;
            title(cond2string)
            
            subplot(2, 2, 3)
            pd = squeeze(mean(cond1.individual, 1));
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['poscluster #' num2str(sig_pos(cl)) ' - p=' num2str(pval)])

            subplot(2, 2, 4)
            pd = squeeze(mean(cond2.individual, 1));
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['poscluster #' num2str(sig_pos(cl)) ' - p=' num2str(pval)])

            saveas(gcf, [PATH_OUTPUT, [titlestring, '_pos_'], num2str(sig_pos(cl)), '.png']);
        end
    end
    if sig_neg
        for cl = 1 : length(sig_neg)

            % Indices of the cluster
            idx = stat.negclusterslabelmat == sig_neg(cl);
            pval = round(stat.negclusters(sig_neg(cl)).prob, 3);

            % Save contour of cluster
            dlmwrite([PATH_VEUSZ, titlestring, '_contour_negcluster_', num2str(sig_neg(cl)), '.csv'], idx);
            
            % Identify significant channels and time points
            chans_sig = find(sum(idx, 2));
            times_sig = find(sum(idx, 1));

            % Plot a topo of effect sizes
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            clim = [-0.5, 0.5];
            figure('Visible', 'off'); clf;
            pd = mean(apes(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            saveas(gcf, [PATH_OUTPUT, [titlestring, '_effectsize_topo_neg_'], num2str(sig_neg(cl)), '.png']);

            % Plot
            markercolor = 'k';
            markersize = 10;
            cmap = 'jet';
            clim = [-5, 5];
            figure('Visible', 'off'); clf;

            subplot(2, 2, 1)
            pd = squeeze(mean(cond1.individual, 1));
            contourf(stat.time, [1 : numel(cond1.label)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(cond1.label)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;
            title(cond1string)

            subplot(2, 2, 2)
            pd = squeeze(mean(cond2.individual, 1));
            contourf(stat.time, [1 : numel(cond1.label)], pd, 40, 'linecolor','none')
            hold on
            contour(stat.time, [1 : numel(cond1.label)], idx, 1, 'linecolor', 'k', 'LineWidth', 2)
            colormap(cmap)
            set(gca, 'clim', [clim])
            colorbar;
            title(cond2string)
            
            subplot(2, 2, 3)
            pd = squeeze(mean(cond1.individual, 1));
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['negcluster #' num2str(sig_neg(cl)) ' - p=' num2str(pval)])

            subplot(2, 2, 4)
            pd = squeeze(mean(cond2.individual, 1));
            pd = mean(pd(:, times_sig), 2);
            topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, markersize, 1});
            colormap(cmap);
            caxis(clim);
            title(['negcluster #' num2str(sig_neg(cl)) ' - p=' num2str(pval)])

            saveas(gcf, [PATH_OUTPUT, [titlestring, '_neg_'], num2str(sig_neg(cl)), '.png']);
        end
    end
end % End function


