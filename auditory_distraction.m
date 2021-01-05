%
% Code written by Stefan Arnau, July 2020
% Email: arnau@ifado.de
%

% Better...
clear all;

% Set paths
PATH_FIELDTRIP       = '/home/plkn/repos/fieldtrip/';
PATH_EEGLAB          = '/home/plkn/repos/eeglab/';
PATH_IN              = '/mnt/projdat/schroeger_2fac/0_data_2fac/';
PATH_SOMETHING_DONE  = '/mnt/projdat/schroeger_2fac/1_something_has_been_done/';
PATH_PLOT            = '/mnt/projdat/schroeger_2fac/2_plots/';
PATH_VEUSZ           = '/mnt/projdat/schroeger_2fac/3_veusz/';

% Part switch
dostuff = {'thing3'};

% THING 1: prepare the data
if ismember('thing1', dostuff)

    % Initialize eeglab for topoplots
    addpath(PATH_EEGLAB);
    eeglab;
        
    % Read data
    performance_data = readtable([PATH_IN, 'vs_schroeger_behavior.csv']);
    n_subjects = size(performance_data, 1);
    erps = [];
    for f = 1 : n_subjects
        load([PATH_IN, performance_data.VPCode{f}, '_dat.mat']);
        erps(f, 1, 1, :, :) = squeeze(mean(ERP_Dat.data(1, :, :), 1)); % std short
        erps(f, 1, 2, :, :) = squeeze(mean(ERP_Dat.data(2, :, :), 1)); % std long
        erps(f, 2, 1, :, :) = squeeze(mean(ERP_Dat.data(3, :, :), 1)); % dev short
        erps(f, 2, 2, :, :) = squeeze(mean(ERP_Dat.data(4, :, :), 1)); % dev long
    end

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
    for x1 = 1 : size(erps, 1)
        for x2 = 1 : size(erps, 2)
            for x3 = 1 : size(erps, 3)
                ref_val = squeeze(mean(erps(x1, x2, x3, [10, 21], :), 4));
                for x4 = 1 : size(erps, 4)
                    erps(x1, x2, x3, x4, :) = squeeze(erps(x1, x2, x3, x4, :)) - ref_val;
                end
            end
        end
    end
    
    % Remove mastoid channels
    erps(:, :, :, [10, 21], :) = [];
    VP_Header.chanlocs([10, 21]) = [];

    % Adjust new order
    new_order(new_order > 21) = new_order(new_order > 21) - 1;
    new_order(new_order > 10) = new_order(new_order > 10) - 1;

    % Re-order channels
    erps = erps(:, :, :, new_order, :);
    chanlocs = VP_Header.chanlocs;
    for ch = 1 : numel(VP_Header.chanlocs)
        chanlocs(ch) = VP_Header.chanlocs(new_order(ch));
    end

    % Save erps
    save([PATH_SOMETHING_DONE 'erps'], 'erps');

    % Save chanloc struct
    save([PATH_SOMETHING_DONE 'chanlocs'], 'chanlocs');

    % Init fieldtrip
    addpath(PATH_FIELDTRIP);
    ft_defaults;

    % Build chanlocs
    chanlabs = {};
    coords = [];
    for c = 1 : numel(chanlocs)
        chanlabs{c} = chanlocs(c).labels;
        coords(c, :) = [chanlocs(c).X, chanlocs(c).Y, chanlocs(c).Z];
    end

    % Save chanlabel struct
    save([PATH_SOMETHING_DONE 'chanlabs'], 'chanlabs');

    % A sensor struct
    sensors = struct();
    sensors.label = chanlabs;
    sensors.chanpos = coords;
    sensors.elecpos = coords;

    % Prepare neighbor struct
    cfg                 = [];
    cfg.elec            = sensors;
    cfg.feedback        = 'no';
    cfg.method          = 'triangulation'; 
    neighbours          = ft_prepare_neighbours(cfg);

    % Save neighbor struct
    save([PATH_SOMETHING_DONE 'neighbours'], 'neighbours');

    % Save erp time
    erp_time = VP_Header.erp_time;
    save([PATH_SOMETHING_DONE 'erp_time'], 'erp_time');

end % End part1

% THING 2: Test conditions
if ismember('thing2', dostuff)

    % Load
    load([PATH_SOMETHING_DONE 'erps']); % erps
    load([PATH_SOMETHING_DONE 'erp_time']); % erp_time
    load([PATH_SOMETHING_DONE 'neighbours']); % neighbours
    load([PATH_SOMETHING_DONE 'chanlabs']); % chanlabs
    load([PATH_SOMETHING_DONE 'chanlocs']); % chanlocs

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
    prune_idx = erp_time >= -200 & erp_time <= 1000;
    d.time = erp_time(prune_idx);

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
	dlmwrite([PATH_VEUSZ, 'erp_time.csv'], erp_time(prune_idx));

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


% THING 3: Correlations
if ismember('thing3', dostuff)

    % Read data
    performance_data = readtable([PATH_IN, 'vs_schroeger_behavior.csv']);
	
	% Load
    load([PATH_SOMETHING_DONE 'erps']); % erps
    load([PATH_SOMETHING_DONE 'erp_time']); % erp_time
    load([PATH_SOMETHING_DONE 'neighbours']); % neighbours
    load([PATH_SOMETHING_DONE 'chanlabs']); % chanlabs
    load([PATH_SOMETHING_DONE 'chanlocs']); % chanlocs

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
    prune_idx = erp_time >= -200 & erp_time <= 1000;
    d.time = erp_time(prune_idx);

    % Test interaction based on length differences
    D = {};
    for s = 1 : size(erps, 1)
        long_minus_short_in_std = squeeze(erps(s, 1, 2, :, prune_idx)) - squeeze(erps(s, 1, 1, :, prune_idx));
        long_minus_short_in_dev = squeeze(erps(s, 2, 2, :, prune_idx)) - squeeze(erps(s, 2, 1, :, prune_idx));
        d.avg = long_minus_short_in_dev - long_minus_short_in_std; % Difference of difference
        D{s} = d;
    end 
    GA_diff_std = ft_timelockgrandaverage(cfg, D{1, :});

    % Cluster statistics configuration
    cfg = [];
    cfg.tail             = 0; 
    cfg.statistic        = 'ft_statfun_correlationT'; 
    cfg.alpha            = 0.025;
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2;
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = 0.01;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = 1000;
    cfg.computecritval   = 'yes'; 
    cfg.ivar             = 1;

    aa=bb;
    % Correlation age
    performance = performance_data
    cfg.design = ages;
    [stat] = ft_timelockstatistics(cfg, GA_diff_std);
    plot_correlation_clusters(stat, 0.025, 'correlation-lengthdiff-in-std-x-age', [-0.5, 0.5], PATH_VEUSZ, PATH_PLOT, chanlabs, chanlocs); 

	
	
	
	
	
	
	
	
	
	
	
	
	

end % End part3















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


