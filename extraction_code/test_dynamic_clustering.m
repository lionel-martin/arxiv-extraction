clear all;
attempt_name = 'attempt6';
%load('../matgraphs/new_graphs/graphs_th20_100_final.mat');
load(sprintf('../matgraphs/%s/authors_main_grand_topic.mat', attempt_name));
load(sprintf('../matgraphs/%s/presence_th20_100_connected.mat', attempt_name));

start_point = 1;
nb_cols = size(pres_final_gt, 2);

k = 5;
p = 0.15;
N = 366572;

redo = false;
plots = true;

eigs_opts = struct('issym', true);
parpool;

%% unique code

if exist('gts') == 0
    classes_full = zeros(size(pres_final_gt, 1), 1);
    classes_full(authors_grand_topic(:, 1)) = authors_grand_topic(:, 2);
    gts = repmat(classes_full, 1, size(pres_final_gt, 2)) .* pres_final_gt;
end

% for ii = 1:nb_cols
%     G = graphz{ii};
%     save(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii), 'G');
%     graphz{ii} = [];
% end

% clear graphz;


for ii = start_point:nb_cols
    clear G;
    load(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii));
    fprintf('\nComputations for graph %d, N = %d\n', ii, G.N);

%     rm_idx = find(G.d>1e3);
%     sub_idx = setdiff(1:G.N, rm_idx);
%     G_tmp = gsp_graph(G.W(sub_idx, sub_idx));

%     comps_id = gsp_components_v2(G_tmp, bfs_par);
%     [sizes, vals] = hist(comps_id, unique(comps_id));
%     [components(ii), comp_arg] = max(sizes);
%     comp_id = vals(comp_arg);
%     
%     if G_tmp.N == sum(pres_final_gt_rem(:, ii))
%         pres_idx = find(pres_final_gt_rem(:, ii) == 1);
%         prez_idx = pres_idx(comps_id == comp_id);
%         pres_final_gt_rem(setdiff(pres_idx, prez_idx), ii) = 0;
%     end
%     G = gsp_graph(G_tmp.W(comps_id == comp_id, comps_id == comp_id));

%     assert(G.N == sum(pres_final_gt_rem(:, ii)));

%     gts(:, ii) = classes_full .* pres_final_gt(:, ii);

    if (strcmp(G.lap_type, 'combinatorial'))
        fprintf('\tCompute new Laplacian...\n');
        G = gsp_create_laplacian(G, 'normalized');
    end

    if (~isfield(G, 'gt') || redo)
        s = gts(pres_final_gt(:, ii) == 1, ii);
        G.gt = s;
    end

    k_old = 0;
    if isfield(G, 'k')
        k_old = G.k;
    end

    G.k = k;
    if G.N > k+5
        if (~isfield(G, 'U') || k ~= k_old || redo)
            fprintf('\tCompute Fourier basis...\n');
            [U, e] = eigs(G.L, G.k+5, 'sm', eigs_opts);
            e = diag(e);
            [G.e, e_idx] = sort(e);
            G.U = U(:, e_idx);
            G.Uk = G.U(:, 1:G.k);

            t_eigs = tic;
            [U, e] = eigs(G.L, G.k, 'sm', eigs_opts);
            sort(diag(e)); % only for timing
            G.time.SC.eigs = toc(t_eigs);

            fprintf('\tCompute Spectral Clustering...\n');
            t_SC = tic;
            [G.SC.ass, ~, kmcost] = kmeans(G.Uk, G.k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
            G.SC.kmcost = sum(kmcost);
            G.time.SC.kmeans = toc(t_SC);
            G.time.SC.total = G.time.SC.kmeans + G.time.SC.eigs;
            
            G.SC.ncut = compute_ncut(G, G.SC.ass);
            G.SC.acc_gt = compute_accuracy(G.gt, G.SC.ass, G.gt);
        end

        if (~isfield(G, 'CSC') || k ~= k_old || redo)
            fprintf('\tCompute CSC...\n');
            params = struct('d_value', round(max(20, 10*log(G.N))), 'features', 1, 'assignment', 1, 'poly_order', 100, ...
                            'sampling', 'VD', 'n_factor', max(min(2, G.N/(G.k*log(G.k))), 0.1*G.N/(G.k*log(G.k))), 'regu', 1e-4);
            [G.CSC.ass, G.CSC.features, G.time.CSC, G.CSC.outputs, G.CSC.params] = my_CSC(G, params);
            
            fprintf('Computation of the full kmeans...');
            t_CSC_km = tic;
            G.CSC.ass_full = kmeans(G.CSC.features, G.k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
            G.time.CSC.kmeans = toc(t_CSC_km);
            if ~isfield(G.time.CSC, 'total')
                G.time.CSC.total = sum(cell2mat(struct2cell(G.time.CSC)));
            end
            fprintf('\t\tDone.\n');

            G.CSC.kmcost = compute_kmeans_cost(G.Uk, G.CSC.ass);
            G.CSC.ncut = compute_ncut(G, G.CSC.ass);
            G.CSC.acc_gt = compute_accuracy(G.gt, G.CSC.ass, G.gt);
            
            G.CSC.kmcost_full = compute_kmeans_cost(G.Uk, G.CSC.ass_full);
            G.CSC.ncut_full = compute_ncut(G, G.CSC.ass_full);
            G.CSC.acc_gt_full = compute_accuracy(G.gt, G.CSC.ass_full, G.gt);
        end
        save(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii), 'G');
    end
end

% clear classes authors_grand_topic classes_full gts;

%% run code
start_point = 2;

load(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, start_point -1));
d = size(G.CSC.features, 2);
G_old = G;

for ii = start_point:nb_cols
    fprintf('\nComputations for graph %d\n', ii);
    clear G;
    load(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii));

    if p > 0 && p < 0.5
        d = size(G_old.CSC.features, 2);
        cols_to_reuse = randperm(round(d*(1-p)), round(p*d)); % only among the lastly created signals
        new_sigs = sparse(N, round(p*d));
        old_idx = find(pres_final_gt(:, ii-1)==1);
        new_idx = find(pres_final_gt(:, ii)==1);
        
        all_idx = union(old_idx, new_idx);
        Ut1 = zeros(numel(all_idx), k);
        Ut1(ismember(all_idx, old_idx), :) = G_old.Uk;
        Ut = zeros(numel(all_idx), k);
        Ut(ismember(all_idx, new_idx), :) = G.Uk;
        G.rho = norm(Ut1*Ut1'- Ut*Ut', 'fro');

        new_sigs(old_idx, :) = G_old.CSC.features(:, cols_to_reuse);
        G.CSC.params.filt_sig = full(new_sigs(new_idx, :));

        new_weights = sparse(N, 1);
        new_weights(old_idx) = G_old.CSC.outputs.weight_VD;
        G.CSC.params.weight_VD = full(new_weights(new_idx));

        G.CSC.params.lk_est = G_old.CSC.outputs.lk_est;
    elseif p == 0
        if isfield(G.CSC.params, 'filt_sig'), G.CSC.params = rmfield(G.CSC.params, 'filt_sig'); end
        if isfield(G.CSC.params, 'lk_est'), G.CSC.params = rmfield(G.CSC.params, 'lk_est'); end
        if isfield(G.CSC.params, 'weight_VD'), G.CSC.params = rmfield(G.CSC.params, 'weight_VD'); end
    else
        error('pi is too high');
    end

%     params = struct('d_value', round(max(20, 10*log(G.N))), 'features', 1, 'assignment', 1, 'poly_order', 100, ...
%                             'sampling', 'VD', 'n_factor', max(min(2, G.N/(G.k*log(G.k))), 0.1*G.N/(G.k*log(G.k))), 'regu', 1e-4);
    [G.dCSC.ass, G.dCSC.features, G.time.dCSC, G.dCSC.outputs, ~] = my_CSC(G, G.CSC.params);

    if ~isfield(G.time.dCSC, 'lk_est'), G.time.dCSC.lk_est = 0; end
    G.dCSC.kmcost = compute_kmeans_cost(G.Uk, G.dCSC.ass);
    G.dCSC.ncut = compute_ncut(G, G.dCSC.ass);
    G.dCSC.acc_gt = compute_accuracy(G.gt, G.dCSC.ass, G.gt);

    fprintf('Computation of the full kmeans...');
    t_dCSC_km = tic;
    G.dCSC.ass_full = kmeans(G.dCSC.features, G.k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
    G.time.dCSC.kmeans = toc(t_dCSC_km);

    fprintf('\t\tDone.\n');
    if ~isfield(G.time.dCSC, 'total')
        G.time.dCSC.total = sum(cell2mat(struct2cell(G.time.dCSC)));
    end

    G.dCSC.kmcost_full = compute_kmeans_cost(G.Uk, G.dCSC.ass_full);
    G.dCSC.ncut_full = compute_ncut(G, G.dCSC.ass_full);
    G.dCSC.acc_gt_full = compute_accuracy(G.gt, G.dCSC.ass_full, G.gt);
    
    save(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii), 'G');
    
    G_old = G;
end

%% plotting
start_point = 2;
nb_cols = 120;

TITLES = cell({'Accuracy', '$k$-means cost', 'NCut', 'Time'});

measures = zeros(nb_cols - start_point + 1, 22); % kmcost (SC, CSC, dCSC), time (SC, CSC, dCSC)
e = zeros(nb_cols - start_point + 1, k+5);

for ii = start_point:nb_cols
    clear G;
    fprintf('\nLoading data for graph %d\n', ii);
    load(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, ii));
    
    measures(ii-start_point+1, 1) = G.SC.acc_gt;
    measures(ii-start_point+1, 2) = G.CSC.acc_gt;
    measures(ii-start_point+1, 3) = G.dCSC.acc_gt;
    measures(ii-start_point+1, 4) = G.CSC.acc_gt_full;
    measures(ii-start_point+1, 5) = G.dCSC.acc_gt_full;
    measures(ii-start_point+1, 6) = G.SC.kmcost;
    measures(ii-start_point+1, 7) = G.CSC.kmcost;
    measures(ii-start_point+1, 8) = G.dCSC.kmcost;
    measures(ii-start_point+1, 9) = G.CSC.kmcost_full;
    measures(ii-start_point+1, 10) = G.dCSC.kmcost_full;
    measures(ii-start_point+1, 11) = G.SC.ncut;
    measures(ii-start_point+1, 12) = G.CSC.ncut;
    measures(ii-start_point+1, 13) = G.dCSC.ncut;
    measures(ii-start_point+1, 14) = G.CSC.ncut_full;
    measures(ii-start_point+1, 15) = G.dCSC.ncut_full;
    measures(ii-start_point+1, 16) = G.time.SC.total;
    measures(ii-start_point+1, 17) = G.time.CSC.total;
    measures(ii-start_point+1, 18) = G.time.dCSC.total;
    measures(ii-start_point+1, 19) = G.time.CSC.lk_est + G.time.CSC.filtering + G.time.CSC.kmeans;
    if ~isfield(G.time.dCSC, 'lk_est'), G.time.dCSC.lk_est = 0; end
    measures(ii-start_point+1, 20) = G.time.dCSC.lk_est + G.time.dCSC.filtering + G.time.dCSC.kmeans;
    measures(ii-start_point+1, 21) = G.rho;
    measures(ii-start_point+1, 22) = G.N;

    e(ii-start_point+1, :) = G.e;
end

if plots
    for ii = 0:3
        figure(ii+1); hold on;
        for meas = [1, 4, 5] % 1:5
            if meas == 1, cod = 'o-'; else, cod = 'square'; end
            plot(start_point:nb_cols, measures(:, meas+5*ii), cod);
        end
        hold off;
        ylabel(TITLES{ii+1}, 'Interpreter', 'LaTex');
        xlabel('Graph index', 'Interpreter', 'LaTex');
        set(gca, 'FontSize', 16);
        set(gcf, 'Color', 'w');
        xlim([0, 120]);
%         legend('SC', 'CSC', 'dCSC', 'CSC full', 'dCSC full');
        legend('SC', 'CSC full', 'dCSC full');
    end


    figure; plot(start_point:nb_cols, measures(:, 21), 'o');
    ylabel('$\rho$', 'Interpreter', 'LaTex');
    xlabel('Graph index', 'Interpreter', 'LaTex');
    set(gca, 'FontSize', 16);
    set(gcf, 'Color', 'w');
    xlim([0, 120]);
    set(gcf, 'Position', [100 100 600 300]);
    
    figure; plot(start_point:nb_cols, measures(:, 22), '-');
    ylabel('Number of nodes', 'Interpreter', 'LaTex');
    xlabel('Graph index', 'Interpreter', 'LaTex');
    set(gca, 'FontSize', 16);
    set(gcf, 'Color', 'w');
    xlim([0, 120]);
    set(gcf, 'Position', [100 100 600 300]);

    figure; scatter(vec(repmat(start_point:nb_cols, 1, k+5)), vec(e), 10, 'filled');
end

save(sprintf('../matgraphs/%s/measures.mat', attempt_name), 'measures', 'e');