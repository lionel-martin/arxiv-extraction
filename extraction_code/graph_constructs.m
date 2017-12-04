%%%%%%% Matlab GraphConstructs
% Illustrate examples of graph constructs from matfiles

%%%%%%%

%% One month graph construction

%load('../mat_citation/2000-09.mat');
W = sparse(double(i+1), double(j+1), double(data));
W(1:N+1:end) = 0;
W = W + W';
G = gsp_graph(W(active_authors, active_authors));
G.k = 50;

params_CSC = struct('features', 1, 'assignment', 0, 'poly_order', 100, 'sampling', 'VD');
[C_est, data.raw, ~, ~] = CSC(G, params_CSC);
[~, data.labels] = max(C_est./repmat(sqrt(sum(C_est.^2, 1)), G.N, 1), [], 2);

dembed = compute_compressive_embedding(G, data);

figure;
scatter(dembed(:,1), dembed(:,2), 5, data.labels, 'filled')
colormap jet

G.coords = dembed(:, 1:2);
figure;
gsp_plot_graph(G);

%% Aggregation of several months
clear all;
N = 366572;
nb_years = 3;
start_year = 1992;
study_end = 2008;

thres = 20;
k_thres = 0;

filter_not_frequent = false;
attempt_name = 'attempt6';

bfs_par = struct('verbose', 0, 'threshold', 0.8);

load(sprintf('../matgraphs/%s/authors_main_grand_topic.mat', attempt_name));
auth = authors_grand_topic(:, 1);

all_W = sparse(N, N);
nb_cols = (study_end - start_year - nb_years + 1) * 12 + 1;
pres_thres20 = sparse(N, nb_cols);
pres_th20_100 = sparse(N, nb_cols);
pres_final_gt = sparse(N, nb_cols);

load(sprintf('/mnt/data/Arxiv/matgraphs/%s/presence_thres20.mat', attempt_name));
load(sprintf('/mnt/data/Arxiv/matgraphs/%s/presence_th20_100times.mat', attempt_name));

components = zeros(1, nb_cols);
% graphz = cell(nb_cols, 1);

STEPS = {'thresholding', 'frequent presence', 'connectivity'};

%%
for step = 1:3

    fprintf('Step %d - %s\n\n', step, STEPS{step});

    if step == 2

        if filter_not_frequent
            pres100 = find(sum(pres_thres20, 2) >= 100);
            pres_th20_100(pres100, :) = pres_thres20(pres100, :);
        else
            pres_th20_100 = pres_thres20;
        end

        fprintf('Done.\n');

    else
        all_W = sparse(N, N);
        for ii = start_year:start_year+nb_years-1
            for jj = 1:12
                mon = num2str(jj);
                if jj < 10
                    mon = strcat('0', mon);
                end
                fprintf('Processing %s-%d\n', mon, ii);
                load(strcat('../mat_citation/', num2str(ii), '-', mon, '.mat'));
                all_W = all_W + sparse(double(i+1), double(j+1), double(data), N, N);
                clear i j data;
            end
        end

        all_W = all_W + all_W';
        all_W(1:size(all_W, 1) + 1:end) = 0;

        if step == 1
            pres_idx = find(sum(all_W(auth, auth), 1) >= thres);
            pres_thres20(auth(pres_idx), 1) = 1;

%             pres_idx = find(sum(all_W, 1) >= thres);
%             pres_thres20(pres_idx, 1) = 1;

            fprintf('\t%d.\n', numel(pres_idx));

        elseif step == 3
            pres_idx = find(pres_th20_100(:, 1) == 1);
            Ni = numel(pres_idx);

            if numel(pres_idx) > 0

                GW = all_W(pres_idx, pres_idx);
                if k_thres > 0

                    Gii = zeros(1, Ni*k_thres);
                    Gjj = zeros(1, Ni*k_thres);
                    Gvv = zeros(1, Ni*k_thres);
                    next_idx = 1;
                    for nr = 1:Ni
                        [~, nc, nv] = find(GW(nr, :));
                        [~, ord_idx] = sort(nv, 'descend');
                        ord_thres = min(k_thres, numel(ord_idx));
                        sub_idx = ord_idx(1:ord_thres);
                        Gii(next_idx:next_idx + ord_thres - 1) = repmat(nr, 1, ord_thres);
                        Gjj(next_idx:next_idx + ord_thres - 1) =  nc(sub_idx);
                        Gvv(next_idx:next_idx + ord_thres - 1) =  nv(sub_idx);
                        next_idx = next_idx + ord_thres;
                    end
                    GW = sparse(Gii(1:next_idx-1), Gjj(1:next_idx-1), Gvv(1:next_idx-1), Ni, Ni);
                    Gdiff = GW - GW';
                    [Gii, Gjj, Gvv] = find(Gdiff);
                    add_idx = find(Gvv<0);
                    GW = GW + sparse(Gii(add_idx), Gjj(add_idx), -Gvv(add_idx), Ni, Ni);
                end

                G_tmp = gsp_graph(GW);
                comps_id = gsp_components_v2(G_tmp, bfs_par);
                [sizes, vals] = hist(comps_id, unique(comps_id));
                [components(1), comp_arg] = max(sizes);
                comp_id = vals(comp_arg);
                prez_idx = pres_idx(comps_id == comp_id);
                pres_final_gt(prez_idx, 1) = 1;
                G = gsp_graph(all_W(prez_idx, prez_idx));
                G.nodes_idx = prez_idx;
                assert(numel(prez_idx) == G.N);
%                 graphz{1} = G;
                save(sprintf('../matgraphs/%s/graphs/G_1.mat', attempt_name), 'G');
            end
            fprintf('\t%d, %d\n', numel(pres_idx), components(1));

        end

        for ii = start_year+nb_years:study_end
            for jj = 1:12
                mon = num2str(jj);
                if jj < 10
                    mon = strcat('0', mon);
                end

                fprintf('Processing %s-%d', mon, ii);

                load(strcat('../mat_citation/', num2str(ii), '-', mon, '.mat'));
                all_W = all_W + sparse(double(i+1), double(j+1), double(data), N, N) + sparse(double(j+1), double(i+1), double(data), N, N);
                clear i j data;

                load(strcat('../mat_citation/', num2str(ii-nb_years), '-', mon, '.mat'));
                all_W = all_W - sparse(double(i+1), double(j+1), double(data), N, N) - sparse(double(j+1), double(i+1), double(data), N, N);
                clear i j data;

                all_W(1:size(all_W, 1) + 1:end) = 0;

                idx = 12*(ii-start_year-nb_years)+jj+1;


                if step == 1
                    pres_idx = find(sum(all_W(auth, auth), 1) >= thres);
                    pres_thres20(auth(pres_idx), idx) = 1;

%                     pres_idx = find(sum(all_W, 1) >= thres);
%                     pres_thres20(pres_idx, idx) = 1;

                    fprintf('\t%d.\n', numel(pres_idx));

                elseif step == 3
                    pres_idx = find(pres_th20_100(:, idx) == 1);

                    if numel(pres_idx) > 0

                        GW = all_W(pres_idx, pres_idx);
                        if k_thres > 0
                            Ni = numel(pres_idx);
                            
                            Gii = zeros(1, Ni*k_thres);
                            Gjj = zeros(1, Ni*k_thres);
                            Gvv = zeros(1, Ni*k_thres);
                            next_idx = 1;
                            for nr = 1:Ni
                                [~, nc, nv] = find(GW(nr, :));
                                [~, ord_idx] = sort(nv, 'descend');
                                ord_thres = min(k_thres, numel(ord_idx));
                                sub_idx = ord_idx(1:ord_thres);
                                Gii(next_idx:next_idx + ord_thres - 1) = repmat(nr, 1, ord_thres);
                                Gjj(next_idx:next_idx + ord_thres - 1) =  nc(sub_idx);
                                Gvv(next_idx:next_idx + ord_thres - 1) =  nv(sub_idx);
                                next_idx = next_idx + ord_thres;
                            end
                            GW = sparse(Gii(1:next_idx-1), Gjj(1:next_idx-1), Gvv(1:next_idx-1), Ni, Ni);
                            Gdiff = GW - GW';
                            [Gii, Gjj, Gvv] = find(Gdiff);
                            add_idx = find(Gvv<0);
                            GW = GW + sparse(Gii(add_idx), Gjj(add_idx), -Gvv(add_idx), Ni, Ni);
                        end

                        G_tmp = gsp_graph(GW);
                        comps_id = gsp_components_v2(G_tmp, bfs_par);
                        [sizes, vals] = hist(comps_id, unique(comps_id));
                        [components(idx), comp_arg] = max(sizes);
                        comp_id = vals(comp_arg);
                        prez_idx = pres_idx(comps_id == comp_id);
                        pres_final_gt(prez_idx, idx) = 1;
                        G = gsp_graph(all_W(prez_idx, prez_idx));
                        G.nodes_idx = prez_idx;
                        assert(numel(prez_idx) == G.N);
%                         graphz{idx} = G;
                        save(sprintf('../matgraphs/%s/graphs/G_%d.mat', attempt_name, idx), 'G');

                    end
                    fprintf('\t%d, %d\n', numel(pres_idx), components(idx));
                end
            end
        end
    end

    switch step
        case 1
        save(sprintf('../matgraphs/%s/presence_thres20.mat', attempt_name), 'pres_thres20');
        case 2
        save(sprintf('../matgraphs/%s/presence_th20_100times.mat', attempt_name), 'pres_th20_100');
        case 3
        save(sprintf('../matgraphs/%s/presence_th20_100_connected.mat', attempt_name), 'pres_final_gt');
%         save(sprintf('../matgraphs/%s/graphs_th20_100_final.mat', attempt_name), 'graphz', '-v7.3');
        clear pres_final_gt pres_thres20 pres_th20_100;
%         clear graphz;
    end
end