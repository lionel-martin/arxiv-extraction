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

N = 366572;
nb_years = 1;
start_year = 2005;
study_end = 2006;

all_W = sparse(N, N);
presence = zeros(N, (study_end-start_year-nb_years+1)*12+1);
all_authors = [];

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
        %all_authors = unique([all_authors, active_authors]);
    end
end

%%
all_W = all_W + all_W';
thres = 20;
comp_params = struct('threshold', thres, 'only_first', true);

all_W(1:size(all_W, 1) + 1:end) = 0;
connected_nodes = gsp_components_v2(gsp_graph(all_W), comp_params);
presence(find(connected_nodes), 1) = 1;

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
        connected_nodes = gsp_components_v2(gsp_graph(all_W), comp_params);
        presence(find(connected_nodes), 12*(ii-start_year-nb_years)+jj+1) = 1;
        fprintf('\t%d\n', sum(connected_nodes));
    end
end

diff_pres = diff(presence, 1, 2);
never_disapearing_nodes = find(sum(abs(diff_pres)-diff_pres, 2)==0);
always_on_nodes = find(presence(never_disapearing_nodes, 1));

[disapearing_node, disapering_month] = ind2sub(size(presence), find(diff_pres<0));

given_year = 2002; given_month = 1;
node_exist_in_month = find(presence(:, 12*(given_year-start_year-nb_years)+given_month+1));


idx_remove = setdiff(1:N, all_authors);
nb_auth = length(all_authors);
all_W(idx_remove, :) = []; all_W(:, idx_remove) = []; all_W(1:nb_auth+1:end) = 0;
G = gsp_graph(all_W);
% TODO : Generate coordinates
plot_params.show_edges = 1;
gsp_plot_graph(G, plot_params);