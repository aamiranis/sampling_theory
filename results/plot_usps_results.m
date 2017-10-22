addpath(genpath('shadedErrorBar'));
addpath(genpath('~/matlab/plotting'));

load('usps_results.mat');

% figure, errorbar(sample_size, 1/1000*mean(error_set_U_SR,2), std(1/1000*error_set_U_SR,0,2), '-r', 'LineWidth', 2);
% hold on;
% errorbar(sample_size, 1/1000*mean(error_set_S_U,2), std(1/1000*error_set_S_U,0,2), 'color', [0 0.5 0], 'LineWidth', 2);
% % errorbar(sample_size, 1/1000*mean(error_set_rand_U_VR,2), std(1/1000*error_set_rand_U_VR,0,2), '-k', 'LineWidth', 2);
% errorbar(sample_size, 1/1000*mean(error_set_L_k,2), std(1/1000*error_set_L_k,0,2), '-b', 'LineWidth', 2);


% figure, plot(sample_size, mean(error_rand_uni,2), '-ko', 'LineWidth', 2, 'LineStyle','-');
% hold on;
figure;
plot(sample_size, mean(error_U_SR,2), '-r^', 'LineWidth', 2, 'LineStyle','-');
hold on;
plot(sample_size, mean(error_S_U,2), 'color', [0 0 1], 'marker', 'v', 'LineWidth', 2, 'LineStyle','-');
plot(sample_size, mean(error_L_k,2), 'color', [0 0 0], 'marker', 'o', 'LineWidth', 2, 'LineStyle','-');



font_size = 14;
set(gca,'FontSize',font_size);

ylim([0.15 0.22]);
xlim([50 160]);

xlabel('Number of labels','FontSize',font_size);
ylabel('Classification error','FontSize',font_size);
% legend('M1', 'M2', 'M3','Proposed method');

% legend1 = legend('Uni. rand.', 'M1', 'M2', 'Proposed');
% set(legend1, 'Position',[0.16 0.18 0.28 0.2],'FontSize',font_size-2);

legend1 = legend('M1', 'M2', 'Proposed');
set(legend1, 'Position',[0.6 0.7 0.28 0.2],'FontSize',font_size-2);

export_fig('compare_methods_usps_norm_exact_BL.pdf', '-transparent');


%%

load('usps_results_compare_bases.mat');

% figure, errorbar(sample_size, 1/1000*mean(error_set_U_SR,2), std(1/1000*error_set_U_SR,0,2), '-r', 'LineWidth', 2);
% hold on;
% errorbar(sample_size, 1/1000*mean(error_set_S_U,2), std(1/1000*error_set_S_U,0,2), 'color', [0 0.5 0], 'LineWidth', 2);
% % errorbar(sample_size, 1/1000*mean(error_set_rand_U_VR,2), std(1/1000*error_set_rand_U_VR,0,2), '-k', 'LineWidth', 2);
% errorbar(sample_size, 1/1000*mean(error_set_L_k,2), std(1/1000*error_set_L_k,0,2), '-b', 'LineWidth', 2);


% figure, plot(sample_size, mean(error_rand_uni,2), '-ko', 'LineWidth', 2, 'LineStyle','-');
% hold on;
figure;
plot(sample_size, mean(error_adj,2), '-bo', 'LineWidth', 2, 'LineStyle','-');
hold on;
plot(sample_size, mean(error_ha,2), 'color', [1 0 0], 'marker', 's', 'LineWidth', 2, 'LineStyle','-');
plot(sample_size, mean(error_rw,2), 'color', [0 0 0], 'marker', '^', 'LineWidth', 2, 'LineStyle','-');



font_size = 14;
set(gca,'FontSize',font_size);

ylim([0.15 0.22]);
xlim([50 160]);

xlabel('Number of labels','FontSize',font_size);
ylabel('Classification error','FontSize',font_size);
% legend('M1', 'M2', 'M3','Proposed method');

legend1 = legend('Adjacency', 'Hub-authority', 'Random walk');
set(legend1, 'Position',[0.58 0.7 0.30 0.2],'FontSize',font_size-2);

export_fig('compare_bases_usps_norm_exact_BL.pdf', '-transparent');

%%

load('usps_results_compare_bases.mat');


figure;
plot(1:1:N, energy_fraction_adj, '-b', 'LineWidth', 2, 'LineStyle','-');
hold on;
plot(1:1:N, energy_fraction_ha, '-r', 'LineWidth', 2, 'LineStyle','-');
plot(1:1:N, energy_fraction_rw, '-k', 'LineWidth', 2, 'LineStyle','-');

font_size = 14;
set(gca,'FontSize',font_size);

ylim([0 1.1]);
xlim([1 1000]);

set(gca,'XTick',0:200:1000);

xlabel('Number of GFT coefficients','FontSize',font_size);
ylabel('Energy fraction','FontSize',font_size);
% legend('M1', 'M2', 'M3','Proposed method');

legend1 = legend('Adjacency', 'Hub-authority', 'Random walk');
set(legend1, 'Position',[0.58 0.25 0.30 0.2],'FontSize',font_size-2);

export_fig('gft_bases_usps_norm.pdf', '-transparent');