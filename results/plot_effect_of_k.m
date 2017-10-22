addpath(genpath('shadedErrorBar'));
addpath(genpath('exportfig'));

clear

filename = 'result1_ER_P10_F2';

load([filename '.mat']);

% figure, errorbar(sample_size, 1/1000*mean(error_set_U_SR,2), std(1/1000*error_set_U_SR,0,2), '-r', 'LineWidth', 2);
% hold on;
% errorbar(sample_size, 1/1000*mean(error_set_S_U,2), std(1/1000*error_set_S_U,0,2), 'color', [0 0.5 0], 'LineWidth', 2);
% % errorbar(sample_size, 1/1000*mean(error_set_rand_U_VR,2), std(1/1000*error_set_rand_U_VR,0,2), '-k', 'LineWidth', 2);
% errorbar(sample_size, 1/1000*mean(error_set_L_k,2), std(1/1000*error_set_L_k,0,2), '-b', 'LineWidth', 2);


figure, plot(sample_size, 1/1000*mean(error_set_L_k1,2), '-ko', 'LineWidth', 2);
hold on;
plot(sample_size, 1/1000*mean(error_set_L_k2,2), '-r^', 'LineWidth', 2);
plot(sample_size, 1/1000*mean(error_set_L_k3,2), 'color', [0 0.5 0], 'marker', 'v', 'LineWidth', 2);
% plot(sample_size, 1/1000*mean(error_set_rand_U_VR,2), '-k', 'LineWidth', 2);
% plot(sample_size, 1/1000*mean(error_set_L_k1,2), '-bs', 'LineWidth', 2);


font_size = 14;
set(gca,'FontSize',font_size);

ylim([0 10e-3]);
xlim([0 200]);

xlabel('Sample size','FontSize',font_size);
ylabel('Reconstruction MSE','FontSize',font_size);
% legend('M1', 'M2', 'M3','Proposed method');
legend1 = legend('k = 2', 'k = 8', 'k = 14');
set(legend1,...
    'Position',[0.63 0.72 0.27 0.16]);

export_fig([filename '.pdf'], '-transparent');

