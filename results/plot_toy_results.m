addpath(genpath('shadedErrorBar'));
addpath(genpath('~/matlab/plotting'));

filename = 'result_G1_F2';

load([filename '.mat']);

% figure, errorbar(sample_size, 1/1000*mean(error_set_U_SR,2), std(1/1000*error_set_U_SR,0,2), '-r', 'LineWidth', 2);
% hold on;
% errorbar(sample_size, 1/1000*mean(error_set_S_U,2), std(1/1000*error_set_S_U,0,2), 'color', [0 0.5 0], 'LineWidth', 2);
% % errorbar(sample_size, 1/1000*mean(error_set_rand_U_VR,2), std(1/1000*error_set_rand_U_VR,0,2), '-k', 'LineWidth', 2);
% errorbar(sample_size, 1/1000*mean(error_set_L_k,2), std(1/1000*error_set_L_k,0,2), '-b', 'LineWidth', 2);


figure, plot(sample_size, 1/1000*mean(error_set_rand_uni,2), '-ko', 'LineWidth', 1, 'LineStyle','-');
hold on;
plot(sample_size, 1/1000*mean(error_set_U_SR,2), '-r^', 'LineWidth', 1, 'LineStyle','-');
plot(sample_size, 1/1000*mean(error_set_S_U,2), 'color', [0 0.5 0], 'marker', 'v', 'LineWidth', 1, 'LineStyle','-');
% plot(sample_size, 1/1000*mean(error_set_rand_U_VR,2), '-k', 'LineWidth', 2);

plot(sample_size, 1/1000*mean(error_set_L_2,2), 'color', [0 0 1], 'marker', 's', 'LineWidth', 1, 'LineStyle','-');
plot(sample_size, 1/1000*mean(error_set_L_8,2), 'color', [1 0 1], 'marker', 's', 'LineWidth', 1, 'LineStyle','-');
plot(sample_size, 1/1000*mean(error_set_L_14,2), 'color', [1 0.7 0], 'marker', 's', 'LineWidth', 1, 'LineStyle','-');



font_size = 14;
set(gca,'FontSize',font_size);

ylim([0 10e-3]);
xlim([0 200]);

xlabel('Sample size','FontSize',font_size);
ylabel('Reconstruction MSE','FontSize',font_size);
% legend('M1', 'M2', 'M3','Proposed method');

legend1 = legend('Uni. rand.', 'M1', 'M2', 'Prop. k = 2', 'Prop. k = 8', 'Prop. k = 14');
set(legend1, 'Position',[0.16 0.18 0.28 0.2],'FontSize',font_size-2);

% legend1 = legend('Uni. rand.', 'M1', 'M2', 'Prop. k = 2', 'Prop. k = 8', 'Prop. k = 14');
% set(legend1, 'Position',[0.6 0.65 0.28 0.2],'FontSize',font_size-2);

export_fig([filename '.pdf'], '-transparent');

