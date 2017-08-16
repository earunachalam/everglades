clear; clc; close all;

a = csvread('pole_ito_k_536_00.csv');
b = csvread('pole_ito_k_516_00.csv');
c = csvread('pole_ito_k_536_04.csv');
d = csvread('pole_ito_k_516_04.csv');

figure('PaperPositionMode', 'auto');
hold on
plot(a(:,1), 0.0*a(:,1))
plot(a(:,1),a(:,2),'LineWidth',2)
plot(a(:,1),b(:,2),'LineWidth',2)
plot(a(:,1),c(:,2),'LineWidth',2)
plot(a(:,1),d(:,2),'LineWidth',2)

set(gca, 'FontSize', 14)
xlabel('$|\mathbf{k}|$', 'Interpreter', 'latex')
ylabel('growth rate $\sigma(\mathbf{k})$', 'Interpreter', 'latex')

print(gcf, '-dpdf', '~/Dropbox/weekly-report/03Aug17/disp_reln.pdf')