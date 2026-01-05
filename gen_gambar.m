clc; clear; close all;

%% ===================== DATA =====================

attacks = {'JPEG QF60','JPEG QF40','Rotasi +5','Rotasi -5', ...
           'Scale 0.8','Scale 1.2','Noise 0.01','Noise 0.02'};

% ===== BER =====
BER_cam = [0.0225 0.0537 0.8457 0.7998 0.0684 0.1533 0.8115 0.8750];
BER_air = [0.0361 0.0381 0.8604 0.8760 0.0498 0.1416 0.8193 0.8701];

%% ===================== FIGURE =====================

figure('Color','w','Position',[200 200 1100 500]);

%% ---------- SUBPLOT : BER ----------
subplot(1,2,1);

dataBER = [BER_cam; BER_air]';   % ukuran [8 x 2]
b = bar(dataBER,'grouped');
grid on;

% ===== WARNA BAR =====
b(1).FaceColor = [0.3 0.3 0.3];   % Cameraman (abu)
b(2).FaceColor = [0.2 0.6 0.9];   % Airplane (biru)

set(gca,'XTickLabel',attacks,'FontSize',10);
xtickangle(90);

ylabel('Bit Error Rate (BER)');
title('BER Comparison of Watermark Extraction');
legend({'Cameraman','Airplane'},'Location','northwest');

ylim([0 1]);

%% ===================== EXPORT (OPTIONAL) =====================
exportgraphics(gcf,'BER_Comparison_Cameraman_vs_Airplane.png','Resolution',300);
