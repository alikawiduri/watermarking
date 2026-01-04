clc; clear; close all;
rng(1);

%% ================= LOAD TRAINED MODELS =================
load bit_corrector_models   % models, mu, sg

%% ================= PARAMETER =================
resolusihost = 256;
blok  = 8;
L     = blok^2;

alfa   = 1;
alfass = 25;

Ldek = 1;
sub  = 2;

%% ================= LOAD HOST =================
host = imread('host_image/cameraman.jpg');
host = im2gray(host);
host = double(imresize(host,[resolusihost resolusihost]));
host1d = reshape(host,[],1);

%% ================= LOAD WATERMARK =================
logo = imread('watermark/logo tel-u.png');
logo = rgb2gray(imresize(logo,[resolusihost/blok resolusihost/blok]));
logobw = imbinarize(logo);

wm = double(logobw(:));
wm(wm==0) = -1;   % {-1,+1}

%% ================= PN SEQUENCE =================
pn = 2*randi([0 1],blok,blok)-1;
pn_block = pn(1:2,1:2);

%% ================= EMBEDDING =================
xw = zeros(size(host1d));
So = cell(length(wm),1);

for i = 1:length(wm)

    seg = alfa * host1d(L*i-L+1:L*i);

    Sw = swt(seg,Ldek,'haar');
    sdst = dst(Sw(sub,:));
    Sdct = dct(sdst);
    Sdct2 = reshape(Sdct,[blok blok]);

    [U,S,V] = svd(Sdct2);
    So{i} = S;

    % Spread Spectrum on dominant singular
    Semb = S;
    Semb(1:2,1:2) = S(1:2,1:2) + alfass * wm(i) * pn_block;

    Srec = U*Semb*V';
    Sw(sub,:) = idst(idct(reshape(Srec,[],1)));
    xw(L*i-L+1:L*i) = iswt(Sw,'haar');
end

xwr = uint8(min(max(reshape(xw,size(host)),0),255));

%% ================= PSNR =================
mse = mean((double(host(:)) - double(xwr(:))).^2);
PSNR = 10*log10(255^2/mse);
fprintf('PSNR Watermarked = %.2f dB\n',PSNR);

%% ================= MULTI ATTACK EVALUATION =================
attack_list = {
    'JPEG_QF60',        [2 60];
    'JPEG_QF40',        [2 40];
    'Rotation_5deg',    [3 5];
    'Rotation_-5deg',   [3 -5];
    'Scaling_08',       [6 0.8];
    'Scaling_12',       [6 1.2];
    'Noise_001',        [17 0.01];
    'Noise_002',        [17 0.02];
};

model_list = fieldnames(models);

base_out = 'results_multi_attack';
if ~exist(base_out,'dir'), mkdir(base_out); end

for a = 1:size(attack_list,1)

    attack_name = attack_list{a,1};
    jenisserangan = attack_list{a,2};

    fprintf('\n=== ATTACK: %s ===\n', attack_name);

    %% ===== Apply attack =====
    attacked_img = stirmark_imagebenchmark(xwr, jenisserangan);
    xwn = reshape(double(attacked_img),[],1);

    %% ===== Output folder =====
    outdir = fullfile(base_out, attack_name);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    %% ===== Init =====
    wm_hat_noML = zeros(length(wm),1);
    wm_hat = struct();
    for m = 1:length(model_list)
        wm_hat.(model_list{m}) = zeros(length(wm),1);
    end

    %% ===== Extraction =====
    for i = 1:length(wm)

        seg = alfa * xwn(L*i-L+1:L*i);
        Sw = swt(seg,Ldek,'haar');
        sdst = dst(Sw(sub,:));
        Sdct = dct(sdst);
        Sdct2 = reshape(Sdct,[blok blok]);
        [~,S,~] = svd(Sdct2);

        D = S(1:2,1:2) - So{i}(1:2,1:2);

        corr_val = mean(D(:) .* pn_block(:));
        energy   = mean(abs(D(:)));
        meanD    = mean(D(:));
        stdD     = std(D(:));

        Xtest   = [corr_val energy meanD stdD];
        Xtest_s = (Xtest - mu) ./ sg;

        % ===== No ML =====
        wm_hat_noML(i) = sign(corr_val);
        if wm_hat_noML(i)==0, wm_hat_noML(i)=1; end

        % ===== ML =====
        for m = 1:length(model_list)
            name = model_list{m};
            pred = predict(models.(name), Xtest_s);
            if iscell(pred), pred = str2double(pred{1}); end
            wm_hat.(name)(i) = pred*2 - 1;
        end
    end

    %% ================= RF POST-PROCESSING (MAJORITY FILTER) =================
    wm_rf_img = reshape(wm_hat.rf==1, size(logobw));
    wm_rf_pp  = medfilt2(wm_rf_img,[3 3]);
    wm_hat_rf_pp = wm_rf_pp(:)*2 - 1;

    %% ================= EVALUATION =================
    BER_noML = mean(wm_hat_noML ~= wm);
    fprintf('No ML | BER = %.4f\n', BER_noML);

    BER_ml = zeros(1,length(model_list));
    for m = 1:length(model_list)
        BER_ml(m) = mean(wm_hat.(model_list{m}) ~= wm);
        fprintf('%-7s | BER = %.4f\n', upper(model_list{m}), BER_ml(m));
    end

    BER_RF_PP = mean(wm_hat_rf_pp ~= wm);
    fprintf('RF+PP | BER = %.4f\n', BER_RF_PP);

    %% ================= SAVE IMAGES =================
    imwrite(attacked_img, fullfile(outdir,'attacked_image.png'));
    imwrite(reshape(wm_hat_noML==1,size(logobw)), fullfile(outdir,'wm_NoML.png'));

    for m = 1:length(model_list)
        imwrite(reshape(wm_hat.(model_list{m})==1,size(logobw)), ...
            fullfile(outdir,['wm_' upper(model_list{m}) '.png']));
    end

    imwrite(wm_rf_pp, fullfile(outdir,'wm_RF_postprocess.png'));

    %% ================= SAVE BER GRAPH =================
    fig = figure('Visible','off');

    BER_plot = [ ...
        BER_noML * ones(1,length(model_list)); ...
        BER_ml; ...
        BER_RF_PP * ones(1,length(model_list)) ...
    ]';

    bar(BER_plot,'grouped'); grid on
    set(gca,'XTickLabel',upper(model_list),'FontSize',11)
    ylabel('Bit Error Rate (BER)')
    xlabel('Model')
    title(['BER Comparison - ' attack_name])

    legend({'No ML','ML','RF + PostProcess'},'Location','northwest')

    saveas(fig, fullfile(outdir,'BER_comparison.png'));
    close(fig)

end

fprintf('\nAll attack results saved in folder: %s\n', base_out);
