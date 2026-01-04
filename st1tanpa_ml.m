clc; clear; close all;
rng(1);   % reproducible

%% ================= PARAMETER =================
resolusihost = 256;
blok  = 8;
L     = blok^2;

alfa   = 1;        % jangan dibesarin
alfass = 25;       % strength SS (20–35 aman)

Ldek = 1;
sub  = 2;          % SWT subband

%% ================= LOAD HOST =================
host = imread('host_image/airplane.png');
host = im2gray(host);
host = double(imresize(host,[resolusihost resolusihost]));
host1d = reshape(host,[],1);

%% ================= LOAD WATERMARK =================
logo = imread('watermark/logo tel-u.png');
logo = rgb2gray(imresize(logo,[resolusihost/blok resolusihost/blok]));
logobw = imbinarize(logo);

% NRZ mapping {-1,+1}
wm = double(logobw(:));
wm(wm==0) = -1;

%% ================= PN SEQUENCE =================
pn = 2*randi([0 1],blok,blok)-1;

%% ================= EMBEDDING =================
xw = zeros(size(host1d));
So = cell(length(wm),1);

for i = 1:length(wm)

    seg = alfa * host1d(L*i-L+1:L*i);

    % SWT
    Sw = swt(seg,Ldek,'haar');

    % DST
    sdst = dst(Sw(sub,:));

    % DCT
    Sdct = dct(sdst);

    % reshape
    Sdct2 = reshape(Sdct,[blok blok]);

    % SVD
    [U,S,V] = svd(Sdct2);
    So{i} = S;

    % ===== SS EMBEDDING (DOMINANT SINGULAR ONLY) =====
    Semb = S;
    Semb(1:2,1:2) = S(1:2,1:2) + ...
        alfass * wm(i) * pn(1:2,1:2);

    % ISVD
    Srec = U * Semb * V';

    % inverse transform
    Srec1 = reshape(Srec,[],1);
    idct1 = idct(Srec1);
    Sw(sub,:) = idst(idct1);
    st = iswt(Sw,'haar');

    xw(L*i-L+1:L*i) = st / alfa;
end

xwr = reshape(xw,size(host));
xwr = uint8(min(max(xwr,0),255));

%% ================= PSNR =================
host_uint8 = uint8(host);        % host asli
mse = mean((double(host_uint8(:)) - double(xwr(:))).^2);
PSNR = 10 * log10(255^2 / mse);

%% ================= EXTRACTION =================
%xwn = reshape(double(xwr),[],1);
jenisserangan = [14 1];     % contoh: histogram modification
% jenisserangan = [1 80];   % JPEG QF = 80
attacked_img = stirmark_imagebenchmark(xwr, jenisserangan);

%% ================= EXTRACTION =================
xwn = reshape(double(attacked_img),[],1);
wm_hat = zeros(length(wm),1);

for i = 1:length(wm)

    seg = alfa * xwn(L*i-L+1:L*i);

    Sw = swt(seg,Ldek,'haar');
    sdst = dst(Sw(sub,:));
    Sdct = dct(sdst);
    Sdct2 = reshape(Sdct,[blok blok]);

    [~,S,~] = svd(Sdct2);

    % ===== SS DETECTION (DOMINANT AREA) =====
    D = S(1:2,1:2) - So{i}(1:2,1:2);
    pn_block = pn(1:2,1:2);
    corr_val = mean( D(:) .* pn_block(:) );


    if corr_val >= 0
        wm_hat(i) = 1;
    else
        wm_hat(i) = -1;
    end
end

%% ================= EVALUATION =================
BER = sum(wm_hat ~= wm) / length(wm);
NC  = corr(double(wm_hat),double(wm));

fprintf('BER  (No Attack, No ML) = %.6f\n', BER);
fprintf('NC   (No Attack, No ML) = %.6f\n', NC);
fprintf('PSNR (Watermarked)     = %.2f dB\n', PSNR);

%% ================= DISPLAY =================
wm_hat_bin = wm_hat;
wm_hat_bin(wm_hat_bin==-1) = 0;

figure
subplot(1,3,1)
imshow(logobw), title('Original Watermark')

subplot(1,3,2)
imshow(reshape(wm_hat_bin,size(logobw)))
title('Extracted Watermark (No Attack)')

subplot(1,3,3)
imshow(xwr), title('Watermarked Image')

save embed_state xwr So pn wm alfa alfass blok L sub Ldek
