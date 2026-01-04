clc; clear; close all;
rng(1);

SAVE_FILE = 'best_parameters_all_attacks.mat';

if exist(SAVE_FILE,'file')
    load(SAVE_FILE,'Results','last_attack_done');
    fprintf('Resuming from attack %d\n', last_attack_done+1);
else
    Results = struct([]);
    last_attack_done = 0;
    fprintf('Starting from attack 1\n');
end

%% =========================================================
% LOAD HOST & WATERMARK
%% =========================================================
host = im2double(imread('host_image/airplane.png'));
if ndims(host)==3
    host = rgb2gray(host);
end
host = imresize(host,[256 256]);
host1d = reshape(host,[],1);

logo = imread('watermark/logo tel-u.png');
logo = imbinarize(rgb2gray(imresize(logo,[32 32])));
wm = double(logo(:));
wm(wm==0) = -1;

%% =========================================================
% PARAMETER SEARCH SPACE
%% =========================================================
sub_list    = [1 2];
blok_list   = [8 16];
alfa_list   = [0.5 1];
alfass_list = [10 20 30];

Ldek = 1;

%% =========================================================
% ATTACK POOL (SEMUA SERANGAN)
%% =========================================================
attack_pool = {};

for qf = [20 30 40 50 60 70 80 90]
    attack_pool{end+1} = [2 qf];            % JPEG
end
for ang = [-7 -5 -3 3 5 7]
    attack_pool{end+1} = [3 ang];           % Rotation
end
for c = [5 10 15]
    attack_pool{end+1} = [5 c];             % Cropping
end
for sc = [0.7 0.8 0.9 1.1 1.2 1.3]
    attack_pool{end+1} = [6 sc];            % Scaling
end
attack_pool{end+1} = [7 0.81];              % Non-uniform scale
attack_pool{end+1} = [7 1.21];

for k = [3 5 7 9]
    attack_pool{end+1} = [12 k];            % LPF
end

attack_pool{end+1} = [14 1];                % Histogram

for g = [0.6 0.8 1.2 1.4]
    attack_pool{end+1} = [15 g];            % Gamma
end

for q = [8 16 32]
    attack_pool{end+1} = [16 q];             % Color quant
end

for n = [0.005 0.01 0.02 0.03]
    attack_pool{end+1} = [17 n];             % Noise
end

fprintf('Total attack variants = %d\n', numel(attack_pool));

%% =========================================================
% RESULT STORAGE
%% =========================================================
Results = struct();
res_idx = 1;

%% =========================================================
% GRID SEARCH
%% =========================================================
for a = last_attack_done+1 : length(attack_pool)

    attack = attack_pool{a};
    fprintf('\nAttack %d / %d\n', a, length(attack_pool));

    bestBER = inf;
    bestParam = struct();

    for sub = sub_list
    for blok = blok_list
    for alfa = alfa_list
    for alfass = alfass_list

        % ================= SAFETY CHECK =================
        L = blok^2;
        if L * length(wm) > numel(host1d)
            continue;   % skip konfigurasi tidak valid
        end

        try
            % ===== EMBEDDING =====
            pn = 2*randi([0 1],blok,blok)-1;
            pn_block = pn(1:2,1:2);

            xw = zeros(size(host1d));
            So = cell(length(wm),1);

            for i = 1:length(wm)
                seg = alfa * host1d(L*i-L+1:L*i);
                Sw = swt(seg,Ldek,'haar');
                Sdct = dct(dst(Sw(sub,:)));
                Sdct2 = reshape(Sdct,[blok blok]);

                [U,S,V] = svd(Sdct2);
                So{i} = S;

                S(1:2,1:2) = S(1:2,1:2) + alfass*wm(i)*pn_block;
                Sw(sub,:) = idst(idct(reshape(U*S*V',[],1)));
                xw(L*i-L+1:L*i) = iswt(Sw,'haar');
            end

            xwr = reshape(xw,size(host));

            % ===== ATTACK =====
            attacked = stirmark_imagebenchmark(xwr, attack);
            xwn = reshape(double(attacked),[],1);

            % ===== EXTRACTION =====
            wm_hat = zeros(length(wm),1);
            for i = 1:length(wm)
                seg = alfa * xwn(L*i-L+1:L*i);
                Sw = swt(seg,Ldek,'haar');
                Sdct = dct(dst(Sw(sub,:)));
                Sdct2 = reshape(Sdct,[blok blok]);
                [~,S2,~] = svd(Sdct2);

                D = S2(1:2,1:2) - So{i}(1:2,1:2);
                wm_hat(i) = sign(mean(D(:).*pn_block(:)));
                if wm_hat(i)==0, wm_hat(i)=1; end
            end

            BER = mean(wm_hat ~= wm);

            if BER < bestBER
                bestBER = BER;
                bestParam.sub = sub;
                bestParam.blok = blok;
                bestParam.alfa = alfa;
                bestParam.alfass = alfass;
            end

        catch ME
            fprintf('⚠️  Skip config: %s\n', ME.message);
        end

    end,end,end,end

    % ===== SAVE RESULT PER ATTACK =====
    Results(a).attack  = attack;
    Results(a).bestBER = bestBER;
    Results(a).param   = bestParam;

    last_attack_done = a;
    save(SAVE_FILE,'Results','last_attack_done');

    fprintf('✔ Attack %d saved (best BER = %.4f)\n', a, bestBER);
end
