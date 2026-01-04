clc; clear; close all;
rng(1);

%% =========================================================
% LOAD EMBED STATE
% =========================================================
% WAJIB ada dari proses embedding baseline
% berisi: xwr, wm, So, pn, alfa, L, blok, Ldek, sub
load embed_state

pn_block = pn(1:2,1:2);

%% =========================================================
% BUILD ATTACK POOL
% =========================================================

attack_pool = {};

% ---------- JPEG ----------
for qf = [20 30 40 50 60 70 80 90]
    attack_pool{end+1} = [2 qf];
end

% ---------- ROTATION ----------
for ang = [-7 -5 -3 3 5 7]
    attack_pool{end+1} = [3 ang];
end

% ---------- RANDOM CROPPING ----------
for c = [5 10 15]
    attack_pool{end+1} = [5 c];
end

% ---------- UNIFORM SCALING ----------
for sc = [0.7 0.8 0.9 1.1 1.2 1.3]
    attack_pool{end+1} = [6 sc];
end

% ---------- NON-UNIFORM SCALING (CELL PARAM) ----------
attack_pool{end+1} = {7 [0.8 1.2]};
attack_pool{end+1} = {7 [1.2 0.8]};

% ---------- LPF ----------
for k = [3 5 7 9]
    attack_pool{end+1} = [12 k];
end

% ---------- HISTOGRAM ----------
attack_pool{end+1} = [14 1];

% ---------- GAMMA ----------
for g = [0.6 0.8 1.2 1.4]
    attack_pool{end+1} = [15 g];
end

% ---------- COLOR QUANT ----------
for q = [8 16 32]
    attack_pool{end+1} = [16 q];
end

% ---------- NOISE ----------
for n = [0.005 0.01 0.02 0.03]
    attack_pool{end+1} = [17 n];
end

fprintf('Total attack variants = %d\n', numel(attack_pool));

%% =========================================================
% DATASET GENERATION
% =========================================================

Xbit = [];
Ybit = [];

for k = 1:numel(attack_pool)

    jenisserangan = attack_pool{k};

    attacked_img = stirmark_imagebenchmark(xwr, jenisserangan);
    xwn = reshape(double(attacked_img),[],1);

    for i = 1:length(wm)

        % ===== EXTRACTION (SAME AS BASELINE) =====
        seg = alfa * xwn(L*i-L+1:L*i);

        Sw = swt(seg, Ldek, 'haar');
        sdst = dst(Sw(sub,:));
        Sdct = dct(sdst);
        Sdct2 = reshape(Sdct,[blok blok]);

        [~,S,~] = svd(Sdct2);

        % ===== RESIDUAL =====
        D = S(1:2,1:2) - So{i}(1:2,1:2);

        % ===== FEATURES (BLIND & OBSERVABLE) =====
        corr_val = mean(D(:) .* pn_block(:));
        energy   = mean(abs(D(:)));
        meanD    = mean(D(:));
        stdD     = std(D(:));

        Xbit = [Xbit;
            corr_val, energy, meanD, stdD
        ];

        Ybit = [Ybit; wm(i)];
    end
end

%% =========================================================
% LABEL CONVERSION
% =========================================================
Ybit(Ybit==-1) = 0;   % {-1,+1} → {0,1}

%% =========================================================
% CLASS BALANCING (WAJIB)
% =========================================================
idx1 = find(Ybit==1);
idx0 = find(Ybit==0);

N = min(length(idx1), length(idx0));

idx_bal = [idx1(1:N); idx0(1:N)];
idx_bal = idx_bal(randperm(length(idx_bal)));

Xbit = Xbit(idx_bal,:);
Ybit = Ybit(idx_bal);

%% =========================================================
% SHUFFLE
% =========================================================
perm = randperm(size(Xbit,1));
Xbit = Xbit(perm,:);
Ybit = Ybit(perm);

%% =========================================================
% SAVE
% =========================================================
save dataset_bits Xbit Ybit

fprintf('Dataset size (balanced) = %d samples\n', length(Ybit));
