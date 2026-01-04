clc; clear; close all;

%% =========================================================
% LOAD HASIL GRID SEARCH
% =========================================================
load best_parameters_all_attacks   % berisi struct: Results

N = length(Results);
fprintf('Total serangan dianalisis = %d\n\n', N);

%% =========================================================
% CEK STRUKTUR (AMAN)
% =========================================================
assert(isfield(Results,'attack'),'Field attack tidak ada');
assert(isfield(Results,'bestBER'),'Field bestBER tidak ada');
assert(isfield(Results,'param'),'Field param tidak ada');

%% =========================================================
% AMBIL SEMUA PARAMETER TERBAIK
% =========================================================
subs    = zeros(N,1);
bloks   = zeros(N,1);
alfas   = zeros(N,1);
alfasss = zeros(N,1);
bers    = zeros(N,1);

for i = 1:N
    subs(i)    = Results(i).param.sub;
    bloks(i)   = Results(i).param.blok;
    alfas(i)   = Results(i).param.alfa;
    alfasss(i) = Results(i).param.alfass;
    bers(i)    = Results(i).bestBER;
end

%% =========================================================
% TAMPILKAN RINGKASAN GLOBAL
% =========================================================
fprintf('================ PARAMETER GLOBAL TERBAIK ================\n');
fprintf('sub     (mode) = %d\n', mode(subs));
fprintf('blok    (mode) = %d\n', mode(bloks));
fprintf('alfa    (mode) = %.3f\n', mode(alfas));
fprintf('alfass  (mode) = %.3f\n', mode(alfasss));
fprintf('==========================================================\n\n');

fprintf('BER minimum = %.4f\n', min(bers));
fprintf('BER maksimum = %.4f\n', max(bers));
fprintf('BER rata-rata = %.4f\n\n', mean(bers));

%% =========================================================
% SERANGAN TERMUDAH & TERSULIT
% =========================================================
[~, idx_best]  = min(bers);
[~, idx_worst] = max(bers);

fprintf('================ SERANGAN TERMUDAH =================\n');
disp(Results(idx_best))

fprintf('================ SERANGAN TERSULIT =================\n');
disp(Results(idx_worst))

%% =========================================================
% TABEL RINGKAS (SIAP EKSPOR)
% =========================================================
AttackID = zeros(N,1);
for i = 1:N
    atk = Results(i).attack;
    if isnumeric(atk)
        AttackID(i) = atk(1);
    else
        AttackID(i) = atk{1};
    end
end

T = table( ...
    (1:N)', AttackID, subs, bloks, alfas, alfasss, bers, ...
    'VariableNames', ...
    {'AttackIndex','AttackID','Subband','Block','Alfa','AlfaSS','BestBER'} ...
);

disp('================ TABEL PARAMETER TERBAIK ================');
disp(T)

%% =========================================================
% SIMPAN KE FILE
% =========================================================
writetable(T,'parameter_terbaik_per_serangan.csv');
save parameter_summary.mat T subs bloks alfas alfasss bers

fprintf('\nFile disimpan:\n');
fprintf('- parameter_terbaik_per_serangan.csv\n');
fprintf('- parameter_summary.mat\n');
