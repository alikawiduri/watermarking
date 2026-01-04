clc; clear; close all;

%% =========================================================
% LOAD DATASET
% =========================================================
load dataset_bits    % Xbit (Nx4), Ybit (Nx1 {0,1})

fprintf('Total samples = %d\n', size(Xbit,1));

%% =========================================================
% SAFETY CHECK
% =========================================================
assert(size(Xbit,2)==5, 'Feature dimension must be 4');
assert(all(ismember(unique(Ybit),[0 1])), 'Labels must be {0,1}');

%% =========================================================
% NORMALIZATION (GLOBAL, SHARED)
% =========================================================
mu = mean(Xbit,1);
sg = std(Xbit,[],1);
sg(sg==0) = 1;

Xn = (Xbit - mu) ./ sg;

%% =========================================================
% TRAIN MODELS
% =========================================================
models = struct();

% ---------------------------------------------------------
% 1. SVM (RBF)
% ---------------------------------------------------------
models.svm = fitcsvm(Xn,Ybit,...
    'KernelFunction','rbf',...
    'KernelScale','auto',...
    'Standardize',false,...
    'ClassNames',[0 1]);

% ---------------------------------------------------------
% 2. Random Forest
% ---------------------------------------------------------
models.rf = TreeBagger(150,Xn,Ybit,...
    'Method','classification',...
    'OOBPrediction','on',...
    'MinLeafSize',5);

% ---------------------------------------------------------
% 3. KNN
% ---------------------------------------------------------
models.knn = fitcknn(Xn,Ybit,...
    'NumNeighbors',5,...
    'Distance','euclidean',...
    'Standardize',false);

% ---------------------------------------------------------
% 4. Decision Tree
% ---------------------------------------------------------
models.tree = fitctree(Xn,Ybit,...
    'MinLeafSize',5);

% ---------------------------------------------------------
% 5. Logistic Regression
% ---------------------------------------------------------
models.logreg = fitclinear(Xn,Ybit,...
    'Learner','logistic',...
    'Regularization','ridge',...
    'Lambda',1e-3);

% ---------------------------------------------------------
% 6. Neural Network
% ---------------------------------------------------------
models.nn = fitcnet(Xn,Ybit,...
    'LayerSizes',10,...
    'Activations','relu');

% ---------------------------------------------------------
% 7. AdaBoost
% ---------------------------------------------------------
models.ada = fitcensemble(Xn,Ybit,...
    'Method','AdaBoostM1',...
    'NumLearningCycles',100);

% ---------------------------------------------------------
% 8. Naive Bayes
% ---------------------------------------------------------
models.nb = fitcnb(Xn,Ybit,...
    'DistributionNames','normal');

%% =========================================================
% SAVE
% =========================================================
save bit_corrector_models models mu sg

disp('====================================');
disp(' Training 8 models DONE successfully ');
disp('====================================');
