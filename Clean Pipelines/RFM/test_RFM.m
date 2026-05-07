function[accuracy_score, d_prime_score, precision_score, recall_score, f1_score,...
    balacc_score, mcc_score, prauc_score, rocauc_score, brier_score] = ...
    test_RFM(p1, y_pred, y_test)

% Calculate accuracy
accuracy_score = sum(y_pred == y_test) / length(y_test);

% Confusion counts
hits  = sum(y_pred == 1 & y_test == 1);
miss  = sum(y_pred == 0 & y_test == 1);
fa    = sum(y_pred == 1 & y_test == 0);
cr    = sum(y_pred == 0 & y_test == 0);

precision = hits / max(hits+fa,1);
recall    = hits / max(hits+miss,1);
spec      = cr / max(cr+fa,1);
f1        = 2*precision*recall / max(precision+recall, eps);
balacc    = 0.5*(recall + spec);
mcc       = (hits*cr - fa*miss) / sqrt(max((hits+fa)*(hits+miss)*(cr+fa)*(cr+miss), eps));
brier = mean((p1 - y_test).^2);

% PR-AUC + ROC-AUC from probabilities
[rec_curve, prec_curve, ~, prAUC] = perfcurve(y_test, p1, 1, 'xCrit','reca','yCrit','prec');
[~,~,~, rocAUC] = perfcurve(y_test, p1, 1);

nSignal = hits + miss;   % number of positives in test
nNoise  = fa + cr;       % number of negatives in test

% Compute hit/false alarm rates with log-linear (Hautus) correction
hit_rate = (hits + 0.5) / (nSignal + 1);
false_alarm_rate = (fa + 0.5) / (nNoise + 1);
   
% Compute d-prime
d_prime_score = (norminv(hit_rate) - norminv(false_alarm_rate));

precision_score = precision;
recall_score    = recall;
f1_score      = f1;
balacc_score    = balacc;
mcc_score     = mcc;
prauc_score     = prAUC;
rocauc_score    = rocAUC;
brier_score     = brier;

end