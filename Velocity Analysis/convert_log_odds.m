function[raw_probability] = convert_log_odds(log_odds)

raw_probability = zeros(size(log_odds));

for idx = 1:numel(log_odds)
    odds_ratio = exp(log_odds(idx));
    raw_probability(idx) = odds_ratio;
end

for idx = 1:numel(raw_probability)
    raw_probability(idx) = raw_probability(idx)/(1 + raw_probability(idx));
end

