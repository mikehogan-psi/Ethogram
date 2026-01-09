function[extracted_features] = extract_behavioural_features(Xfit, b, R, T, missing, behaviour)
%EXTRACT_BEHAVIOURAL_FEATURES: Extracts features of interest from
%SSM-fitted data for training/applying a behavioural predictive model
%
% This currently accepts 'Grooming', 'Rearing' and 'Darting' as inputs to
% the behaviour field. To add additional behaviours, select/add the
% appropriate features and add these as a feature list index under the new
% behaviour name.
%
% Outputs:
%   -extracted_features: a frame x feature array containing values for each
%   feature for every frame

num_frames = size(b, 2);

% List of potential features use for feature extraction
feature_strings = {
 
'mean(b_window(3, :))';                          % Feature 1:  Average "elongation/hunching" (General horizontal body deformation )
'var(b_window(3, :))';                           % Feature 2:  Variability in "elongation/hunching"  ( horizontal deformation dynamism / general fluctuation)
'mean(diff(b_window(3, :)))';                    % Feature 3:  Trend of change "elongation/hunching"  posture (rate and  direction of horizontal deformation dynamism)
'mean(b_window(1, :))';                          % Feature 4:  Average "bending side to side" (general lateral deformation) 
'var(b_window(1, :))';                           % Feature 5:  Variability in "bending side to side" (lateral motion dynamism)
'mean(diff(b_window(1, :)))';                    % Feature 6:  Trend of change in "bending side to side" (rate and  direction of lateral motion dynamism)
'mean(b_window(2, :))';                          % Feature 7:  Average "looking up/down" (General vertical posture trend)
'var(b_window(2, :))';                           % Feature 8:  Variability in "looking up/down" posture (vertical motion dynamism / general fluctuation)
'mean(diff(b_window(2, :)))';                    % Feature 9:  Trend of change in "looking up/down" posture (rate and  direction of vertical motion dynamism)

'mean([0, sqrt(sum(diff(T_window(1:2,:), 1, 2).^2, 1))])';                           % Feature 10: Average 2D velocity - ONLY USING X AND Y COODINATES !!!
'mean([0, abs(diff(R_window(1,:)))])'                                                % Feature 11: Average Angular Velocity / Average yaw (Describes rotational movement around z-axis)
'var([0, diff([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])])'                 % Feature 12: Variability in acceleration of movement / How much acceleration fluctuates (2D)
'sum(abs([0, diff([0, diff([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])])]))' % Feature 13: Cumulative change in acceleration (2D)  
'max([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                            % Feature 14: Peak velocity within the window (2D)
'sum(abs(diff([0, diff(b_window(1, :))])))'                                          % Feature 15: Cumulative change in rate of bending
'sum([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                            % Feature 16: Total distance travelled (2D)
'mean([0, abs(diff(R_window(2,:))) + abs(diff(R_window(3,:)))])';                    % Feature 17: average angular velocity (magnitude only) in pitch and roll directions combined (Describes rotational movement around x- and y-axis)


'Xfit(4,3,frame)-Xfit(11,3,frame)'                                         % Feature 18: z-difference between nose (body-marker 4) and tail-anterior (body-marker 11)
'mean([0, sqrt(sum(diff(T_window(1:2,:), 1, 2).^2, 1))])';                         % Feature 19: Average velocity (overall movement intensity) - ONLY USING X AND Y COODINATES
'sum([0, sqrt(sum(diff(T_window(1:2, :), 1, 2).^2, 1))])'                          % Feature 20: Total distance travelled, using XY coords only

'sum(abs(diff(T_window(3, :))))'                	   % Feature 21: cumulative vertical movement 
'Xfit(7,3,frame)-Xfit(10,3,frame)'             % Feature 22: z-difference between neck base (body-marker 7) and tail-base (body-marker 10)
'norm(Xfit(7,1:2,frame) - Xfit(10,1:2,frame))' % Feature 23: euclidian distance in x-y plane bewteen neck base and tail base 
'Xfit(4,3,frame)'                                  % Feature 24: z-value of nose
'max(T_window(3, :))'                                  % Feature 25: max vertical translation
'Xfit(1,3,frame)'                                  % Feature 26: z-value of cable tip
'Xfit(4,3,frame)-Xfit(7,3,frame)'              % Feature 27: z-difference between nose (body-marker 4) and tail-anterior (body-marker 11) - implanted animals

 };

% Determine which behaviour is being extracted and select appropriate
% features
if strcmp(behaviour, 'Grooming')
    feature_strings_idx = [2 5 9 24 21 27 20];
    window_size = 30;

elseif strcmp(behaviour, 'Rearing')   
    feature_strings_idx = [21, 24, 25, 26, 27];
    window_size = 4;

elseif strcmp(behaviour, 'Darting')
    feature_strings_idx = [11, 12, 14, 15, 16];
    window_size = 26;

else
    warning('This behaviour is not recognised')
    extracted_features = [];
    return
end

feature_selection = feature_strings(feature_strings_idx);
extracted_features = NaN(num_frames, length(feature_selection));   

% Extract features for appropriate window size for behaviour
for frame = 1:num_frames
    
    % Sort data into windows
    start_idx = max(1, frame - (window_size / 2));
    end_idx = min(num_frames, frame + (window_size / 2));
    
    b_window = b(:, start_idx:end_idx);
    T_window = T(:, start_idx:end_idx);
    R_window = R(:, start_idx:end_idx);
    missing_window = missing(:, start_idx:end_idx);

    if sum(missing_window(:)) > 6
        continue % too many missing points → leave as NaN
    end

    valid_cols = ~any(isnan([b_window; T_window; R_window]), 1);
    b_window = b_window(:, valid_cols);
    T_window = T_window(:, valid_cols);
    R_window = R_window(:, valid_cols);


    % Compute all features for this frame  
    for feature = 1:length(feature_selection)
        try
        extracted_features(frame, feature) = eval(feature_selection{feature});
        catch
        extracted_features(frame, feature) = NaN;
        end
    end

end

end