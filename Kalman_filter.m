% Function KALMAN_FILTER uses output from KALMAN_CUSTOM_TRAINING or
% KALMAN_CUSTOM_TRAINING_IG and applies the filter to given model forecasts.
%
% Inputs:
%       Tmod: values forecasted by model
%       dim: rank of polynomial for Kalman correction
%       XX: Kalman state matrix

function PRED = Kalman_filter(Tmod,dim,XX)
% max value of difference between original model value and corrected one
% (suitable for wind speed corrections)
hard_limit = 20; 

pre = XX(1,1);
if dim>=2
    for m=2:dim
        pre = pre + XX(m,1)*Tmod^(m-1);
    end
end


pre = pre + Tmod;

if abs(pre-Tmod)>hard_limit
    pre = Tmod;
end
if pre<0
    pre = 0;
end

% Predictions corrected by Kalman filter
PRED = pre;

