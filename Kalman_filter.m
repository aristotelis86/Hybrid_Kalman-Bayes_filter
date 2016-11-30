% Function KALMAN_FILTER uses output from KALMAN_CUSTOM_TRAINING and
% applies the filter to given model forecasts.

function PRED = Kalman_filter(Tmod1,dim,XX)

hard_limit = 20; % max value of difference between original model value and corrected one

pre = XX(1,1);
if dim>=2
    for m=2:dim
        pre = pre + XX(m,1)*Tmod1^(m-1);
    end
end


pre = pre + Tmod1;

if abs(pre-Tmod1)>hard_limit
    pre = Tmod1;
end
if pre<0
    pre = 0;
end


PRED = pre;

