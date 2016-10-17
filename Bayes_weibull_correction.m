% Function BAYES_WEIBULL_CORRECTION applies Bayesian Inference to input data
% assuming Weibull distribution

function corrected = Bayes_weibull_correction(model,dbase_real,dbase_model,f_pr)

% Essentials to numerical integration 
df = abs(f_pr(2) - f_pr(1));

length_of_pred = length(model);

% Ignore non-positive values
dbase_real_non_zero = dbase_real(dbase_real>0);

WEIBULL = wblfit(dbase_real_non_zero);
sigma_error = var(dbase_model-dbase_real);

Alpha = WEIBULL(1);
Beta = WEIBULL(2);

% Vector initialization
PRED = nan(length_of_pred,1);
for ij=1:length_of_pred
    
    INTEGRAL_up = 0;
    INTEGRAL_down = 0;
    for ll=1:length(f_pr)-1
        
        B1 = ((f_pr(ll))^Beta)*exp((-(model(ij)-f_pr(ll))^2)/(2*sigma_error))*exp(-(f_pr(ll)/Alpha)^Beta);
        B2 = ((f_pr(ll+1))^Beta)*exp((-(model(ij)-f_pr(ll+1))^2)/(2*sigma_error))*exp(-(f_pr(ll+1)/Alpha)^Beta);
        INTEGRAL_up = INTEGRAL_up + (B1+B2)*(df/2);
        
        B1 = ((f_pr(ll))^(Beta-1))*exp((-(model(ij)-f_pr(ll))^2)/(2*sigma_error))*exp(-(f_pr(ll)/Alpha)^Beta);
        B2 = ((f_pr(ll+1))^(Beta-1))*exp((-(model(ij)-f_pr(ll+1))^2)/(2*sigma_error))*exp(-(f_pr(ll+1)/Alpha)^Beta);
        INTEGRAL_down = INTEGRAL_down + (B1+B2)*(df/2);
    end
    PRED(ij) = INTEGRAL_up/INTEGRAL_down;
       
end

corrected = PRED;