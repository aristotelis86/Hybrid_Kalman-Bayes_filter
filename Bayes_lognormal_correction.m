% Function BAYES_LOGNORMAL_CORRECTION applies Bayesian Inference to input data
% assuming Weibull distribution

function corrected = Bayes_lognormal_correction(model,dbase_real,dbase_model,f_pr)

% Essentials to numerical integration 
df = abs(f_pr(2) - f_pr(1));

length_of_pred = length(model);

% Ignore non-positive values
dbase_real_non_zero = dbase_real(dbase_real>0);

LOGNORMAL = lognfit(dbase_real_non_zero);
sigma_error = var(dbase_model-dbase_real);

mu = LOGNORMAL(1);
sigmaf = LOGNORMAL(2);

% Vector initialization
PRED = nan(length_of_pred,1);
for ij=1:length_of_pred
    
    INTEGRAL_up = 0;
    INTEGRAL_down = 0;
    for ll=1:length(f_pr)-1
        
        a1 = ((model(ij)-f_pr(ll))^2)/(2*sigma_error^2);
        b1 = ((log(f_pr(ll))-mu)^2)/(2*sigmaf^2);
        a12 = ((model(ij)-f_pr(ll+1))^2)/(2*sigma_error^2);
        b12 = ((log(f_pr(ll+1))-mu)^2)/(2*sigmaf^2);
        
        B1 = exp(-(a1+b1));
        B2 = exp(-(a12+b12));
        
        INTEGRAL_up = INTEGRAL_up + (B1+B2)*(df/2);
        
        B12 = exp(-(a1+b1))./f_pr(ll);
        B22 = exp(-(a12+b12))./f_pr(ll+1);
        
        INTEGRAL_down = INTEGRAL_down + (B12+B22)*(df/2);

    end
    PRED(ij) = INTEGRAL_up/INTEGRAL_down;
end

corrected = PRED;