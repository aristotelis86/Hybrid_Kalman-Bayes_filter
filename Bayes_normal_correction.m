% Function BAYES_NORMAL_CORRECTION applies Bayesian Inference to input data
% assuming Normal distribution

function corrected = Bayes_normal_correction(model,dbase_real,dbase_record,varargin)

model_length = length(model);

mean_dbase_real = mean(dbase_real);
sigma_dbase_real = var(dbase_real);
% mean_error = mean(dbase_record-dbase_real);
sigma_error = var(dbase_record-dbase_real);



PRED = nan(model_length,1);
% Create corrected predictions
for ij=1:length(model)
    % Define Distribution Parameters
    sigma_predic = ((1/sigma_error)+(1/sigma_dbase_real))^(-1); % Variance of temperature for predicted values
    mean_val_predic = ((1/sigma_error)*model(ij)+(1/sigma_dbase_real)*mean_dbase_real)*sigma_predic;
    
    
    % Corrected values!
%     Interm=zeros(accuracy,1);
%     for t=1:accuracy
%         Interm(t) = random('norm',mean_val_predic,sqrt(sigma_predic));
%     end
%     PRED(ij) = mean(Interm);
    PRED(ij) = mean_val_predic;
end


corrected = PRED;