% Hybrid Kalman-Bayes filter to correct forecasted values of a wind/wave
% parameter from a numerical model based on observation data AND former
% model results.
% 


function varargout = Hybrid_Kal_Bayes_filter(obs,model,varargin)

clear Kalman_custom_training_IG Kalman_custom_training

%% Gather all input arguments
[~, args, ~] = axescheck(varargin{:});

%% Check data size consistency
if size(obs,1)~=size(model,1) && size(obs,2)~=size(model,2)
    
    error(message({'Check data size...'; 'Observations and model predictions do not match in size.'}));
    
end

%% Detect input of further arguments
indx_dim = find(strcmpi(args,'dim'),1);
indx_history = find(strcmpi(args,'history'),1);
indx_KalmanIG = find(strcmpi(args,'kalmanIG'),1);
indx_fcst = find(strcmpi(args,'forecast'),1);
indx_pdf = find(strcmpi(args,'Distribution'),1);
indx_acc = find(strcmpi(args,'IntegSteps'),1);
indx_write = find(strcmpi(args,'File'),1);
indx_frame = find(strcmpi(args,'Limit'),1);

% Degree of polynomial for Kalman filter.
if ~isempty(indx_dim)    
    dim = args{indx_dim+1};
else
    dim = 2;
end

% Length of history from observations (training purposes, number of values)
if ~isempty(indx_history)
    history_index = args{indx_history+1};
else
    history_index = 72;
end

% Decide whether to use Kalman with Information Geometry technique or not.
if ~isempty(indx_KalmanIG)
    var_factor = args{indx_KalmanIG+1};
    use_kalIG_flag = 1;
else
    use_kalIG_flag = 0;
end

% Length of forecast (output, number of values)
if ~isempty(indx_fcst)
    length_of_fcst = args{indx_fcst+1};
else
    length_of_fcst = 24;
end

% Test discovering distridution with best fit or just define one
if ~isempty(indx_pdf)
    Dist_Name = args{indx_pdf+1};
    if strcmpi(Dist_Name,'weibull')
        Dist_name = 'weibull';
    elseif strcmpi(Dist_Name,'lognormal')
        Dist_name = 'lognormal';
    elseif strcmpi(Dist_Name,'normal')
        Dist_name = 'normal';
    end
    A = 1;
else
    A = [];
    try
        A = DiscoverDist(obs);
        Dist_name = A(1).DistName;
    catch
        
        fprintf('No Distribution detected, continuing correction using Kalman filter only... \n\n')
        Dist_name = 'none';
        
    end
    
    if ~isempty(A)
        fprintf('Distribution detected: %s from %d values. \n\n',Dist_name,length(obs))
    end
    
end

% Determine accuracy for the calculation of the integrals (number of
% divisions)
if ~isempty(indx_acc)
    ISteps = args{indx_acc+1};
else
    ISteps = 100;
end

% Create file with outputs, if needed
if ~isempty(indx_write)
    fid = fopen(args{indx_write+1},'w');
else
    fid = fopen('Hybrid_filter_results.txt','w');
end

% Create graphic with the results, if needed
xlim_switch = 0;
if ~isempty(indx_frame)
    XLIM = args{indx_frame+1};
    if length(XLIM)>2
        error(message({'Check figure limit size...'; 'Too many inputs in array, should be only 2...'}));
    else 
        if max(XLIM(:))>length(obs) && min(XLIM(:))<1
            error(message({'Check figure limit indexes...'; 'Values exceed data dimensions...'}));
        end
    end
    xlim_switch = 1;
end
    
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Last size check before calculations 
if length(obs)<(history_index+length_of_fcst)
    
    error(message({'Observation data are not enough to initiate calculations'; 'Check values of parameters.'}))

end



%% Initialize output array

TT = 1:length_of_fcst:length(model)-(history_index+length_of_fcst)+1;
TJ = length(history_index+1:TT(end)+history_index+length_of_fcst-1);
Hyb_Corr = nan(TJ,1);



%% Apply Kalman/KalmanIG Filter

for t=TT
    
    Tobs = obs(t:t+history_index-1);
    Tmod1 = model(t:t+history_index-1);
    Tmod2 = model(t+history_index:t+history_index+length_of_fcst-1);
    
    if use_kalIG_flag
        XX = Kalman_custom_training_IG(Tobs,Tmod1,dim,history_index,t,var_factor); % Advanced Training Algorithm based on IG
    else
        XX = Kalman_custom_training(Tobs,Tmod1,dim,history_index,t); % Training Algorithm
    end

    Tmod2_kal = nan(length_of_fcst,1);
    for i=1:length_of_fcst

        Tmod_in = Tmod2(i);

        PRED = Kalman_filter(Tmod_in,dim,XX);
        Tmod2_kal(i) = PRED;


    end
    

%% Apply Bayesian Inference Scheme
    if ~isempty(A)
        
        % Calculate Integration Intervals for Bayesian Technique
        f_pr = min(obs(:)):((max(obs(:))-min(obs(:)))/ISteps):max(obs(:));
        if f_pr(1)==0
            f_pr = f_pr(2:end);
        end
        
        switch Dist_name

            case 'normal'

                Tmod_hybrid = Bayes_normal_correction(Tmod2_kal,Tobs,Tmod1,10);

                Hyb_Corr(t:t+length_of_fcst-1) = Tmod_hybrid(:);


            case 'weibull'

                Tmod_hybrid = Bayes_weibull_correction(Tmod2_kal,Tobs,Tmod1,f_pr);
                Hyb_Corr(t:t+length_of_fcst-1) = Tmod_hybrid(:);


            case 'lognormal'

                Tmod_hybrid = Bayes_lognormal_correction(Tmod2_kal,Tobs,Tmod1,f_pr);
                Hyb_Corr(t:t+length_of_fcst-1) = Tmod_hybrid(:);


            case 'extreme value'

                Hyb_Corr(t:t+length_of_fcst-1) = Tmod2_kal(:);


        end
    else
        
        % In case NO distribution is detected
        Hyb_Corr(t:t+length_of_fcst-1) = Tmod2_kal(:);
    
    end
    
end

OBS = obs(history_index+1:t+(history_index+length_of_fcst)-1);
ORIG_MODEL = model(history_index+1:t+(history_index+length_of_fcst)-1);


% File output section
fprintf(fid,'OBS \t MODEL_original \t MODEL_corrected \n');
fprintf(fid,'%f \t %f \t %f \n',[OBS, ORIG_MODEL, Hyb_Corr]');

fclose(fid);


% Figure output section
HF = figure;

hold on 
plot(OBS,'--+b')
plot(ORIG_MODEL,'r')
plot(Hyb_Corr,'-og')

ylabel('parameter values')
xlabel('index')
title('Comparison Figure')

if xlim_switch
    set(gca,'XLim',XLIM)
else
    set(gca,'XLim',[1 length(OBS)])
end
legend('obs','original model','corrected model','Location','best');

[bias, rmse, NS] = Datasets_statistics(OBS,Hyb_Corr);
xpos = min(get(gca,'XLim'));
ypos = max(get(gca,'YLim'));
my_text = {sprintf('OBS-Corrected Model');sprintf('BIAS=%.2f',bias); sprintf('RMSE=%.2f',rmse); sprintf('NS=%.2f',NS)};
text(xpos,ypos,my_text,'VerticalAlignment','top')

[bias, rmse, NS] = Datasets_statistics(OBS,ORIG_MODEL);
xpos = max(get(gca,'XLim'));
ypos = max(get(gca,'YLim'));
my_text = {sprintf('OBS-Original Model');sprintf('BIAS=%.2f',bias); sprintf('RMSE=%.2f',rmse); sprintf('NS=%.2f',NS)};
text(xpos,ypos,my_text,'VerticalAlignment','top','HorizontalAlignment','right')

print(HF,'-dpng','-r250','Comparison_figure.png')



varargout{1} = Hyb_Corr;
varargout{2} = ORIG_MODEL;
varargout{3} = OBS;
varargout{4} = Dist_name;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function for statistics

function [BIAS, RMSE, NS] = Datasets_statistics(obs,model)

BIAS = obs - model;

try
    RMSE = rms(BIAS); 
catch
    RMSE = (mean(BIAS.^2))^.5;
end

NS = 1-((sum(BIAS.^2))/(sum((obs-mean(obs)).^2)));

BIAS = mean(BIAS);

