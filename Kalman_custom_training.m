% Function KALMAN_CUSTOM_TRAINING uses obs and mod arrays to train a Kalman
% filter assigning values to state matrices. Classic Kalman filter.
%
% Inputs:
%           Tobs: observations array
%           Tmod: model results array
%           dim: power of the polynomial for Kalman (>1, non-linear filter)
%           history_index: number of observations to be used as history
%           t: variable indicating if it is a first run(=1) or not (>1)


function XX = Kalman_custom_training(Tobs,Tmod,dim,history_index,t)

persistent P  yV x_matrix x first_run

% first_run

if isempty(first_run) || t==1
    
    init_val = 10;
    % Initialize necessary matrices
    yV = init_val*ones(1,history_index);
    x_matrix = init_val*ones(dim+1,history_index+1);
    P = init_val*ones(dim,dim);
    x = init_val*ones(dim,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    first_run = 1;
    
end


% hard_limit = 20;
ID = eye(dim);

for ij=1:length(Tobs)
    
    ss = Tmod(ij);
    z = Tobs(ij);

    y = z - ss;

    H = nan(1,dim);
    for n=1:dim
        H(1,n) = ss^(n-1);
    end

    xV = nan(dim,history_index);
    for i1=1:history_index
        for i2=1:dim

            xV(i2,i1) = x_matrix(i2,i1+1) - x_matrix(i2,i1);

        end
    end
    
    % Variance matrices
    VV = var(yV);

    W = cov(xV');


    % Kalman Gain calculation

    P1 = P + W;

    Q = H*P1*H';
    Q1 = Q(1,1);
    T = H'; 
    PR = P1*T;

    KG = PR./(Q1+VV); % Kalman Gain
    
    PR = H*x;
    D = y - PR(1,1);
    PR = KG*D;

    S = x + PR(1,:);
    
    % Update Coefficients
    PR = KG*H;
    D = ID-PR;
    PR = D*P1;

    % P_matrix update
    P = PR;

    % yV update
    yV_last = Tobs(ij) - Tmod(ij);
    for m=1:dim
        aa = Tmod(ij);
        yV_last = yV_last-x(m,1)*(aa^(m-1));
    end

    yV = yV(2:history_index);
    yV = [yV, yV_last];

    % x_matrix update
    x_matrix = x_matrix(1:dim,1:history_index);
    x_matrix = [x_matrix, x(1:dim,1)];

    % size(x_matrix)
    % x vector update
    x = S;


end
% Kalman state matrix
XX = x;
end
