function [data,truth] = create_piecewise_cst_1D_multivariate(M,N,Ep,sig,Nval)
    
    % Generate noisy piecewise constant multivariate signals
    % The hops are placed randomly.
    % The number of hops is random but likely to be around Ep.
    % On each tray the value of data is randomly selected among 10
    % regularly spaced values.
    % The hops are common to all components but the value of the trays for
    % each components are independant random variables.
    %
    % Inputs:   M number of components
    %           N length of signals
    %           Ep expected number of hops (default 10)
    %           sig standard deviation of the noise (default 1e-1) possibly
    %           one noise level per component
    %           Nval expected number of differents values taken by the
    %           signal (default 10) possibly one per component
    %
    % Outputs:  data noisy multivariate signals
    %           underlying ground truth
    %
    % Implemented by B. PASCAL, ENS de Lyon
    % April 2020
    
    if nargin < 5
        Nval = [10, 10, 10];
        if nargin < 4
            sig = 1e-1;
            if nargin < 3
                Ep = 10;
            end
        else
            if numel(Nval) == 1
                Nval = Nvel*ones(1,M);
            elseif ~numel(Nval) == M
                disp('ERROR: number of different values differs from number of component')
            end
        end
        
    end
    
    
    p = Ep/N;
    indic_hop = rand(1,N) < p;
    N_hop = sum(indic_hop);
    trays = cumsum(indic_hop);
    
    truth = zeros(M,N);
    for nh = 1:N_hop
        rnd_val = zeros(M,1);
        for m = 1:M
            rnd_val(m) = floor(Nval(m)*rand(1,1));
        end
        rnd_kron = kron(rnd_val,ones(1,sum(trays == nh)));
        truth(:,trays == nh) = rnd_kron;
    end
    
    for m = 1:M
        truth(m,:) = (truth(m,:) - min(truth(m,:)))/(max(truth(m,:)) - min(truth(m,:)));
    end
    if numel(sig) == 1 || numel(sig) == M
        data = truth + diag(sig)*randn(M,N);
    else
        error('Number of noise levels different from number of component')
    end
    
end
