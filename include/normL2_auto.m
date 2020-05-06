function nL = normL2_auto(x, filter_def, computation,param, op, lambda)
% define the norm of the linear operator
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019


if strcmp(param.type,'1D')

    if ~exist('computation') %#ok<EXIST>
        computation = 'fourier';
    end

    if ~exist('filter_def') %#ok<EXIST>
        filter_def = 'gradient';
    end

    if ischar(filter_def)
        if strcmp(filter_def,'gradient')
            if numel(lambda) == 1
                if strcmp(computation,'fourier') 
                    nL = lambda^2;
                else
                    nL = lambda^2;
                end
            else
                nL = max(lambda)^2;
            end
        elseif strcmp(filter_def,'laplacian')
            if numel(lambda) == 1
                if strcmp(computation,'fourier')
                    nL = lambda^2; 
                else
                    nL = lambda^2;
                end
            else
                nL = max(lambda)^2;
            end
        end
    else
        nL = lambda^2;
    end
else
    if ~exist('computation') %#ok<EXIST>
        computation = 'fourier';
    end

    if ~exist('filter_def') %#ok<EXIST>
        filter_def = 'gradient';
    end

    if ischar(filter_def)
        if strcmp(filter_def,'gradient')
            if length(lambda)==1
                nL = (lambda*sqrt(2))^2;
            else
                nL = powermethod(x,op);
            end
            

        elseif strcmp(filter_def,'laplacian')
            if strcmp(computation,'fourier')
                nL = lambda^2;
                
            end
        end
    else
         nL = lambda^2; 
    end
end