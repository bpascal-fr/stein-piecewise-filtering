function nL = normL2(x,filter_def, computation,param,op)
% define the norm of the linear operator
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

dim = size(x);
if strcmp(param.type,'1D')

    if ~exist('computation') %#ok<EXIST>
        computation = 'fourier';
    end

    if ~exist('filter_def') %#ok<EXIST>
        filter_def = 'gradient';
    end

    if ischar(filter_def)
        if strcmp(filter_def,'gradient');
            if numel(param.lambda) == 1
                if strcmp(computation,'fourier') 
                    nL = param.lambda^2;
                else
                    nL = param.lambda^2;
                end
            else
                nL = max(param.lambda)^2;
            end
        elseif strcmp(filter_def,'laplacian');
            if numel(param.lambda) == 1
                if strcmp(computation,'fourier')
                    nL = param.lambda^2; 
                else
                    nL = param.lambda^2;
                end
            else
                nL = max(param.lambda)^2;
            end
        end
    else
        nL = param.lambda^2;
    end
else
    if ~exist('computation') %#ok<EXIST>
        computation = 'fourier';
    end

    if ~exist('filter_def') %#ok<EXIST>
        filter_def = 'gradient';
    end

    if ischar(filter_def)
        if strcmp(filter_def,'gradient');
            if length(param.lambda)==1
                nL = (param.lambda*sqrt(2))^2;
            else
                nL = powermethod(x,op);
            end
            

        elseif strcmp(filter_def,'laplacian');
            if strcmp(computation,'fourier')
                nL = param.lambda^2;
                
            end
        end
    else
         nL = param.lambda^2; 
    end
end