function nL = normL(x,filter_def, computation,param,op)
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
                    nL = param.lambda;
                else
                    nL = param.lambda;
                end
            else
                nL = max(param.lambda);
            end
        elseif strcmp(filter_def,'laplacian');
            if numel(param.lambda) == 1
                if strcmp(computation,'fourier')
                    nL = param.lambda; 
                else
                    nL = param.lambda;
                end
            else
                nL = max(param.lambda);
            end
        end
    else
        nL = param.lambda;
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
                nL = param.lambda*sqrt(2);
            else
                nL = powermethod(x,op);
            end
            

        elseif strcmp(filter_def,'laplacian');
            if strcmp(computation,'fourier')
                nL = param.lambda;
                
            end
        end
    else
        if size(filter_def,2) == dim(2)
            H = filter_def;
        else
            filter_def = param.lambda*filter_def./sum(abs(filter_def(:)));
            H = psf2otf(filter_def,[dim(1) dim(2)]); 
        end
        xt = real(ifft2(fft2(x).*H));  
    end
end