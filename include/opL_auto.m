function xt = opL_auto(x,filter_def, computation,param, lambda)
% Definition of the operator used in the penalization
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
        if strcmp(filter_def,'gradient')
            if numel(lambda) == 1
                if strcmp(computation,'fourier') 
                    h = lambda*[1/2, -1/2];
                    H = psf2otf(h,dim);
                    xt = real(ifft(fft(x,[],2).*H,[],2)); 
                else
                    xt = lambda*[x(:,2:dim(2))/2-x(:,1:dim(2)-1)/2,zeros(dim(1),1)];
                end
            else
                xt = diag(lambda)*[x(:,2:dim(2))/2-x(:,1:dim(2)-1)/2,zeros(dim(1),1)];
            end
        elseif strcmp(filter_def,'laplacian')
            if numel(lambda) == 1
                if strcmp(computation,'fourier')
                    h = lambda*[1/4, -1/2, 1/4];
                    H = psf2otf(h,size(x));
                    xt = real(ifft(fft(x,[],2).*H,[],2));   
                else
                    xt = lambda*[x(:,3:dim(2))/4 - x(:,2:dim(2)-1)/2 + x(:,1:dim(2)-2)/4,zeros(dim(1),2)];
                end
            else
                xt = diag(lambda)*[x(:,3:dim(2))/4 - x(:,2:dim(2)-1)/2 + x(:,1:dim(2)-2)/4,zeros(dim(1),2)];
            end
        end
    else
        if numel(lambda) == 1
            filter_def = lambda*filter_def./sum(abs(filter_def));
            H = psf2otf(filter_def,dim); 
            xt = real(ifft(fft(x,[],2).*H,[],2));  
        else
            disp('The reglarization parameter \lambda should be of length 1\n');
        end
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
                lambda = [lambda,lambda];
            end
            xt = zeros(dim(1), dim(2),2);
            if strcmp(computation,'fourier')
                    h = lambda(1)*[1/2, -1/2];
                    H = psf2otf(h,[dim(1) dim(2)]);
                    v = lambda(2)*[1/2, -1/2]';
                    V = psf2otf(v,[dim(1) dim(2)]);                    
                    xt(:,:,1) = real(ifft2(fft2(x).*H)); 
                    xt(:,:,2) = real(ifft2(fft2(x).*V));                     
            else
                    tau = 1;
                    [n,m]=size(x);
                    xt(:,:,1) = lambda(1)*[x(:,(1+tau):m)/2-x(:,1:m-tau)/2,zeros(n,tau)];
                    x = x';
                    [n,m]=size(x);
                    xt(:,:,2) = lambda(2)*[x(:,(1+tau):m)/2-x(:,1:m-tau)/2,zeros(n,tau)]';

            end

        elseif strcmp(filter_def,'laplacian')
            if strcmp(computation,'fourier')
                h = lambda*[1/4, -1/2, 1/4]'*[1/4, -1/2, 1/4];
                H = psf2otf(h,size(x));
                xt = real(ifft2(fft2(x).*H));   
            %else
            %    xt = lambda[x(:,3:dim(2))/4 - x(:,2:dim(2)-1)/2 + x(:,1:dim(2)-2)/4,zeros(dim(1),2)];
            end
        end
    else
        filter_def = lambda*filter_def./sum(abs(filter_def(:)));
        H = psf2otf(filter_def,[dim(1) dim(2)]); 
        xt = real(ifft2(fft2(x).*H));  
    end
end
