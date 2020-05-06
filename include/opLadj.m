function x = opLadj(y,filter_def, computation,param)
% Adjoint of the operator used in the penalization
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

if strcmp(param.type,'1D')
    dim = size(y);

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
                        h = param.lambda*[1/2, -1/2];
                        H = psf2otf(h,dim);
                        x = real(ifft(fft(y,[],2).*conj(H),[],2)); 
                else
                    x = -param.lambda*([y(:,1:1)/2,y(:,2:dim(2)-1)/2- y(:,1:dim(2)-2)/2,-y(:,dim(2)-2+1:dim(2)-1)/2]);
                end
            else
                    x = -diag(param.lambda)*([y(:,1:1)/2,y(:,2:dim(2)-1)/2- y(:,1:dim(2)-2)/2,-y(:,dim(2)-2+1:dim(2)-1)/2]);
            end
        elseif strcmp(filter_def,'laplacian');
            if numel(param.lambda) == 1
                if strcmp(computation,'fourier')
                    h = param.lambda*[1/4, -1/2, 1/4];
                    H = psf2otf(h,size(y));
                    x = real(ifft(fft(y,[],2).*conj(H),[],2));   
                else
                    x = param.lambda*([y(:,1:1)/4,-y(:,1)/2 + y(:,2)/4, y(:,3:dim(2)-2)/4 - y(:,2:dim(2)-3)/2 + y(:,1:dim(2)-4)/4, y(:,dim(2)-3)/4 - y(:,dim(2)-2)/2, y(:,dim(2)-2)/4]);
                end
            else
                x = diag(param.lambda)*([y(:,1:1)/4,-y(:,1)/2 + y(:,2)/4, y(:,3:dim(2)-2)/4 - y(:,2:dim(2)-3)/2 + y(:,1:dim(2)-4)/4, y(:,dim(2)-3)/4 - y(:,dim(2)-2)/2, y(:,dim(2)-2)/4]);

            end
        end
    else
        filter_def = param.lambda*filter_def./sum(abs(filter_def));
        H = psf2otf(filter_def,dim); 
        x = real(ifft(fft(y,[],2).*conj(H),[],2));  
    end
else
    dim = size(y);

    if ~exist('computation') %#ok<EXIST>
        computation = 'fourier';
    end

    if ~exist('filter_def') %#ok<EXIST>
        filter_def = 'gradient';
    end

    if ischar(filter_def)
        if strcmp(filter_def,'gradient');
            if length(param.lambda)==1
                param.lambda = [param.lambda,param.lambda];
            end
            if strcmp(computation,'fourier')
                    h = param.lambda(1)*[1/2, -1/2];
                    H = psf2otf(h,[dim(1) dim(2)]);
                    v = param.lambda(2)*[1/2, -1/2]';
                    V = psf2otf(v,[dim(1) dim(2)]);                       
                    x = real(ifft2(fft2(y(:,:,1)).*conj(H))) + real(ifft2(fft2(y(:,:,2)).*conj(V))); 
            else
                    tau = 1;
                    yt = y(:,:,1);
                    [~,m] = size(yt);
                    x = -param.lambda(1)*([yt(:,1:tau)/2,yt(:,1+tau:m-tau)/2- yt(:,1:m-2*tau)/2,-yt(:,m-2*tau+1:m-tau)/2]);
                    yt = y(:,:,2)';
                    [~,m] = size(yt);
                    x = x -param.lambda(2)*([yt(:,1:tau)/2,yt(:,1+tau:m-tau)/2- yt(:,1:m-2*tau)/2,-yt(:,m-2*tau+1:m-tau)/2])';
            end

        elseif strcmp(filter_def,'laplacian');
            if strcmp(computation,'fourier')
                h = param.lambda*[1/4, -1/2, 1/4]'*[1/4, -1/2, 1/4];
                H = psf2otf(h,[dim(1), dim(2)]);
                x = real(ifft2(fft2(y).*conj(H)));   
            %else
            %    x = ([y(:,1:1)/4,-y(:,1)/2 + y(:,2)/4, y(:,3:dim(2)-2)/4 - y(:,2:dim(2)-3)/2 + y(:,1:dim(2)-4)/4, y(:,dim(2)-3)/4 - y(:,dim(2)-2)/2, y(:,dim(2)-2)/4]);
            end
        end
    else
        filter_def = param.lambda*filter_def./sum(abs(filter_def(:)));
        H = psf2otf(filter_def,[dim(1), dim(2)]); 
        x = real(ifft2(fft2(y,[],2).*conj(H)));  
    end
end