function [xfilter,lambda, objective] = linear_filtering(data,filter_def,lambda,bound)
% Perform linear filtering
%
% [xfilter,lambda,objective] = linear_filtering(data,filter_def,bound, lambda)
% data        : data to be filtered
% filter_def  : filter can either be defined as 'gradient' or 'laplacian'
%               or on a form of the PSF (i.e. [1/2, -1/2]) or the OTF
% lambda      : regularisation parameter (the larger, the smoother)
%               default value N^(3/2)*sigma/4
% bound       : (=1) to deal with boundary effects 
% xfilter     : denoised signal
% objective   : objective function value
%
% minimization problem:
% --------------------
%       min_x   1/2|| data - x ||_2^2 + lambda ||h*x||_2^2
%
% exemple: 
% --------
% >>  data = create_piecewise_cst_1D(1000,0.005,0.07);plot(data);
% >>  [xfilter,lambda] = linear_filtering(data,[1/2,-1/2],[],1);
% >>  hold on; plot(xfilter,'linewidth',2);
%
% Implementation N. PUSTELNIK, ENS Lyon
% June 2019

if ~exist('lambda')  || isempty(lambda)%#ok<EXIST>
    [~,cH] = dwt(data,'db1');
    C = abs(cH);
    sigma = median(C)/0.6745;        
    lambda = length(data)^(1+1/2) * sigma/4;
end


if ~exist('bound')  || isempty(bound)%#ok<EXIST>
    bound = 0;
end


if bound==1  %#ok<EXIST>
    data = [data,fliplr(data)];
end

dim = size(data);
if dim(1)>dim(2);
    data = data';
    dim = size(data);
end


if ~exist('filter_def') %#ok<EXIST>
    filter_def = 'gradient';
end



    
if ischar(filter_def)
    if strcmp(filter_def,'gradient');
            h = [1/2, -1/2];
            H = psf2otf(h,dim);
    elseif strcmp(filter_def,'laplacian');
            h = [1/4, -1/2, 1/4];
            H = psf2otf(h,size(data));
    else
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('!!! Change second argument in ''linear_filtering''  !!!\n');
        fprintf('!!! Possible entries are ''gradient'' or ''laplacian''!!!\n');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n');
        return;
    end
else
    if size(filter_def,2) == dim(2)
        H = filter_def;
    else
        filter_def = filter_def./sum(abs(filter_def));
        H = psf2otf(filter_def,dim); 
    end
end


xfilter0 = real(ifft( fft(data,[],2)./(1+lambda*conj(H).*H),[],2));
if bound==1  %#ok<EXIST>
    xfilter = xfilter0(1:dim(2)/2);
else
    xfilter = xfilter0;
end


objective = 1/2*sum(abs(xfilter(:)-data(:)).^2) + lambda/2*sum(abs(real(ifft(H.*fft(xfilter,[],2),[],2))).^2);

