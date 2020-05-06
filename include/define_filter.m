% N. Pustelnik, B. Pascal, ENS de Lyon, France
% April 2020
%---------------------------------------------

function define_filter(filter_def, dim, computation)


if ~exist('computation') %#ok<EXIST>
    computation = 'fourier';
end

if ischar(filter_def)
    if strcmp(filter_def,'gradient');
        if strcmp(computation,'fourier')
                h = [1/2, -1/2];
                H = psf2otf(h,dim);
                xt = real(ifft(fft(x).*H)); 
        elseif strcmp(opt.computation,'direct')
            xt = [x(2:dim(2))/2-x(1:dim(2)-1)/2,zeros(1,1)];
        end

    elseif strcmp(opt.filter,'laplacian');
        if strcmp(opt.computation,'fourier')
            h = [1/4, -1/2, 1/4];
            H = psf2otf(h,size(x));
            xt = real(ifft(fft(x).*H));   
        elseif strcmp(opt.computation,'direct')
            xt = [x(3:l)/4 - x(2:l-1)/2 + x(1:l-2)/4,zeros(1,2)];
        end
    end
end