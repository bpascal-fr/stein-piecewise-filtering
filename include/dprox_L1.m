function dwp = dprox_L1(wx,dwx,gamma)
% Differential of the proximity operator of the l1 norm with respect to its
% first argument wx, applied to dwx
% dwp = (d/dwx) prox_{gamma || .||_1}(wx)[dwx]
% 
% Implementation B. PASCAL, ENS Lyon
% April 2020

        
        ny = abs(wx);        
        dwp = zeros(size(dwx));

        dproj = (dwx(ny >gamma).*wx(ny >gamma))./ny(ny >gamma).^2;
        projy = dwx(ny >gamma) -   dproj.*wx(ny >gamma);
        
        dwp(ny >gamma) = dwx(ny >gamma) - gamma./ny(ny >gamma).*projy;
        
end