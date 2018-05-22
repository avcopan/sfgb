% returns value of all the basis functions (in a row) at coordinate vector R
function [value]=basisGaussianH2CO(G, R) % R is Npts x Ncoords, G are Gaussian centers
[M ndim] = size(R);
[Ng ndim] = size(G);
beta = ones(M,1)*( [0.8 1.4 1.4 1.1 1.1 1.8]/(1.35*Ng^(1/6)) );

value = zeros(M,Ng);
for k=1:Ng,
    value1 = prod( exp(-0.5*((R-ones(M,1)*G(k,:))./beta).^2)./beta, 2);   
    value(:,k) = value1;
end;

end

