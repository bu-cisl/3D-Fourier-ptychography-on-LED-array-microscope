function [ O,phi ] = Proj_OslicePhi(O0,phi0,psi,psi0,alpha,beta,iters)
%GDUPDATE_MULTIPLICATION update estimate of O and P according to gradient
%descent method, where psi = O*P
%   Inputs:
%   O0: object estimate, n1xn2
%   P0: pupil function estimate: m1xm2
%   psi: update estimate field estimate
%   psi0: previous field estimate
%   alpha: gradient descent step size for O
%   betta: gradient descent step size for P
%   Ps: support constraint for P0, e.g. spatially confined probe or
%   objective with known NA
%   iters: # of iterations to run on updates
%
% last modified by Lei Tian, lei_tian@alum.mit.edu, 5/27/2014

% init guess
O = O0;
phi = phi0;
it = 0;

dpsi = psi-psi0;
while (it<iters)    
    O = O+abs(phi).*conj(phi)./(abs(phi).^2+alpha).*dpsi/max(abs(phi(:)));
    phi = phi+abs(O).*conj(O)./(abs(O).^2+beta).*dpsi/max(abs(O(:)));
    it = it+1;
end

end

