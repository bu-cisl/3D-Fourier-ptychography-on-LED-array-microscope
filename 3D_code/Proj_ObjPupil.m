function [ O,P ] = Proj_ObjPupil(O0,P0,G,G0,Ps,alpha,beta,iters)
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

% size of O
No = size(O0); No = No(:);
% size of P, Np<=No
Np = size(P0); Np = Np(:);
cen0 = round((No+1)/2);

% init guess
O = O0;
P = P0;
it = 0;

dG = G-G0;

while (it<iters)
    % operator to put P at proper location at the O plane
    upsamp = @(x) padarray(x,(No-Np)/2);
    % operator to crop region of O from proper location at the O plane
    if mod(Np(1),2) == 1
        downsamp = @(x) x(cen0(1,m)-(Np(1)-1)/2:cen0(1,m)+(Np(1)-1)/2,...
            cen0(2,m)-(Np(2)-1)/2:cen0(2,m)+(Np(2)-1)/2);
    else
        downsamp = @(x) x(cen0(1)-Np(1)/2:cen0(1)+Np(1)/2-1,...
            cen0(2)-Np(2)/2:cen0(2)+Np(2)/2-1);
    end
    
    O = O+upsamp(abs(P).*conj(P)./(abs(P).^2+alpha).*dG)/max(abs(P(:)));
    P = P+downsamp(abs(O).*conj(O)./(abs(O).^2+beta)).*dG/max(abs(O(:))).*Ps;
    it = it+1;
end

end

