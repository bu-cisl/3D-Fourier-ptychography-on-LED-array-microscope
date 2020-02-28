function [ o_slice ] = Back_Prop_MultiSlice_v2( O, k2, dz, o_slice0, phi0, psi0, ...
                                            i0, alpha, beta, iters)
%FWD_PROP_MULTISLICE computes the field using multislice approach, with
%propagator H
% Inputs:
%   O: total object field (from multi-slice propagation) at the pupil plane
%   H: fwd propagator between slices
%   H0: fwd propagator from Nth slice to focal plane of objective
%   o_slice0: current estimate of multi-slice object
%   phi0: current estimate of incident field at each slice
%   psi0: current estimate of output field at each slice
%   i0: KNOWN illumination

% by Lei Tian (lei_tian@alum.mit.edu)
% last modified on 5/28/2014

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% define propagation operator, f: input field, h: propagation transfer
% function
Prop = @(f,h) Ft(F(f).*h);

%% backpropagate means take conj of transfer function
% H = conj(H);
% H0 = conj(H0);

% N: lateral dimension, Nslice: # of total z-slices
[N,~,Nslice] = size(o_slice0);
o_slice = zeros(N,N,Nslice);
%% backprop starts from here
% 1) the last slice is special from Pupil plane to Obj plan
% psi = Ft(O.*H0);
% [o_slice(:,:,Nslice), phi] = Proj_OslicePhi(O0,phi0(:,:,Nslice),psi,psi0(:,:,Nslice),alpha,beta,iters)
psi = Ft(O);
% 2) from Nslice-1 to 2nd slice 
for m = Nslice:-1:2
    % update o_slice at current slice
    [o_slice(:,:,m), phi] = Proj_OslicePhi(o_slice0(:,:,m),phi0(:,:,m),...
        psi,psi0(:,:,m),alpha,beta,iters);
    H = exp(-1i*k2*dz(m-1));
    % propagate next plane
    psi = Prop(phi, H);
end
% 3) the 1st slice also is different
o_slice(:,:,1) = psi./i0;
% dpsi = psi-psi0;
%o_slice(:,:,1) = o_slice(:,:,1) ...
%    +abs(i0).*conj(i0)./(abs(i0).^2+alpha).*(psi-psi0(:,:,1))/max(abs(i0(:)));
% o_slice(:,:,1) = o_slice(:,:,1)+(psi-psi0(:,:,1))./i0/(1+alpha);

end

