function I_est = Fwd_Prop_MultiSlice_Intensity( i0, o_slice, k2, dz, P, H0)
%FWD_PROP_MULTISLICE computes the field using multislice approach, with
%propagator H
% Inputs:
%   H: fwd propagator between slices
%   H0: fwd propagator from Nth slice to focal plane of objective
%   o_slice0: current estimate of multi-slice object
%   i0: KNOWN illumination

% by Lei Tian (lei_tian@alum.mit.edu)
% last modified on 5/28/2014

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% define propagation operator, f: input field, h: propagation transfer
% function
Prop = @(f,h) Ft(F(f).*h);

% N: lateral dimension, Nslice: # of total z-slices
[N,~,Nslice] = size(o_slice); 
Np = size(P,1);

cen0 = round((N+1)/2);
downsamp = @(x) x(cen0-Np/2:cen0+Np/2-1,cen0-Np/2:cen0+Np/2-1);

% initialize incident field at each slice
% phi = zeros(N,N,Nslice);
% phi = i0; % incident field of 1st slice is illumination
% initialize output field at each slice
% psi = zeros(N,N,Nslice);
psi = i0.*o_slice(:,:,1);
for m = 2:Nslice
    H = exp(1i*k2*dz(m-1));
    % propagate from neiboring slices
    phi = Prop(psi,H);
    % output field = incidence * object
    psi = phi.*o_slice(:,:,m);
end

% estimated intensity w/o correction term
I_est = abs(Ft(downsamp(F(psi)).*P.*H0)).^2;



end

