% function [ varargout ] = SystemSetup( varargin )
%SYSTEMSETUP initilize general system parameters for LED array microscope

% Last modofied on 4/22/2014 
% Lei Tian (lei_tian@berkeley.edu)

% addpath(['..\..\Source_coding']);
% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% F = @(x) fftshift(fft2(x));
% Ft = @(x) ifft2(ifftshift(x));

% Define RMSE function operator
RMSE = @(x,x0) sqrt(sum(abs(x(:)-x0(:)).^2)/sum(x0(:).^2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelength of illumination, assume monochromatic
% R: 643nm +- 50nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.643;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical aperture of the objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA = 0.25; 
% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA
dx0 = 1/um_m/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnification of the system,
% need to calibrate with calibration slides
% on 2x objective, front port with an extra 2x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mag = 8/4*10;

dpix_c = 6.5; %6.5um pixel size on the sensor plane
% effective image pixel size on the object plane
dpix_m = dpix_c/mag; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # of pixels at the output image patch
% each patch will assign a single k-vector, the image patch size cannot be
% too large to keep the single-k assumption holds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Np = 100;

% FoV in the object space
FoV = Np*dpix_m;
% sampling size at Fourier plane set by the image size (FoV)
% sampling size at Fourier plane is always = 1/FoV
if mod(Np,2) == 1
    du = 1/dpix_m/(Np-1);
else
    du = 1/FoV;
end

% low-pass filter diameter set by the NA = bandwidth of a single measurment
% in index
% N_NA = round(2*um_m/du_m);
% generate cutoff window by NA
m = 1:Np;
[mm,nn] = meshgrid(m-round((Np+1)/2));
ridx = sqrt(mm.^2+nn.^2);
um_idx = um_m/du;
% assume a circular pupil function, lpf due to finite NA
w_NA = double(ridx<um_idx);
% h = fspecial('gaussian',10,5);
% w_NA = imfilter(w_NA,h);

% support of OTF is 2x of ATF(NA)
Ps_otf = double(ridx<2*um_idx);

phC = ones(Np);
% phC(ridx<0.8*um_idx&ridx>0.7*um_idx) = 0.5;
% aberration modeled by a phase function
aberration = ones(Np);
% aberration = exp(pi/2*1i*(exp(-(mm-20).^2/50^2-(nn+40).^2/150^2))-...
%     pi/8*1i*(exp(-(mm+40).^2/100^2-(nn-80).^2/80^2))+...
%     pi/3*1i*(exp(-(mm).^2/60^2-(nn-10).^2/30^2)));


%aberration = ones(N_m);
pupil = w_NA.*phC.*aberration;

clear m mm nn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up image corrdinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original image size: 2160x2560
% can calibrate the center of the illumination with respect to the image by
% looking at the data from the dark/bright field image transitions
ncent = [1080,1280];
% nstart = ncent-Np/2;
% start pixel of the image patch
% nstart = [981,1181];
% center, start & end of the image patch
img_center = (nstart-ncent+Np/2)*dpix_m;
img_start = nstart*dpix_m;
img_end = (nstart+Np)*dpix_m;
%img_center = [0,0];
%% LED array geometries and derived quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacing between neighboring LEDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_led = 4e3; %4mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from the LED to the object
% experientally determined by placing a grating object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_led = 74e3; %mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diameter of # of LEDs used in the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dia_led = 19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up LED coordinates
% h: horizontal, v: vertical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lit_cenv = 13;
lit_cenh = 13;
vled = [0:31]-lit_cenv;
hled = [0:31]-lit_cenh;

[hhled,vvled] = meshgrid(hled,vled);
rrled = sqrt(hhled.^2+vvled.^2);
LitCoord = rrled<dia_led/2;
% total number of LEDs used in the experiment
Nled = sum(LitCoord(:));
% index of LEDs used in the experiment
Litidx = find(LitCoord);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to flip x&y(spatial frequency u&v) based on expt geometry
% also with minus (-) on both v & h based on scanning order
% also need to know the relative orientation between led array and the
% image, otherwise get garbge near the edge of the FoV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corresponding angles for each LEDs
dd = sqrt((-hhled*ds_led-img_center(1)).^2+(-vvled*ds_led-img_center(2)).^2+z_led.^2);
sin_thetav = (-hhled*ds_led-img_center(1))./dd;
sin_thetah = (-vvled*ds_led-img_center(2))./dd;

Tanv = -(-hhled*ds_led)/z_led;
Tanh = -(-vvled*ds_led)/z_led;
Tanh_lit = Tanh(LitCoord);
Tanv_lit = Tanv(LitCoord);

illumination_na = sqrt(sin_thetav.^2+sin_thetah.^2);
% corresponding spatial freq for each LEDs
%
vled = sin_thetav/lambda;
uled = sin_thetah/lambda;
% spatial freq index for each plane wave relative to the center
idx_u = round(uled/du);
idx_v = round(vled/du);

illumination_na_used = illumination_na(LitCoord);

% number of brightfield image
NBF = sum(illumination_na_used<NA);
idx_BF = find(illumination_na_used<NA);
% darkfield images indices
NDF = sum(illumination_na_used>NA);
idx_DF = find(illumination_na_used>NA);
% maxium spatial frequency achievable based on the maximum illumination
% angle from the LED array and NA of the objective
um_p = max(illumination_na_used)/lambda+um_m;
% resolution achieved after freq post-processing
dx0_p = 1/um_p/2;

NA_s = um_p*lambda;
disp(['synthetic NA is ',num2str(NA_s)]);

% assume the max spatial freq of the original object
% um_obj>um_p
% assume the # of pixels of the original object
% upsample by 2 for intensity
N_obj = round(um_p/du)*2;%*2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to enforce N_obj/Np = integer to ensure no FT artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_obj = ceil(N_obj/Np)*Np;
% max spatial freq of the original object
um_obj = du*N_obj/2;

% sampling size of the object (=pixel size of the test image)
dx_obj = 1/um_obj/2;
% spatial coordiates
[x,y] = meshgrid([-N_obj/2:N_obj/2-1]*dx_obj);

[xn,yn] = meshgrid(-1/2:1/N_obj:1/2-1/N_obj);


%% spatial frequency coordinates
[u,v] = meshgrid(-um_obj:du:um_obj-du);
[u0,v0] = meshgrid(-du*Np/2:du:du*Np/2-du);


%% compute depth sectioning capability
DOF0 = lambda/NA^2/2;
DOFs = lambda/NA_s^2/2;

%% resolution
res0 = lambda/NA/2;
res_s = lambda/NA_s/2;

%% define propagation transfer function
eva = double(u0.^2+v0.^2<1/lambda^2);
H0 = exp(-1i*2*pi*sqrt(1/lambda^2-u0.^2-v0.^2)*z0).*eva;%.*mask_evanescent;
% H0 = exp(1i*2*pi*sqrt(1/lambda^2-u.^2-v.^2)*z0).*mask_evanescent;
k2 = pi*lambda*(u.^2+v.^2);
% k02 = pi*lambda*(u0.^2+v0.^2);
