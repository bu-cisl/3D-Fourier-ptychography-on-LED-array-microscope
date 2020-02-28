% By Lei Tian, lei_tian@berkeley.edu
% last modified 5/27/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

addpath(['..\3D_code']);

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% F = @(x) fftshift(fft2(x));
% Ft = @(x) ifft2(ifftshift(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse problem use alternating projection
% 2/28/2014
% experiments, 4/1/2014
% account for geometry WITHOUT condenser, 3/22/2014

% By Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 LED expt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numlit = 1;

filedir = ['G:\Project_Backup\LED_Array_Microscopy\Expt\NoCondenser\TE300\2014-5-20\CardiacTissue\1LED\'];

imglist = dir([filedir,'ILED*.tif']);
% out_dir = ['.\Res',num2str(numlit),'LED-Result'];
% mkdir(out_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the current processing patch starting coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nstart = [784,1375];
% nstart = [981,1181];
% nstart = [1444,701];
Np = 140;

ns1 = 1:Np-Np/10:2160; ns1 = ns1(1:end-1);
ns2 = 11:Np-Np/10:2560; ns2 = ns2(1:end-1);
[ns2,ns1] = meshgrid(ns2,ns1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for multi-slice model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance between slices, unit: um
dz = 10;
% distance from the Nth slice to the focal plane
z0 = 0;

Nslice = 2;


%%
for l = 100%1:length(ns1(:))
nstart = [ns1(l),ns2(l)];
fn = [filedir,'Iled_0147.tif'];
I = imread(fn);
figure(30); imagesc(I(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1));  
axis image; colormap gray; axis off; drawnow;
% setup output folder for each patch
out_dir = ['.\Res-patch-',num2str(nstart(1)),'-',num2str(nstart(2)),'-',...
    num2str(numlit),'LED-Result'];
mkdir(out_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in general system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SystemSetup4x_Multislice();
%% load in data
% LED intensity normalization from the calibration data
% load('..\Intensity-LED-calibrate\ILEDMean40x');
Nimg = Nled;
Imea = zeros(Np,Np,Nimg);
Ibk = zeros(Nimg,1);
for m = 1:Nimg
    fn = [filedir,imglist(m).name];
    I = imread(fn);
    Imea(:,:,m) = double(I(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1));
    bk1 = mean2(double(I(2101:end,1:400)));
    bk2 = mean2(double(I(2101:end,1501:1800)));
%     bk3 = mean2(double(I(650:700,1100:1500)));
    Ibk(m) = min([bk1,bk2]);
    %     Inorm(:,:,m) = Imea(:,:,m)/ILEDMean40x(m);
    if Ibk(m)>300
        Ibk(m) = Ibk(m-1);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% Denoise I. remove high freq noise beyond support of OTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ilpf = zeros(Np,Np,Nimg);
for m = 1:Nimg
    % filter out the high freq noise
    Ilpf(:,:,m) = Ft(F(Imea(:,:,m)).*Ps_otf);
end

%% corresponding LED locations
% find the on-LED indices
ledidx = 1:Nled;
ledidx = reshape(ledidx,numlit,Nimg);
lit = Litidx(ledidx);
lit = reshape(lit,numlit,Nimg);
[litv,lith] = ind2sub([32,32],lit);
% find the index to reorder the measurements so that the image contains the
% central LEDs will be used first during the updates
dis_lit = sqrt((litv-lit_cenv-1).^2+(lith-lit_cenh-1).^2);
[dis_lit2,idx_led] = sort(min(dis_lit,[],1));

Nsh_lit = zeros(numlit,Nimg);
Nsv_lit = zeros(numlit,Nimg);

for m = 1:Nimg
    % should make sure it always covers all the leds
    % index of LEDs are lit for each pattern
    %lit = condenseridx(ceil(rand(numlit,1)*Nled));
    % corresponding index of spatial freq for the LEDs are lit
    lit0 = lit(:,m);
    Nsh_lit(:,m) = idx_u(lit0);
    Nsv_lit(:,m) = idx_v(lit0);
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit
Ns = [];
Ns(:,:,1) = Nsv_lit;
Ns(:,:,2) = Nsh_lit;

Imea_reorder = Imea(:,:,idx_led);
Ilpf_reorder = Ilpf(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% Denoise II. background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ithresh_reorder = Ilpf_reorder;
Ithresh_reorder = Ilpf_reorder;
for m = 1:Nimg
    Itmp = Ithresh_reorder(:,:,m);
    %     Itmp(Itmp<mean2(Itmp(1:10,:))) = 0;
    %     Itmp = Itmp-mean2(Itmp(1:6,:));
    %     Itmp(Itmp<0) = 0;
    %     Ithresh_reorder(:,:,m) = Itmp;
    Itmp = Itmp-Ibk_reorder(m);
    Ithresh_reorder(:,:,m) = Itmp;
    
    %     Ithresh_reorder(:,:,m) = Itmp-min(Itmp(:));
    Ithresh_reorder(:,:,m) = Ft(F(Ithresh_reorder(:,:,m)).*Ps_otf);
    Ithresh_reorder(Ithresh_reorder<0) = 0;
end

% Imea_norm_reorder = Imea_norm(:,:,idx_led);
Ns_reorder = Ns(:,idx_led,:);

%% this part check if the calculation of brightfield and darkfield  matches the experiments
illumination_na_reorder = illumination_na_used(idx_led);

for m = 1:Nimg
    Imean(m) = mean2(Ithresh_reorder(:,:,m));
end

snr = Imean(:)./Ibk_reorder(:);

%% reconstruction algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use only a sub-set number of the measurements to reconstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nused_vec = 10:10:293;

% for qq = 1:length(Nused_vec)
%     Nused = Nused_vec(qq);
% Nused = 293;
Nused = 293;
% NBF = 37;
% Ntrans = 46;
% if Nused<=NBF
%     I = Imea_reorder(:,:,1:Nused);
%     Ns2 = Ns_reorder(:,1:Nused,:);
% else
%     I = Imea_reorder(:,:,[1:NBF,Ntrans:Ntrans+Nused-NBF-1]);
%     Ns2 = Ns_reorder(:,[1:NBF,Ntrans:Ntrans+Nused-NBF-1],:);
% end
% I = Ithresh_reorder(:,:,1:Nused);
idx_used = [1:9,find(Imean(10:Nused)<Imean(1)/5)+9];

idx_err = find(Imean(10:Nused)>Imean(1)/5)+9;
disp(['problematic frames are ',num2str(idx_err),' and are discarded']);


% idx_used = [1:11,13:16,18:Nused];
I = Ithresh_reorder(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);

% reconstruction algorithm
opts.tol = 1;
opts.maxIter = 300;
opts.minIter = 150;
opts.monotone = 1;
% 'full', display every subroutin,
% 'iter', display only results from outer loop
% 0, no display
opts.display = 'iter';
% opts.saveIterResult = 0;
% opts.out_dir = ['.\tmp2'];
% mkdir(opts.out_dir);
% upsample the intensity
% I0interp = real(Ft(padarray(F(I(:,:,1)),[(N_obj-Np)/2,(N_obj-Np)/2])));
% opts.O0 = F(sqrt(I0interp));
% this method does not work for global minimization method
opts.Ps = w_NA;
opts.iters = 1;
% index of led used in the experiment
% opts.ledidx = ledidx(:,idx_led);
opts.OP_alpha = 1e-2;
opts.OP_beta = 1e-2;
opts.BP_alpha = 1e-2;
opts.BP_beta = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 0: initializatin
% idea: initialize with lightfield refocused intensities at different
% slices with only bright field data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
O0 = F(sqrt(I(:,:,1))).*w_NA;
O0 = padarray(O0,[(N_obj-Np)/2,(N_obj-Np)/2]);
% estimated object field at focal plane
o0 = Ft(O0);
% define propagation operator, f: input field, h: propagation transfer
% function
Prop = @(f,h) Ft(F(f).*h);
% propagate to the last slice
o0 = Prop(o0,conj(H0));
% one reasonable (or trial solution) is the last slice = o0 and the rest
% slices = 1:
o_slice0 = ones(N_obj,N_obj,Nslice);
% o_slice0(:,:,Nslice) = o0;
for l = 1:Nslice
    o_slice0(:,:,l) = o0;
end
opts.O0 = o_slice0;
opts.P0 = w_NA;


%%
[O,P,err] = AlterMin_MultiSlice(I,[N_obj,N_obj],Nslice, round(Ns2),H, H0, opts);

f3 = figure(88);
for m = 1:Nslice
subplot(2,2,2*m-1); imagesc(abs(O(:,:,m))); axis image; colormap gray; colorbar;
title('ampl(o)');
subplot(2,2,2*m); imagesc(angle(O(:,:,m))); axis image; colormap gray; colorbar;
title('phase(o)');
% subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
% title('ampl(P)');
% subplot(224); imagesc(angle(P)); axis image; colormap gray; colorbar;
% title('phase(P)');
end

% f4 = figure(79); plot(c(1:Nused));
% title('adaptive intensity correction factor');

fn = ['RandLit-',num2str(numlit),'-',num2str(Nused)];
save([out_dir,'\',fn],'O','P','err_pc','c','idx_err','snr');

saveas(f3,[out_dir,'\R-',fn,'.png']);
% saveas(f2,[out_dir,'\err-',fn,'.png']);
fprintf([fn,' saved\n']);

end










