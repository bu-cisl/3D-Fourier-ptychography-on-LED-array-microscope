% By Lei Tian, lei_tian@berkeley.edu
% last modified 5/27/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

%addpath(['C:\Users\Lei\Dropbox\Berkeley\LEDArray\MatlabCodes\Coded_Illumination\Source_coding']);

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

filedir = ['G:\Project_Backup\LED_Array_Microscopy\Expt\NoCondenser\NewArray\3D-2\Spiral2\10x_3\'];

imglist = dir([filedir,'ILED*.tif']);
% out_dir = ['.\Res',num2str(numlit),'LED-Result'];
% mkdir(out_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the current processing patch starting coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np = 1024;
% original image size
n1 = 2160;
n2 = 2560;

% ns1 = 1:Np-Np/10:2160; ns1 = ns1(1:end-1);
% ns2 = 11:Np-Np/10:2560; ns2 = ns2(1:end-1);
% [ns2,ns1] = meshgrid(ns2,ns1);
nstart = [1,1];

%%
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for multi-slice model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vector to define distance between slices, unit: um

dz = [20*ones(1,5)];
%  dz = [dz0(ll),40*ones(1,1)];
% distance from the Nth slice to the focal plane
z0 = 10;
% number of z-slices
Nslice = length(dz)+1;
% z vector
% z = (-z0-[Nslice-1:-1:0]*dz);
z = fliplr([-z0,-z0-(cumsum(fliplr(dz)))]);
m0 = find(abs(z)==min(abs(z)));
m0 = m0(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess the data by shift and add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SystemSetupV2Array10x_Multislice_v4();
% I0 = zeros(2160,2560,Nslice);
if Np<n1
    Ibf = zeros(n1,n2,Nslice);
    for m = 1:length(z)
        for n = 1:NBF % BF only
            nn = idx_BF(n);
            fn = [filedir,imglist(nn).name];
            tp = double(imread(fn));
            shift_x = z(m) * Tanh_lit(nn);
            shift_y = z(m) * Tanv_lit(nn);
            shift_xp = round(shift_x/dpix_m);
            shift_yp = round(shift_y/dpix_m);
            % shift
            tp = circshift(tp,[shift_yp,shift_xp]);
            Ibf(:,:,m) = Ibf(:,:,m)+tp;
        end
    end
    Ibf = (Ibf/NBF);
    % Ibf = (Ibf/NBF).^(1/Nslice);
    
    a = min(min(Ibf(:,:,m0)));
    b = max(max(Ibf(:,:,m0)));
    
%     for m = 1:Nslice
%         figure; imagesc(Ibf(:,:,m),[a,b]); axis image; colormap gray;
%         title(['z=',num2str(z(m))])
%     end
else
    Ibf = zeros(Np,Np,Nslice);
    for m = 1:length(z)
        for n = 1:NBF % BF only
            nn = idx_BF(n);
            fn = [filedir,imglist(nn).name];
            tp = double(imread(fn));
            bk1 = mean2(tp(100:200,400:600));
            bk2 = mean2(tp(1800:1900,500:600));
            %     bk3 = mean2(double(I(650:700,1100:1500)));
            bk = min([bk1,bk2]);
            %     Inorm(:,:,m) = Imea(:,:,m)/ILEDMean40x(m);
            tp = padarray(tp,[(Np-n1)/2,(Np-n2)/2],bk);
            
            shift_x = z(m) * Tanh_lit(nn);
            shift_y = z(m) * Tanv_lit(nn);
            shift_xp = round(shift_x/dpix_m);
            shift_yp = round(shift_y/dpix_m);
            % shift
            tp = circshift(tp,[shift_yp,shift_xp]);
            Ibf(:,:,m) = Ibf(:,:,m)+tp;
        end
    end
    Ibf = (Ibf/NBF);
    % Ibf = (Ibf/NBF).^(1/Nslice);
    
    a = min(min(Ibf(:,:,m0)));
    b = max(max(Ibf(:,:,m0)));
    
    for m = 1:Nslice
        figure; imagesc(Ibf(:,:,m),[a,b]); axis image; colormap gray;
        title(['z=',num2str(z(m))])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear up some memory space
clear tp u u0 v v0 yn y xn x phC pupil ridx aberration dd hhled vvled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% for l = 1%1:length(ns1(:))
nstart = [700,500];
% nstart = [1600,800]; % check z0
% nstart = [800,1320]; % check z0
% fn = [filedir,'Iled_0147.tif'];
% I = imread(fn);
if Np<n1
    figure(30); imagesc(Ibf(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1,2));
    axis off; axis image; colormap gray; drawnow;
end
% setup output folder for each patch
out_dir = ['.\Res-patch-',num2str(nstart(1)),'-',num2str(nstart(2)),'-',...
    num2str(numlit),'LED-Result'];
mkdir(out_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in general system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SystemSetupV2Array10x_Multislice_v4();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear up some memory space
clear tp u u0 v v0 yn y xn x phC pupil ridx aberration dd hhled vvled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in data
% LED intensity normalization from the calibration data
% load('..\Intensity-LED-calibrate\ILEDMean40x');
Nimg = Nled;
Imea = zeros(Np,Np,Nimg);
Ibk = zeros(Nimg,1);
for m = 1:Nimg
    fn = [filedir,imglist(m).name];
    I = double(imread(fn));
    bk1 = mean2(I(200:400,200:400));
    bk2 = mean2(I(1800:2000,2000:2200));
    %     bk3 = mean2(double(I(650:700,1100:1500)));
    Ibk(m) = min([bk1,bk2]);
    %     Inorm(:,:,m) = Imea(:,:,m)/ILEDMean40x(m);
    if Np<n1
        Imea(:,:,m) = I(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1);
    else
        Imea(:,:,m) = padarray(I,[(Np-n1)/2,(Np-n2)/2],Ibk(m));
    end
    
    if Ibk(m)>600
        Ibk(m) = Ibk(m-1);
    end
    %     Ibk(m) = 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear up some memory space
clear I 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % pre-processing the data to DENOISING is IMPORTANT
% % Denoise I. remove high freq noise beyond support of OTF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ilpf = zeros(Np,Np,Nimg);
% for m = 1:Nimg
%     % filter out the high freq noise
%     Ilpf(:,:,m) = Ft(F(Imea(:,:,m)).*Ps_otf);
%     m
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % clear up some memory space
% clear Imea I
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Imea_reorder = Imea(:,:,idx_led);
% Ilpf_reorder = Ilpf(:,:,idx_led);
% Ilpf_reorder = Imea(:,:,idx_led);
% Ithresh_reorder = Imea(:,:,idx_led);
Imea = Imea(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear up some memory space
clear Itmp Ibk tp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% Denoise II. background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ithresh_reorder = Ilpf_reorder;
for m = 1:Nimg
%     Itmp = Ithresh_reorder(:,:,m);
    Itmp = Imea(:,:,m);
    %     Itmp(Itmp<mean2(Itmp(1:10,:))) = 0;
    %     Itmp = Itmp-mean2(Itmp(1:6,:));
    %     Itmp(Itmp<0) = 0;
    %     Ithresh_reorder(:,:,m) = Itmp;
    Itmp = Itmp-Ibk_reorder(m);
%     Ithresh_reorder(:,:,m) = Itmp;
    Imea(:,:,m) = Itmp;
    %     Ithresh_reorder(:,:,m) = Itmp-min(Itmp(:));
    %     Ithresh_reorder(:,:,m) = Ft(F(Ithresh_reorder(:,:,m)).*Ps_otf);
end
Imea(Imea<0) = 0;

% for m = NBF+1:Nimg
% %     Itmp = Ithresh_reorder(:,:,m);
%     Itmp = Imea(:,:,m);
%     %     Itmp(Itmp<mean2(Itmp(1:10,:))) = 0;
%     %     Itmp = Itmp-mean2(Itmp(1:6,:));
%     %     Itmp(Itmp<0) = 0;
%     %     Ithresh_reorder(:,:,m) = Itmp;
%     min(Itmp(:))
%     Itmp = Itmp-min(Itmp(:));
%     min(Itmp(:))
% %     Ithresh_reorder(:,:,m) = Itmp;
%     Imea(:,:,m) = Itmp;
%     %     Ithresh_reorder(:,:,m) = Itmp-min(Itmp(:));
%     %     Ithresh_reorder(:,:,m) = Ft(F(Ithresh_reorder(:,:,m)).*Ps_otf);
% end
% 
% % Ithresh_reorder(Ithresh_reorder<0) = 0;
% Imea(Imea<0) = 0;
% Imea_norm_reorder = Imea_norm(:,:,idx_led);
Ns_reorder = Ns(:,idx_led,:);
Tanh_reorder = Tanh(idx_led);
Tanv_reorder = Tanv(idx_led);

%% this part check if the calculation of brightfield and darkfield  matches the experiments
illumination_na_reorder = illumination_na_used(idx_led);

for m = 1:Nimg
    Imean(m) = mean2(Imea(:,:,m));
end

snr = Imean(:)./Ibk_reorder(:);

%% reconstruction algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use only a sub-set number of the measurements to reconstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 9;
switch d
    case 3
        Nused = 9;
    case 5
        Nused = 21;
    case 7
        Nused = 37;
    case 9
        Nused = 69;
    case 11
        Nused = 97;
    case 13
        Nused = 137;
    case 15
        Nused = 177;
    case 17
        Nused = 225;
    case 19
        Nused = 293;
end


if Nused>NBF
    idx_used = [find(Imean(1:NBF)>Imean(1)*2/3),find(Imean(NBF+1:Nused)<Imean(1)/4)+NBF];
    idx_err = [find(Imean(1:NBF)<Imean(1)*2/3),find(Imean(NBF+1:Nused)>Imean(1)/4)+NBF];
else
    idx_used = find(Imean(1:Nused)>Imean(1)*2/3);
    idx_err = [find(Imean(1:Nused)<Imean(1)*2/3)];
    
end

disp(['problematic frames are ',num2str(idx_err),' and are discarded']);

% I = Ithresh_reorder(:,:,idx_used);
I = Imea(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);
Tanh2 = Tanh_reorder(idx_used);
Tanv2 = Tanh_reorder(idx_used);

% reconstruction algorithm
opts.tol = 1;
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
opts.BP_alpha = 1e-8;
opts.BP_beta = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 0: initializatin
% idea: initialize with lightfield refocused intensities at different
% slices with only bright field data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one reasonable (or trial solution) is the last slice = o0 and the rest
% slices = 1:
o_slice0 = zeros(N_obj,N_obj,Nslice);
o_slice = Ibf(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1,:);
for n = 1:Nslice
    O0 = F(sqrt(o_slice(:,:,n)));
    O0 = padarray(O0,[(N_obj-Np)/2,(N_obj-Np)/2]);
    % estimated object field at focal plane
    o_slice0(:,:,n) = abs(Ft(O0));
end
% o_slice0(:,:,Nslice) = o0;
% for l = 1:Nslice
o_slice0 = o_slice0.^(1/Nslice);
% end
% O0 = F(sqrt(I(:,:,1))).*w_NA;
% O0 = padarray(O0,[(N_obj-Np)/2,(N_obj-Np)/2]);
% % estimated object field at focal plane
% o0 = Ft(O0);
% opts.O0 = o0;
% % define propagation operator, f: input field, h: propagation transfer
% % function
% Prop = @(f,h) Ft(F(f).*h);
% % propagate to the last slice
% o0 = Prop(o0,conj(H0));

opts.O0 = o_slice0;
opts.P0 = w_NA;

opts.maxIter = 1;
opts.minIter = 1;

%     title(['slice',num2str(m0)],'fontsize',10); axis off;
if Nslice>1
    f4 = figure(87);
    switch Nslice
        case 2
            p = 2; q = 1;
        case {3,4}
            p = 2; q = 2;
        case {5,6}
            p = 3; q = 2;
        otherwise
            p = 3; q = 3;
    end
    for m = 1:min(Nslice,9)
        subplot(p,q,m); imagesc(abs(opts.O0(:,:,m))); axis image; colormap gray;
        h= colorbar; set(h,'fontsize',8); axis off;
        title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
    end
    
end
%     title('phase(P)'); axis off;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear up some memory space
clear pupil phC o_slice0 o_slice aberration Ps_otf O0
clear Itmp Imea_reorder Ilpf_reorder Ilpf Ibk_reorder
clear Ibk tp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start optimization
[O,P,err] = AlterMin_MultiSlice_v2(I,[N_obj,N_obj], round(Ns2),k2, dz, H0, opts);

f3 = figure(88);
%     title(['slice',num2str(m0)],'fontsize',10); axis off;
if Nslice>1
    switch Nslice
        case 2
            p = 2; q = 1;
        case {3,4}
            p = 2; q = 2;
        case {5,6}
            p = 3; q = 2;
        otherwise
            p = 3; q = 3;
    end
    for m = 1:min(Nslice,9)
        subplot(p,q,m); imagesc(abs(O(:,:,m))); axis image; colormap gray;
        h= colorbar; set(h,'fontsize',8); axis off;
        title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
    end
else
    subplot(221); imagesc(abs(O(:,:,Nslice))); axis image; colormap gray;
    h= colorbar; set(h,'fontsize',8); axis off;
    subplot(222); imagesc(angle(O(:,:,Nslice))); axis image; colormap gray; colorbar;
    h= colorbar; set(h,'fontsize',8); axis off;
    %         title('phase(o)','fontsize',10); axis off;
    subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
    h= colorbar; set(h,'fontsize',8); axis off;
    %     title('ampl(P)'); axis off;
    subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
    h= colorbar; set(h,'fontsize',8); axis off;
end
%     title('phase(P)'); axis off;
drawnow;

% f4 = figure(79); plot(c(1:Nused));
% title('adaptive intensity correction factor');

fn = ['3D-',num2str(Np),'-',num2str(Nused),'-ns-',num2str(Nslice),'-dz-',num2str(dz),'-z0-',num2str(z0)];
saveas(f3,[out_dir,'\R-',fn,'.png']);
if Nslice>1
    saveas(f4,[out_dir,'\LF-',fn,'.png']);
end
save([out_dir,'\',fn],'O','P','err');

% saveas(f2,[out_dir,'\err-',fn,'.png']);
fprintf([fn,' saved\n']);

%% equalizing background
bk0 = mean2(abs(O(770:860,560:690,3)));
mb = 3;
ma = 0;
s0 = 50;
for m = 1:Nslice
    oo = abs(O(s0+1:end-s0,s0+1:end-s0,m));
    bkt = mean2(abs(oo(720:810,510:640)));
    oo = oo-(bkt-bk0);
    f1 = figure('Colormap',...
    [0 0 0;0.0011834615143016 0.0011834615143016 0.0011834615143016;0.0023669230286032 0.0023669230286032 0.0023669230286032;0.00355038465932012 0.00355038465932012 0.00355038465932012;0.00473384605720639 0.00473384605720639 0.00473384605720639;0.00591730745509267 0.00591730745509267 0.00591730745509267;0.00710076931864023 0.00710076931864023 0.00710076931864023;0.00828423071652651 0.00828423071652651 0.00828423071652651;0.00946769211441278 0.00946769211441278 0.00946769211441278;0.0106511535122991 0.0106511535122991 0.0106511535122991;0.0118346149101853 0.0118346149101853 0.0118346149101853;0.0130180763080716 0.0130180763080716 0.0130180763080716;0.0142015386372805 0.0142015386372805 0.0142015386372805;0.0153850000351667 0.0153850000351667 0.0153850000351667;0.016568461433053 0.016568461433053 0.016568461433053;0.0177519228309393 0.0177519228309393 0.0177519228309393;0.0189353842288256 0.0189353842288256 0.0189353842288256;0.0476649329066277 0.0476649329066277 0.0476649329066277;0.0763944834470749 0.0763944834470749 0.0763944834470749;0.105124033987522 0.105124033987522 0.105124033987522;0.133853584527969 0.133853584527969 0.133853584527969;0.162583127617836 0.162583127617836 0.162583127617836;0.191312685608864 0.191312685608864 0.191312685608864;0.22004222869873 0.22004222869873 0.22004222869873;0.248771786689758 0.248771786689758 0.248771786689758;0.277501344680786 0.277501344680786 0.277501344680786;0.306230872869492 0.306230872869492 0.306230872869492;0.334960430860519 0.334960430860519 0.334960430860519;0.363689988851547 0.363689988851547 0.363689988851547;0.392419517040253 0.392419517040253 0.392419517040253;0.421149075031281 0.421149075031281 0.421149075031281;0.449878633022308 0.449878633022308 0.449878633022308;0.478608191013336 0.478608191013336 0.478608191013336;0.507337749004364 0.507337749004364 0.507337749004364;0.536067306995392 0.536067306995392 0.536067306995392;0.564796805381775 0.564796805381775 0.564796805381775;0.593526363372803 0.593526363372803 0.593526363372803;0.622255921363831 0.622255921363831 0.622255921363831;0.650985479354858 0.650985479354858 0.650985479354858;0.679715037345886 0.679715037345886 0.679715037345886;0.708444595336914 0.708444595336914 0.708444595336914;0.737174153327942 0.737174153327942 0.737174153327942;0.765903651714325 0.765903651714325 0.765903651714325;0.794633209705353 0.794633209705353 0.794633209705353;0.823362767696381 0.823362767696381 0.823362767696381;0.852092325687408 0.852092325687408 0.852092325687408;0.874098122119904 0.874098122119904 0.874098122119904;0.896103858947754 0.896103858947754 0.896103858947754;0.918109655380249 0.918109655380249 0.918109655380249;0.940115451812744 0.940115451812744 0.940115451812744;0.962121188640594 0.962121188640594 0.962121188640594;0.98412698507309 0.98412698507309 0.98412698507309;0.985449731349945 0.985449731349945 0.985449731349945;0.986772477626801 0.986772477626801 0.986772477626801;0.988095223903656 0.988095223903656 0.988095223903656;0.989417970180511 0.989417970180511 0.989417970180511;0.990740716457367 0.990740716457367 0.990740716457367;0.992063522338867 0.992063522338867 0.992063522338867;0.993386268615723 0.993386268615723 0.993386268615723;0.994709014892578 0.994709014892578 0.994709014892578;0.996031761169434 0.996031761169434 0.996031761169434;0.997354507446289 0.997354507446289 0.997354507446289;0.998677253723145 0.998677253723145 0.998677253723145;1 1 1]);
    imagesc(oo,[ma,mb]); axis image; %colormap gray
    axis off;
    fn = ['ampl_z=',num2str(z(m)),'um.tif'];
    export_fig(f1, fn, '-m2');

%     title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
end

%%
bk0 = mean2(angle(O(770:860,560:690,1)));
% cm = ColormapRGW( -.6, 6 ,1);

for m = 1:Nslice
    oo = angle(O(s0+1:end-s0,s0+1:end-s0,m));
    bkt = mean2(oo(720:810,510:640));
    oo = oo-(bkt-bk0);
    f1 = figure('Colormap',...
    [0 0.498039215803146 0;0.00230654748156667 0.499197006225586 0.00230654748156667;0.00461309496313334 0.500354826450348 0.00461309496313334;0.0069196424447 0.501512587070465 0.0069196424447;0.00922618992626667 0.502670407295227 0.00922618992626667;0.0115327378734946 0.503828227519989 0.0115327378734946;0.0138392848894 0.504985988140106 0.0138392848894;0.016145832836628 0.506143808364868 0.016145832836628;0.0184523798525333 0.507301568984985 0.0184523798525333;0.0207589287310839 0.508459389209747 0.0207589287310839;0.0230654757469893 0.509617209434509 0.0230654757469893;0.0253720227628946 0.510774970054626 0.0253720227628946;0.0276785697788 0.511932790279388 0.0276785697788;0.0299851186573505 0.513090550899506 0.0299851186573505;0.0322916656732559 0.514248371124268 0.0322916656732559;0.0908203125 0.543627440929413 0.0908203125;0.149348959326744 0.573006510734558 0.149348959326744;0.207877606153488 0.602385640144348 0.207877606153488;0.266406238079071 0.631764709949493 0.266406238079071;0.324934899806976 0.661143779754639 0.324934899806976;0.383463531732559 0.690522909164429 0.383463531732559;0.441992193460464 0.719901978969574 0.441992193460464;0.500520825386047 0.749281048774719 0.500520825386047;0.559049487113953 0.778660118579865 0.559049487113953;0.617578148841858 0.80803918838501 0.617578148841858;0.676106750965118 0.8374183177948 0.676106750965118;0.734635412693024 0.866797387599945 0.734635412693024;0.793164074420929 0.89617645740509 0.793164074420929;0.851692736148834 0.92555558681488 0.851692736148834;0.910221338272095 0.954934656620026 0.910221338272095;0.96875 0.984313726425171 0.96875;0.984375 0.992156863212585 0.984375;1 1 1;1 0.983870983123779 0.983870983123779;1 0.967741906642914 0.967741906642914;1 0.912778854370117 0.912778854370117;1 0.857815861701965 0.857815861701965;1 0.802852809429169 0.802852809429169;1 0.747889816761017 0.747889816761017;1 0.69292676448822 0.69292676448822;1 0.637963712215424 0.637963712215424;1 0.583000719547272 0.583000719547272;1 0.528037667274475 0.528037667274475;1 0.473074644804001 0.473074644804001;1 0.418111622333527 0.418111622333527;1 0.36314857006073 0.36314857006073;1 0.308185547590256 0.308185547590256;1 0.253222525119781 0.253222525119781;1 0.198259502649307 0.198259502649307;1 0.143296465277672 0.143296465277672;1 0.088333435356617 0.088333435356617;1 0.0333704091608524 0.0333704091608524;1 0.0305895414203405 0.0305895414203405;1 0.0278086736798286 0.0278086736798286;1 0.0250278078019619 0.0250278078019619;1 0.02224694006145 0.02224694006145;1 0.0194660723209381 0.0194660723209381;1 0.0166852045804262 0.0166852045804262;1 0.0139043368399143 0.0139043368399143;1 0.011123470030725 0.011123470030725;1 0.00834260229021311 0.00834260229021311;1 0.0055617350153625 0.0055617350153625;1 0.00278086750768125 0.00278086750768125;1 0 0]);
    imagesc(oo,[-.6,.6]); axis image; axis off
    fn = ['ph_z=',num2str(z(m)),'um.tif'];
    export_fig(f1, fn, '-m2');
%     title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
end

% %% light field refocusing
% bk0 = mean2(abs(opts.O0(770:860,560:690,1)).^2);
% ma = min(min(abs(opts.O0(s0+1:end-s0,s0+1:end-s0,1)).^2));
% mb = max(max(abs(opts.O0(s0+1:end-s0,s0+1:end-s0,1)).^2));
% for m = 1:Nslice
%     oo = abs(opts.O0(s0+1:end-s0,s0+1:end-s0,m)).^2;
%     bkt = mean2(oo(720:810,510:640));
%     oo = oo-(bkt-bk0);
%     f1 = figure; imagesc(oo,[ma,mb]); axis image; colormap gray; axis off;
%     fn = ['lf_z=',num2str(z(m)),'um.tif'];
%     export_fig(f1, fn, '-m2');
% 
% %     title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
% end
% 
% for m = 1:Nslice
%     f1 = figure; imagesc(Ibf(:,:,m)); axis image; colormap gray; axis off;
%     fn = ['lf-bf_z=',num2str(z(m)),'um'];
%     export_fig(f1, fn, '-m1');
% 
% %     title(['z=',num2str(z(m)),'um'],'fontsize',10);axis off;
% end

% %% physical actual focus
% fdir = ['G:\Project_Backup\LED_Array_Microscopy\Expt\NoCondenser\NewArray\3D-2\Spiral2\10x_m3\'];
% ilist = dir([fdir,'*.tif']);
% Nc = Np-2*s0;
% nstartc = nstart+s0;
% for n = 1:length(ilist) % BF only
%     fn = [fdir,ilist(n).name];
%     Iz = double(imread(fn));
%     Iz = Iz(nstartc(1):nstartc(1)+Nc-1,nstartc(2):nstartc(2)+Nc-1);
%     if n ==1
%         bk0 = mean2(Iz(720:810,510:640));
%         ma = 0;
%         mb = max(Iz(:));
%     end
%     bkt = mean2(Iz(720:810,510:640));
%     Iz = Iz-(bkt-bk0);
%     f1 = figure; imagesc(Iz,[ma,mb]); axis image; colormap gray; axis off;
%     fn = ['mf_z=',ilist(n).name];
%     export_fig(f1, fn, '-m2');
% end

