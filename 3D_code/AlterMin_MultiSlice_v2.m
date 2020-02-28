function [o_slice, P, err] = AlterMin_MultiSlice_v2( I, No, Ns, k2, dz, H0, opts)
%AlterMinGlobal_Adaptive Implement alternative minimization sequentially on a stack of
%measurement I (n1 x n2 x nz). It consists of 2 loop. The main loop update
%the reconstruction results r. the inner loop applies projectors/minimizers
%P1 and P2 on each image I and steps through the entire dataset
%   Outputs:
%   r: reconsturcted high-res image
%   err: errors at each iteration
%
%   Inputs:
% Measurements data
%   I: intensity measurements by different LEDs
%   du: sampling pixel size in spatial freq domain
%   um: Max spatial freq of I set by NA
% Reconstruction parameters
%   No = [Ny_obj,Nx_obj]: size of the reconstructed image
%   k2 = pi*lambda*(u^2+v^2), used in calculating Fresnel propagation
%   dz: % vector to define distance between slices, unit: um
%   z0: % distance from the Nth slice to the focal plane
%
% Illumination coding parameters
%   Ns = [Nsy,Nsx]: centers of corresponding lpf regions for
%   the illumination pattern

% Iteration parameters: opts
%   tol: maximum change of error allowed in two consecutive iterations
%   maxIter: maximum iterations allowed
%   minIter: minimum iterations
%   monotone (1, default): if monotone, error has to monotonically dropping
%   when iters>minIter
%   display: display results (0: no (default) 1: yes)
%   saveIterResult: save results at each step as images (0: no (default) 1: yes)
%   R0: initial guess of R
%
%% reconstruction algorithm: partial coherence effect in both spatial and
% Fourier domain
% spatial updating method:
% ref [1] C. Rydberg, J. Bengtsson, Opt. Express, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the idea of the inverse algorithm is that each point on LED is a coherent
% mode. The total intensity is the sum of the intensities due to all the
% modes.
% although we do not know the actual intensity for each modes, as show in
% [1], an update rule for amplitude can be found by properly scale the
% estimated amplitude for each modes. Again, the phase part is left
% unchanged.
%
% within each LED, different sub-images' spectrum are replaced sequentially
% last modified 08/19/2013
%
% within each LED, weighted average of the spectrum overlapping patches of
% all the sub-images
% last modified 08/19/2013
%
% last modified 3/1/2014
%
% implementing the adaptive routine described in
% Z. Bian, S. Dong & G. Zheng, OE, 2013
% Last modified on 03/25/2014

% Last modified 4/24/2014
% by Lei Tian, lei_tian@alum.mit.edu

% Last modified on 9/7/2014
% by Lei Tian, lei_tian@alum.mit.edu

%% derived constants
% size of measurement
[Nmy,Nmx,Nimg] = size(I);
Np = [Nmy,Nmx];
% r = # of LEDs lit up in each pattern, # of coherent modes
[r,~,~] = size(Ns);
cen0 = round((No+1)/2);

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% F = @(x) fftshift(fft2(x));
% Ft = @(x) ifft2(ifftshift(x));

% col = @(x) x(:);
% row = @(x) x(:).';


%% options for the method
if nargin<4
    % default values
    opts.tol = 1;
    opts.maxIter = 50;
    opts.minIter = 3;
    opts.monotone = 1;
    opts.display = 0;
    opts.saveIterResult = 0;
    opts.out_dir = [];
    opts.O0 = Ft(sqrt(I(:,:,1)))/r;
    opts.O0 = padarray(opts.O0,(No-Np)/2);
    opts.P0 = ones(Np);
    opts.OP_alpha = 1;
    opts.OP_beta = 1;
    opts.BP_alpha = 1;
    opts.BP_beta = 1;
    % index of led used in the experiment
    %     opts.ledidx = 1:Nled;
    %     opts.scale_tau = 1e-6;
    %     opts.min_mode = 'seq';
    %     opts.fourier_mode = 'projection';
else
    if ~isfield(opts,'tol')
        opts.tol = 1;
    end
    if ~isfield(opts,'maxIter')
        opts.maxIter = 50;
    end
    if ~isfield(opts,'minIter')
        opts.minIter = 3;
    end
    if ~isfield(opts,'monotone')
        opts.monotone = 1;
    end
    if ~isfield(opts,'display')
        opts.display = 0;
    end
    if ~isfield(opts,'saveIterResult')
        opts.saveIterResult = 0;
    end
    if ~isfield(opts,'out_dir')
        opts.out_dir = ['IterResults'];
        if opts.saveIterResult
            mkdir(opts.out_dir);
        end
    end
    if ~isfield(opts,'O0')
        opts.O0 = Ft(sqrt(I(:,:,1)))/r;
        opts.O0 = padarray(opts.O0,(No-Np)/2);
    end
    if ~isfield(opts,'P0')
        opts.P0 = ones(Np);
    end
    if ~isfield(opts,'OP_alpha')
        opts.OP_alpha = 1;
    end
    if ~isfield(opts,'OP_beta')
        opts.OP_beta = 1;
    end
    if ~isfield(opts,'BP_alpha')
        opts.BP_alpha = 1;
    end
    if ~isfield(opts,'BP_beta')
        opts.BP_beta = 1;
    end
    if ~isfield(opts,'Ps')
        opts.Ps = 1;
    end
    if ~isfield(opts,'iters')
        opts.iters = 1;
    end
    %     if ~isfield(opts,'scalecorrect')
    %         opts.scalecorrect = 0;
    %     end
    %     if ~isfield(opts,'scale')
    %         opts.scale = ones(Nled,1);
    %     end
    %     if ~isfield(opts,'ledidx')
    %         opts.ledidx = 1:Nled;
    %     end
    %     if ~isfield(opts,'scale_tau')
    %         opts.scale_tau = 1e-6;
    %     end
    %     if ~isfield(opts,'min_mode')
    %         opts.min_mode = 'seq';
    %     end
    %     if ~isfield(opts,'fourier_mode')
    %         opts.fourier_mode = 'projection';
    %     end
end



T0 = clock;

fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 0: initializatin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = opts.P0; opts.P0 = 0;
o_slice = opts.O0; opts.O0 = 0;
err1 = inf;
err2 = 50;
iter = 0;
Ps = opts.Ps;
iters = opts.iters;


% Num of slices
m0 = length(dz)+1;


if opts.display
    f1 = figure(88);
%     title(['slice',num2str(m0)],'fontsize',10); axis off;
    if m0>1
        switch m0
            case 2
                p = 2; q = 1;
            case {3,4}
                p = 2; q = 2;
            case {5,6}
                p = 3; q = 2;
            otherwise
                p = 3; q = 3;
        end
        for m = 1:min(m0,9)
            subplot(p,q,m); imagesc(abs(o_slice(:,:,m))); axis image; colormap gray; 
            h= colorbar; set(h,'fontsize',8); axis off;
%             title(['slice',num2str(m0-m)],'fontsize',10);axis off;
        end
    else
        subplot(221); imagesc(abs(o_slice(:,:,m0))); axis image; colormap gray; 
        h= colorbar; set(h,'fontsize',8); axis off;
        subplot(222); imagesc(angle(o_slice(:,:,m0))); axis image; colormap gray; colorbar;
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
end

if opts.saveIterResult
    saveas(f1,[opts.out_dir,'\R_',num2str(iter),'.png']);
    %     saveas(f2,[opts.out_dir,'\Ph_',num2str(iter),'.png']);
end


%% main algorithm starts here
% stopping criteria: when relative change in error falls below some value,
% can change this value to speed up the process by using a larger value but
% will trading off the reconstruction accuracy
% error is defined by the difference b/w the measurement and the estimated
% images in each iteration
% fprintf('| %2d   | %.2e |\n',iter,err1);


downsamp = @(x) x(cen0(1)-Np(1)/2:cen0(1)+Np(1)/2-1,...
    cen0(2)-Np(2)/2:cen0(2)+Np(2)/2-1);

% project onto the intensity measurement space
Proj_FT = @(f,I) F(sqrt(I).*exp(1i*angle(f)));

if mod(No(1),2) == 1
    x = [-1/2:1/(No(2)-1):1/2];
    y = [-1/2:1/(No(1)-1):1/2];
else
    x = [-1/2:1/No(2):1/2-1/No(2)];
    y = [-1/2:1/No(2):1/2-1/No(2)];
end
[x,y] = meshgrid(x,y);

% H0 = exp(1i*k02*z0);
while abs(err1-err2)>opts.tol&&iter<opts.maxIter
    err1 = err2;
%     err2 = 0;
    err = zeros(1,Nimg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % incremetal updating method, process image by image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m = 1:Nimg
        % compute illumination
        i0 = exp(1i*2*pi*(x*Ns(:,m,2)+y*Ns(:,m,1)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 1.1: Fwd propagation:
        %           compute incident (phi) and output (psi) field for each slice
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [phi_fwd, psi_fwd] = Fwd_Prop_MultiSlice_v2( i0, o_slice, k2, dz);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 1.2: compute object field & estimating intensity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        O = F(psi_fwd(:,:,m0));
        G_fwd = downsamp(O).*P;
        %G_fwd2 = G_fwd.*H0;
        % compute field in real/measurement space
        %g = Ft(G_fwd2);
        g = Ft(G_fwd.*H0);
        % estimated intensity w/o correction term
        I_est = abs(g).^2;
        % measured intensity
        I_mea = I(:,:,m);
        %G_new2 = Proj_FT(g, I_mea);
        %G_new = G_new2./H0;
        G_new = Proj_FT(g, I_mea)./H0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stopping criteria: residual between estimated & measured intensities
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        err(m) = sum((I_mea(:)-I_est(:)).^2);
%         err2 = err2+sum((I_mea(:)-I_est(:)).^2);
        %         if opts.scalecorrect
        %             if iter>1
        %                 %% Step 1a) Update scaling factor due to LED non-uniformity
        %                 scale(m) = scaling_update(scale(m), I_mea, I_est, ...
        %                     opts.scale_tau, opts.scale_mode, opts.scale_alpha);
        %             end
        %         end
        %         I_est = scale(m)*I_est;
        %         Psi = Proj_Fourier(G, I_mea, scale(m));
        %             Psi = Psi.*Ps;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 2: update O & P from product constraint G = O*P
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [O, P] = Proj_ObjPupil(O, P, G_new, G_fwd, Ps,...
            opts.OP_alpha, opts.OP_beta, iters);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 3: back-propgate field based on multi-slice approach & update
        % estimationo of o_slice
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        o_slice = Back_Prop_MultiSlice_v2(O, k2, dz, o_slice, phi_fwd, psi_fwd, ...
            i0, opts.BP_alpha, opts.BP_beta, iters);
        
        
    end
    err2 = sqrt(sum(err));
    fprintf('| %2d   | %.2e |\n',iter,err2);
    
    iter = iter+1;
    if strcmp(opts.display,'full')
        f1 = figure(88);
        %     title(['slice',num2str(m0)],'fontsize',10); axis off;
        if m0>1
            switch m0
                case 2
                    p = 2; q = 1;
                case {3,4}
                    p = 2; q = 2;
                case {5,6}
                    p = 3; q = 2;
                otherwise
                    p = 3; q = 3;
            end
            for m = 1:min(m0,9)
                subplot(p,q,m); imagesc(abs(o_slice(:,:,m))); axis image; colormap gray;
                h= colorbar; set(h,'fontsize',8); axis off;
                %             title(['slice',num2str(m0-m)],'fontsize',10);axis off;
            end
        else
            subplot(221); imagesc(abs(o_slice(:,:,m0))); axis image; colormap gray;
            h= colorbar; set(h,'fontsize',8); axis off;
            subplot(222); imagesc(angle(o_slice(:,:,m0))); axis image; colormap gray; colorbar;
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
        
    end
    %end
    
%     %% compute error
%     % record the error and can check the convergence later.
%     err = [err,err2];
    
    if strcmp(opts.display,'iter')
        f1 = figure(88);
        %     title(['slice',num2str(m0)],'fontsize',10); axis off;
        if m0>1
            switch m0
                case 2
                    p = 2; q = 1;
                case {3,4}
                    p = 2; q = 2;
                case {5,6}
                    p = 3; q = 2;
                otherwise
                    p = 3; q = 3;
            end
            for m = 1:min(m0,9)
                subplot(p,q,m); imagesc(abs(o_slice(:,:,m))); axis image; colormap gray;
                h= colorbar; set(h,'fontsize',8); axis off;
                %             title(['slice',num2str(m0-m)],'fontsize',10);axis off;
            end
        else
            subplot(221); imagesc(abs(o_slice(:,:,m0))); axis image; colormap gray;
            h= colorbar; set(h,'fontsize',8); axis off;
            subplot(222); imagesc(angle(o_slice(:,:,m0))); axis image; colormap gray; colorbar;
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
    end
    
    if opts.saveIterResult
        saveas(f1,[opts.out_dir,'\R_',num2str(iter),'.png']);
    end
    
    if opts.monotone&&iter>opts.minIter
        if err2>err1
            break;
        end
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare to end the process, calculate the error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err = zeros(1,Nimg);
for m = 1:Nimg
    % compute illumination
    i0 = exp(1i*2*pi*(x*Ns(:,m,2)+y*Ns(:,m,1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 1.1: Fwd propagation:
    %           compute incident (phi) and output (psi) field for each slice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_est = Fwd_Prop_MultiSlice_Intensity(i0, o_slice, k2, dz, P, H0);
   
    I_mea = I(:,:,m);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stopping criteria: residual between estimated & measured intensities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    err(m) = sum((I_mea(:)-I_est(:)).^2);
end

fprintf('| %2d   | %.2e |\n',iter,sqrt(sum(err)));

fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));

end

