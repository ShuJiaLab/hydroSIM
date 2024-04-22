%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    This algorithm is for the paper:                     %
%    "hydroSIM: Super-resolution speckle illumination microscopy with a   %
%                             hydrogel diffuser"                          %
%            Please cite our paper on Biomedical Optics Express           %
%-------------------------------------------------------------------------%
%    This algorithm is written based on L.-H. Yeh et al work with major   %
%             modification, please also cite their work:                  %
%       "Structured illumination microscopy with unknown patterns and a   %
%            statistical prior," Biomed. Opt. Express 8695-711 (2017)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Path setup
clear all; close all; clc
img_dir = 'The dictionary containing your image';
cd(img_dir)
addpath('The dictionary containing the functions');
F  = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

%% Parameter setting
I_image_sum = tiffreadVolume("Stack.tif"); % The microtubule image of Fig. 3 (a) and (b), all parameter in the example code is consistent with the paper.
I_image_sum = double(I_image_sum);

lambda    = 0.515; % wavelength in microns
k         = 2*pi/lambda;
mag       = 100; % Magnification of objective lens
PixelSize = 6.5; % Physical pixel size of camera sensor (microns)
pscrop    = PixelSize/mag; % Effective pixels size (microns)
NA_obj    = 1.45;% NA of objective lens

%% Upsampling data

[Ncrop,Mcrop,Nimg] = size(I_image_sum);


N = Ncrop*2; M = Mcrop*2; ps = pscrop/2;

xh = (-M/2:(M/2-1)).*ps; 
yh = (-N/2:(N/2-1)).*ps; 
fx = (-M/2:(M/2-1))./(ps*M); 
fy = (-N/2:(N/2-1))./(ps*N); 
NAx = fx*lambda; 
NAy = fy*lambda;
[xhh,yhh] = meshgrid(xh,yh);
[fxx,fyy] = meshgrid(fx,fy);

dark_current = 100;
photon_count = 5000;

I_image_sum=I_image_sum/max(I_image_sum(:))*photon_count;



I_image_up = zeros(N,M,Nimg);
tic;
for i = 1:Nimg
    temp = max(0,I_image_sum(:,:,i) - 100);
    I_image_up(:,:,i) = abs(iF(padarray(F(temp),[(N-Ncrop)/2,(M-Mcrop)/2])));
end


%% PSF calculation/simulation
Pupil_obj = zeros(N,M);
r_obj     = (fxx.^2+fyy.^2).^(1/2); % spatial frequency of obj
Pupil_obj(r_obj<NA_obj/(lambda))=1;

T_incoherent = abs(F(abs(iF(Pupil_obj)).^2));
T_incoherent = T_incoherent/max(T_incoherent(:));

H_eff_final = padarray(interp2(T_incoherent.^3),[1,1],'post');
H_eff_final = H_eff_final(N/2+1:N*3/2,M/2+1:M*3/2);
H_eff_final = abs(H_eff_final/max(abs(H_eff_final(:))));

Mean_PSF  = abs(iF(T_incoherent));
Final_PSF = abs(iF(H_eff_final));

%% Deconvolution of Wide-Field Image
I_mean       = mean(I_image_up,3);
dec_itr_mean = 10;        % iteration for deconv.
dec_reg_mean = 1e-5;     % small constant to reduce artifact
max_itr      = 50;
I_mean_dec_t = deconvlucy(I_mean,Mean_PSF,dec_itr_mean,dec_reg_mean);

figure;
subplot 121
imshow(I_mean,[]); title 'Wide-Field Image'
subplot 122
imshow(I_mean_dec_t,[]); title 'Deconv. Wide-Field Image'
% use this step to find optimize the deconvolution iteration;

%% Pattern Estimation
reg_delta                  = 5*1e-3;
[Ip_est,I_mean,I_mean_dec] = PatternEstimation(I_image_up,T_incoherent,max_itr,reg_delta,I_mean_dec_t);

%% Intermediate Image Formation
range_ps = ceil(lambda/2/NA_obj/ps)+(~mod(ceil(lambda/2/NA_obj/ps),2));
[~,I_INT] = IntImage(Ip_est,I_image_up,range_ps);
%% Final Deconvolution
dec_itr_cov_rec = 10;
dec_reg_rec     = 1e-5;
reg_shading_pr  = 1e-3;
[SRImage,SRImage_Shading] = INT_DEC(Ip_est,gather(I_INT),reg_shading_pr,dec_itr_cov_rec,Final_PSF,dec_reg_rec);

figure;
subplot 131
imshow(I_INT,[]); title 'INT image'
subplot 132
imshow(SRImage,[]); title 'Super-resolution Image'
subplot 133
imshow(SRImage_Shading,[]); title 'Shaded Super-resolution Image'