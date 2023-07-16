clc
clear all
close all
%%

% Tent weighting is used, according to the paper by Debevec et al.
% Uniform (=1), Gaussian(=3) or Photon(=3) weighting may also be used,
% however the results are not that good in these cases.
weightfun_used = 2;
namefun = {"Uniform","Tent","Gaussian","Photon"};
namechan = {"Red","Green","Blue"};
exposuretimes1 =  [1/2500, 1/1000, 1/500, 1/250, 1/125, 1/60, 1/30, 1/15, 1/8, 1/4, 1/2, 1, 2, 4, 8, 15];
exposuretimes2 =  [1/400, 1/250, 1/100, 1/40, 1/25, 1/8, 1/3];
exposuretimes = {exposuretimes1,exposuretimes2};
Image2_names = {'sample2-00.jpg','sample2-01.jpg','sample2-02.jpg','sample2-03.jpg','sample2-04.jpg','sample2-05_rotated.jpg','sample2-06.jpg'};
numimgs = [length(exposuretimes1),length(exposuretimes2)];
numsets = 2;

%Initialize output structures
responseCurve = zeros(256,3,numsets);
radiancemaps_cell = cell(1,numsets);
output_imgs = cell(1,2);

%Parameters for each set
resize_factor = [1/64, 1/32];
Zmin = round(0.05*255);
Zmax = round(0.99*255);
lamdas = [300,100];
gammas = [1,0.8];

%%

for set = 1:numsets
    
    % Read first image from set to get some parameters (sizes, channels)
    if(set==1)
        im1 = imread(sprintf('%s%s%s','exposure',num2str(1),'.jpg'));
    else
        im1 = imread(Image2_names{1});
    end
    M = size(im1,1);
    N = size(im1,2);
    K = length(exposuretimes{set});
    chans = size(im1,3);
    %%
    
    % Parse image set into Q matrix. 
    % Images are assumed to be directly accessible, either being in the 
    % same folder, or in PATH.
    Q = zeros(M,N,chans,numimgs(set),'like',im1);
    Q(:,:,:,1) = im1;
    if (set == 1)
        for i = 2:numimgs(set)
            Q(:,:,:,i) = imread(sprintf('%s%s%s','exposure',num2str(i),'.jpg'));
        end
    else
        for i = 2:numimgs(set)
            Q(:,:,:,i) = imread(Image2_names{i});
        end
    end
    %%
    % Initialize radiance map 
    radiancemap = zeros(M,N,chans);
    
    % Estimate response curve for each channel separately
    figure("windowstate","maximized");
    for c = 1:chans
        % Stack of images of the same channel
        imgStack = double(squeeze(Q(:,:,c,:)));
        
        % Estimate G(Z) for the current channel of the image set
        responseCurve(:,c,set) = estimateResponseCurve(Q(:,:,c,:),exposuretimes{set},lamdas(set),weightfun_used,resize_factor(set),Zmin,Zmax);
        
        % Plot response function
        subplot(2,2,c);
        mat = imgStack(imgStack<=Zmax & imgStack>=Zmin);
        [MM,NN,KK] = size(mat);
        p1 = plot(reshape(responseCurve(mat+1,c,set),[MM*NN*KK, 1]),reshape(mat,[MM*NN*KK, 1]),'o','LineWidth',0.01);
        hold on
        plot(responseCurve(Zmin+1:Zmax+1,c,set),Zmin:Zmax,'k','LineWidth',2);
        title({
            ["Plot of the estimation of the response curve"]
            ["for the "+namechan{c}+" channel for the image set \#"+num2str(set)]
            },"interpreter","latex");
        ylabel("pixel value $Z$","Interpreter","Latex");
        xlabel("log exposure $X$","Interpreter","Latex");

        %%
        % Repeat mergeLDRstack routine in a (slightly...) more compact and suitable way       
        WZ = zeros(M,N,K);
        [r,cc,v] =ind2sub(size(imgStack), find(imgStack<=Zmax & imgStack>=Zmin));
        switch weightfun_used
            case 1
                w = @(z) (1);
                for i = 1:length(r)
                     WZ(r(i),cc(i),v(i)) =  w(imgStack(r(i),cc(i),v(i)));
                end
            case 2
                w = @(z)(min(z,255-z));
                for i = 1:length(r)
                     WZ(r(i),cc(i),v(i)) =  w(imgStack(r(i),cc(i),v(i)));
                end
            case 3
                w = @(z)(exp(-4*(z-255/2).^2/(255/2)^2) );
                for i = 1:length(r)
                     WZ(r(i),cc(i),v(i)) =  w(imgStack(r(i),cc(i),v(i)));
                end
            case 4
                w = @(k)(exposuretimes{set}(k));                
                for i = 1:length(r)
                     WZ(r(i),cc(i),v(i)) =  w(v(i));
                end
        end
        T = zeros(M,N,numimgs(set));
        for i = 1:K
            T(:,:,i) = repmat(log(exposuretimes{set}(i)),M,N);
        end
        G = reshape( responseCurve(Q(:,:,c,:)+1,c,set) ,[M,N,K]);     
        radiancemap(:,:,c) = rescale(sum((WZ.*(G-T)),3)./sum(WZ,3));
        
    end
    
    % Plot for the response curves for each channel in the same axes
    subplot(2,2,4);
    plot(responseCurve(Zmin+1:Zmax+1,1,set),Zmin:Zmax,'r');
    hold on
    plot(responseCurve(Zmin+1:Zmax+1,2,set),Zmin:Zmax,'g--');
    plot(responseCurve(Zmin+1:Zmax+1,3,set),Zmin:Zmax,'b-.');
    ylabel("pixel value $Z$","Interpreter","Latex");
    xlabel("log exposure $X$","Interpreter","Latex");
    title({
            ["Plot of the estimation of the response curve"]
            ["for the all chanels on the same axes for the image set \#"+num2str(set)]
            },"interpreter","latex");   
    
    % Put the resulting radiance map in a cell matrix
    radiancemaps_cell{set} = radiancemap;
    % Use tone mapping
    output_imgs{set} = toneMapping(radiancemap,gammas(set));
    if(set == 1)  
        
        % Check palette pixels for the 1st image
        y_check = [235,290,345,400,445,510];
        x_check = 1330*ones(size(y_check));
        
        figure();
        for j = 1:length(y_check)
            plot(j,rgb2gray(output_imgs{set}(y_check(j),x_check(j),:)),'x');
            hold on
        end
        line([1,length(y_check)],[rgb2gray(output_imgs{set}(y_check(1),x_check(1),:)), rgb2gray(output_imgs{set}(y_check(end),x_check(end),:))]);
        title({
            ["Intensities of the palette pixels for $\gamma$ = 1"]
            ["and Tent weighting"]
            },"interpreter","latex");
    end
    
    % Display and save the final image
    figure();
    imshow(output_imgs{set});
    title("HDR image created from imageset \#"+num2str(set),"interpreter","latex");
    imwrite(output_imgs{set},sprintf('%s%s%s','demo3_output_',num2str(set),'.bmp'));
end
% save("responseCurves.mat","responseCurve");
% save("radianceMaps.mat","radiancemaps_cell");