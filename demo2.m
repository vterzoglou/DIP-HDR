%demo2 is supposed to be run after demo1, so the variable maps will
%probably have been calculated.
%If for whatever reason this is not the case, maps is calculated repeating
%most steps from demo1.

clearvars -except maps
clc
close all

%if variable maps not existing
if( exist('maps','var') == 0)
    %if the maps file exists load from there
    if( exist('maps','file') == 2)
        data = load('maps.mat');
        maps = data.maps;
        clearvars -except maps
    %else repeat demo1 steps to calculate maps
    else
        numimgs = 16;
        numfun = 4;
        exposuretimes =  [1/2500, 1/1000, 1/500, 1/250, 1/125, 1/60, 1/30, 1/15, 1/8, 1/4, 1/2, 1, 2, 4, 8, 15];
        im1 = imread(sprintf('%s%s%s','exposure',num2str(1),'.jpg'));
        M = size(im1,1);
        N = size(im1,2);
        chans = size(im1,3);
        Q = zeros(M,N,chans,numimgs);
        Q(:,:,:,1) = im1;
        for i = 2:numimgs
            Q(:,:,:,i) = imread(sprintf('%s%s%s','exposure',num2str(i),'.jpg'));
        end
        Q = rescale(Q);
        maps = zeros(M,N,chans,numfun);
        for c = 1:chans
            for f = 1:numfun
                maps(:,:,c,f) = mergeLDRstack(squeeze(Q(:,:,c,:)),exposuretimes,f);
            end
        end
        clearvars -except maps
    end
end

%points on palette to check for
y_check = [235,290,345,400,445,510];
x_check = 1330*ones(size(y_check));

namefun = {"Uniform","Tent","Gaussian","Photon"};
g = 1.1;

%output images
hdr_images = toneMapping(maps,g);
for i = 1:size(namefun,2)
    %show toned image
    figure("windowstate","maximized");
    imagesc(hdr_images(:,:,:,i));
    title({"After gamma correction, using $\gamma$ = "+num2str(g)+" and "+namefun{i}+" weighting"},"interpreter","latex");
    
    %show brightness curve on palette pixels
    figure();
    for j = 1:length(y_check)
        plot(j,rgb2gray(hdr_images(y_check(j),x_check(j),:,i)),'x');
        hold on
    end
    line([1,length(y_check)],[rgb2gray(hdr_images(y_check(1),x_check(1),:,i)), rgb2gray(hdr_images(y_check(end),x_check(end),:,i))]);
    title({
            ["Intensities of the palette pixels for $\gamma$ = "+num2str(g)]
            ["and "+namefun{i}+" weighting"]},"interpreter","latex");
end

% pretty hefty file, better delete it :P
delete maps.mat