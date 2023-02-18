clc
clear all
close all

namefun = {"Uniform","Tent","Gaussian","Photon"};
namechannels = {"Red","Green","Blue"};
exposuretimes =  [1/2500, 1/1000, 1/500, 1/250, 1/125, 1/60, 1/30, 1/15, 1/8, 1/4, 1/2, 1, 2, 4, 8, 15];
numfun = length(namefun);

%used to test how many (and which) images are needed to extract a good a
%estimation for the final HDR image.
imgsUsed = [1:16];
numimgs = length(exposuretimes(imgsUsed));

im1 = imread(sprintf('%s%s%s','exposure',num2str(imgsUsed(1)),'.jpg'));
M = size(im1,1);
N = size(im1,2);
chans = size(im1,3);

%Fill 4D matrix containing each image(3D) for each exposure(+1D)
Q = zeros(M,N,chans,numimgs);
Q(:,:,:,1) = im1;
for i = 2:numimgs
    Q(:,:,:,i) = imread(sprintf('%s%s%s','exposure',num2str(imgsUsed(i)),'.jpg'));
end
Q = rescale(Q);

%output initialization
maps = zeros(M,N,chans,numfun);

for c = 1:chans
    figure("windowstate","maximized");
    suptitle("Channel: "+namechannels{c});
    figure("windowstate","maximized");
    suptitle("Channel: "+namechannels{c});
    for f = 1:numfun
        
        %For each channel and function get the output map and rescale the
        %values to [0,1]
        maps(:,:,c,f) = rescale(mergeLDRstack(squeeze(Q(:,:,c,:)),exposuretimes(imgsUsed),f));
        
        figure(2*c-1);
        subplot(2,2,f);
        imagesc(maps(:,:,c,f))
        colorbar
        title("Weight Function: "+namefun(f));
        
        figure(2*c);
        subplot(2,2,f);
        histogram(maps(:,:,c,f))
        title("Weight Function: "+namefun(f));
    end
end
save("maps.mat","maps");
