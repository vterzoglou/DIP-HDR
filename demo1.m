% Clear console, turn off some warnigns
clc
clear all
close all
warning('off','imageio:tifftagsread:expectedTagDataFormatMultiple');
warning('off','imageio:tifftagsread:expectedTagDataFormat');
%%

% User input of folder with images
fprintf("Please provide folder of images.\n");
imgFolderInfo = dir(uigetdir);
imgFolderInfo = imgFolderInfo(~ismember({imgFolderInfo.name},{'.','..'}));
imgFilePaths = strcat(imgFolderInfo(1).folder,'\',{imgFolderInfo(:).name}');
clear imgFolderInfo

% Extract some info for images (sizes, number of channels)
numimgs = size(imgFilePaths,1);
exposures = zeros(numimgs,1);
for i = 1:numimgs
    exposures(i) = imfinfo(imgFilePaths{i}).DigitalCamera.ExposureTime;
end

M = imfinfo(imgFilePaths{1}).Height;
N = imfinfo(imgFilePaths{1}).Width;
chans = imfinfo(imgFilePaths{1}).NumberOfSamples;

%Fill 4D matrix containing each image(3D) for each exposure(4th Dimension)
Q = zeros(M,N,chans,numimgs);
for i = 1:numimgs
    Q(:,:,:,i) = imread(imgFilePaths{i});
end
Q = rescale(Q);

%%
namefun = {"Uniform","Tent","Gaussian","Photon"};
namechannels = {"Red","Green","Blue"};
numfun = length(namefun);

%output initialization
maps = zeros(M,N,chans,numfun);

for c = 1:chans
    % Set up two figures per channel; one for irradiance visualization and
    % one for plotting histogram.
    figure("windowstate","maximized");
    tcl1 = tiledlayout(2,2);
    title(tcl1,"Channel: "+namechannels{c});
    figure("windowstate","maximized");
    tcl2 = tiledlayout(2,2);
    title(tcl2,"Channel: "+namechannels{c});
    for f = 1:numfun
        
        %For each channel and function get the output map and rescale the
        %values to [0,1]
        maps(:,:,c,f) = rescale(mergeLDRstack(squeeze(Q(:,:,c,:)),exposures,f));
        
        % Plot the visualization of the irradiance for each weighting
        % function
        figure(2*c-1);
        nexttile
        imagesc(maps(:,:,c,f))
        colorbar
        title("Weight Function: "+namefun(f));
        
        % Plot the histogram of values for each weighting function
        figure(2*c);
        nexttile
        histogram(maps(:,:,c,f))
        title("Weight Function: "+namefun(f));
    end
    
end
save("maps.mat","maps");
