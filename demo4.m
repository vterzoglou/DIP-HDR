clc
close all
clear all
%%
% Fix rotated image based on img #6
img1 = imread('sample2-06.jpg');
img2 = imread('sample2-05_rotated.jpg');

% Grayscale images from original ones
img1g = im2double(rgb2gray(img1));
img2g = im2double(rgb2gray(img2));

% Detect salient points
points1 = detectSURFFeatures(img1g,'MetricThreshold',750,'NumOctaves',3,'NumScaleLevels',5);
points2 = detectSURFFeatures(img2g,'MetricThreshold',750,'NumOctaves',3,'NumScaleLevels',5);

% Feature extraction and matching
[feat1,vp1] = extractFeatures(img1g,points1,'Upright',false);
[feat2,vp2] = extractFeatures(img2g,points2,'Upright',false);
idxp = matchFeatures(feat1,feat2,'MatchThreshold',10.5422,'MaxRatio',0.1054);
mp1 = vp1(idxp(:,1),:);
mp2 = vp2(idxp(:,2),:);

% Estimate transformation
tform = estimateGeometricTransform(mp2,mp1,'affine');

% Display  and save recovered image
outputView = imref2d(size(img1));
Ir = imwarp(img2,tform,'OutputView',outputView);
figure
imshow(Ir);
title('Recovered Image');

imwrite(Ir,'sample2-05.jpg');

%%
% Repeat demo3 algorithm...

weightfun_used = 2;
exposuretimes =  [1/400, 1/250, 1/100, 1/40, 1/25, 1/8, 1/3];
Image2_names = {'sample2-00.jpg','sample2-01.jpg','sample2-02.jpg','sample2-03.jpg','sample2-04.jpg','sample2-05.jpg','sample2-06.jpg'};
numimgs = length(exposuretimes);
resize_factor = 1/32;
responseCurve = zeros(256,3);

im1 = imread(Image2_names{1});
M = size(im1,1);
N = size(im1,2);
K = length(exposuretimes);
chans = size(im1,3);

Q = zeros(M,N,chans,numimgs,'like',im1);
Q(:,:,:,1) = im1;
for i = 2:numimgs
    Q(:,:,:,i) = imread(Image2_names{i});
end
calibrated = zeros(size(Q));
radiancemap = zeros(M,N,chans);

% Estimate response curve for each channel separately
lamda = 100;
Zmin = round(0.05*255);
Zmax = round(0.99*255);
for c = 1:chans
    responseCurve(:,c) = estimateResponseCurve(Q(:,:,c,:),exposuretimes,lamda,weightfun_used,resize_factor,Zmin,Zmax);

    % Repeat mergeLDRstack routine in a (slightly...) more compact and suitable way 
    imgStack = squeeze(Q(:,:,c,:));
    WZ = zeros(M,N,K);
    
    w = @(z)(min(z,255-z));
    [r,cc,v] =ind2sub(size(imgStack), find(imgStack<=Zmax & imgStack>=Zmin));
    for i = 1:length(r)
         WZ(r(i),cc(i),v(i)) =  w(imgStack(r(i),cc(i),v(i)));
    end
    T = zeros(M,N,numimgs);
    for i = 1:numimgs
        T(:,:,i) = repmat(log(exposuretimes(i)),M,N);
    end
    G = reshape( responseCurve(Q(:,:,c,:)+1,c) ,[M,N,K]);
    radiancemap(:,:,c) = rescale(sum((WZ.*(G-T)),3)./sum(WZ,3));
end

% Use tone mapping, display and save the final image
output_img = toneMapping(radiancemap,0.8);

figure();
imshow(output_img);
title("HDR image created from imageset \#2","interpreter","latex");
imwrite(output_img,'demo4_output.bmp');
