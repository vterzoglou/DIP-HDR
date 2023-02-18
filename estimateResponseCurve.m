
function responseCurve = estimateResponseCurve(imgStack,exposureTimes,smoothingLamda,weightingFcn,resize_factor,Zmin,Zmax)
    %imgstack:MxNxnumimgs
    M = size(imgStack,1);
    N = size(imgStack,2);
    numimgs = length(exposureTimes);
    
    B = log(exposureTimes);
    
    %Resize images and convert subscripts for spatial dimensions to a
    %linear index
    newM = ceil(M*resize_factor);
    newN = ceil(N*resize_factor);
    newQ = zeros(newM*newN, numimgs);
    for i = 1:numimgs
       newQ(:,i) = double( reshape( imresize(imgStack(:,:,i),resize_factor,'method','nearest'), [newM*newN, 1] ) ); 
    end
  
    
    zvals = 0:255;
    %flag is passed to gsolve to notify if photon weighting is used, in
    %order to make some changes...
    flag = 0;
    switch weightingFcn
        case 1
            W = ones(1,256);
        case 2
            W = min(zvals,255-zvals);
        case 3
            W = (exp(-4*(zvals-mean(zvals)).^2/mean(zvals)^2) );
        case 4
            W = exposureTimes;
            flag = 1;
    end
    %Set W values to Z outside of the interval [Zmin,Zmax]
    W([1:Zmin+1,Zmax+1:256]) = 0;
    [responseCurve,~] = gsolve(newQ,B,smoothingLamda,W,flag);

end