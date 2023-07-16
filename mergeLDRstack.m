function radiancemap = mergeLDRstack(imgStack, exposureTimes, weightingFcn)

    assert(size(imgStack,3) == length(exposureTimes),'Number of images must match number of exposure times');
    K = length(exposureTimes);
    M = size(imgStack,1);
    N = size(imgStack,2);
    
    %function may be used with weightFcn being a function handle
    %alternatively
    if(class(weightingFcn) == "double")
        Zmin = 0.05;
        Zmax = 0.99;
        
        %Set up WZ matrix containing the values of the weighting function
        WZ = zeros(M,N,K);
        % Find the indices (r,c,v) of the 3D array imgStack where its
        % values lie between Zmin and Zmax
        [r,c,v] =ind2sub(size(imgStack), find(imgStack<=Zmax & imgStack>=Zmin));
        
        % For the indices found above, fill in the weight matrix WZ,
        % according to the weight function. The remaining indices are left
        % unchanged (equal to zero).
        switch weightingFcn
            case 1
    %             W = @Wuniform;
                W =@(z)1;
                for i = 1:length(r)
                     WZ(r(i),c(i),v(i)) =  W(imgStack(r(i),c(i),v(i)));
                end
            case 2
    %             W = @Wtent;
                W = @(z)(min(z,255-z));
                for i = 1:length(r)
                     WZ(r(i),c(i),v(i)) =  W(imgStack(r(i),c(i),v(i)));
                end
            case 3
    %             W = @Wgaussian;
                W = @(z)(exp(-4*(z-0.5).^2/0.5^2));
                for i = 1:length(r)
                     WZ(r(i),c(i),v(i)) =  W(imgStack(r(i),c(i),v(i)));
                end
            case 4
    %             W = @Wphoton;
                W = @(k)(exposureTimes(k));
                for i = 1:length(r)
                     WZ(r(i),c(i),v(i)) =  W(v(i));
                end
        end
    
    elseif(class(weightingFcn) == "function_handle")
        WZ = weightingFcn(imgStack);
    end
        
%     radiancemap = zeros(M,N);

    % Set up matrix T(3D) containing log of exposure times for each pixel
    % and each exposure
    T = zeros(M,N,K);
    for i = 1:K
        T(:,:,i) = repmat(log(exposureTimes(i)),M,N);
    end
    
    
    % Set up the matrix G containing the log of the pixels' values
    
    % For the values Z == 0, it holds that log(Z) == infty, but also,
    % W(Z)=0. We don't want these values to be used in the calculation (and
    % result in infty), therefore we manually set their value of G(Z) equal
    % to 0.
    G = log(imgStack);
    idx = isinf(G);
    G(idx) = 0;
    
    % Calculate the transformation to radiance map
    radiancemap = sum((WZ.*(G-T)),3)./sum(WZ,3);
    
    
    % Find any pixels for which the weights were all zero, so consequently
    % the denominator would be 0.
    [r,c] = find(sum(WZ,3) == 0);
    
    % (Forunately) Matlab excludes infinite values from min/max, then we
    % can find the min and max values of the rest of the pixels.
    m = min(min(radiancemap));
    M = max(max(radiancemap));
    mid = (m+M)/2;
    
    % For each pixel corresponding to weights = 0 for all exposures
    for i = 1:length(r)
        
        % If the pixel is always underexposed, set its output value to the
        % min of the outputs.
        if( all(imgStack(r(i),c(i),:)<Zmin) )
            radiancemap(r(i),c(i)) = m;
            
        % If the pixel is always overexposed, set its output value to the
        % max of the outputs.
        elseif( all(imgStack(r(i),c(i),:)>Zmax) )
            radiancemap(r(i),c(i)) = M;
            
        % If for some (strange) reason the pixel doesnt correspond to the
        % cases above (ie the pixel sometimes is underexposed and some
        % other times is overexposed, but never well exposed...) set its
        % output value to the mean of max and min... 
        else
            radiancemap(r(i),c(i)) = mid;
        end
    end
        

%     radiancemap = exp(radiancemap);
end