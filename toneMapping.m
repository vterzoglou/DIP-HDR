function tonedImage = toneMapping(radianceMap,gamma)
    tonedImage = uint8(rescale(radianceMap.^gamma,0,255));
    
end