function cleanup()
    % Delete all stored Cantera objects and reclaim memory.
    %
    checklib;
    calllib(ct, 'ct_clearOneDim');
    calllib(ct, 'ct_clearMix');
    calllib(ct, 'ct_clearFunc');
    calllib(ct, 'ct_clearStorage');
    calllib(ct, 'ct_clearReactors');
    calllib(ct, 'ct_clearReactionPath');
    clear all
end
