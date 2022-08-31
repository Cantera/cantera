function cleanup()
    % CLEANUP
    % Delete all stored Cantera objects and reclaim memory.
    %
    checklib;
    callct('ct_clearOneDim');
    callct('ct_clearMix');
    callct('ct_clearFunc');
    callct('ct_clearStorage');
    callct('ct_clearReactors');
    callct('ct_clearReactionPath');
    clear all
end
