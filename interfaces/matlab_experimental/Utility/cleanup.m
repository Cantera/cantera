function cleanup()
% Delete all stored Cantera objects and reclaim memory
    checklib;
    calllib('cantera_shared', 'ct_clearOneDim');
    calllib('cantera_shared', 'ct_clearMix');
    calllib('cantera_shared', 'ct_clearXML');
    calllib('cantera_shared', 'ct_clearFunc');
    calllib('cantera_shared', 'ct_clearStorage');
    calllib('cantera_shared', 'ct_clearReactors');
    calllib('cantera_shared', 'ct_clearReactionPath');
    clear all
end
