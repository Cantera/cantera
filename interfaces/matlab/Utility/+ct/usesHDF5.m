function v = usesHDF5()
    % Returns true if Cantera was compiled with HDF5 support. ::
    %
    %     >> ct.usesHDF5()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera.

    ct.isLoaded(true);
    v = logical(ct.impl.call('mCt_usesHDF5'));
end
