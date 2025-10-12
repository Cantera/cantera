function v = ctUsesHDF5()
    % Returns true if Cantera was compiled with HDF5 support. ::
    %
    %     >> ctUsesHDF5()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera.

    ctIsLoaded;
    v = logical(ctFunc('ct_usesHDF5'));
end
