function v = ctMakeDeprecationWarningsFatal()
    % Turns deprecation warnings into exceptions. ::
    %
    %     >> ctMakeDeprecationWarningsFatal()

    ctIsLoaded;
    ctFunc('mCt_makeDeprecationWarningsFatal');
end
