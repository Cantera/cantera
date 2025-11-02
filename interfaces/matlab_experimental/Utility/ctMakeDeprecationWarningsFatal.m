function v = ctMakeDeprecationWarningsFatal()
    % Turns deprecation warnings into exceptions. ::
    %
    %     >> ctMakeDeprecationWarningsFatal()

    ctIsLoaded(true);
    ctFunc('mCt_makeDeprecationWarningsFatal');
end
