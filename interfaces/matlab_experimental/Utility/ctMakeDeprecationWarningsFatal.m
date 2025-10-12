function v = ctMakeDeprecationWarningsFatal()
    % Turns deprecation warnings into exceptions. ::
    %
    %     >> ctMakeDeprecationWarningsFatal()

    ctIsLoaded;
    ctFunc('ct_makeDeprecationWarningsFatal');
end
