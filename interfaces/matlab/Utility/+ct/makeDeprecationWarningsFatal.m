function v = makeDeprecationWarningsFatal()
    % Turns deprecation warnings into exceptions. ::
    %
    %     >> ct.makeDeprecationWarningsFatal()

    isLoaded(true);
    ctFunc('mCt_makeDeprecationWarningsFatal');
end
