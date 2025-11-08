function v = makeDeprecationWarningsFatal()
    % Turns deprecation warnings into exceptions. ::
    %
    %     >> ct.makeDeprecationWarningsFatal()

    ct.isLoaded(true);
    ct.impl.call('mCt_makeDeprecationWarningsFatal');
end
