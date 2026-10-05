function release(funcName, id, session)
    % Delete a library object from a class destructor.
    %
    % Objects issued before the library was last unloaded were freed by the unload,
    % and their handle values may since have been reused for unrelated objects, so
    % they are skipped. Checking `ct.isLoaded` also avoids restarting an unloaded
    % out-of-process library just to delete an object.

    if id < 0 || session ~= ct.impl.session() || ~ct.isLoaded()
        return
    end
    ct.impl.call(funcName, id);
end
