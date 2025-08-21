function success = test_examples()
    ctLoad

    examples = {
        'equil', 'isentropic', 'reactor1', 'reactor2', 'surf_reactor', ...
        'periodic_cstr', 'plug_flow_reactor', 'lithium_ion_battery', ...
        'rankine', 'prandtl1', 'prandtl2', 'catcomb','ignite', 'diamond_cvd', ...
        % 'flame1', 'flame2', 'diff_flame', 'ignite_hp', 'ignite_uv'
    };

    success = true;
    passed = {};
    failed = {};

    for idx = 1:length(examples)
        scriptName = examples{idx};
        try
            evalin('base', sprintf('run(''%s.m'')', scriptName));
            passed{end+1} = scriptName;
        catch ME
            fprintf('An error occurred while running %s: %s\n', scriptName, ME.message);
            fprintf('Identifier: %s\n', ME.identifier);
            if strcmp(ME.identifier, 'Cantera:ctError')
                disp('Caught a CanteraError. Continuing execution...\n');
            end
            failed{end+1} = scriptName;
            success = false;
        end
    end


    % Summary report
    disp(' ');
    disp('============================');
    disp('Summary of example runs');
    disp('============================');

    if ~isempty(passed)
        fprintf('✅ Passed: %s\n', strjoin(passed, ', '));
    else
        disp('✅ Passed: (none)');
    end

    if ~isempty(failed)
        fprintf('❌ Failed: %s\n', strjoin(failed, ', '));
    else
        disp('❌ Failed: (none)');
    end
    ctUnload

    disp('Test example run complete.');
end