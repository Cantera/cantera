function s = GRI30(tr)
% GRI30  Create an object with the GRI-Mech 3.0 reaction mechanism.
% s = GRI30(tr)
% Create a Solution instance representing
% reaction mechanism GRI-Mech 3.0.
%
% GRI-Mech 3.0 is a widely-used reaction mechanism for natural gas
% combustion. It contains 53 species composed of the elements H,
% C, O, N, and/or Ar, and 325 reactions, most of which are
% reversible. GRI-Mech 3.0, like most combustion mechanisms, is
% designed for use at pressures where the ideal gas law holds.
% GRI-Mech 3.0 is available from http://www.me.berkeley.edu/gri_mech/
%
% Function :mat:func:`GRI30` creates the solution according to the
% specifications in file gri30.yaml. The ideal gas equation of
% state is used. Transport property evaluation is mixture-averaged by
% default. To change or disable transport properties, supply the name of
% the transport model to use.
%
% .. code-block:: matlab
%
%     g1 = GRI30           % mixture-averaged transport properties
%     g2 = GRI30('Mix')    % mixture-averaged transport properties
%     g3 = GRI30('Multi')  % miulticomponent transport properties
%     g4 = GRI30('None')   % no transport properties
%
% :param tr:
%     Transport modeling, ``'None'``, ``'Mix'``, or ``'Multi'``
% :return:
%     Instance of class :mat:func:`Solution`
%

if nargin == 0
    s = Solution('gri30.yaml', 'gri30');
elseif nargin == 1
    if strcmp(tr, 'None')
        s = Solution('gri30.yaml', 'gri30', 'None');
    elseif strcmp(tr, 'Mix')
        s = Solution('gri30.yaml', 'gri30', 'Mix');
    elseif strcmp(tr, 'Multi')
        s = Solution('gri30.yaml', 'gri30', 'Multi')
    else
        error('Unknown transport specified. "None", "Mix", or "Multi" are supported.')
    end
else
    error('Wrong number of arguments.');
end
