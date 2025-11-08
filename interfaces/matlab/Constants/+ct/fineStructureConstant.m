function r = fineStructureConstant
    % Get the fine-structure constant (dimensionless). ::
    %
    %     >> r = ct.fineStructureConstant
    %
    % :return:
    %     The fine-structure constant (dimensionless).

    r = ct.impl.call('mCt_fineStructureConstant');
end
