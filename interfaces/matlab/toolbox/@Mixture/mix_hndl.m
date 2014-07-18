function i = mix_hndl(self)
% MIX_HNDL  Get the integer used to access the kernel object.
% i = mix_hndl(self)
% :param self:
%     Instance of :mat:func:`Mixture`
% :return:
%     Integer used to access the kernel object
%

i = self.mixindex;
