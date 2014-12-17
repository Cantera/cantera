function n = elementIndex(self, name)
% ELEMENTINDEX  Get the index of an element.
% n = elementIndex(self, name)
% :param self:
%     Instance of class :mat:func:`Mixture`
% :param name:
%     Name of the element whose index is desired
% :return:
%     Index of element with name ``name``
%

n = mixturemethods(22, mix_hndl(self), name);
