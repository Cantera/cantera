function i = xml_hndl(x)
% XML_HNDL  Get the integer used to access the kernel object.
% i = xml_hndl(x)
%
% :param x:
%     Instance of class :mat:func:`XML_Node`
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object
%

i = x.id;
