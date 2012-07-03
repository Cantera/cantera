function v = isIdealGas(a)
% ISIDEALGAS - True if the phase is an ideal gas or ideal gas
% mixture, and false otherwise.
%
if eosType(a) == 1
    v = 1;
else
    v = 0;
end
