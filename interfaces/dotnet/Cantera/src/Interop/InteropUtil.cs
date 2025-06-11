// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera.Interop;

static class InteropUtil
{
    /// <summary>
    /// Checks the return code of the call into the native Cantera library
    /// and throws a <see cref="CanteraException"/> if necessary
    /// </summary>
    /// <remarks>
    /// This method is automatically called by P/Invoke return type marshallers.
    /// You should not need to call it elsewhere.
    /// </remarks>
    public static int CheckReturn(int code)
    {
        // Cantera returns this value when the function resulted in an
        // error internal to Cantera
        const int Error1 = -1;
        // Sometimes Cantera returns this value when the function resulted in an
        // error internal to Cantera and related to an object handle
        const int Error2 = -2;
        // Cantera returns this value when the function resulted in an external error
        // Some functions also return a negative value as the amount of space they need
        // to fill a buffer with a string. There is no way to account for the ambiguity
        // that arises when such a function returns -999!
        const int Error999 = -999;

        CallbackException.ThrowIfAny();

        if (code is Error1 or Error2 or Error999)
        {
            CanteraException.ThrowLatest();
        }

        return code;
    }

    public static int GetInteropBool(bool value) =>
        value ? InteropConsts.True : InteropConsts.False;
}
