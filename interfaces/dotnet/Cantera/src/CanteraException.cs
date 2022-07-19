// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Represents a error that occurred in the native Cantera library.
/// </summary>
public class CanteraException : ExternalException
{
    private CanteraException(string message) : base(message) { }

    unsafe internal static void ThrowLatest()
    {
        var errorMessage = InteropUtil.GetString(500, LibCantera.ct_getCanteraError);
        throw new CanteraException(errorMessage);
    }
}
