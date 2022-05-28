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
