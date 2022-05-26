using System.Runtime.InteropServices;

namespace Cantera;

/// <summary>
/// Represents a error that occurred in the native Cantera library.
/// </summary>
public class CanteraException : ExternalException
{
    static class LibCantera
    {
        [DllImport(Interop.LIBCANTERA)]
        unsafe public static extern int ct_getCanteraError(int buflen, byte* buf);
    }

    private CanteraException(string message) : base(message) { }

    unsafe internal static void ThrowLatest()
    {
        var errorMessage = Interop.GetString(500, LibCantera.ct_getCanteraError);
        throw new CanteraException(errorMessage);
    }
}