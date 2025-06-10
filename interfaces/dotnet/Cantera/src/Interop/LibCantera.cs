// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;
using System.Runtime.InteropServices.Marshalling;

namespace Cantera.Interop;

/// <summary>
/// Contains P/Invoke signatures for calling the native Cantera library.
/// </summary>
/// <remarks>
/// Each P/Invoke signature represents a single function exported by the
/// native Cantera library C API. Signatures are marked with LibraryImportAttribute
/// to facilitate generation of marshalling code at build-time. In addition,
/// functions that get strings or get or set arrays are marked private and
/// exposed via wrappers which simplify calling these functions.
/// </remarks>
static partial class LibCantera
{
    const string LibFile = "cantera_shared";
    const int BufferSize = 360;

    static void ThrowOnBadString() =>
        throw new InvalidOperationException("Could not retrieve a string value from Cantera!");

    public delegate void LogCallback(LogLevel logLevel,
                                     [MarshalAs(UnmanagedType.LPUTF8Str)] string category,
                                     [MarshalAs(UnmanagedType.LPUTF8Str)] string message);

    [CustomMarshaller(typeof(int), MarshalMode.ManagedToUnmanagedOut, typeof(ReturnCodeChecker))]
    [CustomMarshaller(typeof(double), MarshalMode.ManagedToUnmanagedOut, typeof(ReturnCodeChecker))]
    static class ReturnCodeChecker
    {
        public static int ConvertToManaged(int value) =>
            InteropUtil.CheckReturn(value);

        public static double ConvertToManaged(double value)
        {
            // Cantera returns this value when the function resulted in error
            const double Error = -999.999;

            CallbackException.ThrowIfAny();

            if (value == Error)
            {
                CanteraException.ThrowLatest();
            }

            return value;
        }
    }
}
