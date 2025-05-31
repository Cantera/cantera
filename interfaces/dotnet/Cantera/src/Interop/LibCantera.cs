// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;
using System.Runtime.InteropServices.Marshalling;

namespace Cantera.Interop;

static partial class LibCantera
{
    const string LibFile = "cantera_shared";

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
