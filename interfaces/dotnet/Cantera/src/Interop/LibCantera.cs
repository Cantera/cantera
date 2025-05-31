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
}
