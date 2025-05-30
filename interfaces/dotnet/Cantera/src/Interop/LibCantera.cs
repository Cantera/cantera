// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;

namespace Cantera.Interop;

static partial class LibCantera
{
    const string LibFile = "cantera_shared";

    public unsafe delegate void LogCallback(LogLevel logLevel, byte* category, byte* message);

    [LibraryImport(LibFile, StringMarshalling = StringMarshalling.Utf8)]
    private static partial string yo();
}
