// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera.Interop;

static partial class LibCantera
{
#if IS_WINDOWS
    const string LibFile = "cantera.3.0.0.dll";
#elif IS_MACOS
    const string LibFile = "cantera.3.0.0.dylib";
#elif IS_LINUX
    const string LibFile = "cantera.3.0.0.so";
#endif

    public delegate void LogCallback(LogLevel logLevel, string category, string message);
}
