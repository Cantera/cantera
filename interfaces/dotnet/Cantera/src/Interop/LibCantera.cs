// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera.Interop;

static partial class LibCantera
{
    const string LibFile = "cantera_shared";

    public delegate void LogCallback(LogLevel logLevel, string category, string message);
}
