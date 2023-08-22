// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera.Interop;

static partial class LibCantera
{
    public delegate void LogCallback(LogLevel logLevel, string category, string message);

    const string LibFile = "cantera_shared";

    /// <summary>
    /// Represents the delegate that is marshalled to LibCantera as a function pointer.
    /// </summary>
    /// <remarks>
    /// ct_setLogWriter() needs a delegate which is marshalled as a function pointer to
    /// the C++ Cantera library. We could create one implicitly when calling
    /// <c>LibCantera.ct_setLogWriter(LogMessage)</c>, but the garbage collector would
    /// not know the native method is using it and could reclaim it. By explicitly
    /// storing it as a class member, we ensure it is not collected until the class is.
    /// </remarks>
    static readonly LogCallback s_invokeMessageLoggedDelegate;

    static LibCantera()
    {
        s_invokeMessageLoggedDelegate = Application.RaiseMessageLogged;
        InteropUtil.CheckReturn(ct_setLogCallback(s_invokeMessageLoggedDelegate));
    }
}
