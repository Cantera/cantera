// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.CompilerServices;
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
        var errorMessage = InteropUtil.GetString(500, LibCantera.ct3_getCanteraError);
        throw new CanteraException(errorMessage);
    }
}

/// <summary>
/// Wraps one or more exceptions that are thrown when a callback written in
/// managed code is invoked by the native Cantera library.
/// </summary>
public class CallbackException : AggregateException
{
    // Cannot use an initializer because it only fires on the first thread.
    [ThreadStatic]
    static Stack<Exception>? t_exceptions;

    private CallbackException(string message, Exception[] innerExceptions)
        : base(message, innerExceptions) { }

    internal static void Register(Exception ex)
    {
        if (t_exceptions == null)
        {
            t_exceptions = new();
        }

        t_exceptions.Push(ex);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static void ThrowIfAny()
    {
        if (t_exceptions == null)
        {
            return;
        }

        DoThrow();

        // Separate this bit out so the most likely case can be inlined.
        static void DoThrow()
        {
            var exceptionsToThrow = t_exceptions!.ToArray();

            // Typically I would prefer mark a field as readonly and not reset it,
            // but in this case it’s ThreadStatic, and thus cannot be reliably
            // eager-initialized. Plus, we’re in a hot path, so setting and comparing to
            // null is the fastest thing to do.
            t_exceptions = null;

            throw new CallbackException("An exception occurred while executing a callback",
                exceptionsToThrow);
        }
    }
}
