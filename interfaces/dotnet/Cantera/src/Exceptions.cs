// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.CompilerServices;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Represents a error that occurred in the native Cantera library.
/// </summary>
public class CanteraException : Exception
{
    private CanteraException(string message) : base(message) { }

    internal static void ThrowLatest() =>
        throw new CanteraException(LibCantera.ct_getCanteraError());
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

    private CallbackException(string message, IEnumerable<Exception> innerExceptions)
        : base(message, innerExceptions) { }

    internal static void Register(Exception ex) =>
        (t_exceptions ??= new()).Push(ex);

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
            CallbackException exceptionToThrow =
                new("An exception occurred while executing a callback", t_exceptions!);

            t_exceptions = null;

            throw exceptionToThrow;
        }
    }
}
