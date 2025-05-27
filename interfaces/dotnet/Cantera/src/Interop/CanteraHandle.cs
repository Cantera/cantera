// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;
using System.Runtime.InteropServices.Marshalling;

namespace Cantera.Interop;

/// <summary>
/// The base class for a handle to a native Cantera object.
/// </summary>
/// <remarks>
/// <c>CanteraHandle</c> is a thin-wrapper around the integer handle to a native Cantera object.
/// It provides cleanup of the native object via the dispose pattern and custom marshalling
/// for use in P/Invoke signatures with <see cref="LibraryImportAttribute"/>.
/// This allows the P/Invoke signatures to have a touch of object-orientedness instead of using
/// plain <c>int</c>s to pass the handles.
/// </remarks>
abstract class CanteraHandle : IDisposable
{
    int _value;
    // I would prefer an enum or a bool, but int is needed for Interlocked.Exchange call
    // 0 = uninitialized or disposed, 1 = active
    // Only the marshaller should set the state to active
    int _state;

    abstract protected void Close();

    ~CanteraHandle()
    {
        if (_state == 1)
        {
            Close();
        }
    }

    public void Dispose()
    {
        if (Interlocked.Exchange(ref _state, 0) == 1)
        {
            Close();
            GC.SuppressFinalize(this);
        }
    }

    public sealed override string ToString() =>
        $"{GetType().Name} {{{_value}}}";

    [CustomMarshaller(typeof(CustomMarshallerAttribute.GenericPlaceholder), MarshalMode.ManagedToUnmanagedIn,
        typeof(Marshaller<>))]
    [CustomMarshaller(typeof(CustomMarshallerAttribute.GenericPlaceholder), MarshalMode.ManagedToUnmanagedOut,
        typeof(Marshaller<>))]
    public static class Marshaller<T> where T : CanteraHandle, new()
    {
        public static int ConvertToUnmanaged(T handle)
            => handle._value;

        public static T ConvertToManaged(int value)
        {
            InteropUtil.CheckReturn(value);
            T handle = new();
            // Set the field value separately because the generic constraint only allows parameterless constructors.
            handle._value = value;
            handle._state = 1;

            return handle;
        }
    }
}
