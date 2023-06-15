// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Runtime.InteropServices;

namespace Cantera.Interop;

/// <summary>
/// The base class for a handle to a Cantera object.
/// </summary>
/// <remarks>
/// We use the SafeHandle class, which has low-level support in the runtime to ensure
/// proper reference counting and cleanup and is thread-safe. This allows us to use
/// the dispose pattern easily and safely and without having to write our own finalizer.
/// Cantera uses signed 32-bit ints for handles, yet SafeHandle uses
/// "native int" IntPtr. The Value and IsValid properties are designed to account for
/// this and only consider the lower 32-bit of the IntPtr on 64-bit systems.
/// </remarks>
abstract class CanteraHandle : SafeHandle
{
    static readonly IntPtr s_invalid = IntPtr.Size == 4
        ? new IntPtr(-1)                    // 32-bit IntPtr: 0xFFFFFFFF
        : new IntPtr(unchecked((uint) -1)); // 64-bit IntPtr: 0x00000000FFFFFFFF

    public CanteraHandle() : base(s_invalid, true) { }

    protected int Value =>
        IntPtr.Size == 4
            ? (int) handle
            : (int)(long) handle; // removes any leading bits

    public override bool IsInvalid =>
        Value < 0;

    public void EnsureValid()
    {
        if (IsInvalid)
        {
            CanteraException.ThrowLatest();
        }
    }
}
