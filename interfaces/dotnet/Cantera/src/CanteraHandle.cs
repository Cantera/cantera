using System.Runtime.InteropServices;

namespace Cantera;

/// <summary>
/// The base class for a handle to a cantera object.
/// </summary>
/// </remarks>
/// Cantera uses signed 32-bit ints for handles, yet SafeHandle uses "native int" IntPtr.
/// The constuctor and the IsInvalid property are designed to account for this.
/// <remarks>
internal abstract class CanteraHandle : SafeHandle
{
    static IntPtr INVALID = IntPtr.Size == 4
        ? new IntPtr(-1)                           // 32-bit IntPtr: 0xFFFFFFFF
        : new IntPtr((long) unchecked((uint) -1)); // 64-bit IntPtr: 0x00000000FFFFFFFF

    public CanteraHandle() : base(INVALID, true) { }

    public override bool IsInvalid =>
         (int)(long)handle < 0; // removes any leading bits

    public void EnsureValidHandleReceived()
    {
        if (IsInvalid)
            CanteraException.ThrowLatest();
    }

}