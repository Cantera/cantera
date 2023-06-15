// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Buffers;

namespace Cantera;

static class MemoryExtensions
{
    /// <summary>
    /// Provides a convenience method to rent a block of memory and get a reference
    /// to the <see cref="Span{T}" /> it encapsulates simultaneously.
    /// </summary>
    /// <remarks>
    /// The <see cref="IMemoryOwner{T}" /> returned by this method can be used with
    /// an anonymous using block to avoid declaring an unnecessary local variable,
    /// for example: <c>using (myMemoryPool.Rent(200, out var mySpan)){ ... }</c>.
    /// The returned <see cref="Span{T}" /> must only be used within that using block.
    /// </remarks>
    public static IMemoryOwner<T> Rent<T>(this MemoryPool<T> pool, int minimumSize,
                                          out Span<T> span)
    {
        var owner = pool.Rent(minimumSize);
        span = owner.Memory.Span;

        return owner;
    }
}
