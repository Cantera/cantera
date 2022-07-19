// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Buffers;

namespace Cantera;

static class MemoryExtensions
{
    public static IMemoryOwner<T> Rent<T>(this MemoryPool<T> pool, int minimumSize,
                                          out Span<T> span)
    {
        var owner = pool.Rent(minimumSize);
        span = owner.Memory.Span;

        return owner;
    }
}
