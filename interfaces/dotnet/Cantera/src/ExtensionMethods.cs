using System.Buffers;

namespace Cantera;

static class MemoryExtensions
{
    public static IMemoryOwner<T> Rent<T>(this MemoryPool<T> pool, int minmumSize, out Span<T> span)
    {
        var owner = pool.Rent(minmumSize);
        span = owner.Memory.Span;

        return owner;
    }
}
