using System.Buffers;
using System.Diagnostics.CodeAnalysis;

namespace Cantera.Interop;

static class InteropUtil
{
    /// <summary>
    /// Represents a function that gets the length of an array for a particular property
    /// of a Cantera object represented by <c>handle</c>.
    /// </summary>
    public delegate nuint GetSizeFunc<THandle>(THandle handle) where THandle : CanteraHandle;

    /// <summary>
    /// Represents a function that fills an array pointed to by <c>buffer</c> for particular property
    /// of a Cantera object represented by <c>handle</c>.
    /// </summary>
    public unsafe delegate int FillDoubleBufferFunc<THandle>(THandle handle, nuint size, double* buffer)
        where THandle : CanteraHandle;

    /// <summary>
    /// Represents a function that fills a byte buffer representing a native string.
    /// </summary>
    /// <remarks>
    /// The Cantera C API is not very consistent on whether the size should be specifed as an int
    /// or a nuint (size_t), so users of this delegate may need to wrap their LibCantera call in
    /// a lambda to perform the appropriate conversions, and/or pass in the Cantera handle
    /// as a capture. For example:
    /// <code>
    /// (size, buffer) => kin_getType(_handle, (nuint) size, buffer)
    /// </code>
    /// </remarks>
    public unsafe delegate int FillStringBufferFunc(int size, byte* buffer);

    public static double CheckReturn(double code)
    {
        // Cantera returns this value when the function resulted in error
        const double Error = -999.999;

        if (code == Error)
            CanteraException.ThrowLatest();

        return code;
    }

    public static nuint CheckReturn(nuint code)
    {
        var error = unchecked((nuint) (-1));

        if (code == error)
            CanteraException.ThrowLatest();

        return code;
    }

    public static int CheckReturn(int code)
    {
        // Cantera returns this value when the function resulted in an error internal to Cantera
        const int Error1 = -1;
        // Cantera returns this value when the function resulted in an external error
        // Some functions also return a negative value as the amount of space they need to
        // fill a buffer with a string. There is no way to account for the ambiguity that arises
        // when such a function returns -999!
        const int Error999 = -999;

        if (code == Error1 || code == Error999)
            CanteraException.ThrowLatest();

        // some functions return negative when they want more chars, others positive!
        return Math.Abs(code);
    }

    public static double[] GetDoubles<T>(T handle, GetSizeFunc<T> getSizefunc, FillDoubleBufferFunc<T> fillBufferfunc)
        where T : CanteraHandle
    {
        var size = getSizefunc(handle);
        var array = new double[size];

        GetDoubles(handle, array, fillBufferfunc);

        return array;
    }

    public static double[] GetDoubles<T>(T handle, int count, FillDoubleBufferFunc<T> fillBufferfunc)
        where T : CanteraHandle
    {
        var array = new double[count];

        GetDoubles(handle, array, fillBufferfunc);

        return array;
    }

    public static unsafe void GetDoubles<T>(T handle, Span<Double> span, FillDoubleBufferFunc<T> fillBufferfunc)
        where T : CanteraHandle
    {
        fixed(double* buffer = span)
        {
            CheckReturn(fillBufferfunc(handle, (nuint) span.Length, buffer));
        }
    }

    [SuppressMessage("Reliability", "CA2014:NoStackallocInLoops", Justification = "Loop is executed at most twice.")]
    public static string GetString(int initialSize, FillStringBufferFunc func)
    {
        // take up to two tries
        // 1) use the initial size
        //    if the initial size was large enough, return the string
        //    if the initial size was not large enough ...
        // 2) try again with the needed size
        //    if the needed size was large enough, return the string
        //    otherwise, catastrophe, throw!
        for(var i = 0; i < 2; i++)
        {
            int neededSize;
            
            if (initialSize <= 120)
            {
                Span<byte> span = stackalloc byte[initialSize];
                if (TryGetString(span, func, out var value, out neededSize))
                    return value;
            }
            else 
            {
                using (MemoryPool<byte>.Shared.Rent(initialSize, out var span))
                {
                    if (TryGetString(span, func, out var value, out neededSize))
                        return value;
                }
            }

            initialSize = neededSize;
        }

        throw new InvalidOperationException("Could not retrieve a string value from Cantera!");

        static unsafe bool TryGetString(Span<byte> span, FillStringBufferFunc func,
            [NotNullWhen(true)] out string? value, out int neededSize)
        {
            var initialSize = span.Length;

            fixed(byte* buffer = span)
            {
                neededSize = CheckReturn(func(initialSize, buffer));

                if(initialSize >= neededSize)
                {
                    value = new String((sbyte*) buffer);
                    return true;
                }
            }

            value = null;
            return false;
        }
    }
}
