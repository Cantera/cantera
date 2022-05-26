using System.Buffers;
using System.Diagnostics.CodeAnalysis;

namespace Cantera;

internal static class Interop
{
    public unsafe delegate int GetStringFunc(int initialSize, byte* buffer);

    static ArrayPool<byte> Pool = ArrayPool<byte>.Shared;

    public const string LIBCANTERA = "cantera.2.6.0";

    public static double CheckReturn(double code)
    {
        // Cantera returns this value when the function resulted in error
        const double Error = -999.999;

        if (code == Error)
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

    [SuppressMessage("Reliability", "CA2014:NoStackallocInLoops", Justification = "Loop is executed at most twice.")]
    public static string GetString(int initialSize, GetStringFunc func)
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
                var array = Pool.Rent(initialSize);
                try
                {
                    if (TryGetString(array, func, out var value, out neededSize))
                        return value;
                }
                finally
                {
                    Pool.Return(array);
                }
            }

            initialSize = neededSize;
        }

        throw new InvalidOperationException("Could not retrieve a string value from Cantera!");

        unsafe static bool TryGetString(Span<byte> span, GetStringFunc func, [NotNullWhen(true)] out string? value, out int neededSize)
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