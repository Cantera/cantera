using System.Runtime.CompilerServices;
using Xunit;

namespace Cantera.Tests;

public class EnumsTests
{
    [Fact]
    public void LogLevel_MapsCorrectly() =>
        TestInteropEnumInvariants<LogLevel>(true, 0, 2);

    static void TestInteropEnumInvariants<TEnum>(bool contiguous = false, int? min = null, int? max = null)
        where TEnum : struct, Enum
    {
        Assert.Equal(typeof(int), Enum.GetUnderlyingType(typeof(TEnum)));

        var values = Enum
            .GetValues<TEnum>()
            .Select(v => Unsafe.As<TEnum, int>(ref v))
            .OrderBy(v => v)
            .ToList();

        if (contiguous) for (var i = 1; i < values.Count; i++)
        {
            Assert.Equal(1, values[i] - values[i - 1]);
        }

        if (min is not null)
            Assert.Equal(min, values.Min());

        if (max is not null)
            Assert.Equal(max, values.Max());
    }
}