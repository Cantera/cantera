// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Reflection;
using System.Runtime.CompilerServices;
using Xunit;

namespace Cantera.Tests;

public class EnumsTests
{
    [Fact]
    public void ThermoPair_ToStringsCorrectly()
    {
        var fields = typeof(ThermoPair)
            .GetFields(BindingFlags.Public | BindingFlags.Static)
            .Select(f =>
            (
                name: f.Name,
                interopString: ((ThermoPair) f.GetValue(null)!).ToInteropString()
            ));

        var withShortNames = fields
            .Where(f => f.name.Length == 2);

        var withLongNames = fields
            .Except(withShortNames);

        // Ensure we have both long and short versions
        Assert.Equal(withShortNames.Count(), withLongNames.Count());

        foreach (var (name, interopString) in withShortNames)
        {
            Assert.Equal(name, interopString);

            Assert.Single(withLongNames, f => f.interopString == name);
        }
    }

    [Fact]
    public void LogLevel_MapsCorrectly() =>
        TestInteropEnumInvariants<LogLevel>(true, 0, 2);

    static void TestInteropEnumInvariants<TEnum>(bool contiguous = false,
                                                 int? min = null, int? max = null)
        where TEnum : struct, Enum
    {
        Assert.Equal(typeof(int), Enum.GetUnderlyingType(typeof(TEnum)));

        var values = Enum
            .GetValues<TEnum>()
            .Select(v => Unsafe.As<TEnum, int>(ref v))
            .OrderBy(v => v)
            .ToList();

        if (contiguous)
        {
            for (var i = 1; i < values.Count; i++)
            {
                Assert.Equal(1, values[i] - values[i - 1]);
            }
        }

        if (min is not null)
        {
            Assert.Equal(min, values.Min());
        }

        if (max is not null)
        {
            Assert.Equal(max, values.Max());
        }
    }
}
