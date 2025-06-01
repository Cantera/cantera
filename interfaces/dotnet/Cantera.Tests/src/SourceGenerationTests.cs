using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.Marshalling;
using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

public static class SourceGenerationTests
{
    private static IReadOnlyList<MethodInfo> s_interopMethods =
        typeof(LibCantera).GetMethods(
            BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static)
        .Where(m => m.IsDefined(typeof(LibraryImportAttribute)))
        .ToList();

    /// <remarks>
    /// Get name rather than object for better compatibility with test runners
    /// and IDE displays which may struggle with non-primitive theory data.
    /// </remarks>
    public static TheoryData<string> InteropMethodNames =>
        new(s_interopMethods.Select(m => m.Name));

    [Theory]
    [MemberData(nameof(InteropMethodNames))]
    public static void LibCantera_PInvokeReturnCodesChecked(string methodName)
    {
        var method = s_interopMethods.Single(m => m.Name == methodName);

        var returnCodeCheckerType = typeof(LibCantera).GetNestedType(
            "ReturnCodeChecker", BindingFlags.NonPublic);
        Assert.NotNull(returnCodeCheckerType);

        var returnType = method.ReturnType;
        var returnMarshaller = method.ReturnParameter
            .GetCustomAttribute<MarshalUsingAttribute>()?.NativeType;

        // XOR - only one of these should be true
        Assert.True(returnType.IsSubclassOf(typeof(CanteraHandle))
            ^ returnMarshaller == returnCodeCheckerType
            // del functions should be unchecked
            ^ method.Name.EndsWith("_del", StringComparison.Ordinal));
    }
}
