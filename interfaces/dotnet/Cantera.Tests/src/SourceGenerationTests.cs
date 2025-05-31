using System.Reflection;
using System.Runtime.InteropServices.Marshalling;
using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

public static class SourceGenerationTests
{
    public static TheoryData<string> InteropMethodNames =>
        new(typeof(LibCantera).GetMethods(
                BindingFlags.Public | BindingFlags.Static)
            // get name rather than object for better compatibility
            // with test runners and IDE displays which may struggle with
            // non-primitive theory data
            .Select(m => m.Name));

    [Theory]
    [MemberData(nameof(InteropMethodNames))]
    public static void LibCantera_PInvokeReturnCodesChecked(string methodName)
    {
        var method = typeof(LibCantera).GetMethod(methodName);
        Assert.NotNull(method);

        var ReturnCodeCheckerType = typeof(LibCantera).GetNestedType(
            "ReturnCodeChecker", BindingFlags.NonPublic);
        Assert.NotNull(ReturnCodeCheckerType);

        var returnType = method.ReturnType;
        var returnMarshaller = method.ReturnParameter
            .GetCustomAttribute<MarshalUsingAttribute>()?.NativeType;

        // XOR - only one of these should be true
        Assert.True(returnType.IsSubclassOf(typeof(CanteraHandle))
            ^ returnMarshaller == ReturnCodeCheckerType
            // del functions should be unchecked
            ^ method.Name.EndsWith("_del", StringComparison.Ordinal));
    }
}