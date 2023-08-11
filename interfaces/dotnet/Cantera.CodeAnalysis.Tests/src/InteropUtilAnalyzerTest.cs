using Microsoft.CodeAnalysis.Testing;
using Xunit;
using Verifier =
    Microsoft.CodeAnalysis.CSharp.Testing.XUnit.AnalyzerVerifier<Cantera.CodeAnalysis.InteropUtilAnalyzer>;

namespace Cantera.CodeAnalysis.Tests;

public static class InteropUtilAnalyzerTest
{
    [Fact]
    public static async Task InteropUtilAnalyzer_CreatesDiagnostics()
    {
        var expectedDiagnostics = new[]
        {
            DiagnosticResult.CompilerError("CS8820")
                .WithLocation(16, 70),

            Verifier.Diagnostic("CT0001")
                .WithSpan(19, 39, 19, 41)
                .WithArguments("ExampleUtilMethod"),

            Verifier.Diagnostic("CT0001")
                .WithSpan(22, 47, 22, 50)
                .WithArguments("ExampleUtilMethodOneParam"),

            Verifier.Diagnostic("CT0002")
                .WithSpan(25, 39, 25, 53)
                .WithArguments("ExampleUtilMethod")
        };

        await Verifier.VerifyAnalyzerAsync(await Helpers.GetTestSource(), expectedDiagnostics);
    }
}
