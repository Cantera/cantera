// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Text.RegularExpressions;
using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

public class ExceptionTest
{
    class FooException : Exception { }

    [Fact]
    public void CanteraException_Thrown()
    {
        var handle = LibCantera.soln_newSolution(".yaml", "", "");

        Assert.True(handle.IsInvalid);

        // test that an error message is gathered from the native library
        Assert.Throws<CanteraException>(() => CanteraException.ThrowLatest());
    }

    [Fact]
    public void CallbackException_Thrown()
    {
        CallbackException.Register(new FooException());

        var thrown =
            Assert.Throws<CallbackException>(() => CallbackException.ThrowIfAny());

        Assert.NotNull(thrown.InnerException);
        Assert.IsType<FooException>(thrown.InnerException);
    }

    [Fact]
    public void CallbackException_ThrowIfAnyInlined()
    {
        CallbackException.Register(new FooException());

        var thrown =
            Assert.ThrowsAny<Exception>(() => CallbackException.ThrowIfAny());

        var methodCalls = thrown.StackTrace!
            .Split('\n')
            .Select(f => Regex.Match(f, "^   at (.*)\\(").Groups[1].Value);

        var methodName = typeof(CallbackException).FullName + '.'
            + nameof(CallbackException.ThrowIfAny);

        // Test that the method was inlined,
        // meaning it does not appear in the stack trace.
        Assert.DoesNotContain(methodName, methodCalls);
    }

    [Fact]
    public void CallbackException_NotThrown()
    {
        CallbackException.ThrowIfAny();
    }
}
