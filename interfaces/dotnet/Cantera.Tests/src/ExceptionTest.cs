// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Diagnostics;
using System.Reflection;
using Cantera.Interop;
using Xunit;
using Xunit.Sdk;

namespace Cantera.Tests;

[Collection(LibCanteraFixture.Collection)]
public class ExceptionTest
{
    class FooException : Exception { }

    [Fact]
    public void CanteraException_Thrown()
    {
        var handle = LibCantera.thermo_newFromFile(".yaml", "");

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

        var methodCalls = new StackTrace(thrown).GetFrames()
            .Select(f=> f.GetMethod());

        var method = typeof(CallbackException).GetMethod(nameof(CallbackException.ThrowIfAny),
            BindingFlags.NonPublic | BindingFlags.Static);

        Assert.NotNull(method);

        try
        {
            // Test that the method was inlined,
            // meaning it does not appear in the stack trace.
            Assert.DoesNotContain(method, methodCalls);
        }
        catch (DoesNotContainException)
        {
            // The method may have not been inlined because we’re running in Debug mode,
            // or tiered compilation hasn’t inlined it yet.
            // So, check that the appropriate flag is set so that it _could_ be inlined.
            Assert.False((method.MethodImplementationFlags & MethodImplAttributes.AggressiveInlining) == 0);
        }
    }

    [Fact]
    public void CallbackException_NotThrown()
    {
        CallbackException.ThrowIfAny();
    }
}
