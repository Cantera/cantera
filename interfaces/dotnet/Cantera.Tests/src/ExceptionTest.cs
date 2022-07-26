// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

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
    public void CallbackException_NotThrown()
    {
        CallbackException.ThrowIfAny();
    }
}
