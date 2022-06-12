using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

public class CanteraExceptionTest
{
    [Fact]
    public void CanteraException_Thrown()
    {
        Assert.Throws<CanteraException>(() =>
        {
            LibCantera.thermo_newFromFile(".yaml", "");

            // test that an error message is gathered from the native library
            CanteraException.ThrowLatest();
        });
    }
}
