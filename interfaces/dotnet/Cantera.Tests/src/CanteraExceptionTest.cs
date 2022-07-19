// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
