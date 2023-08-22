// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Xunit;

namespace Cantera.Tests;

public class LibCanteraFixture
{
    internal const string Collection = "LibCantera";

    public LibCanteraFixture()
    {
        Application.DataDirectories.AddAssemblyDirectory();
    }
}

[CollectionDefinition(LibCanteraFixture.Collection)]
public class LibCanteraCollection : ICollectionFixture<LibCanteraFixture> { }
