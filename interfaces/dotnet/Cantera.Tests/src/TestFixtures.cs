// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Xunit;

namespace Cantera.Tests;

public class ApplicationFixture
{
    public ApplicationFixture()
    {
        Application.DataDirectories.AddAssemblyDirectory();
    }
}

[CollectionDefinition("Application")]
public class DatabaseCollection : ICollectionFixture<ApplicationFixture> { }
