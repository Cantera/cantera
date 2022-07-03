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
