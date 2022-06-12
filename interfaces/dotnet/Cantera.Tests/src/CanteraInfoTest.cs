using Xunit;

namespace Cantera.Tests;

public class CanteraInfoTest
{
    [Fact]
    public void CanteraInfo_VersionRetrieved()
    {
        // version string should start with a number
        Assert.Matches("^[1-9]", CanteraInfo.Version);
    }

    [Fact]
    public void CanteraInfo_GitCommitRetrieved()
    {
        Assert.NotEmpty(CanteraInfo.GitCommit);
    }

    [Fact]
    public void CanteraInfo_DataDirectoryAdded()
    {
        var dirs = CanteraInfo.DataDirectories;
        var originalCount = dirs.Count;
        var longestDir = dirs.MaxBy(d => d.FullName.Length);

        Assert.NotNull(longestDir);

        dirs.Add(Path.Join(longestDir!.FullName, "garbazh"));

        Assert.Equal(originalCount + 1, dirs.Count);
    }
}