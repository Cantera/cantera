using System.Reflection;
using System.Text.RegularExpressions;
using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

[Collection("Application")]
public class ApplicationTest
{
    [Fact]
    public void CanteraInfo_VersionRetrieved()
    {
        // version string should start with a number
        Assert.Matches("^[1-9]", Application.Version);
    }

    [Fact]
    public void CanteraInfo_GitCommitRetrieved()
    {
        Assert.NotEmpty(Application.GitCommit);
    }

    [Fact]
    public void CanteraInfo_DataDirectoryAdded()
    {
        var dirs = Application.DataDirectories;
        var originalCount = dirs.Count;
        var longestDir = dirs.MaxBy(d => d.FullName.Length);

        Assert.NotNull(longestDir);

        dirs.Add(Path.Join(longestDir!.FullName, "garbazh"));

        Assert.Equal(originalCount + 1, dirs.Count);
    }

    [Fact]
    public void LogWriter_MessageLogged()
    {
        LogLevel? logLevel = null;
        string? report = null;

        Application.MessageLogged += (sender, e) =>
        {
            logLevel = e.LogLevel;
            report = e.Message;
        };

        ProduceLogOutput();

        Assert.NotNull(logLevel);
        Assert.NotNull(report);

        Assert.Equal(LogLevel.Info, logLevel);
        Assert.NotEmpty(report);
    }

    [Fact]
    public void LogWriter_ConsoleLogged()
    {
        var stdOut = Console.Out;
        var consoleOut = new StringWriter();

        try
        {
            Console.SetOut(consoleOut);
            Application.AddConsoleLogging();
            ProduceLogOutput();
            
            var output = consoleOut.ToString();

            const string prefix = "INFO (Info) ";
            var iso8601FormatString = "yyyy-MM-ddThh:mm:ss.fffzzz";
            var lengthOfIso8601FormattedString = 29;

            Assert.Matches('^' + Regex.Escape(prefix), output);

            var nowString = output.Substring(prefix.Length, lengthOfIso8601FormattedString);

            Assert.True(DateTimeOffset.TryParseExact(nowString, iso8601FormatString, null, default, out _));
        }
        finally
        {
            Application.RemoveConsoleLogging();
            Console.SetOut(stdOut);
        }
    }

    static void ProduceLogOutput()
    {
        using var thermo = Application.CreateThermoPhase("gri30.yaml");

        var handle = (ThermoPhaseHandle) typeof(ThermoPhase)
            .GetField("_handle", BindingFlags.NonPublic | BindingFlags.Instance)!
            .GetValue(thermo)!;

        InteropUtil.CheckReturn(LibCantera.thermo_print(handle, InteropConsts.True, 0));
    }
}
