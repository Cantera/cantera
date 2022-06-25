using System.Reflection;
using System.Text.RegularExpressions;
using Cantera.Interop;
using Xunit;

namespace Cantera.Tests;

public class LoggerTest
{
    [Fact]
    public void LogWriter_MessageLogged()
    {
        LogLevel? logLevel = null;
        string? report = null;

        CanteraLogger.MessageLogged += (sender, e) =>
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
            CanteraLogger.AddConsoleLogging();
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
            CanteraLogger.RemoveConsoleLogging();
            Console.SetOut(stdOut);
        }
    }

    static void ProduceLogOutput()
    {
        CanteraInfo.DataDirectories.AddAssemblyDirectory();

        using var thermo = new ThermoPhase("gri30.yaml");

        var handle = (ThermoPhaseHandle) typeof(ThermoPhase)
            .GetField("_handle", BindingFlags.NonPublic | BindingFlags.Instance)!
            .GetValue(thermo)!;

        InteropUtil.CheckReturn(LibCantera.thermo_print(handle, InteropConsts.True, 0));
    }
}