// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;

namespace Cantera;

public class LogMessageEventArgs
{
    public LogLevel LogLevel { get; }
    public string Category { get; }
    public string Message { get; }

    internal LogMessageEventArgs(LogLevel logLevel, string category, string message)
    {
        LogLevel = logLevel;
        Category = category;
        Message = message;
    }
}

/// <summary>
/// The primary API for accessing the Cantera library.
/// </summary>
/// </remarks>
/// All access to Cantera should funnel through this class.
/// This ensures that any necessary initialization can be run in
/// the static constructor.
/// <remarks>
public static class Application
{
    static Application()
    {
        _invokeMessageLoggedDelegate = (level, category, message) =>
            MessageLogged
                ?.Invoke(null, new LogMessageEventArgs(level, category, message));

        InteropUtil.CheckReturn(
            LibCantera.ct_setLogWriter(_invokeMessageLoggedDelegate));
    }

    /// <summary>
    /// Represents the delegate that is marshalled to LibCantera as a function pointer.
    /// </summary>
    /// <remarks>
    /// ct_setLogWriter() needs a delegate which is marshalled as a function pointer to
    /// the C++ Cantera library. We could create one implicitly when calling
    /// <c>LibCantera.ct_setLogWriter(LogMessage)</c>, but the garbage collector would
    /// not know the native method is using it and could reclaim it. By explicitly
    /// storing it as
    /// a class member, we ensure it is not collected until the class is.
    /// </remarks>
    static LibCantera.Writer? _invokeMessageLoggedDelegate;

    unsafe static readonly Lazy<string> _version =
        new(() => InteropUtil.GetString(10, LibCantera.ct_getCanteraVersion));

    unsafe static readonly Lazy<string> _gitCommit =
        new(() => InteropUtil.GetString(10, LibCantera.ct_getGitCommit));

    unsafe static readonly Lazy<DataDirectoryCollection> _dataDirectories =
        new(() => new DataDirectoryCollection());

    public static event EventHandler<LogMessageEventArgs>? MessageLogged;

    public static string Version =>
        _version.Value;

    public static string GitCommit =>
        _gitCommit.Value;

    unsafe public static DataDirectoryCollection DataDirectories =>
        _dataDirectories.Value;

    /// <summary>
    /// Convenience method to add logging to the console.
    /// </summary>
    public static void AddConsoleLogging() =>
        MessageLogged += LogToConsole;

    /// <summary>
    /// Convenience method to remove logging to the console.
    /// </summary>
    public static void RemoveConsoleLogging() =>
        MessageLogged -= LogToConsole;

    static void LogToConsole(object? sender, LogMessageEventArgs e)
    {
        var logLevel = e.LogLevel.ToString().ToUpperInvariant();

        var nowString = DateTimeOffset.Now.ToString("yyyy-MM-ddThh:mm:ss.fffzzz");

        var message = $"{logLevel} ({e.Category}) {nowString}: {e.Message}";

        if (e.LogLevel == LogLevel.Error)
        {
            Console.Error.WriteLine(message);
        }
        else
        {
            Console.WriteLine(message);
        }
    }

    public static ThermoPhase CreateThermoPhase(string filename,
                                                string? phasename = null) =>
        new ThermoPhase(filename, phasename);
}
