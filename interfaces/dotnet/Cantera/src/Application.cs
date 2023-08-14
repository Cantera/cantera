// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Contains information about a log message.
/// </summary>
public class LogMessageEventArgs
{
#pragma warning disable CS1591
    public LogLevel LogLevel { get; }
    public string Category { get; }
    public string Message { get; }
#pragma warning restore CS1591

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
/// <remarks>
/// All access to Cantera should funnel through this class.
/// This ensures that any necessary initialization can be run in
/// the static constructor.
/// </remarks>
public static class Application
{
    unsafe static readonly Lazy<string> s_version =
        new(() => InteropUtil.GetString(10, LibCantera.ct_getCanteraVersion));

    unsafe static readonly Lazy<string> s_gitCommit =
        new(() => InteropUtil.GetString(10, LibCantera.ct_getGitCommit));

    static readonly Lazy<DataDirectoryCollection> s_dataDirectories =
        new(() => new DataDirectoryCollection());

    /// <summary>
    /// Occurs when the Cantera native library attempts to log a message.
    /// </summary>
    public static event EventHandler<LogMessageEventArgs>? MessageLogged;

#pragma warning disable CS1591
    public static string Version =>
        s_version.Value;

    public static string GitCommit =>
        s_gitCommit.Value;

    public static DataDirectoryCollection DataDirectories =>
        s_dataDirectories.Value;
#pragma warning restore CS1591

    internal static void RaiseMessageLogged(LogLevel logLevel, string category, string message)
    {
        try
        {
            MessageLogged?.Invoke(null, new(logLevel, category, message));
        }
        catch (Exception ex)
        {
            CallbackException.Register(ex);
        }
    }

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
}
