// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Runtime.CompilerServices;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Contains information about a log message.
/// </summary>
public class LogMessageEventArgs : EventArgs
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
/// Provides APIs for accessing the Cantera library.
/// </summary>
public static class Application
{
    /// <remarks>
    /// This method is automatically called by the runtime to initialize
    /// the native Cantera library. You should not need to call it elsewhere.
    /// </remarks>
    [ModuleInitializer]
    [SuppressMessage("Usage", "CA2255: No ModuleInitializerAttribute in library code.",
        Justification = "Initialization code is essential.")]
    internal static void Initialize()
    {
        LibCantera.ct_setLogCallback(s_invokeMessageLoggedDelegate);
    }

    /// <summary>
    /// Represents the delegate that is marshalled to LibCantera as a function pointer.
    /// </summary>
    /// <remarks>
    /// ct_setLogWriter() needs a delegate which is marshalled as a function pointer to
    /// the C++ Cantera library. We could create one implicitly when calling
    /// <c>LibCantera.ct_setLogWriter(LogMessageEventArgs)</c>, but the garbage collector would
    /// not know the native method is using it and could reclaim it. By explicitly
    /// storing it as
    /// a class member, we ensure it is not collected until the class is.
    /// </remarks>
    static readonly LibCantera.LogCallback s_invokeMessageLoggedDelegate =
        (level, categoryStr, messageStr) =>
    {
        try
        {
            Application.
            MessageLogged
                ?.Invoke(null, new LogMessageEventArgs(level, categoryStr, messageStr));
        }
        catch (Exception ex)
        {
            CallbackException.Register(ex);
        }
    };

    static readonly Lazy<string> s_version =
        new(LibCantera.ct_version);

    static readonly Lazy<string> s_gitCommit =
        new(LibCantera.ct_gitCommit);

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

        var nowString = DateTimeOffset.Now.ToString(
            "yyyy-MM-ddThh:mm:ss.fffzzz", CultureInfo.InvariantCulture);

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

    /// <summary>
    /// Returns a new <see cref="ThermoPhase" /> object by loading and parsing the
    /// given configuration file. Optionally chooses the phase to load by
    /// looking up the given name.
    /// </summary>
    public static ThermoPhase CreateThermoPhase(string filename,
                                                string? phaseName = null) =>
        new ThermoPhase(filename, phaseName);
}
