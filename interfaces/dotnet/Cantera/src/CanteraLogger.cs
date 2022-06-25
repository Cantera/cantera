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

public static class CanteraLogger
{
    /// <summary>
    /// Represents the delegate that is marshalled to LibCantera as a function pointer.
    /// </summary>
    /// <remarks>
    /// ct_setLogWriter() needs a delagate which is marshalled as a function pointer to
    /// the C++ Cantera library. We could create one implicitly when calling
    /// <c>LibCantera.ct_setLogWriter(LogMessage)</c>, but the garbage collector would not know
    /// the native method is using it and could reclaim it. By explicitly storing it as
    /// a class member, we ensure it is not collected until the class is.
    /// </remarks>
    static LibCantera.Writer? _invokeMessageLoggedDelegate;

    // 0 => not hooked up to Cantera via CLIB
    // 1 => hooked up to Cantera
    static int _state;

    static event EventHandler<LogMessageEventArgs>? _messageLogged;

    public static event EventHandler<LogMessageEventArgs> MessageLogged
    {
        add
        {
            if (Interlocked.CompareExchange(ref _state, 1, 0) == 0) try
            {
                _invokeMessageLoggedDelegate = (level, category, message) =>
                    _messageLogged?.Invoke(null, new LogMessageEventArgs(level, category, message));

                InteropUtil.CheckReturn(LibCantera.ct_setLogWriter(_invokeMessageLoggedDelegate));
            }
            catch
            {
                _state = 0;
                throw;
            }

            _messageLogged += value;
        }

        remove
        {
            _messageLogged -= value;
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
            Console.Error.WriteLine(message);
        else
            Console.WriteLine(message);
    }
}