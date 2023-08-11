using System;

public class Example
{
    public void ExampleMethod()
    {
        var year = DateTime.Now.Year;

        // OK – Lambda expression marked static
        InteropUtil.ExampleUtilMethod(static () => Console.WriteLine("foo"));

        // OK – Referenced method is static
        InteropUtil.ExampleUtilMethod(static () => Console.WriteLine("foo"));

        // Compiler Error – Static lambda contains captures
        InteropUtil.ExampleUtilMethod(static () => Console.WriteLine(year));

        // Analyzer Error – Lambda not marked static (explicit param list)
        InteropUtil.ExampleUtilMethod(() => Console.WriteLine(year));

        // Analyzer Error – Lambda not marked static (1 naked param)
        InteropUtil.ExampleUtilMethodOneParam(obj => Console.WriteLine(obj));

        // Analyzer Error – Referenced method is not static
        InteropUtil.ExampleUtilMethod(InstanceMethod);
    }

    static void StaticMethod() { }

    void InstanceMethod() { }
}

static class InteropUtil
{
    public static void ExampleUtilMethod(Action action) { }

    public static void ExampleUtilMethodOneParam(Action<object> action) { }
}
