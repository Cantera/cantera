using System.Diagnostics;
using System.Reflection;
using System.Runtime.CompilerServices;

namespace Cantera.CodeAnalysis.Tests;

static class Helpers
{
    [MethodImpl(MethodImplOptions.NoInlining)]
    public static Task<string> GetTestSource()
    {
        int i = 1;
        Type testClass;
        do
        {
            testClass = new StackTrace().GetFrame(i++)!.GetMethod()!.DeclaringType!;
        } while (testClass.Assembly != typeof(Helpers).Assembly
            || testClass.GetCustomAttribute<CompilerGeneratedAttribute>() is not null);

        return File.ReadAllTextAsync("testSources/" + testClass.Name + ".testSource.cs");
    }
}
