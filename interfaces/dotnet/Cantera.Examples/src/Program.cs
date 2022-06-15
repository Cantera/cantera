using Cantera.Examples;

var iExampleType = typeof(IExample);

var examples = iExampleType.Assembly
    .GetTypes()
    .Except(new[] {iExampleType})
    .Where(t => t.IsAssignableTo(iExampleType))
    .OrderBy(t => t.Name)
    .Select((t, i) => (index: i + 1, type: t))
    .ToList();

if (GetExample() is not Type exampleType)
{
    Console.WriteLine("\nUsage: Cantera.Examples [#/ExampleClassName]");
    Console.WriteLine("Available examples:\n");
    foreach (var example in examples)
    {
        Console.WriteLine($"    {example.index}: {example.type.Name}");
    }
    Console.WriteLine();

    return 1;
}

iExampleType
    .GetMethod(nameof(IExample.Run))!
    .Invoke(Activator.CreateInstance(exampleType), null);

return 0;

Type? GetExample()
{
    if (args.Length == 0)
        return null;

    if (examples!.SingleOrDefault(e => e.index.ToString() == args[0]).type is Type byIndex)
        return byIndex;

    if (examples!.SingleOrDefault(e => e.type.Name == args[0]).type is Type byName)
        return byName;

    return null;
}
