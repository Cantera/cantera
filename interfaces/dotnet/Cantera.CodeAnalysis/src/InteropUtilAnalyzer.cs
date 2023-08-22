using System.Collections.Immutable;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.Diagnostics;
using Microsoft.CodeAnalysis.Operations;

namespace Cantera.CodeAnalysis;

/// <summary>
/// An analyzer that finds improper usages of the InteropUtil class.
/// </summary>
[DiagnosticAnalyzer(LanguageNames.CSharp)]
public class InteropUtilAnalyzer : DiagnosticAnalyzer
{
    static readonly DiagnosticDescriptor s_NonStaticLambdaExpression = new
    (
        id: "CT0001",
        category: "Interop",
        title: "All lambda expressions must be static in InteropUtil method invocations",
        description: "To avoid delegate and closure allocations in interop code, "
                       + "lambda expressions in InteropUtil method invocations must be marked static.",
        messageFormat: "Lambda expression in InteropUtil.{0} call must be marked static",
        defaultSeverity: DiagnosticSeverity.Error,
        isEnabledByDefault: true
    );

    static readonly DiagnosticDescriptor s_NonStaticMethodReference = new
    (
        id: "CT0002",
        category: "Interop",
        title: "All referenced method must be static in InteropUtil method invocations",
        description: "To avoid delegate and closure allocations in interop code, "
                     + "methods referenced in InteropUtil method invocations must be static methods.",
        messageFormat: "Reference to non-static method in InteropUtil.{0} call is not allowed",
        defaultSeverity: DiagnosticSeverity.Error,
        isEnabledByDefault: true
    );

    // It’s mandatory to report what diagnostics may be reported
    public override ImmutableArray<DiagnosticDescriptor> SupportedDiagnostics { get; } =
        ImmutableArray.Create(s_NonStaticLambdaExpression, s_NonStaticMethodReference);

    public override void Initialize(AnalysisContext context)
    {
        context.ConfigureGeneratedCodeAnalysis(
            GeneratedCodeAnalysisFlags.Analyze | GeneratedCodeAnalysisFlags.ReportDiagnostics);
        context.EnableConcurrentExecution();
        context.RegisterOperationAction(AnalyzeOperation, OperationKind.Invocation);
    }

    /// <summary>
    /// Executed on the completion of the semantic analysis associated with the Invocation operation.
    /// </summary>
    static void AnalyzeOperation(OperationAnalysisContext context)
    {
        if (context.Operation is not IInvocationOperation invocationOperation)
        {
            return;
        }

        var methodSymbol = invocationOperation.TargetMethod;

        // Check for invocations that include lambda/delegate arguments
        if (methodSymbol.MethodKind != MethodKind.Ordinary
            || methodSymbol.ReceiverType?.Name != "InteropUtil"
            || !methodSymbol.Parameters.Any(p => p.Type.TypeKind == TypeKind.Delegate))
        {
            return;
        }

        var delegateTarget = invocationOperation.Arguments
            .Select(a => a.Value)
            .OfType<IDelegateCreationOperation>()
            .SingleOrDefault()?.Target;

        switch (delegateTarget)
        {
            case IAnonymousFunctionOperation { Symbol.IsStatic: false }:
                var diagnostic1 = Diagnostic.Create(s_NonStaticLambdaExpression,
                    delegateTarget.Syntax.ChildNodes().First().GetLocation(), methodSymbol.Name);

                context.ReportDiagnostic(diagnostic1);
                break;

            case IMethodReferenceOperation { Method.IsStatic: false }:
                var diagnostic2 = Diagnostic.Create(s_NonStaticMethodReference,
                    delegateTarget.Syntax.GetLocation(), methodSymbol.Name);

                context.ReportDiagnostic(diagnostic2);
                break;
        }
    }
}
