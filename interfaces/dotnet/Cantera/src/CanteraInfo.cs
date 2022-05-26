using System.Runtime.InteropServices;

namespace Cantera;

public static class CanteraInfo
{
    static class LibCantera
    {
        [DllImport(Interop.LIBCANTERA)]
        unsafe public static extern int ct_getCanteraVersion(int buflen, byte* buf);

        [DllImport(Interop.LIBCANTERA)]
        unsafe public static extern int ct_getGitCommit(int buflen, byte* buf);

    }

    unsafe static readonly Lazy<string> _version =
        new Lazy<string>(() => Interop.GetString(10, LibCantera.ct_getCanteraVersion));

    unsafe static readonly Lazy<string> _gitCommit =
        new Lazy<string>(() => Interop.GetString(10, LibCantera.ct_getGitCommit));

    public static string Version =>
        _version.Value;

    public static string GitCommit =>
        _gitCommit.Value;
}