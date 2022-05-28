using System.Collections;
using Cantera.Interop;

namespace Cantera;

public static class CanteraInfo
{
    public class DataDirectoryCollection : IReadOnlyList<DirectoryInfo>
    {
        unsafe static IEnumerable<DirectoryInfo> GetDirs()
        {
            const char sep = ';';

            return InteropUtil
                .GetString(500, (length, buffer) =>
                    LibCantera.ct_getDataDirectories(length, buffer, sep.ToString()))
                .Split(sep)
                .Select(d => new DirectoryInfo(d));
        }

        readonly List<DirectoryInfo> _dirs;

        public DirectoryInfo this[int index] => _dirs[index];

        public int Count => _dirs.Count;

        internal DataDirectoryCollection()
        {
            _dirs = GetDirs().ToList();
        }

        public void Add(string dir) =>
            Add(new DirectoryInfo(dir));

        public void Add(DirectoryInfo dir)
        {
            InteropUtil.CheckReturn(LibCantera.ct_addCanteraDirectory((nuint) dir.FullName.Length, dir.FullName));

            _dirs.Clear();
            _dirs.AddRange(GetDirs());
        }

        public IEnumerator<DirectoryInfo> GetEnumerator() =>
            _dirs.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() =>
            _dirs.GetEnumerator();
    }

    unsafe static readonly Lazy<string> _version =
        new Lazy<string>(() => InteropUtil.GetString(10, LibCantera.ct_getCanteraVersion));

    unsafe static readonly Lazy<string> _gitCommit =
        new Lazy<string>(() => InteropUtil.GetString(10, LibCantera.ct_getGitCommit));

    unsafe static readonly Lazy<DataDirectoryCollection> _dataDirectories =
        new Lazy<DataDirectoryCollection>(() => new DataDirectoryCollection());

    public static string Version =>
        _version.Value;

    public static string GitCommit =>
        _gitCommit.Value;

    unsafe public static DataDirectoryCollection DataDirectories =>
        _dataDirectories.Value;
}