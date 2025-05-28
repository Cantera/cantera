// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Collections;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// A collection of directories used by Cantera to locate data files.
/// </summary>
public class DataDirectoryCollection : IReadOnlyList<DirectoryInfo>
{
    static IEnumerable<DirectoryInfo> GetDirs()
    {
        const char sep = ';';

        return InteropUtil
            .GetString(500, (size, buffer) =>
                LibCantera.ct3_getDataDirectories(sep.ToString(), size, buffer))
            .Split(sep)
            .Select(d => new DirectoryInfo(d));
    }

    readonly List<DirectoryInfo> _dirs;

    /// <inheritdoc />
    public DirectoryInfo this[int index] => _dirs[index];

    /// <inheritdoc />
    public int Count => _dirs.Count;

    internal DataDirectoryCollection()
    {
        _dirs = GetDirs().ToList();
    }

    /// <summary>
    /// Adds the given directory to beginning of the collection. If the directory
    /// is already present in the collection, removes the existing entry.
    /// </summary>
    public void Add(string dir) =>
        Add(new DirectoryInfo(dir));

    /// <inheritdoc cref="Add(string)"/>
    public void Add(DirectoryInfo dir)
    {
        InteropUtil.CheckReturn(LibCantera.ct3_addCanteraDirectory(dir.FullName));

        _dirs.Clear();
        _dirs.AddRange(GetDirs());
    }

    /// <inheritdoc />
    public IEnumerator<DirectoryInfo> GetEnumerator() =>
        _dirs.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() =>
        _dirs.GetEnumerator();

    /// <summary>
    /// Adds the directory named “data” at the same location as the
    /// primary Cantera assembly.
    /// </summary>
    public void AddAssemblyDirectory() =>
        Add(Path.Combine(Path.GetDirectoryName(GetType().Assembly.Location)!, "data"));
}
