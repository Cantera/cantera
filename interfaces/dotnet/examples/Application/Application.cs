/*
 * C# Application Example
 * ======================
 *
 * Demonstrates how to retrieve and modify global properties of the Application class
 * using Cantera's .NET API.
 *
 * .. tags:: .NET
 */

/// This file is part of Cantera. See License.txt in the top-level directory or
/// at https://cantera.org/license.txt for license and copyright information.

using Cantera;

// read the current Version and the git commit
Console.WriteLine("Version: " + Application.Version);
Console.WriteLine("Git Commit: " + Application.GitCommit);

// modify and read the Cantera data directories
var homeDir = Environment.GetFolderPath(Environment.SpecialFolder.UserProfile);
Application.DataDirectories.Add(homeDir);
Application.DataDirectories.AddAssemblyDirectory();

Console.WriteLine("Data Directories:");

foreach(var dir in Application.DataDirectories)
{
    Console.WriteLine("    " + dir);
}
