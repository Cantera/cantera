# Configuring Visual Studio Code for Cantera

This page offers an opinionated take on setting up Visual Studio Code (VS Code or
vscode) for Cantera development. Feel free to take or leave any of these tips as you
wish. What works for you may be different.

```{seealso}
For information on running Cantera's test suites and debugging with VS Code,
see the sections [](sec-pytest-vscode) and [](sec-gtest-vscode).
```

## Settings

VS Code has two places it can store settings:

1. A global directory specific to your operating system, called "User Configuration" in
   VS Code
2. Within your local working directory, in a folder called `.vscode`, called "Workspace
   Configuration" in VS Code

In both places, settings are stored in a JSON file so you can edit them manually if you
want, or you can use the point-and-click interface. You can open the settings window by
typing `CTRL+,` (Windows, Linux) or `COMMAND+,` (macOS). You'll see the User and
Workspace settings tabs right below the search box. You can open the JSON settings file
by clicking the "Open Settings (JSON)" button in the top right of the editor pane; the
icon looks like a sheet of paper with an arrow wrapping around it.

Settings defined in your local configuration always override the global configuration.
Many commonly used settings are related to an Extension, so we'll describe those below.
The settings we recommend changing that come by default with VS Code are:

```json
{
    "files.insertFinalNewline": true,
    "files.trimFinalNewlines": true,
    "editor.renderWhitespace": "boundary",
    "files.trimTrailingWhitespace": true,
    "editor.formatOnSave": true,
    "editor.rulers": [
        88
    ]
}
```

You can change these in the User Configuration so they apply to all of your projects.
These five settings make sure that:

1. A file always ends with a newline character
2. A file always ends with a single newline character
3. Whitespace (spaces, tabs) are rendered when they are not between words, for instance
   at the beginning or ending of a line
4. Trailing whitespace is removed from each line
5. This automatic formatting is applied every time you save the file
6. Show the location of the 88th column, which is the preferred maximum line length in
   Cantera's source code.

VS Code also allows you to specify language-specific settings by putting the name of the
language as the mapping key inside square brackets. When VS Code detects a particular
language is used for a file, it applies the settings in that block. The following
settings may be useful for YAML files:

```json
{
    "[yaml]": {
        "editor.tabSize": 2,
    }
}
```

## Extensions

The following extensions may be useful for Cantera work:

- C/C++ (Microsoft)
- Cython (Thomas Walther)
- GitHub Pull Requests (GitHub)
- GitLens (GitKraken)
- Matlab (Xavier Hahn)
- Python (Microsoft)
- Visual Studio IntelliCode (Microsoft)

We'll talk a little more about some of them below.

### GitHub Pull Requests

This extension makes it much easier to do code reviews on big pull requests from within
VS Code, rather than trying to work in the browser.

### GitLens

GitLens allows you to see the commit history for a line in a file from within VS Code. Super useful when you want to start yelling at the terrible developer who last touched this line, because you don't have to go to your browser to find out it was yourself.

### Python

The Python extension is really excellent because it includes the ability to run linters and formatters on Python code as you write it. Unfortunately, it's not super helpful for Cantera development because we write most code in Cython. This is most useful for working on examples.

## Other General Tips/Keyboard Shortcuts

- You can use `F12` (you might need to press `Fn` for your keyboard to actually send
  `F12`) to jump to the definition of the function or class under your cursor. Really,
  really useful in C++ land where everything seems like it's a virtual function defined
  on the parent class.
- You can press `CTRL+k` (Windows, Linux) or `COMMAND+k` (macOS) and then `s` to save a
  file without auto-formatting it, if you have format-on-save enabled by default.
- You can press `CTRL+G` (macOS, yes that means Control not Command, not sure about the
  other OS) to jump to a specific line in the open file.
- `CTRL+SHIFT+T` (Windows, Linux) and `COMMAND+SHIFT+T` (macOS) reopens a tab you just
  closed. This works in Chrome and Firefox too.