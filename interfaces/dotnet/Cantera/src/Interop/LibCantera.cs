namespace Cantera.Interop;

static partial class LibCantera
{
#if IS_WINDOWS
    const string LibFile = "cantera.2.6.0.dll";
#elif IS_MACOS
    const string LibFile = "cantera.2.6.0.dylib";
#elif IS_LINUX
    const string LibFile = "cantera.2.6.0.so";
#endif
}