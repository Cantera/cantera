This Matlab Toolbox for Cantera serves as a proof of concept for direct C-library integration into Matlab.

Move the toolbox_C folder into <installation directory to Cantera>\matlab.

For first time users, launch matlab, navigate to the toolbox folder, then run the 'SetCanteraPath(<installation directory to Cantera>)' command.

To start using the toolbox, use the LoadCantera command to load cantera_shared.dll into the memory.

To stop using Cantera and purge all objects from memory, first run cleanup command, then run UnloadCantera command.
