cd ..\..\..\Cantera\python
%PYTHON_CMD% winsetup.py build
copy ..\..\build\lib\i686-pc-win32\clib.dll build\lib.win32-2.6\Cantera
copy ..\..\build\lib\i686-pc-win32\clib.exp build\lib.win32-2.6\Cantera
copy ..\..\build\lib\i686-pc-win32\clib.lib build\lib.win32-2.6\Cantera
%PYTHON_CMD% winsetup.py bdist_wininst
%PYTHON_CMD% winsetup.py install
echo 'ok' > status
