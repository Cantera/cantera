cd ..\..\..\Cantera\python
%PYTHON_CMD% winsetup.py build
copy ..\..\build\lib\i686-pc-win32\clib_d.dll build\lib.win32-2.6\Cantera
copy ..\..\build\lib\i686-pc-win32\clib_d.exp build\lib.win32-2.6\Cantera
copy ..\..\build\lib\i686-pc-win32\clib_d.lib build\lib.win32-2.6\Cantera
%PYTHON_CMD% winsetup_d.py bdist_wininst
%PYTHON_CMD% winsetup_d.py install
echo 'ok' > status
