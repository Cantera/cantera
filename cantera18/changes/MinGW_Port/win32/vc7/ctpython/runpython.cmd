cd ..\..\..\Cantera\python
%PYTHON_CMD% winsetup.py build
%PYTHON_CMD% winsetup.py bdist_wininst
%PYTHON_CMD% winsetup.py install
echo 'ok' > status
