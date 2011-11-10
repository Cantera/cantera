cd ..\..\..\Cantera\matlab\cantera
echo 'delete me!' > ctmethods.mexw32
echo 'delete me!' > ctmethods.dll
%MATLAB_CMD% -nodisplay -nosplash -nojvm -r buildwin
echo 'ok' > status
