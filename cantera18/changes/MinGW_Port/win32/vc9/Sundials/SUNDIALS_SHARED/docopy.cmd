echo on

set ccwwdd=%CD%
cd ..\..\..\..\..

cd sundials\include\sundials

if exist sundials_config.h exit 0

cd %ccwwdd%

copy sundials_config.h              ..\..\..\..\..\sundials\include\sundials

echo ok
echo off
