echo on
cd ..\..\..
copy winconfig.h config.h
cd ext\f2c_libs
copy arith.hwin32 arith.h
copy sysdep1.h0 sysdep1.h
copy signal1.h0 signal1.h
cd ..\..\win32\vc8\config_h
echo ok
echo off
