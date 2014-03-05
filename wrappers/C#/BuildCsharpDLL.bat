REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

erase *_wrap.cpp

swig.exe -csharp -dllimport "CoolProp" -c++ -outcurrentdir ../../CoolProp/CoolProp.i
cl /c /I../../CoolProp /EHsc CoolProp_wrap.cxx

REM ~ REM ******* compile all the sources ***************
cl /c  /MP3 /I../../CoolProp /EHsc ../../CoolProp/*.cpp
link /DLL CoolProp_wrap.obj *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase CoolProp_wrap.cxx
erase CoolProp.lib
erase CoolProp.exp

rem **** Make a zip file using 7-zip ***
7z a -r Csharp.7z *.cs CoolProp.dll