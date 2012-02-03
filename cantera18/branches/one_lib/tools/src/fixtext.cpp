#include <string>
#include <iostream>
#include <stdio.h>
using namespace std;

////////////////////////////////////////////////////////////////////////
//
//  Convert a text file created on any system (unix, Mac, or Windows)
//  to the line-ending convention used on this machine.  The input is
//  read from the standard input, and the output is written to the
//  standard output.
//
//  > convertfile < input_file.txt > output_file.txt
//
///////////////////////////////////////////////////////////////////////

int main() {
    char ch;
    char last_eol=' ';
    const char char10 = char(10);
    const char char13 = char(13);
    string line;
    while (1 > 0) {
        line = "";
        while (1 > 0) {
            cin.get(ch);
            if (cin.eof()) break;
            if (ch == char13 || (ch == char10 
                    && (last_eol != char13)))  {
                last_eol = ch;
                break;
            }
            if (isprint(ch)) line += ch;
            if (ch == '\t') {line += ' '; line += ' '; line += ' '; line += ' ';}
        }
        cout << line << endl;
        if (cin.eof()) break;
    }
}
