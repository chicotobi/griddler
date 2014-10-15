#include <clocale>
#include <iostream>

using namespace std;

int main() {
    setlocale(LC_ALL, "");
    wcout << "hi" << endl;
    wcout << L'\u25A0';
    wcout << L'\u25A0';
    wcout << L'\u25A1';
    wcout << L'\u25A0';
    wcout << L'\u25A0';
    wcout << L'\u25A1';
    wcout << L'\u25A0';
    wcout << L'\u25A0';
    wcout << L'\u25A1';
    wcout << L'\u25A1';
    wcout << L'\u25A1';
    wcout << L'\u25A0';
    wcout << L'\u25A0';
    wcout << endl;
    wcout << "hu" << endl;
    return 0;
}
