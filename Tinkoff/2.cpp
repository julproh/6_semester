#include "iostream"

using namespace std;

int main(){
    int abc[3] = {0,0,0};
    for (int i = 0; i<3; i++) {
        char symbol;
        cin >> symbol;
        switch (symbol){
            case '>': 
            {
                if(i != 2) abc[0]++;
                else abc[1]++;
                break;
            }
            case '<':
            {
                cout << '<' << endl;
                if(i != 2) abc[i+1]++;
                else abc[2]++;
                break;
            }
        }
    }

        if(abc[0] >abc[1]) {
        if (abc[1] > abc[2]) cout << "cba" << endl;
        else if (abc[1] == abc[2]) cout << "bca" << endl << "cba" << endl;
        else if (abc[1] < abc[2]) {
            if (abc[0] > abc[2]) cout << "bca" << endl;
            else if (abc[0] < abc[2]) cout << "bac" << endl;
            else cout << "bac" << endl << "bca" << endl;
        }
    }
    else if(abc[0] == abc[1]){
        if (abc[1] == abc[2]) cout << "abc" << endl << "acb" 
                                << endl << "bac" << endl << "bca" << endl << "cab" << "cba" << endl; 
        else if (abc[1] < abc[2]) cout << "abc" << endl << "bac" << endl;
        else if (abc[1] > abc[2]) cout << "cab" << endl << "cba" << endl;
    }
    else if (abc[0]<abc[1]){
        if (abc[1] == abc[2]) cout << "abc" << endl << "acb"<< endl;
        else if (abc[1] > abc[2]) {
            if (abc[0] == abc[2]) cout << "acb" << endl << "cab" << endl;
            else if (abc[0] > abc[2]) cout << "cab" << endl;
            else cout << "acb" << endl; 
        }
        else cout << "abc" << endl;
    }

    return 0;
}