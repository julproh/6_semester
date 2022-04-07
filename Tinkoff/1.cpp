#include "iostream"

using namespace std;

int main() {
    int N, A, B;
    cin >> A >> B >> N;

    if(A > B){
        cout << "Yes" << endl;
    } else if (A > (B/N + 1)) {
        cout << "Yes" << endl;
    }
    else{
        cout << "No" << endl;
    }
    return 0;
}