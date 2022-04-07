#include "iostream"
using namespace std;

int main(){
    unsigned int N;
    cin >> N;
    while(N%10 == 0) {
        N /= 10; 
    }
    int c = 0;
    while(N > 9) {
        if (N%10 == 0) c++;
        N/=10;
    }
    cout << c << endl;
    return 0;
}