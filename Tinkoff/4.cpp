#include "iostream"

using namespace std;

int main() {
    int n;
    cin >> n;
    int c = 97;
    char array[n];
    for (int i=0; i<(n+1)/2; i++){
        for (int j = 0; j<=i; j++) {
            array[i - j] = (char)((c+j-97)%26+97);
            array[n-i-1+j] = (char)((c+j-97)%26+97);  
        }
        for (int j = 0; j<(n+1)/2-i; j++) {
            array[i + j] = (char)((c+j-97)%26+97);
            array[n-1-i-j] = (char)((c+j-97)%26+97);  
        }
        for (int j = 0; j<n; j++) { cout << array[j]; }
        cout << endl;
    }

    for (int i=(n+1)/2; i<n; i++){
        for (int j = 0; j<=(n-i); j++) {
            array[i + j] = (char)((c+j-97)%26+97);
            array[n-i+j-1] = (char)((c+j-97)%26+97);
        }
        for (int j = 0; j<=(i-(n+1)/2); j++) { 
            array[i-j] = (char)((c+j-97)%26+97);
            array[n-i+j-1] = (char)((c+j-97)%26+97);  
        }
        for (int j = 0; j<n; j++) { cout << array[j]; }
        cout << endl;
    } 
    return 0;
}