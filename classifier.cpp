#include<bits/stdc++.h>
using namespace::std;

struct zero{
	string real;
	string imag;
	string error;
	zero(){}
	zero(string real, string imag, string error) :
		real(real), imag(imag), error(error) {}
};

vector<zero> v;

int get(string &s){
	int ans = 0;
	for(int i = 0; i < s.size(); i++){
		if(s[i] == '.') break;
		ans = 10 * ans + s[i] - '0';
	}
	return ans;
}

int main(){
	string real, imag, error;
	while(cin >> real){
		cin >> imag >> error;
		v.emplace_back(zero(real, imag, error));
	}
	sort(v.begin(), v.end(), [&] (zero a, zero b){
		int imag_a = get(a.imag);
		int imag_b = get(b.imag);
		if(imag_a / 10000 == imag_b / 10000) return false;
		return (imag_a / 10000) < (imag_b / 10000);
	});
	for(auto e : v){
		cout << e.real << " " << e.imag << " " << e.error << endl;
	}
	map<int, int> frec;
	for(auto e : v){
		frec[get(e.imag) / 10000]++;
	}
	for(auto e : frec){
		cerr << "Im range [" << e.first * 10000 << ", " << e.first * 10000 + 10000 - 1 << "] -> " << e.second << " zeros" << endl;
	}
	return 0;
}
