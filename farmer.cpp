#include<bits/stdc++.h>
using namespace::std;

struct number{
	int signo;
	string entero;
	string decimal;
	number(){};
	number(string &s){
		signo = 1;
		if(s[0] == '-'){
			signo = -signo;
			s = s.substr(1);
		}
		int pos = 0;
		while(pos < s.size() and s[pos] != '.') pos += 1;
		decimal = "0";
		entero = "0";
		entero = s.substr(0,pos);
		if(pos+1 < s.size()) decimal = s.substr(pos+1);
	}
};

int T = 10000;

long double N(int n){
	return T * log(n) / (2.0 * acos(-1));
}

int menor(string a, string b){
	if(a.size() < b.size()) return 1;
	if(a.size() > b.size()) return 0;
	if(a < b) return 1;
	if(b < a) return 0;
	return 2;
}

bool menor2(string a, string b){
	int pos = 0;
	while(pos < min(a.size(),b.size()) and a[pos] == b[pos]){
		pos += 1;
	}
	if(pos < min(a.size(), b.size())) return a[pos] < b[pos];
	if(pos < a.size()) return false;
	return true;
}

bool comp(number a, number b){
	if(a.signo == -1 and b.signo == 1) return true;
	if(a.signo == 1 and b.signo == -1) return false;
	if(a.signo == 1){
		if(menor(a.entero,b.entero) == 1) return true;
		if(menor(b.entero,a.entero) == 1) return false;
		return menor2(a.decimal,b.decimal);
	}
	else{
		if(menor(a.entero,b.entero) == 1) return false;
		if(menor(b.entero,a.entero) == 1) return true;
		return menor2(b.decimal,a.decimal);
	}
}

void print(number &x){
	if(x.signo == -1) putchar('-');
	printf("%s.%s",x.entero.c_str(),x.decimal.c_str());
}

int main(int argc, char* argv[]){
	int n = atoi(argv[1]);
	vector<number> Re;
	vector<number> Im;
	string s;
	while(cin >> s){
		Re.emplace_back(number(s));
		cin >> s;
		Im.emplace_back(number(s));
	}
	sort(Re.begin(),Re.end(),comp);
	sort(Im.begin(),Im.end(),comp);
	printf("%d &%d &%d &%.2Lf &",n,T,(int)Re.size(),N(n));
	print(Re[0]);
	printf(" &");
	print(Re.back());
	printf(" \\\\\n");
	puts("Imaginary part:");
	print(Im[0]);
	puts("");
	print(Im.back());
	return 0;
}
