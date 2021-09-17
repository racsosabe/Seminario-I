#include<bits/stdc++.h>
using namespace::std;

set<string> used;

int main(){
	ios_base::sync_with_stdio(false);
	cin.tie(0);
	string s;
	int cnt = 0;
	while(getline(cin, s)){
		if(s.empty()) continue;
		if(used.count(s)) continue;
		cnt += 1;
		used.emplace(s);
		cout << s << endl;
	}
	cerr << cnt << endl;
	return 0;
}
