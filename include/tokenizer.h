#pragma once

#include <utility>
#include <string>

namespace mrpgo{;

template<class T>
struct basic_string_iterator_pair : std::pair< typename std::basic_string<T>::const_iterator, typename std::basic_string<T>::const_iterator >{
	typedef typename std::basic_string<T>::const_iterator iterator_type;
	typedef std::pair<iterator_type, iterator_type>       base_type;
	
	iterator_type begin()const{ return this->first;  }
	iterator_type end  ()const{ return this->second; }
	bool          empty()const{ return this->second == this->first; }
	size_t        size ()const{ return this->second - this->first;  }

	basic_string_iterator_pair(){}
	basic_string_iterator_pair(iterator_type f, iterator_type s):base_type(f, s){}
	basic_string_iterator_pair(const std::basic_string<T>& c):base_type(c.begin(), c.end()){}
};

typedef basic_string_iterator_pair<char   >  string_iterator_pair;
typedef basic_string_iterator_pair<wchar_t> wstring_iterator_pair;

inline std::string to_string(string_iterator_pair str){
	return std::string(str.begin(), str.end());
}

inline std::wstring to_string(wstring_iterator_pair str){
	return std::wstring(str.begin(), str.end());
}

template<class T>
basic_string_iterator_pair<T> eat_white(basic_string_iterator_pair<T> str){
	if(str.empty())
		return str;

	std::locale loc;

	typename basic_string_iterator_pair<T>::iterator_type i0 = str.begin();
	typename basic_string_iterator_pair<T>::iterator_type i1 = str.end();
	while(i0 != i1 && std::isspace(*(i1-1), loc))
		i1--;
	while(i0 != i1 && std::isspace(*i0, loc))
		i0++;
	return basic_string_iterator_pair<T>(i0, i1);
}

template<typename T>
class TokenizerBase{
public:
	basic_string_iterator_pair<T> str;
	basic_string_iterator_pair<T> token;
	std::basic_string<T>          delim;
	bool                          compress;
	bool                          skip_par;

	bool IsDelim(T c){
		for(int i = 0; i < delim.size(); i++)
			if(delim[i] == c)
				return true;
		return false;
	}

	void SkipParenthesis(){
		int depth = 0;
		while(token.second != str.second){
			if(*token.second == (T)('('))
				depth++;
			if(*token.second == (T)(')'))
				depth--;
			token.second++;
			if(!skip_par || depth == 0)
				break;
		}
	}

public:
	static void Split(std::vector< std::basic_string<T> >& tokens, basic_string_iterator_pair<T> line, const std::basic_string<T>& d, bool comp, bool skip = false){
		tokens.clear();
		TokenizerBase<T> tok(line, d, comp, skip);
		while(!tok.IsEnd()){
			tokens.push_back(to_string(tok.GetToken()));
			tok.Next();
		}
	}

	void Set(basic_string_iterator_pair<T> s, const std::basic_string<T>& d, bool comp, bool skip = false){
		str      = s;
		delim    = d;
		compress = comp;
		skip_par = skip;
		token.first = token.second = str.first;

		while((token.second != str.second) && !IsDelim(*token.second))
			SkipParenthesis();
	}

	basic_string_iterator_pair<T> GetToken(){
		return token;
	}

	void Next(){
		token.first = token.second;
		if(IsEnd())
			return;

		do{
			token.first++;
		}
		while(compress && (token.first != str.second) && IsDelim(*token.first));
	
		token.second = token.first;
		while((token.second != str.second) && !IsDelim(*token.second))
			SkipParenthesis();
	}

	bool IsEnd(){
		return (token.first == str.second);
	}

	TokenizerBase(){
		delim.resize(1);
		delim[0] = (T)(' ');
		compress = false;
		skip_par = false;
	}

	TokenizerBase(basic_string_iterator_pair<T> s, const std::basic_string<T>& d, bool comp, bool skip = false){
		Set(s, d, comp, skip);
	}
};

typedef TokenizerBase<char>    Tokenizer;
typedef TokenizerBase<wchar_t> TokenizerW;

}
