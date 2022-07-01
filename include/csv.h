#pragma once

#include <tokenizer.h>

#include <fstream>

namespace mrpgo{;

class CsvReader{
public:
	struct Row{
		std::string line;
		std::vector<string_iterator_pair>	col;
	};
	std::vector<Row>	row;

	std::string contents;			///< contents string

	Tokenizer rowIt;
	Tokenizer colIt;

public:
	void Clear(){
		contents.clear();
		row.clear();
	}

	bool Read(const std::string& filename, const std::string& delim, bool compress = false){
		std::ifstream ifs;
		ifs.open(filename, std::ifstream::in | std::ifstream::binary);

		if(!ifs.is_open())
			return false;

		ifs.seekg(0, ifs.end);
		int len = (int)ifs.tellg();
		ifs.seekg(0, ifs.beg);

		int nrow = 0;
		std::string line;
		while(!ifs.eof()){
			std::getline(ifs, line);
			nrow++;
		}

		row.resize(nrow);

		ifs.clear();
		ifs.seekg(0, ifs.beg);
		for(int i = 0; !ifs.eof(); i++){
			Row& r = row[i];
			
			std::getline(ifs, r.line);
			if(r.line.empty())
				continue;

			colIt.Set(r.line, delim, compress);
			int ncol = 0;
			for(; !colIt.IsEnd(); colIt.Next())
				ncol++;
			r.col.resize(ncol);
			
			colIt.Set(r.line, delim, compress);
			for(int j = 0; !colIt.IsEnd(); colIt.Next(), j++){
				r.col[j] = colIt.GetToken();
			}
		}

		while(!row.empty() && row.back().col.empty())
			row.pop_back();
			
		return true;
	}

	int	NumRow(){
		return row.size();
	}

	int	NumCol(int i){
		return row[i].col.size();
	}

	int FindColumn(std::string label){
		if(NumRow() == 0)
			return -1;
		for(int c = 0; c < NumCol(0); c++)
			if(Get<std::string>(0, c) == label)
				return c;
		return -1;
	}

	template<class T>
	T Get(string_iterator_pair str){
		T val;
		return val;
	}

	template<class T>
	T Get(int i, int j){
		if(i >= row.size())
			return T();
		if(j >= row[i].col.size())
			return T();

		return Get<T>(row[i].col[j]);
	}

};

template<>
inline int CsvReader::Get<int>(string_iterator_pair str){
	return atoi(to_string(str).c_str());
}
template<>
inline float CsvReader::Get<float>(string_iterator_pair str){
	return (float)atof(to_string(str).c_str());
}
template<>
inline double CsvReader::Get<double>(string_iterator_pair str){
	return (double)atof(to_string(str).c_str());
}
template<>
inline std::string CsvReader::Get<std::string>(string_iterator_pair str){
	return to_string(str);
}

}
