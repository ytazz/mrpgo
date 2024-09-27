#pragma once

#include <vector>
#include <string>

namespace mrpgo{;

class Node;
class Link;

typedef std::vector< Node* >  Nodes;
typedef std::vector< Link* >  Links;

class Posegraph{
public:
    struct Space{
		enum{
			SE2,
			SE3,
		};
	};

    int       space;
 
    Nodes     nodes;
	Links     links;

    std::vector<Links>  linksDisjoint;

    struct Depth{
        Nodes         nodes;
    };
    
    std::vector<Depth>  depths;

    int       dim;
	int       dim_trn;
	int       dim_rot;

public:
    void  SetSpace(std::string str);
	void  Analyze ();
    void  PrintCsv();

    Posegraph();
};

}
