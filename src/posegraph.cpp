#include <posegraph.h>
#include <node.h>
#include <link.h>

using namespace std;

namespace mrpgo{;

Posegraph::Posegraph(){
	space     = Space::SE3;
}

void Posegraph::SetSpace(string str){
	if(str == "se2") space = Space::SE2;
	if(str == "se3") space = Space::SE3;
}

void Posegraph::Analyze(){
	for(int i = 0; i < (int)nodes.size(); i++){
		Node* n = nodes[i];
		n->index     = i;
        n->depth     = -1;
        n->linkUpper = 0;
	}
    for(Link* l : links){
        l->node[0]->links.push_back(l);
        l->node[1]->links.push_back(l);
    }

	// disjoint subsets of links
	vector<bool> connected(nodes.size());
	Links links_cur, links_next;
	links_cur.insert(links_cur.end(), links.begin(), links.end());
	while(!links_cur.empty()){
		linksDisjoint.push_back(Links());
		links_next.clear();
		fill(connected.begin(), connected.end(), false);

		for(Link* l : links_cur){
			if( !connected[l->node[0]->index] &&
				!connected[l->node[1]->index] ){
				linksDisjoint.back().push_back(l);
				connected[l->node[0]->index] = true;
				connected[l->node[1]->index] = true;
			}
			else{
				links_next.push_back(l);
			}
		}
		
		links_next.swap(links_cur);
	}
    
    // create spanning-tree hierarchy
    int d = 0;
    nodes[0]->depth = d;
    depths.resize(1);
    depths[0].nodes.push_back(nodes[0]);
    while(!depths[d].nodes.empty()){
        depths.push_back(Depth());
        
        for(Node* n : depths[d].nodes){
            for(Link* l : n->links){
                Node* n1 = l->Opposite(n);
                if(n1->depth == -1){
                    n1->linkUpper = l;
                    n1->depth = d+1;
                    depths[d+1].nodes.push_back(n1);
                }
            }
        }
        d++;
    }
    depths.pop_back();
   
    for(int d = 0; d < depths.size(); d++){
        //DSTR << "depth " << d << " # nodes: " << depths[d].nodes.size() << endl;
    }

    if(space == Space::SE2){
		dim_trn = 2;
		dim_rot = 1;
	}
	else{
		dim_trn = 3;
		dim_rot = 3;
	}
	dim = dim_trn + dim_rot;

}

}
