#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <map>
#include "bwt.h"

class SuffixNode {
    friend class SuffixTree;
    private:
        SuffixNode* children[5];
        SuffixNode* suffixLink;

        bool internalNode;
        int start;
        int* end;

    public:
        SuffixNode(int s, int* e): start(s), end(e), children() {internalNode = false;};

};

class SuffixTree {
    private:
        std::string str;
        SuffixNode* root;
        int remaining;
        int* current_end;

        SuffixNode* activeNode;
        int activeEdge;
        int activeLength;

        void addSuffix(int);
        void traverse(int);
        SuffixNode* getNode();
        SuffixNode* getNode(int);
        int edgeLength(SuffixNode* node) { return *(node->end) - node->start; }
        int tryNextCharacter(int i);

        void bwt(SuffixNode*, int, int, unsigned long*, char*, int*);
    public:
        SuffixTree(std::string s): str(s) {};
        std::string bwt(int*);

        void build();
};

#endif
