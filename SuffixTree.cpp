#include <iostream>
#include <string>
#include <map>

#include "SuffixTree.h"
#include "bwt.h"

using namespace std;


void SuffixTree::build(){
    // initializations
    current_end = new int[1];
    *(current_end) = 0;
    root = new SuffixNode(1, current_end);
    root->internalNode = true;
    
    remaining = 0;
    activeNode = root;
    activeLength = 0;
    activeEdge = -1;

    for(unsigned long i = 0; i < str.length(); ++i){
        addSuffix(i);
    }
}

int SuffixTree::tryNextCharacter(int i){
    // returns the next character if exists in the tree or -1 for the creation of a new edge
    SuffixNode* node = getNode();

    int lengthOfEdge = edgeLength(node);
    
    if (lengthOfEdge >= activeLength){
        return str[activeNode->children[mapOfGenes[str[activeEdge]]]->start + activeLength];
    }

    if (lengthOfEdge + 1 == activeLength){
        // node exists at this specific point
        if (node->children[mapOfGenes[str[i]]] != nullptr){
            return str[i];
        }
        return -1;
    } else{
        // traverse over node
        activeNode = node;
        activeLength = activeLength - lengthOfEdge - 1;
        activeEdge = activeEdge + lengthOfEdge + 1;

        return tryNextCharacter(i);
    }
}

void SuffixTree::traverse(int i){
    SuffixNode* node = getNode();

    if (edgeLength(node) < activeLength){
        activeNode = node;
        activeLength = activeLength - edgeLength(node);
        activeEdge = node->children[mapOfGenes[str[i]]]->start;
    } else{
        activeLength ++;
    }
}


SuffixNode* SuffixTree::getNode() {
    return activeNode->children[mapOfGenes[str[activeEdge]]]; 
}

SuffixNode* SuffixTree::getNode(int index) { 
    return activeNode->children[mapOfGenes[str[index]]]; 
}

void SuffixTree::addSuffix(int i){
    // cout << "i: " << i << " activeLength: " << activeLength << endl;
    SuffixNode* lastCreated = nullptr;

    *(current_end) = *(current_end) + 1;

    remaining ++;

    while (remaining > 0){
        // cout << "remaining: " << remaining << endl;
        if (activeLength == 0){
            // we are at root
            SuffixNode* node = getNode(i);
            if (node != nullptr){
                // RULE 3 extension
                activeEdge = node->start;
                activeLength ++;
                break;
            } else {
                // cout << "RULE2 upper" << endl;
                // RULE 2
                root->children[mapOfGenes[str[i]]] = new SuffixNode(i, current_end);
                remaining --;
            }
        } else {
            int nextChar = tryNextCharacter(i);
            // cout << "nextChar: " << nextChar << endl;

            if (char(nextChar) == str[i]){
                // cout << "Found it" << endl;
                if (lastCreated != nullptr)
                    lastCreated->suffixLink = getNode();

                traverse(i);
                break;
            } else if (nextChar < 0){
                SuffixNode* node = getNode();
                node->children[mapOfGenes[str[i]]] = new SuffixNode(i, current_end);

                if (lastCreated != nullptr)
                    lastCreated->suffixLink = node;

                lastCreated = node;

                if (activeNode != root)
                    activeNode = activeNode->suffixLink;
                else {
                    activeEdge ++;
                    activeLength --;
                }

                remaining --;
            } else{
                // cout << "Adding internal node.." << endl;
                // RULE 2
                SuffixNode* node = getNode();
                int oldStart = node->start;
                node->start = node->start + activeLength;

                int* thisEnd = new int[1];
                *(thisEnd) = oldStart + activeLength - 1;
                SuffixNode* newInternalNode = new SuffixNode(oldStart, thisEnd);

                SuffixNode* newLeafNode = new SuffixNode(i, current_end);

                newInternalNode->children[mapOfGenes[str[newInternalNode->start + activeLength]]] = node;
                newInternalNode->children[mapOfGenes[str[i]]] = newLeafNode;
                newInternalNode->internalNode = true;
                activeNode->children[mapOfGenes[str[newInternalNode->start]]] = newInternalNode;

                if (lastCreated != nullptr)
                    lastCreated->suffixLink = newInternalNode;

                lastCreated = newInternalNode;
                newInternalNode->suffixLink = root;

                if (activeNode != root)
                    activeNode = activeNode->suffixLink;
                else{
                    activeEdge ++;
                    activeLength --;
                }
                remaining --;
            }
        }
    }
}

// creates the bwt string and also creates the lcoation array needed
string SuffixTree::bwt(int* location_array){
    unsigned long str_len = str.length();

    char* transformation = new char[str_len + 1];
    unsigned long* pos = new unsigned long;
    *pos = 0;

    bwt(root, 0, str_len, pos, transformation, location_array);

    return string(transformation);
}

/*** recurrently creates bwt string and location array 
 * @param node current node
 * @param val current distance from root for node
 * @param size lenfth of the string
 * @param *pos current index to set for bwt and location array (this is based on the sorted suffixes)
*/
void SuffixTree::bwt(SuffixNode* node, int val, int size, unsigned long* pos, char* transformation, int* location_array){
    if (node == nullptr)
        return;

    if (!node->internalNode) {
        // cout << "Added: " << transformation[*pos] << " for position " << (size - (*(node->end) - node->start + val) - 1 + size) % size << endl;
        transformation[*pos] = str[(size - (*(node->end) - node->start + val) - 1 + size) % size];
        location_array[*pos] = node->start - val;

        (*pos) = (*pos) + 1;
        return;
    } else
        if (node != root)
            val += *(node->end) - node->start + 1;

    for(int i = 0; i < 5; ++i)
        bwt(node->children[i], val, size, pos, transformation, location_array);
}