#ifndef QUAD_TREE_H
#define QUAD_TREE_H

#include "common.h"
#include <memory>

const int QuadTreeLeafSize = 8;
// NOTE: Do not remove or edit funcations and variables in this class definition
class QuadTreeNode {
public:
  bool isLeaf = 0;

  // four child nodes are stored in following order:
  //  x0, y0 --------------- x1, y0
  //    |           |           |
  //    |children[0]|children[1]|
  //    | ----------+---------  |
  //    |children[2]|children[3]|
  //    |           |           |
  //  x0, y1 ----------------- x1, y1
  // where x0 < x1 and y0 < y1.

  std::unique_ptr<QuadTreeNode> children[4];

  std::vector<Particle> particles;
};

inline float boxPointDistance(Vec2 bmin, Vec2 bmax, Vec2 p) {
  float dx = fmaxf(fmaxf(bmin.x - p.x, p.x - bmax.x), 0.0f);
  float dy = fmaxf(fmaxf(bmin.y - p.y, p.y - bmax.y), 0.0f);
  return sqrt(dx * dx + dy * dy);
}

// NOTE: Do not remove or edit funcations and variables in this class definition
class QuadTree {
public:
  std::unique_ptr<QuadTreeNode> root = nullptr;
  // the bounds of all particles
  Vec2 bmin, bmax;

  void getParticles(std::vector<Particle> &particles, Vec2 position,
                    float radius) const {
    particles.clear();
    getParticlesImpl(particles, root.get(), bmin, bmax, position, radius);
  }

  static void buildQuadTree(const std::vector<Particle> &particles,
                            QuadTree &tree) {
    // find bounds
    Vec2 bmin(1e30f, 1e30f);
    Vec2 bmax(-1e30f, -1e30f);

    for (auto &p : particles) {
      bmin.x = fminf(bmin.x, p.position.x);
      bmin.y = fminf(bmin.y, p.position.y);
      bmax.x = fmaxf(bmax.x, p.position.x);
      bmax.y = fmaxf(bmax.y, p.position.y);
    }

    // build nodes
    tree.bmin = bmin;
    tree.bmax = bmax;

    tree.root = buildQuadTreeImpl(particles, bmin, bmax);
  }

private:
  static std::unique_ptr<QuadTreeNode>
  buildQuadTreeImpl(const std::vector<Particle> &particles, Vec2 bmin,
                    Vec2 bmax) {
    // TODO: paste your sequential implementation in Assignment 3 here.
    // (or you may also rewrite a new version)
    std::unique_ptr<QuadTreeNode> quadTree = std::make_unique<QuadTreeNode>();
    
    if(particles.size() <= QuadTreeLeafSize){
      quadTree->isLeaf = 1;
      if(particles.size() != 0){
        quadTree->particles = particles; 
      }
      return quadTree;
    } 
    quadTree->isLeaf = 0;
    Vec2 pivot = Vec2((bmin.x + bmax.x) * 0.5f, (bmin.y + bmax.y) * 0.5f);
    std::vector<Particle> topLeft;
    std::vector<Particle> topRight;
    std::vector<Particle> bottomLeft;
    std::vector<Particle> bottomRight;

    for (size_t i = 0; i < particles.size(); i++){
      Particle p = particles.at(i);
      Vec2 position = p.position;

      if(position.x < pivot.x && position.y < pivot.y){
        topLeft.push_back(p);
      } else if(position.x < pivot.x && position.y >= pivot.y){
        bottomLeft.push_back(p);
      } else if(position.x >= pivot.x && position.y < pivot.y){
        topRight.push_back(p);
      } else{
        bottomRight.push_back(p);
      }
    }
    //printf("here6\n");
    quadTree->children[0] = buildQuadTreeImpl(topLeft, bmin, pivot);
    quadTree->children[1] = buildQuadTreeImpl(topRight, Vec2(pivot.x, bmin.y), Vec2(bmax.x, pivot.y));
    quadTree->children[2] = buildQuadTreeImpl(bottomLeft, Vec2(bmin.x, pivot.y), Vec2(pivot.x, bmax.y));
    quadTree->children[3] = buildQuadTreeImpl(bottomRight, pivot, bmax);

    return quadTree;
  }

  static void getParticlesImpl(std::vector<Particle> &particles,
                               QuadTreeNode *node, Vec2 bmin, Vec2 bmax,
                               Vec2 position, float radius) {
    if (node->isLeaf) {
      for (auto &p : node->particles)
        if ((position - p.position).length() < radius)
          particles.push_back(p);
      return;
    }
    Vec2 pivot = (bmin + bmax) * 0.5f;
    Vec2 size = (bmax - bmin) * 0.5f;
    for (int i = 0; i < 4; i++) {
      Vec2 childBMin;
      childBMin.x = (i & 1) ? pivot.x : bmin.x;
      childBMin.y = ((i >> 1) & 1) ? pivot.y : bmin.y;
      Vec2 childBMax = childBMin + size;
      if (boxPointDistance(childBMin, childBMax, position) <= radius)
        getParticlesImpl(particles, node->children[i].get(), childBMin,
                         childBMax, position, radius);
    }
  }
};

#endif
