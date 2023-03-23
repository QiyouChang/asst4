#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"

void simulateStep(const QuadTree &quadTree,
                  std::vector<Particle> &newParticles, StepParameters params, int pid, int nproc, int chunksize) {
  // TODO: paste your sequential implementation in Assignment 3 here.
  // (or you may also rewrite a new version)
  for (int i = 0; i < newParticles.size(); i++) {
    auto pi = particles[i+pid*chunksize];
    Vec2 force = Vec2(0.0f, 0.0f);
    std::vector<Particle> temp; 
    quadTree.getParticles(temp, pi.position, params.cullRadius); 
    // accumulate attractive forces to apply to particle i
    for (size_t j = 0; j < temp.size(); j++) {
      if ((pi.position - temp[j].position).length() < params.cullRadius)
        force += computeForce(pi, temp[j], params.cullRadius);
    }
    // update particle state using the computed force
    newParticles[i] = updateParticle(pi, force, params.deltaTime);
  }
}

// This is the function in which we find the min and max of the particle list newParticles to 
// find the exact rectangular boundary of the list
// the result should be store in pmin and pmax.
void findBoundary(std::vector<Particle> &newParticles, Vec2 pmin, Vec2 pmax){
  
  pmin.x = 1e30f;
  pmin.y = 1e30f;
  pmax.x = -1e30f;
  pmax.y = -1e30f;

  for (auto &p : newParticles) {
    pmin.x = fminf(pmin.x, p.position.x);
    pmin.y = fminf(pmin.y, p.position.y);
    pmax.x = fmaxf(pmax.x, p.position.x);
    pmax.y = fmaxf(pmax.y, p.position.y);
  }
}

// This helper function will assign the totalParticles to each grid. The assignment will be stored in myParticles. 
// pmin, pmax are result from findBoundary(totalParticles)
void assignGrid(std::vector<Particle> &totalParticles, Vec2 pmin, Vec2 pmax, std::vector<Particle> &myParticles, int id, int gridSize){

  int row = id / gridSize;
  int col = id % gridSize;

  float topLeftx = pmin.x + (pmax.x - pmin.x) / gridSize * row;
  float topLefty = pmin.y + (pmax.y - pmin.y) / gridSize * col;
  float bottomRightx = topLeftx + (pmax.x - pmin.x) / gridSize;
  float bottomRighty = topLefty + (pmax.y - pmin.y) / gridSize;

  for (auto &p : totalParticles){
    float x = p.position.x;
    float y = p.position.y;
    if(x >= toLeftx && x < bottomRightx && y >= topLefty && y < bottomRighty){
      myParticles.push_back(p);
    }
  }
}

void findRelatedParticles()


int main(int argc, char *argv[]) {
  int pid;
  int nproc;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes specificed at start of run
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  StartupOptions options = parseOptions(argc, argv);

  std::vector<Particle> totalParticles, relativeParticles, newParticles;
  if (pid == 0) { //load all particles to master 
    loadFromFile(options.inputFile, totalParticles);
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);

  

  relativeParticles.resize(gridPsize);
  newParticles.resize(gridPsize);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Scatterv(totalParticles.data(),(int)totalParticles.size()* sizeof(Particle), displs, MPI_BYTE,
       relativeParticles.data(), recvcounts, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  
  // Don't change the timeing for totalSimulationTime.
  
  Timer totalSimulationTimer;

  for (int i = 0; i < options.numIterations; i++) {
    if(i % RESIZEFACTOR == 0){
      relativeParticles.resize(gridPsize);
      newParticles.resize(gridPsize);

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Scatterv(totalParticles.data(),(int)totalParticles.size()* sizeof(Particle), displs, MPI_BYTE,
          relativeParticles.data(), recvcounts, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
   

    for (auto &p : relativeParticles) {
      pmin.x = fminf(pmin.x, p.position.x);
      pmin.y = fminf(pmin.y, p.position.y);
      pmax.x = fmaxf(pmax.x, p.position.x);
      pmax.y = fmaxf(pmax.y, p.position.y);
    }
    // The following code is just a demonstration.
    QuadTree tree;
    QuadTree::buildQuadTree(particles, tree);
    simulateStep(tree, particles, newParticles, stepParams);
    particles.swap(newParticles);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double totalSimulationTime = totalSimulationTimer.elapsed();

  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }

  MPI_Finalize();
}
