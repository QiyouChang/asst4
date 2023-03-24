#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"
#include <math.h>

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


//this function aims to find possibly overlapping region : L represents min; R represents max
bool findRelatedParticles(Vec2 target_L, Vec2 target_R, Vec2 comp_L, Vec2 comp_R){
  //4 cases left, right, top, bottom; otherwise, overlapping
  return (!((target_L.x > comp_R.x) || (target_R.x < comp_L.x) || (target_L.y > comp_R.y) || (target_R.y < comp_L.y)));
}


int main(int argc, char *argv[]) {
  int pid;
  int nproc;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes specificed at start of run
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int gridSize = sqrt(nproc);

  StartupOptions options = parseOptions(argc, argv);
  int Psize = sizeof(Particle);

  std::vector<Particle> totalParticles, relativeParticles, newParticles, myParticles;
  if (pid == 0) { //load all particles to master 
    loadFromFile(options.inputFile, totalParticles);
    size = (int)particles.size();
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);

  //send size
  MPI_Bcast(&size, 1, MPI_INT, 0, comm);
  particles.resize(size);
  //send particles to all processes
  MPI_Bcast(particles.data(), Psize * size, MPI_BYTE, MASTER, MPI_COMM_WORLD);
  //find boundary and assign particles to different grid
  Vec2 pmin(1e30f, 1e30f);
  Vec2 pmax(-1e30f, -1e30f);
  findBoundary(&particles, pmin, pmax);
  assignGrid(&particles, pmin, pmax, &myParticles, pid, gridSize);
  int gridPsize = myParticles.size();
  myParticles.resize(gridPsize); //??issue with resize
  
  // Don't change the timeing for totalSimulationTime.
  
  Timer totalSimulationTimer;

  for (int i = 0; i < options.numIterations; i++) {
    //check reassigning
    if(i % RESIZEFACTOR == 0){
      relativeParticles.resize(gridPsize);
      newParticles.resize(gridPsize);
      //Gathering first and thenn redo the process before the loop
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Scatterv(totalParticles.data(),(int)totalParticles.size()* sizeof(Particle), displs, MPI_BYTE,
          relativeParticles.data(), recvcounts, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }



    //check overlapping
    for (int index; index < nproc; index++){
      //??how to obtain the boundary of other processes instead of using communications 
      if ((index != pid)&&(findRelatedParticles())){
        send its own particles 
        recv the others particle
      }
    }
    combine related particles together //??how

    //simulate step
    QuadTree tree;
    QuadTree::buildQuadTree(related_particles, tree);
    simulateStep(tree, particles, newParticles, stepParams);

    
    //update boundary with current particles after simulating step
    for (auto &p : relativeParticles) {
      pmin.x = fminf(pmin.x, p.position.x);
      pmin.y = fminf(pmin.y, p.position.y);
      pmax.x = fmaxf(pmax.x, p.position.x);
      pmax.y = fmaxf(pmax.y, p.position.y);
    }
   
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
