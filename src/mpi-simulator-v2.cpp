#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"

int resizeIter = 4;

//quadTree is built from relativeParticles 
//newParticles is particles in the grid ==== myParticles
void simulateStep(const QuadTree &quadTree, std::vector<Particle> &myParticles,
                  std::vector<Particle> &newParticles, StepParameters params) {
  // TODO: paste your sequential implementation in Assignment 3 here.
  // (or you may also rewrite a new version)
  for (int i = 0; i < (int)myParticles.size(); i++) {
    auto pi = myParticles[i];
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
// gridSize is num of grids in one row/col
void assignGrid(std::vector<Particle> &totalParticles, Vec2 pmin, Vec2 pmax, 
    std::vector<Particle> &myParticles, int id, int gridSize){

  int row = id / gridSize;
  int col = id % gridSize;

  float topLeftx = pmin.x + (pmax.x - pmin.x) / gridSize * row;
  float topLefty = pmin.y + (pmax.y - pmin.y) / gridSize * col;
  float bottomRightx = topLeftx + (pmax.x - pmin.x) / gridSize;
  float bottomRighty = topLefty + (pmax.y - pmin.y) / gridSize;

  for (auto &p : totalParticles){
    float x = p.position.x;
    float y = p.position.y;
    if(x >= topLeftx && x < bottomRightx && y >= topLefty && y < bottomRighty){
      myParticles.push_back(p);
    }
  }

}


void assignAll(std::vector<Particle> &totalParticles, 
    std::vector<Particle> &myParticles, Vec2 pmin, 
    Vec2 pmax, int nproc, int pid, int gridSize){

  int sendSize = 0;
  MPI_Request request;
  if(pid == 0){
    for(int receiver = 1; receiver < nproc; receiver++){
      assignGrid(totalParticles, pmin, pmax, myParticles, pid, gridSize);
      sendSize = (int)myParticles.size();
      MPI_Isend(&sendSize, 1, MPI_INT, pid, 1, MPI_COMM_WORLD, &request);
      MPI_Isend(&myParticles, sendSize * sizeof(Particle), MPI_BYTE, pid, 1, MPI_COMM_WORLD, &request);
    }
    assignGrid(totalParticles, pmin, pmax, myParticles, 0, gridSize);
  }
  else{
    MPI_Recv(&sendSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    myParticles.resize(sendSize);
    MPI_Recv(&myParticles, sendSize, MPI_BYTE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}


// check if two grid has force
// l1, r1 are min, max of the grid we are currently on 
// l2, r2 are min, max of the grid we want to check if there is a force on current grid.
// radius the force affecting distance
bool forceCheck(Vec2 l1, Vec2 r1, Vec2 l2, Vec2 r2, float radius){

  // compute the actual affect grid
  l1.x += radius;
  l1.y += radius;
  r1.x += radius;
  r1.y += radius;

  // If one rectangle is on left side of other
  if (l1.x > r2.x || l2.x > r1.x)
      return false;

  // If one rectangle is above other
  if (r1.y > l2.y || r2.y > l1.y)
      return false;

  return true;
}

void constructRelatedP(std::vector<Particle> &myParticles, std::vector<Particle> &relativeParticles, 
        std::vector<float> &allLimits, Vec2 target_L, Vec2 target_R, int nproc, float radius){
  
  for(int i = 0; i < myParticles.size(); i++){
    relativeParticles[i] = myParticles[i];
  }

  Vec2 comp_L, comp_R;
  std::vector<bool> overlap;
  overlap.resize(nproc);
  for(int j = 0; j < nproc; j ++){
    comp_L.x = allLimits[j * 4];
    comp_L.y = allLimits[j * 4 + 1];
    comp_R.x = allLimits[j * 4 + 2];
    comp_R.y = allLimits[j * 4 + 3];
    if (forceCheck(target_L, target_R, comp_L, comp_R, radius)){
      overlap[j] = true; 
    } else{
      overlap[j] = false;
    }
  }
  // for(int k = 0; k < nproc; k++){
  //   if(overlap[k]){
  //     MPI_Isend()
  //   }
  // } 
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

  StartupOptions options = parseOptions(argc, argv);
  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);
  
  int gridSize = (int)sqrt(nproc);
  
  std::vector<Particle> totalParticles, relativeParticles, myParticles, newParticles;
  std::vector<float> allLimits;
  allLimits.resize(nproc * 4);

  if (pid == 0) { //load all particles to master 
    loadFromFile(options.inputFile, totalParticles); 
  }

  // Don't change the timing for totalSimulationTime.
  int *displs = (int *)malloc(nproc*sizeof(int));
  int *recvCounts = (int *)malloc(nproc*sizeof(int));
  int *particledispls = (int *)malloc(nproc*sizeof(int));
  int *particlerecvCounts = (int *)malloc(nproc*sizeof(int));
  Vec2 pmin(1e30f, 1e30f);
  Vec2 pmax(-1e30f, -1e30f);
  if(pid == 0){
    findBoundary(totalParticles, pmin, pmax);
  }
  assignAll(totalParticles,  myParticles, pmin, pmax, nproc, pid, gridSize);
  newParticles.resize(myParticles.size());
  recvCounts[pid] = myParticles.size() * sizeof(Particle);
  MPI_Barrier(MPI_COMM_WORLD);

  int sum = 0;
  for(int p = 0; p < nproc; p++){
    displs[p] = recvCounts[p] + sum;
    sum += recvCounts[p];
    particlerecvCounts[p] = 4 * sizeof(float);
    particledispls[p] = 4 * sizeof(float) * p;
  }

  Timer totalSimulationTimer;

  for (int i = 0; i < options.numIterations; i++) {
    if(i % resizeIter == 0){
      MPI_Allgatherv(myParticles.data(), myParticles.size() * sizeof(Particle), MPI_BYTE,
           totalParticles.data(), recvCounts, displs, MPI_BYTE, MPI_COMM_WORLD);
      findBoundary(totalParticles, pmin, pmax);
      assignAll(totalParticles, myParticles, pmin, pmax, nproc, pid, gridSize);
      newParticles.resize(myParticles.size());
    } 
    Vec2 myMin, myMax;
    findBoundary(myParticles, myMin, myMax);
    std::vector<float> boundary{myMin.x, myMin.y, myMax.x, myMax.y};
    MPI_Allgatherv(boundary.data(), sizeof(float) * 4, MPI_FLOAT, allLimits.data(), recvCounts,
      displs, MPI_FLOAT, MPI_COMM_WORLD);
    constructRelatedP(myParticles, relativeParticles, allLimits, myMin, myMax, nproc, stepParams.cullRadius);
    QuadTree tree;
    QuadTree::buildQuadTree(relativeParticles, tree);
    simulateStep(tree, myParticles, newParticles, stepParams);
  }
  MPI_Allgatherv(myParticles.data(), myParticles.size() * sizeof(Particle), MPI_BYTE,
           totalParticles.data(), recvCounts, displs, MPI_BYTE, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  double totalSimulationTime = totalSimulationTimer.elapsed();

  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, totalParticles);
  }

  MPI_Finalize();
}
