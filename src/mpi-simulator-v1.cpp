
#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"
#include <iostream>

#define MASTER 0

void simulateStep(const QuadTree &quadTree,
                  const std::vector<Particle> &particles,
                  std::vector<Particle> &newParticles, StepParameters params, int pid, int nproc, int chunksize) {
  int displacement = pid*chunksize;
  std::vector<Particle> temp; 
  for (int i = 0; i < newParticles.size(); i++) {
    auto pi = particles[i+displacement];
    Vec2 force = Vec2(0.0f, 0.0f);
    temp.clear();
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

int main(int argc, char *argv[]) {
  int pid;
  int nproc;
  MPI_Status status;
  MPI_Comm comm = MPI_COMM_WORLD;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes specificed at start of run
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int displs[nproc];
  int recvcounts[nproc];

  StartupOptions options = parseOptions(argc, argv);
  
  int size = 0;
  int Psize = sizeof(Particle);
  std::vector<Particle> particles, newParticles;
  if (pid == 0) {
    loadFromFile(options.inputFile, particles);
    size = (int)particles.size();
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);
  MPI_Bcast(&size, 1, MPI_INT, 0, comm);

  int chunksize = size/nproc;
  std::fill_n(recvcounts, nproc, (chunksize)*Psize);
  recvcounts[nproc-1] = (size-(nproc-1)*chunksize)*Psize;
  for (int i = 0; i < nproc; i++){
    displs[i] = (i == 0) ? 0 : recvcounts[i-1]+displs[i-1];
  }

  
  particles.resize(size);
  int childsize = (pid < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
  newParticles.resize(childsize); 
  

  MPI_Bcast(particles.data(), sizeof(Particle) * size, MPI_BYTE, MASTER, MPI_COMM_WORLD);
 
  QuadTree tree;
  // Don't change the timeing for totalSimulationTime.
  Timer totalSimulationTimer;
  
  for (int i = 0; i < options.numIterations; i++) {

    QuadTree::buildQuadTree(particles, tree);
    simulateStep(tree, particles, newParticles, stepParams, pid, nproc, chunksize);
    MPI_Allgatherv(newParticles.data(), childsize*Psize, MPI_BYTE, particles.data(), 
        recvcounts, displs, MPI_BYTE, comm);
  }
  if (pid == 0) {
    double totalSimulationTime = totalSimulationTimer.elapsed();
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }
  
  MPI_Finalize();
}
