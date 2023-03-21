#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"
#include <iostream>

#define MASTER 0

void simulateStep(const QuadTree &quadTree,
                  const std::vector<Particle> &particles,
                  std::vector<Particle> &newParticles, StepParameters params, int pid, int nproc, int chunksize) {

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

  int *displs;
  int *recvcounts;
  int Psize = sizeof(Particle);

  displs = (int *)malloc(nproc*sizeof(int));
  recvcounts = (int *)malloc(nproc*sizeof(int));

  StartupOptions options = parseOptions(argc, argv);
  
  int size = 0;
  std::vector<Particle> particles, newParticles;
  if (pid == 0) {
    loadFromFile(options.inputFile, particles);
    //loadFromFile(options.inputFile, newParticles); //********
    size = (int)particles.size();
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);
  MPI_Bcast(&size, 1, MPI_INT, 0, comm);


  int chunksize = size/nproc;
  int sum = 0;
  for (int i = 0; i < nproc; i++){
    recvcounts[i] = (i < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
    recvcounts[i] *= Psize;
    displs[i] = sum;
    sum += recvcounts[i];
  }

  // Don't change the timeing for totalSimulationTime.
  particles.resize(size);
  newParticles.resize(chunksize); //********
  Timer totalSimulationTimer;

  MPI_Bcast(particles.data(), Psize * size, MPI_BYTE, MASTER, MPI_COMM_WORLD);
  //MPI_Bcast(newParticles.data(), Psize * chunksize, MPI_BYTE, MASTER, MPI_COMM_WORLD);
  
  for (int i = 0; i < options.numIterations; i++) {
    QuadTree tree;
    QuadTree::buildQuadTree(particles, tree);
    simulateStep(tree, particles, newParticles, stepParams, pid, nproc, chunksize);
    MPI_Allgatherv(newParticles.data(), chunksize*Psize, MPI_BYTE, particles.data(), 
        recvcounts, displs, MPI_BYTE, comm);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  double totalSimulationTime = totalSimulationTimer.elapsed();

  // std::cerr<<"done!!!!" << pid << "\n";
  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }
  
  MPI_Finalize();
}

