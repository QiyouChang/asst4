#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"

#define MASTER 0

void simulateStep(const QuadTree &quadTree,
                  const std::vector<Particle> &particles,
                  std::vector<Particle> &newParticles, StepParameters params, int pid, int nproc, int chunksize) {
  if (pid<nproc-1){
    //normal cases
    for (int i = pid*chunksize; i < ((pid+1)*chunksize); i++) {
      auto pi = newParticles[i];
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
  }else{
    //Upper limit to the total size
    int leftover = (int)newParticles.size() - pid*chunksize;
    for (int i = pid*chunksize; i < (int)newParticles.size(); i++) {
      auto pi = newParticles[i];
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
}

int main(int argc, char *argv[]) {
  int pid;
  int nproc;
  MPI_Status status;
  MPI_Comm comm;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes specificed at start of run
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int tag1 = 1; //MASTER sends to Others
  int tag2 = 2; //Others sends to MASTER
  int *displs;
  int *recvcounts;

  displs = (int *)malloc(nproc*sizeof(int));
  recvcounts = (int *)malloc(nproc*sizeof(int));


  StartupOptions options = parseOptions(argc, argv);

  std::vector<Particle> particles, newParticles;
  if (pid == 0) {
    loadFromFile(options.inputFile, particles);
    loadFromFile(options.inputFile, newParticles);
  }

  int size = (int)newParticles.size();
  int chunksize = size/nproc;
  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);

  int sum = 0;
  for (int i = 0; i < nproc; i++){
    recvcounts[i] = (i < nproc-1) ? (chunksize) : ((int)newParticles.size()-(nproc-1)*chunksize);
    recvcounts[i] *= sizeof(Particle);
    displs[i] = sum;
    sum += recvcounts[i];
  }

  // Don't change the timeing for totalSimulationTime.
  MPI_Barrier(MPI_COMM_WORLD);
  Timer totalSimulationTimer;
  for (int i = 0; i < options.numIterations; i++) {
    //check if it's master thread
    //if MASTER: Broadcast + simulate step (self portion) + Gather to MASTER
    //if not MASTER: RECV from MASTER newest particles + simulate step (self portion) + Gather to MASTER
    //MPI_Barrier(MPI_COMM_WORLD);
    if (pid == MASTER){
    //   for (int t = 1; t<nproc; t++){
    //     MPI_Send(&particles, particles.size(), MPI_PACKED, t, tag1, MPI_COMM_WORLD); //update new particles to others
    //   }
      MPI_Bcast(&particles, sizeof(Particle)*particles.size(), MPI_BYTE, MASTER, comm);
      QuadTree tree;
      QuadTree::buildQuadTree(particles, tree);
      simulateStep(tree, particles, newParticles, stepParams, pid, nproc, chunksize);


      // for (int k = 1; k<nproc; i++){
      //   MPI_Recv(&newParticles, newParticles.size(), MPI_PACKED, k, tag2, MPI_COMM_WORLD, &status); //update new particles to others
      //   int limit = (k < nproc-1) ? ((k+1)*chunksize) : (int)newParticles.size();
      //   for (int n = k*chunksize; n<limit ; n++){
      //     particles[n] = newParticles[n];
      //   }
      // }

      MPI_Gatherv(&newParticles, chunksize*sizeof(Particle), MPI_BYTE, &particles, recvcounts, displs, MPI_BYTE, MASTER, comm);
    }else{
      //MPI_Recv(&particles, particles.size(), MPI_PACKED, MASTER, tag1, MPI_COMM_WORLD, &status);
      printf("pid:%d\n", pid);
      int size = (pid < nproc-1) ? (chunksize) : ((int)newParticles.size()-(nproc-1)*chunksize);
      QuadTree tree;
      QuadTree::buildQuadTree(particles, tree);
      simulateStep(tree, particles, newParticles, stepParams, pid, nproc, chunksize);
      //MPI_Send(&newParticles, newParticles.size(), MPI_PACKED, pid, tag2, MPI_COMM_WORLD);

      MPI_Gatherv(&newParticles, size*sizeof(Particle), MPI_BYTE, &particles, recvcounts, displs, MPI_BYTE, MASTER, comm);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double totalSimulationTime = totalSimulationTimer.elapsed();

  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }

  MPI_Finalize();
}
