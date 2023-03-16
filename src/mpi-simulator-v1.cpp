#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"

void simulateStep(const QuadTree &quadTree,
                  const std::vector<Particle> &particles,
                  std::vector<Particle> &newParticles, StepParameters params) {
  // TODO: paste your sequential implementation in Assignment 3 here.
  // (or you may also rewrite a new version)

  for (int i = 0; i < (int)newParticles.size(); i++) {
      auto pi = newParticles[i];
      Vec2 force = Vec2(0.0f, 0.0f);
      std::vector<Particle> temp; 
      quadTree.getParticles(temp, pi.position, params.cullRadius); 
      // accumulate attractive forces to apply to particle i
      for (size_t j = 0; j < temp.size(); j++) {
        if ((pi.position - temp[j].position).length() < params.cullRadius)
        //   force += computeForce(pi, particles[j], params.cullRadius);
          force += computeForce(pi, temp[j], params.cullRadius);
      }
      // update particle state using the computed force
      newParticles[i] = updateParticle(pi, force, params.deltaTime);
  }

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

  std::vector<Particle> particles, newParticles;
  if (pid == 0) {
    loadFromFile(options.inputFile, particles);
    loadFromFile(options.inputFile, newParticles);
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);

  // Don't change the timeing for totalSimulationTime.
  MPI_Barrier(MPI_COMM_WORLD);
  Timer totalSimulationTimer;
  for (int i = 0; i < options.numIterations; i++) {
    // The following code is just a demonstration.
    QuadTree tree;
    QuadTree::buildQuadTree(particles, tree);
    MPI_Comm com;
    MPI_Status stat;
    int particleIndex;
    simulateStep(tree, particles, newParticles, stepParams);
    printf ("Hello from task %d!\n", pid);
    for(int j = 0; j < nproc; j++){
      if(j != pid){
        MPI_Send(&newParticles, particles.size(), MPI_PACKED, j, pid, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int p = 0; p < nproc; p ++){
      std::vector<Particle> inMsg;
      if(p != pid){
        MPI_Recv(&inMsg, particles.size(), MPI_PACKED, p, pid, MPI_COMM_WORLD, &stat);
        for(int num = 0; num < inMsg.size()/nproc; num ++){
          particleIndex = num * nproc + p;
          newParticles[particleIndex] = inMsg[particleIndex];
        }
      }
    }
    particles.swap(newParticles);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double totalSimulationTime = totalSimulationTimer.elapsed();

  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }

  MPI_Finalize();
}
