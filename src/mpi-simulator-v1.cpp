#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"

#define  MASTER		0

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
  MPI_Status status;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes specificed at start of run
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  tag1 = 1; //MASTER sends to Others
  tag2 = 2; //Others sends to MASTER

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
    //check if it's master thread
    //if MASTER: SEND particles to others + simulate step + RECV from others
    //if not MASTER: RECV from MASTER newest particles + simulate step + SEND to MASTER
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == MASTER){
      for (i = 1; i<nproc; i++){
        MPI_Send(&particles, particles.size(), MPI_PACKED, i, tag1, MPI_COMM_WORLD); //update new particles to others
      }
      QuadTree tree;
      QuadTree::buildQuadTree(particles, tree);
      simulateStep(tree, particles, newParticles, stepParams);

      for (i = 1; i<nproc; i++){
        MPI_Recv(&newParticles, newParticles.size(), MPI_PACKED, i, tag2, MPI_COMM_WORLD); //update new particles to others
      }
    }else{
      MPI_Recv(&particles, particles.size(), MPI_PACKED, MASTER, tag1, MPI_COMM_WORLD, &status);
      QuadTree tree;
      QuadTree::buildQuadTree(particles, tree);
      simulateStep(tree, particles, newParticles, stepParams);
      MPI_Send(&newParticles, newParticles.size(), MPI_PACKED, i, tag1, MPI_COMM_WORLD);
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
