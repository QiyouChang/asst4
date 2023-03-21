// #include "common.h"
// #include "mpi.h"
// #include "quad-tree.h"
// #include "timing.h"
// #include <iostream>

// #define MASTER 0

// void simulateStep(const QuadTree &quadTree,
//                   const std::vector<Particle> &particles,
//                   std::vector<Particle> &newParticles, StepParameters params, int pid, int nproc, int chunksize) {
//   //if (pid<nproc-1){
//     //normal cases
//     for (int i = pid*chunksize; i < ((pid+1)*chunksize); i++) {
//       auto pi = newParticles[i];
//       Vec2 force = Vec2(0.0f, 0.0f);
//       std::vector<Particle> temp; 
//       quadTree.getParticles(temp, pi.position, params.cullRadius); 
//       // accumulate attractive forces to apply to particle i
//       for (size_t j = 0; j < temp.size(); j++) {
//         if ((pi.position - temp[j].position).length() < params.cullRadius)
//           force += computeForce(pi, temp[j], params.cullRadius);
//       }
//       // update particle state using the computed force
//       newParticles[i] = updateParticle(pi, force, params.deltaTime);
//     }
//   // }else{
//   //   //Upper limit to the total size
//   //   //int leftover = (int)newParticles.size() - pid*chunksize;
//   //   for (int i = pid*chunksize; i < (int)newParticles.size(); i++) {
//   //     auto pi = newParticles[i];
//   //     Vec2 force = Vec2(0.0f, 0.0f);
//   //     std::vector<Particle> temp; 
//   //     quadTree.getParticles(temp, pi.position, params.cullRadius); 
//   //     // accumulate attractive forces to apply to particle i
//   //     for (size_t j = 0; j < temp.size(); j++) {
//   //       if ((pi.position - temp[j].position).length() < params.cullRadius)
//   //         force += computeForce(pi, temp[j], params.cullRadius);
//   //     }
//   //     // update particle state using the computed force
//   //     newParticles[i] = updateParticle(pi, force, params.deltaTime);
//   //   }
//   // }
// }

// int main(int argc, char *argv[]) {
//   int pid;
//   int nproc;
//   MPI_Status status;
//   MPI_Comm comm = MPI_COMM_WORLD;

//   // Initialize MPI
//   MPI_Init(&argc, &argv);
//   // Get process rank
//   MPI_Comm_rank(MPI_COMM_WORLD, &pid);
//   // Get total number of processes specificed at start of run
//   MPI_Comm_size(MPI_COMM_WORLD, &nproc);

//   // int tag1 = 1; //MASTER sends to Others
//   // int tag2 = 2; //Others sends to MASTER
//   int *displs;
//   int *recvcounts;

//   displs = (int *)malloc(nproc*sizeof(int));
//   recvcounts = (int *)malloc(nproc*sizeof(int));


//   StartupOptions options = parseOptions(argc, argv);

  
//   int size = 0;
//   std::vector<Particle> particles, newParticles;
//   if (pid == 0) {
//     loadFromFile(options.inputFile, particles);
//     loadFromFile(options.inputFile, newParticles);
//     size = (int)particles.size();
//   }

//   StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);
//   //MPI_Barrier(MPI_COMM_WORLD);
//   //std::cerr<<"before broadcast size:" << pid << "\n";
//   MPI_Bcast(&size, 1, MPI_INT, 0, comm);
//   //std::cerr<<"mid broadcast size:" << pid << "\n";
//   //MPI_Barrier(MPI_COMM_WORLD);
//   //std::cerr<<"after broadcast size:" << pid << "\n";

//   int chunksize = size/nproc;
//   int sum = 0;
//   for (int i = 0; i < nproc; i++){
    
//     recvcounts[i] = (i < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
//     recvcounts[i] *= sizeof(Particle);
//     displs[i] = sum;
//     // std::cerr<<"displs: " << displs[i] << "; with iteration: " << i << "\n";
//     // std::cerr<<"recvcounts: " << recvcounts[i] << "; with iteration: " << i << "\n";
//     sum += recvcounts[i];
//   }

  
  
//   // Don't change the timeing for totalSimulationTime.
//   MPI_Barrier(MPI_COMM_WORLD);
//   particles.resize(size);
//   newParticles.resize(size);
//   Timer totalSimulationTimer;

//   int childsize = chunksize;
//   for (int i = 0; i < options.numIterations; i++) {
//     // std::cerr<<"num iters:" << options.numIterations << "\n";
//     // std::cerr<<"pid:" << pid << "\n";
//     // std::cerr<<"master prev_bcast:" << pid << "\n";
//     // std::cerr<<"size of newPaticles:" << sizeof(newParticles)<< "\n";
//     // std::cerr<<"size of Particle:" << sizeof(Particle)<< "\n";
//     // std::cerr<<"number of particles: " << size << "\n";
//     MPI_Barrier(MPI_COMM_WORLD);
//     MPI_Bcast(particles.data(), sizeof(Particle) * size, MPI_BYTE, MASTER, comm);
//     //MPI_Barrier(MPI_COMM_WORLD);
//     //std::cerr<<"master post_bcast:" << pid << "\n";
//     // int childsize = (pid < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
//       QuadTree tree;
//       std::cerr<<"pre quadtree:" << pid << "\n";
//       QuadTree::buildQuadTree(particles, tree);
//       std::cerr<<"post quadtree:" << pid << "\n";
//       simulateStep(tree, particles, newParticles, stepParams, pid, nproc, chunksize);
//       std::cerr<<"child place2:" << pid << "\n";
//       MPI_Gatherv(newParticles.data(), childsize*sizeof(Particle), MPI_BYTE, particles.data(), recvcounts, displs, MPI_BYTE, MASTER, comm);
//       std::cerr<<"child place3:" << pid << "\n";
//     // }
//     std::cerr<<"finish iter :" << particles.size() << "\n";
//     // MPI_Barrier(MPI_COMM_WORLD);
//   }
 
//   MPI_Barrier(MPI_COMM_WORLD);
  

//   double totalSimulationTime = totalSimulationTimer.elapsed();

//   std::cerr<<"done!!!!" << pid << "\n";
//   if (pid == 0) {
//     printf("total simulation time: %.6fs\n", totalSimulationTime);
//     saveToFile(options.outputFile, particles);
//   }
  
//   MPI_Finalize();
// }






#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"
#include <iostream>

#define MASTER 0

void simulateStep(const QuadTree &quadTree,
                  const std::vector<Particle> &particles,
                  std::vector<Particle> &newParticles, StepParameters params) {

  for (int i = 0; i < newParticles.size(); i++) {
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

  // int tag1 = 1; //MASTER sends to Others
  // int tag2 = 2; //Others sends to MASTER
  int *displs;
  int *recvcounts;

  displs = (int *)malloc(nproc*sizeof(int));
  recvcounts = (int *)malloc(nproc*sizeof(int));

  StartupOptions options = parseOptions(argc, argv);
  
  int size = 0;
  std::vector<Particle> particles, newParticles;
  if (pid == 0) {
    loadFromFile(options.inputFile, particles);
    loadFromFile(options.inputFile, newParticles);
    size = (int)particles.size();
  }

  StepParameters stepParams = getBenchmarkStepParams(options.spaceSize);
  //MPI_Barrier(MPI_COMM_WORLD);
  //std::cerr<<"before broadcast size:" << pid << "\n";
  MPI_Bcast(&size, 1, MPI_INT, 0, comm);
  //std::cerr<<"mid broadcast size:" << pid << "\n";
  //MPI_Barrier(MPI_COMM_WORLD);
  //std::cerr<<"after broadcast size:" << pid << "\n";

  int chunksize = size/nproc;
  int sum = 0;
  for (int i = 0; i < nproc; i++){
    recvcounts[i] = (i < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
    recvcounts[i] *= sizeof(Particle);
    displs[i] = sum;
    sum += recvcounts[i];
    std::cerr<< "recvcounts: " << recvcounts[i] << "i is: "<< i  << "\n";
    std::cerr<< "displs: " << displs[i] << "i is: "<< i  << "\n";
  }

  int childsize = (pid < nproc-1) ? (chunksize) : (size-(nproc-1)*chunksize);
  // Don't change the timeing for totalSimulationTime.
  MPI_Barrier(MPI_COMM_WORLD);
  particles.resize(size);
  newParticles.resize(size);
  Timer totalSimulationTimer;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(particles.data(), sizeof(Particle) * size, MPI_BYTE, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(newParticles.data(), sizeof(Particle) * size, MPI_BYTE, MASTER, MPI_COMM_WORLD);

  for (int i = 0; i < options.numIterations; i++) {
    // std::cerr<<"num iters:" << options.numIterations << "\n";
    // std::cerr<<"pid:" << pid << "\n";
    // std::cerr<<"master prev_bcast:" << pid << "\n";
    // std::cerr<<"size of newPaticles:" << sizeof(newParticles)<< "\n";
    // std::cerr<<"size of Particle:" << sizeof(Particle)<< "\n";
    // std::cerr<<"number of particles: " << size << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
    // std::cerr<<"master post_bcast:" << pid << "size of particle" << sizeof(Particle) * size << "\n";
    QuadTree tree;
    QuadTree::buildQuadTree(particles, tree);
    
    // std::cerr<<"child place1:" << pid << "\n";
    simulateStep(tree, particles, newParticles, stepParams);
    // std::cerr<<"child place2:" << pid << "\n";
    MPI_Allgatherv(newParticles.data(), childsize*sizeof(Particle), MPI_BYTE, particles.data(), 
        recvcounts, displs, MPI_BYTE, comm);
    // std::cerr<<"finish iter :" << i << "\n";
    // std::cerr<<"size of newPaticles:" << sizeof(*newParticles.data())<< "\n";
    //MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  double totalSimulationTime = totalSimulationTimer.elapsed();

  std::cerr<<"done!!!!" << pid << "\n";
  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    saveToFile(options.outputFile, particles);
  }
  
  MPI_Finalize();
}
