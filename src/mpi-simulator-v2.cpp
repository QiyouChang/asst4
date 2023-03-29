


#include "common.h"
#include "mpi.h"
#include "quad-tree.h"
#include "timing.h"
#include <algorithm>

int resizeIter = 100; 

//quadTree is built from relativeParticles 
//newParticles is particles in the grid ==== myParticles
void simulateStep(const QuadTree &quadTree, std::vector<Particle> &myParticles,
                  std::vector<Particle> &newParticles, StepParameters params) {
  // TODO: paste your sequential implementation in Assignment 3 here.
  // (or you may also rewrite a new version)
  //// std::cerr<< newParticles.size() << " " << myParticles.size() << "\n";
  std::vector<Particle> temp; 
  for (int i = 0; i < (int)newParticles.size(); i++) {
    auto pi = myParticles[i];
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


// This helper function will assign the totalParticles to each grid. The assignment will be stored in myParticles. 
// pmin, pmax are result from findBoundary(totalParticles)
// gridSize is dimension of grids in one row/col
void assignGrid(std::vector<Particle> &totalParticles, Vec2 pmin, Vec2 pmax, 
    std::vector<Particle> &myParticles, int id, int gridSize){

  int row = id / gridSize;
  int col = id % gridSize;
  
  float topLeftx = pmin.x + (pmax.x - pmin.x) * row / gridSize ;
  float topLefty = pmin.y + (pmax.y - pmin.y) * col / gridSize ;
  float bottomRightx = pmin.x + (pmax.x - pmin.x) * (row+1.0f) / gridSize ;
  float bottomRighty = pmin.y + (pmax.y - pmin.y) * (col+1.0f) / gridSize ;
  
  myParticles.resize(totalParticles.size());
  int count = 0;
  
  for (auto &p : totalParticles){
    
    float x = p.position.x;
    float y = p.position.y;
    if(x >= topLeftx && x < (bottomRightx) && y >= topLefty && y < bottomRighty){
     
      myParticles[count] = p;
      count+= 1;
    }
  }
  myParticles.resize(count);
  
}

//this is the function where all the particles are assigned to each individual grid
void assignAll(std::vector<Particle> &totalParticles, 
    std::vector<Particle> &myParticles, Vec2 pmin, 
    Vec2 pmax, int nproc, int pid, int gridSize){

  int sendSize = 0;
  std::vector<MPI_Request> request;
  MPI_Status status[nproc * 2 - 2];
  if(pid == 0){
    for(int receiver = 1; receiver < nproc; receiver++){
      //first compute particles for the receiver which stored in myParticels
      
      std::vector<Particle> thisParticles;
      assignGrid(totalParticles, pmin, pmax, thisParticles, receiver, gridSize);
      //// std::cerr<<"my particles in assign all size is:  %d" << (int)myParticles.size() << "\n";
      sendSize = (int)thisParticles.size();
      
      //then send the size of particles and later the exact array of particles.
      MPI_Send(&sendSize, 1, MPI_INT, receiver, 1, MPI_COMM_WORLD);
      MPI_Send(thisParticles.data(), sendSize * sizeof(Particle), MPI_BYTE, receiver, 2, MPI_COMM_WORLD);
    }
    //assign master itself particles
    assignGrid(totalParticles, pmin, pmax, myParticles, 0, gridSize);
    //MPI_Waitall(nproc * 2 - 2, request.data(), status);
   
  }
  else{
    //if not the master, receive the particle size first and then receive the entire particle list
    MPI_Recv(&sendSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    myParticles.resize(sendSize);
    
    MPI_Recv(myParticles.data(), sendSize * sizeof(Particle), MPI_BYTE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

// check if two grid has force
// l1, r1 are min, max of the grid we are currently on 
// l2, r2 are min, max of the grid we want to check if there is a force on current grid.
// radius the force affecting distance
bool forceCheck(Vec2 l1, Vec2 r1, Vec2 l2, Vec2 r2, float radius){

  // compute the actual affect grid
  l1.x -= radius;
  l1.y -= radius;
  r1.x += radius;
  r1.y += radius;

  // If one rectangle is on left side of other
  if (l1.x > r2.x || l2.x > r1.x)
      return false;

  // If one rectangle is above other
  if (l1.y > r2.y || l2.y > r1.y)
      return false;

  return true;
}

//this is the function that would compute all the relative particles of myparticles and store in rP vector
void constructRelatedP(std::vector<Particle> &myParticles, std::vector<Particle> &relativeParticles, 
        float* allLimits, Vec2 target_L, Vec2 target_R, int nproc, float radius, int pid, int *recvCounts){
  
  // copy over myParticles to the relativeParticles 

  relativeParticles.resize(myParticles.size());
  for(int i = 0; i < (int)myParticles.size(); i++){
    relativeParticles[i] = myParticles[i];
  }
  int mySize = (int)myParticles.size();
  Vec2 comp_L, comp_R;
  bool overlap[nproc];

  MPI_Request request[nproc];

  int msgcount = 0;

  
  for(int j = 0; j < nproc; j ++){
    if(j != pid){
      //copy over the limits from the allLimits array
      comp_L.x = allLimits[j * 4];
      comp_L.y = allLimits[j * 4 + 1];
      comp_R.x = allLimits[j * 4 + 2];
      comp_R.y = allLimits[j * 4 + 3];
      
      if (forceCheck(target_L, target_R, comp_L, comp_R, radius)){
        //send the particles 
        
        MPI_Isend(myParticles.data(), mySize * sizeof(Particle), MPI_BYTE, j, pid, MPI_COMM_WORLD, &request[msgcount]);
        
        overlap[j] = true; 
        msgcount += 1;
      } else{
        overlap[j] = false;
      }
    }
    else{
      overlap[j] = false;
    }
  }
  
  int offset2 = relativeParticles.size(); //number of particles
  for(int k = 0; k < nproc; k++){
    if(overlap[k]){
        int sendSize = recvCounts[k]; //total bytes in process k
        
        relativeParticles.resize(relativeParticles.size() + sendSize / sizeof(Particle));
        
        MPI_Recv(relativeParticles.data() + offset2, sendSize, MPI_BYTE, k, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //store the particles to myparticles with the offset
        offset2 += sendSize/sizeof(Particle);
        
        //msgcount += 1;
    }
  } 
  MPI_Status status[msgcount];
  MPI_Waitall(msgcount, request, status);
 
}

bool sortParticleId(Particle a, Particle b){
  return a.id < b.id;
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
  float *allLimits = (float*)malloc(nproc * 4*sizeof(float));
  //std::vector<float> allLimits(4*nproc);

  //allLimits.resize(nproc * 4);
  // std::cerr<<__LINE__<<"      "<<pid<<std::endl;
  if (pid == 0) { //load all particles to master 
    loadFromFile(options.inputFile, totalParticles); 
  }

  int displs[nproc];
  int recvCounts[nproc];
  int particledispls[nproc];
  int particlerecvCounts[nproc];


  Vec2 pmin(1e30f, 1e30f);
  Vec2 pmax(-1e30f, -1e30f);

  if(pid == 0){
    for (auto &p : totalParticles) {
      pmin.x = fminf(pmin.x, p.position.x);
      pmin.y = fminf(pmin.y, p.position.y);
      pmax.x = fmaxf(pmax.x, p.position.x);
      pmax.y = fmaxf(pmax.y, p.position.y);
    }
    pmin.x = (pmin.x < 0) ? pmin.x -1 : pmin.x + 1;
    pmin.y = (pmin.y < 0) ? pmin.y -1 : pmin.y + 1;
    pmax.x = (pmax.x < 0) ? pmax.x -1 : pmax.x + 1;
    pmax.y = (pmax.y < 0) ? pmax.y -1 : pmax.y + 1;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  assignAll(totalParticles, myParticles, pmin, pmax, nproc, pid, gridSize);
  newParticles.resize(myParticles.size());

  recvCounts[pid] = myParticles.size() * sizeof(Particle);
  int local_recv = recvCounts[pid];
  MPI_Allgather(&local_recv, 1, MPI_INT, &recvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  for(int p = 0; p < nproc; p++){
    displs[p] = (p == 0) ? 0 : recvCounts[p-1] + displs[p-1];
    //particlerecvCounts are boundaries minx, miny, maxx, maxy
    particlerecvCounts[p] = 4 * sizeof(float);
    particledispls[p] = 4 * sizeof(float) * p;
  }

  Timer totalSimulationTimer;

  for (int i = 0; i < options.numIterations; i++) {

    if(i % resizeIter == 0 && i != 0){
 
        recvCounts[pid] = myParticles.size() * sizeof(Particle);
        local_recv = recvCounts[pid];
       
        MPI_Allgather(&local_recv, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);
        
        for(int p = 0; p < nproc; p++){
            displs[p] = (p == 0) ? 0 : recvCounts[p-1] + displs[p-1];
        }//modify the revcount and displs again
        
        totalParticles.resize((displs[nproc-1] + recvCounts[nproc -1])/sizeof(Particle)); 
        MPI_Gatherv(myParticles.data(), myParticles.size() * sizeof(Particle), MPI_BYTE,
            totalParticles.data(), recvCounts, displs, MPI_BYTE,0, MPI_COMM_WORLD);
        std::vector<Particle> myParticles;
        
        for (auto &p : totalParticles) {
          pmin.x = fminf(pmin.x, p.position.x);
          pmin.y = fminf(pmin.y, p.position.y);
          pmax.x = fmaxf(pmax.x, p.position.x);
          pmax.y = fmaxf(pmax.y, p.position.y);
          
        }
        pmin.x = (pmin.x < 0) ? pmin.x -1 : pmin.x + 1;
        pmin.y = (pmin.y < 0) ? pmin.y -1 : pmin.y + 1;
        pmax.x = (pmax.x < 0) ? pmax.x -1 : pmax.x + 1;
        pmax.y = (pmax.y < 0) ? pmax.y -1 : pmax.y + 1;

        assignAll(totalParticles, myParticles, pmin, pmax, nproc, pid, gridSize);
        newParticles.resize(myParticles.size());

        
    } 
    Vec2 myMin, myMax;
    

    for (auto &p : myParticles) {
      myMin.x = fminf(myMin.x, p.position.x);
      myMin.y = fminf(myMin.y, p.position.y);
      myMax.x = fmaxf(myMax.x, p.position.x);
      myMax.y = fmaxf(myMax.y, p.position.y);
    }

    std::vector<float> boundary{myMin.x, myMin.y, myMax.x, myMax.y};
    MPI_Allgather(boundary.data(), 4, MPI_FLOAT, allLimits, 4*nproc,MPI_FLOAT, MPI_COMM_WORLD);

    std::vector<Particle> relativeParticles;
    constructRelatedP(myParticles, relativeParticles, allLimits, myMin, myMax, nproc, stepParams.cullRadius, pid, recvCounts);

    QuadTree tree;
    QuadTree::buildQuadTree(relativeParticles, tree);
    
    newParticles.resize(myParticles.size());
    simulateStep(tree, myParticles, newParticles, stepParams);
    
    myParticles.swap(newParticles);
    
  }
 
  recvCounts[pid] = myParticles.size() * sizeof(Particle);
  local_recv = recvCounts[pid];
 
  MPI_Allgather(&local_recv, 1, MPI_INT, &recvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  for(int p = 0; p < nproc; p++){
    displs[p] = (p == 0) ? 0 : recvCounts[p-1] + displs[p-1];
  }//modify the revcount and displs again

  totalParticles.resize((displs[nproc-1] + recvCounts[nproc-1])/sizeof(Particle));
  
  MPI_Gatherv(myParticles.data(), myParticles.size() * sizeof(Particle), MPI_BYTE,
            totalParticles.data(), recvCounts, displs, MPI_BYTE,0, MPI_COMM_WORLD);
  double totalSimulationTime = totalSimulationTimer.elapsed();
 
  if (pid == 0) {
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    std::sort(totalParticles.begin(), totalParticles.end(), sortParticleId);
    saveToFile(options.outputFile, totalParticles);
  }
  
  MPI_Finalize();
 
}

