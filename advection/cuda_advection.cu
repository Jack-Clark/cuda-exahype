/* 
	@author Jack Clark

	Simple program to simulate 2D advection using the finite volume approach, with naive averaging at cell boundaries. 
	
	Compile with nvcc -O3 advection.cu -o gpu_advection
*/


#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <cuda.h>


#define NUM_CELLS_X 40
#define NUM_CELLS_Y 40
#define TIMESTEP 0.001
#define NUM_TIMESTEPS 100000
#define DELTA_X 1
#define DELTA_Y 1
#define PLOT_FREQUENCY 100

// CPU data
double h_flux_x[NUM_CELLS_X*NUM_CELLS_Y];
double h_flux_y[NUM_CELLS_X*NUM_CELLS_Y];
double h_q[NUM_CELLS_X*NUM_CELLS_Y];
double h_velocities[2];
double h_max_velocity;

// GPU data
double * d_flux_x;
double * d_flux_y;
double * d_q;
double * d_max_velocity;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int getIndex(int row, int col) {
	return row * NUM_CELLS_X + col;
}

void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<NUM_CELLS_Y; i++) {
  	for(int j=0; j<NUM_CELLS_X; j++) {
	    out << i
	        << ","
	        << j
	 	    << ","
	        << h_q[i*NUM_CELLS_X+j]
	        << std::endl;
	}
  }
}

void printResult(int timestep) {

	printf("\nTimestep %d \n", timestep);
	printf("Q: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", h_q[getIndex(i,j)]);
  		}
  	}
	printf("\n");
	printf("Flux_X: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", h_flux_x[getIndex(i,j)]);
  		}
  	}
	printf("\n");
	printf("Flux_Y: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", h_flux_y[getIndex(i,j)]);
  		}
  	}
	printf("\n");
}

void setup() {
	for (int i=1; i<4; i++) {
  		for(int j=1; j<4; j++) {
  			int index = getIndex(i,j);
  			h_q[index] = 5.0;
  		}
  	}
	h_velocities[0] = 0.5;
	h_velocities[1] = 0.5;
	h_max_velocity = 0.0;
	int i;
	for(i = 0; i<sizeof(h_velocities)/sizeof(double); i++) {
		h_max_velocity += h_velocities[i] * h_velocities[i];
	}
	h_max_velocity = sqrt(h_max_velocity);
}

void reconstruction() {
	for(int i=0; i<NUM_CELLS_Y; i++) {
		for(int j=0; j<NUM_CELLS_X; j++) {
			int index = getIndex(i,j);
			if(j != 0) {
				h_flux_x[index] = ((h_q[index] + h_q[index-1]) / 2) - (h_max_velocity/2 * (h_q[index] - h_q[index-1]));
			}
			if(i != 0) {
				h_flux_y[index] = ((h_q[index] + h_q[index-NUM_CELLS_X]) / 2) - (h_max_velocity/2 * (h_q[index] - h_q[index-NUM_CELLS_X]));
			}
		}
	}
}

__global__ void reconstruction_kernel(const double * const d_q, double * d_flux_x, double * d_flux_y, const double * const d_max_velocity) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= NUM_CELLS_Y || col >= NUM_CELLS_X)
		return;

	int self = row * NUM_CELLS_X + col;

	if(col != 0) {
		d_flux_x[self] = ((d_q[self] + d_q[self-1]) / 2) - (*d_max_velocity/2 * (d_q[self] - d_q[self-1]));
	}

	if(row != 0) {
		d_flux_y[self] = ((d_q[self] + d_q[self-NUM_CELLS_X]) / 2) - (*d_max_velocity/2 * (d_q[self] - d_q[self-NUM_CELLS_X]));
	}
}

__global__ void update_kernel(double * d_q, const double * const d_flux_x, const double * const d_flux_y) {  
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= NUM_CELLS_Y || col >= NUM_CELLS_X)
		return;

	int self = row * NUM_CELLS_X + col;

	if(self > ((NUM_CELLS_X-1) * (NUM_CELLS_Y-1) - 1))
		return;

	double temp = d_q[self] + ((TIMESTEP/DELTA_X) * (d_flux_x[self] - d_flux_x[self+1]));

	if(row < NUM_CELLS_X-1)
		temp += (TIMESTEP/DELTA_Y) * (d_flux_y[self] - d_flux_y[self+NUM_CELLS_X]);

	// 0 limiter to remove numerical artefacts 
	if(temp <= 0.0) {
		d_q[self] = 0.0;
	} else {
		d_q[self] = temp;
	}

}

void update_cells() {
	double temp;
	for(int i=0; i<NUM_CELLS_Y-1; i++) {
		for(int j=0; j<NUM_CELLS_X-1; j++) {
			int index = getIndex(i,j);
			assert(index >= 0);
			assert(index < NUM_CELLS_X * NUM_CELLS_Y);
			temp = h_q[index] + ((TIMESTEP/DELTA_X) * (h_flux_x[index] - h_flux_x[index+1])) + ((TIMESTEP/DELTA_Y) * (h_flux_y[index] - h_flux_y[index+NUM_CELLS_X]));
			if(temp < 0) {
				h_q[index] = 0;
			} else {
				h_q[index] = temp;
			}
		}
	}
}


int main() {
	setup();
	printCSVFile(0);
	//printResult(0);

	// Allocate memory for GPU data - Naive assumption that our GPU can fit all our data in memory at once
	gpuErrchk(cudaMalloc((void **) &d_flux_x, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y));
	gpuErrchk(cudaMalloc((void **) &d_flux_y, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y));
	gpuErrchk(cudaMalloc((void **) &d_q, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y));
	gpuErrchk(cudaMalloc((void **) &d_max_velocity, sizeof(double)));

	gpuErrchk(cudaMemset(d_flux_x, 0, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y));
	gpuErrchk(cudaMemset(d_flux_y, 0, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y));
	gpuErrchk(cudaMemcpy(d_max_velocity, &h_max_velocity, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_q, h_q, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y, cudaMemcpyHostToDevice));
	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	dim3 block(32, 32);

	for(int i=0; i<NUM_TIMESTEPS; i++) {
		reconstruction_kernel<<<3, block>>>(d_q, d_flux_x, d_flux_y, d_max_velocity);
		gpuErrchk( cudaPeekAtLastError() );
		update_kernel<<<3, block>>>(d_q, d_flux_x, d_flux_y);
		gpuErrchk( cudaPeekAtLastError() );

		// Copy GPU data back to GPU for recording - TODO: Improve this so that we only copy data back when we need it.
		gpuErrchk(cudaMemcpy(h_q, d_q, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(h_flux_x, d_flux_x, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(h_flux_y, d_flux_y, sizeof(double) * NUM_CELLS_X * NUM_CELLS_Y, cudaMemcpyDeviceToHost));

		if (i%PLOT_FREQUENCY==0) {
      		printCSVFile(i/PLOT_FREQUENCY+1); // Please switch off all IO if you do performance tests.
      		//printResult(i+1);
    	}
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	printf ("Time for the program: %f ms\n", time);

	cudaFree(d_flux_x);
	cudaFree(d_flux_y);
	cudaFree(d_q);
	cudaFree(d_max_velocity);

	return 0;
}