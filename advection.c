/* 
	@author Jack Clark

	Simple program to simulate 1D advection using the finite volume approach, with naive averaging at cell boundaries. 
	
	Compile with g++ -O3 advection.c -o advection
*/


#include <fstream>
#include <sstream>
#include <math.h>

#define NUM_CELLS_X 100
#define X_VELOCITY 1
#define TIMESTEP 0.1
#define NUM_TIMESTEPS 500
#define DELTA_X 1
#define PLOT_FREQUENCY 1

double flux[NUM_CELLS_X];
double q[NUM_CELLS_X];

void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<NUM_CELLS_X; i++) {
    out << i
        << ","
        << q[i]
        << std::endl;
  }
}

void printResult(int timestep) {

	printf("Timestep %d \n", timestep);
	printf("Q: ");
	for(int j=0; j<NUM_CELLS_X; j++) {
		printf("%f ", q[j]);
	}
	printf("\n");
	printf("Flux: ");
	for(int j=0; j<NUM_CELLS_X; j++) {
		printf("%f ", flux[j]);
	}
	printf("\n");
}

void setup() {
	for(int i=0; i<NUM_CELLS_X; i++) {
		flux[i] = 0.0;
		if(i > 0 && i < 6) {
			q[i] = 5.0;
		} else {
			q[i] = 0.0;
		}
	}
}

void reconstruction() {
	for(int i=1; i<NUM_CELLS_X; i++) { // Note that flux is not computed for the boundaries
		flux[i] = (q[i] + q[i-1]) / 2 * X_VELOCITY;
	}
}

void update_cells() {
	double temp;
	for(int i=0; i<NUM_CELLS_X-1; i++) {
		temp = q[i] + (TIMESTEP/DELTA_X) * (flux[i] - flux[i+1]);
		if(temp < 0) {
			q[i] = 0;
		} else {
			q[i] = temp;
		}
	}
}


int main() {
	setup();
	printCSVFile(0);
	printResult(0);

	for(int i=0; i<NUM_TIMESTEPS; i++) {
		reconstruction();
		update_cells();

		printResult(i+1);

		if (i%PLOT_FREQUENCY==0) {
      		printCSVFile(i/PLOT_FREQUENCY+1); // Please switch off all IO if you do performance tests.
    	}
	}

	return 0;
}