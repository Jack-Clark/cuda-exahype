/* 
	@author Jack Clark

	Simple program to simulate 2D advection using the finite volume approach, with naive averaging at cell boundaries. 
	
	Compile with g++ -O3 advection.c -o advection
*/


#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>

#define NUM_CELLS_X 48
#define NUM_CELLS_Y 48
#define TIMESTEP 0.001
#define NUM_TIMESTEPS 100000
#define DELTA_X 1
#define DELTA_Y 1
#define PLOT_FREQUENCY 100

double flux_x[NUM_CELLS_X*NUM_CELLS_Y];
double flux_y[NUM_CELLS_X*NUM_CELLS_Y];
double q[NUM_CELLS_X*NUM_CELLS_Y];
double velocities[2];
double max_velocity;

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
	        << q[i*NUM_CELLS_X+j]
	        << std::endl;
	}
  }
}

void printResult(int timestep) {

	printf("Timestep %d \n", timestep);
	printf("Q: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", q[getIndex(i,j)]);
  		}
  	}
	printf("\n");
	printf("Flux_X: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", flux_x[getIndex(i,j)]);
  		}
  	}
	printf("\n");
	printf("Flux_Y: ");
	for (int i=0; i<NUM_CELLS_Y; i++) {
		printf("\n");
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			printf("%f ", flux_y[getIndex(i,j)]);
  		}
  	}
	printf("\n");
}

void setup() {
	for (int i=0; i<NUM_CELLS_Y; i++) {
  		for(int j=0; j<NUM_CELLS_X; j++) {
  			int index = getIndex(i,j);
  			assert(index >= 0);
			assert(index < NUM_CELLS_X * NUM_CELLS_Y);
  			if(index == 2*NUM_CELLS_X+2 || index == 2*NUM_CELLS_X+3 || index == 3*NUM_CELLS_X+2 || index == 3*NUM_CELLS_X+3) {
  				q[index] = 2.0;
  			} else {
  				q[index] = 0.0;
  			}
  		}
  	}
	velocities[0] = 0.5;
	velocities[1] = 0.5;
	max_velocity = 0.0;
	int i;
	for(i = 0; i<sizeof(velocities)/sizeof(double); i++) {
		max_velocity += velocities[i] * velocities[i];
	}
	max_velocity = sqrt(max_velocity);
}

void reconstruction() {
	for(int i=0; i<NUM_CELLS_Y; i++) {
		for(int j=0; j<NUM_CELLS_X; j++) {
			int index = getIndex(i,j);
			assert(index >= 0);
			assert(index < NUM_CELLS_X * NUM_CELLS_Y);
			if(j != 0) {
				flux_x[index] = ((q[index] + q[index-1]) / 2) - (max_velocity/2 * (q[index] - q[index-1]));
			}
			if(i != 0) {
				flux_y[index] = ((q[index] + q[index-NUM_CELLS_X]) / 2) - (max_velocity/2 * (q[index] - q[index-NUM_CELLS_X]));
			}
		}
	}
}

void update_cells() {
	double temp;
	for(int i=0; i<NUM_CELLS_Y-1; i++) {
		for(int j=0; j<NUM_CELLS_X-1; j++) {
			int index = getIndex(i,j);
			assert(index >= 0);
			assert(index < NUM_CELLS_X * NUM_CELLS_Y);
			temp = q[index] + ((TIMESTEP/DELTA_X) * (flux_x[index] - flux_x[index+1])) + ((TIMESTEP/DELTA_Y) * (flux_y[index] - flux_y[index+NUM_CELLS_X]));
			if(temp < 0) {
				q[index] = 0;
			} else {
				q[index] = temp;
			}
		}
	}
}


int main() {
	setup();
	printCSVFile(0);
	//printResult(0);
	clock_t begin = clock();

	for(int i=0; i<NUM_TIMESTEPS; i++) {
		reconstruction();
		update_cells();

		//printResult(i+1);

		if (i%PLOT_FREQUENCY==0) {
      		printCSVFile(i/PLOT_FREQUENCY+1); // Please switch off all IO if you do performance tests.
    	}
	}
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("CPU took %f \n", time_spent);

	return 0;
}