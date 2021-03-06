There are many different versions of the code included in this folder. 

advection:
This folder contains both the CPU sequential and GPU versions of the 2D advection simulations.
Instructions for compiling and running these programs are included in their respective source files.

exahype:
This contains the version of the ExaHyPE code that my project worked with. It is important to use this version instead of other versions of ExaHyPE.

seq_cuda, stream_cuda, stream_cuda_refactored:
These folders contain the different CUDA implementations of the patch level godunov solver. seq_cuda contains a non-stream version of the CUDA code which is roughly ~8% slower than the stream based versions. stream_cuda contains the stream version of the CUDA code. stream_cuda_refactored contains stream based cuda code that has been refactored to make it more user friendly, by separating out functions that require user definitions into a separate file, as well as separating out utility functions. This is the "final" version of the code.

To use one of these versions of the code, please replace all the contents of the exahype/ExaHyPE/kernels/finitevolumes/godunov/c/2d folder with the version of the code that you wish to use, with the exception of the makefile. The makefile should replace the makefile in the exahype/ExaHyPE directory. You should then be able to build the project.

NOTES:

1. A lot of this code can be simplified in newer versions of ExaHyPE by accessing solverType class information rather than passing information by parameter.

2. For some reason, with certain patch sizes the code gets a Finite volume solver time step size harmed CFL condition warning. I'm not sure exactly what the reason for this is.

3. The main findings from the project were that the current ExaHyPE threading model made it difficult to provide the GPU with enough work to do in order for the GPU version to be efficient. There was a suggestion that in order to make it work, a GPU manager could be created to intercept GPU offload requests and then forward them to the GPU. This has the advantage of abstracting away the details of GPU programming from most of the CPU code and also providing the possibility to make the offloading process efficient through keeping state on the GPUs across requests. This seemed particularly appealing for the multi-GPU case.


