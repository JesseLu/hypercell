all:
	nvcc -I/home/hypercell/NVIDIA_GPU_Computing_SDK/C/common/inc/ -lhdf5 -lhdf5_hl fdtd.cu -o hypercell.out

