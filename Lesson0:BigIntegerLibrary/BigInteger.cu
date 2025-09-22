#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>
#define BASE ((uint64_t)1 << 32)
#define MAX_LIMBS_256 8  // 256 bits
#define MAX_LIMBS_512 16 // 512 bits

struct kibicho_big_256_int_struct
{
	uint32_t digits[MAX_LIMBS_256];
	int sign; //0 is positive, 1 is negative
};
struct kibicho_big_512_int_struct
{
	uint32_t digits[MAX_LIMBS_512];
	int sign;
};

typedef struct kibicho_big_512_int_struct Kib512[1];
typedef struct kibicho_big_256_int_struct Kib256[1];

//Device functions
//memset is not on device. Replace with for loop
__device__ void d_Kib512Zero(Kib512 x)
{
	for(int i = 0; i < MAX_LIMBS_512; ++i)
	{
		x[0].digits[i] = 0;
	}
	x[0].sign = 0;
}

__device__ void d_Kib256Zero(Kib256 x)
{
	for(int i = 0; i < MAX_LIMBS_256; ++i)
	{
		x[0].digits[i] = 0;
	}
	x[0].sign = 0;
}

__device__ int d_Kib256CompareMagnitude(Kib256 a, Kib256 b)
{
	for(int i = MAX_LIMBS_256 - 1; i >= 0; i--)
	{
		if(a[0].digits[i] < b[0].digits[i]) return -1;
		if(a[0].digits[i] > b[0].digits[i]) return 1;
	}
	return 0;
}

__device__ int d_Kib512CompareMagnitude(Kib512 a, Kib512 b)
{
	for(int i = MAX_LIMBS_512 - 1; i >= 0; i--)
	{
		if(a[0].digits[i] < b[0].digits[i]) return -1;
		if(a[0].digits[i] > b[0].digits[i]) return 1;
	}
	return 0;
}

__device__ void d_Kib512AddMagnitudes(Kib512 res, Kib512 a, Kib512 b) 
{
	uint64_t carry = 0;
	for(int i = 0; i < MAX_LIMBS_512; i++)
	{
		uint64_t sum = (uint64_t)a[0].digits[i] + b[0].digits[i] + carry;
		res[0].digits[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
		carry = sum >> 32;
	}
}

__device__ void d_Kib256AddMagnitudes(Kib256 res, Kib256 a, Kib256 b) 
{
	uint64_t carry = 0;
	for(int i = 0; i < MAX_LIMBS_256; i++)
	{
		uint64_t sum = (uint64_t)a[0].digits[i] + b[0].digits[i] + carry;
		res[0].digits[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
		carry = sum >> 32;
	}
}

// Subtract magnitudes (assumes a >= b)
__device__ void d_Kib512SubtractMagnitudes(Kib512 res, Kib512 a, Kib512 b) 
{
	int64_t borrow = 0;
	for(int i = 0; i < MAX_LIMBS_512; i++)
	{
		int64_t diff = (int64_t)a[0].digits[i] - b[0].digits[i] - borrow;
		if(diff < 0)
		{
			diff += BASE;
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		res[0].digits[i] = (uint32_t)diff;
	}
}

// Subtract magnitudes (assumes a >= b)
__device__ void d_Kib256SubtractMagnitudes(Kib256 res, Kib256 a, Kib256 b) 
{
	int64_t borrow = 0;
	for(int i = 0; i < MAX_LIMBS_256; i++)
	{
		int64_t diff = (int64_t)a[0].digits[i] - b[0].digits[i] - borrow;
		if(diff < 0)
		{
			diff += BASE;
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		res[0].digits[i] = (uint32_t)diff;
	}
}

__device__ void d_Kib512Add(Kib512 res, Kib512 a, Kib512 b) 
{
	if(a[0].sign == b[0].sign)
	{
		d_Kib512AddMagnitudes(res, a, b);
		res[0].sign = a[0].sign;
	}
	else
	{
		int cmp = d_Kib512CompareMagnitude(a, b);
		if(cmp == 0)
		{
			d_Kib512Zero(res);
		}
		else if(cmp > 0)
		{
			d_Kib512SubtractMagnitudes(res, a, b);
			res[0].sign = a[0].sign;
		}
		else
		{
			d_Kib512SubtractMagnitudes(res, b, a);
			res[0].sign = b[0].sign;
		}
	}
}

__device__ void d_Kib256Add(Kib256 res, Kib256 a, Kib256 b) 
{
	if(a[0].sign == b[0].sign)
	{
		d_Kib256AddMagnitudes(res, a, b);
		res[0].sign = a[0].sign;
	}
	else
	{
		int cmp = d_Kib256CompareMagnitude(a, b);
		if(cmp == 0)
		{
			d_Kib256Zero(res);
		}
		else if(cmp > 0)
		{
			d_Kib256SubtractMagnitudes(res, a, b);
			res[0].sign = a[0].sign;
		}
		else
		{
			d_Kib256SubtractMagnitudes(res, b, a);
			res[0].sign = b[0].sign;
		}
	}
}

__device__ void d_Kib256Mul(Kib512 res, Kib256 a, Kib256 b)
{
	d_Kib512Zero(res);
	for(int i = 0; i < MAX_LIMBS_256; i++)
	{
		uint64_t carry = 0;
		for(int j = 0; j < MAX_LIMBS_256; j++)
		{
			uint64_t sum = (uint64_t)a[0].digits[i] * b[0].digits[j] + res[0].digits[i + j] + carry;
			res[0].digits[i + j] = (uint32_t)(sum & 0xFFFFFFFF);
			carry = sum >> 32;
		}
		res[0].digits[i + MAX_LIMBS_256] = (uint32_t)carry;
	}
	//Added sign by XOR
	res[0].sign = a[0].sign ^ b[0].sign;
}

//Kernel
__global__ void kib256MulKernel(Kib512* results, Kib256* inputsA, Kib256* inputsB, int count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < count)
	{
		d_Kib256Mul(results[idx], inputsA[idx], inputsB[idx]);
	}
}

//Host functions
cudaError_t kib256MulBatch(Kib512* d_results, Kib256* d_inputsA, Kib256* d_inputsB, int count)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (count + threadsPerBlock - 1) / threadsPerBlock;

	kib256MulKernel<<<blocksPerGrid, threadsPerBlock>>>(d_results, d_inputsA, d_inputsB, count);
	return cudaGetLastError();
}
cudaError_t allocateKib256Array(Kib256** d_array, int count)
{
	return cudaMalloc(d_array, count * sizeof(Kib256));
}

cudaError_t allocateKib512Array(Kib512** d_array, int count)
{
	return cudaMalloc(d_array, count * sizeof(Kib512));
}

cudaError_t copyToDeviceKib256(Kib256* d_dest, Kib256* h_src, int count)
{
	return cudaMemcpy(d_dest, h_src, count * sizeof(Kib256), cudaMemcpyHostToDevice);
}

cudaError_t copyFromDeviceKib512(Kib512* h_dest, Kib512* d_src, int count)
{
	return cudaMemcpy(h_dest, d_src, count * sizeof(Kib512), cudaMemcpyDeviceToHost);
}

//CPU functions
//Rand DNE on CUDA
void Kib512Zero(Kib512 x)
{
	memset(x[0].digits, 0, sizeof(x[0].digits));
	x[0].sign = 0;
}

void Kib256Zero(Kib256 x)
{
	memset(x[0].digits, 0, sizeof(x[0].digits));
	x[0].sign = 0;
}

void Kib512PrintHex(Kib512 a)
{
	int started = 0;
	for (int i = MAX_LIMBS_512 - 1; i >= 0; i--)
	{
		if (!started && a[0].digits[i] == 0) continue;
		if (!started){printf("%x", a[0].digits[i]);started = 1;}
		else{printf("%08x", a[0].digits[i]);}
	}
	if(!started)printf("0");
}

void Kib256PrintHex(Kib256 a)
{
	int started = 0;
	for (int i = MAX_LIMBS_256 - 1; i >= 0; i--)
	{
		if (!started && a[0].digits[i] == 0) continue;
		if (!started){printf("%x", a[0].digits[i]);started = 1;}
		else{printf("%08x", a[0].digits[i]);}
	}
	if(!started)printf("0");
}


void Kib512Random(Kib512 a)
{
	for(int i = 0; i < MAX_LIMBS_512; i++)
	{
		a[0].digits[i] = (uint32_t)rand();
	}
	a[0].sign = rand() % 2;
}

void Kib256Random(Kib256 a)
{
	for(int i = 0; i < MAX_LIMBS_256; i++)
	{
		a[0].digits[i] = (uint32_t)rand();
	}
	a[0].sign = 0;
}

int main()
{
	const int NUM_OPERATIONS = 1000000;
	cudaEvent_t start, stop;
	float milliseconds = 0;

	// Create CUDA events for timing
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// Allocate host memory
	Kib256* h_inputsA = (Kib256*)malloc(NUM_OPERATIONS * sizeof(Kib256));
	Kib256* h_inputsB = (Kib256*)malloc(NUM_OPERATIONS * sizeof(Kib256));
	Kib512* h_results = (Kib512*)malloc(NUM_OPERATIONS * sizeof(Kib512));

	// Initialize with random data
	for(int i = 0; i < NUM_OPERATIONS; i++)
	{
		Kib256Random(h_inputsA[i]);
		Kib256Random(h_inputsB[i]);
	}

	// Allocate device memory
	Kib256* d_inputsA, *d_inputsB;
	Kib512* d_results;

	allocateKib256Array(&d_inputsA, NUM_OPERATIONS);
	allocateKib256Array(&d_inputsB, NUM_OPERATIONS);
	allocateKib512Array(&d_results, NUM_OPERATIONS);

	// Copy data to device
	copyToDeviceKib256(d_inputsA, h_inputsA, NUM_OPERATIONS);
	copyToDeviceKib256(d_inputsB, h_inputsB, NUM_OPERATIONS);

	// Start timing (includes kernel launch + GPU execution)
	cudaEventRecord(start);

	// Launch kernel
	kib256MulBatch(d_results, d_inputsA, d_inputsB, NUM_OPERATIONS);

	// Wait for kernel to complete and stop timing
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	// Calculate elapsed time
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("GPU execution time: %.3f ms\n", milliseconds);
	printf("Operations per second: %.0f\n", NUM_OPERATIONS / (milliseconds / 1000.0f));

	// Copy results back
	copyFromDeviceKib512(h_results, d_results, NUM_OPERATIONS);

	// Cleanup
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(d_inputsA);
	cudaFree(d_inputsB);
	cudaFree(d_results);
	free(h_inputsA);
	free(h_inputsB);
	free(h_results);

	return 0;
}


