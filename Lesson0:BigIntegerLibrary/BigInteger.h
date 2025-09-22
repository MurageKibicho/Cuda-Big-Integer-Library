#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define BASE ((uint64_t)1 << 32)
#define MAX_LIMBS_256 8  // 256 bits
#define MAX_LIMBS_512 16 // 512 bits

struct kibicho_big_256_int_struct
{
	uint32_t digits[MAX_LIMBS_256];
	int sign; //0 is positive, anything else is negative
};
struct kibicho_big_512_int_struct
{
	uint32_t digits[MAX_LIMBS_512];
	int sign;
};

typedef struct kibicho_big_512_int_struct Kib512[1];
typedef struct kibicho_big_256_int_struct Kib256[1];


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

void Kib512AllOnes(Kib512 x)
{
	for(int i = 0; i < MAX_LIMBS_512; ++i){x[0].digits[i] = 0xFFFFFFFFu;}
	x[0].sign = 0;
}

void Kib256AllOnes(Kib256 x)
{
	for(int i = 0; i < MAX_LIMBS_256; ++i){x[0].digits[i] = 0xFFFFFFFFu;}
	x[0].sign = 0;
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

int Kib256CompareMagnitude(Kib256 a, Kib256 b)
{
	for(int i = MAX_LIMBS_256 - 1; i >= 0; i--)
	{
		if(a[0].digits[i] < b[0].digits[i]) return -1;
		if(a[0].digits[i] > b[0].digits[i]) return 1;
	}
	return 0;
}

int Kib512CompareMagnitude(Kib512 a, Kib512 b)
{
	for(int i = MAX_LIMBS_512 - 1; i >= 0; i--)
	{
		if(a[0].digits[i] < b[0].digits[i]) return -1;
		if(a[0].digits[i] > b[0].digits[i]) return 1;
	}
	return 0;
}

void Kib512AddMagnitudes(Kib512 res, Kib512 a, Kib512 b) 
{
	uint64_t carry = 0;
	for(int i = 0; i < MAX_LIMBS_512; i++)
	{
		uint64_t sum = (uint64_t)a[0].digits[i] + b[0].digits[i] + carry;
		res[0].digits[i] = (uint32_t)(sum & 0xFFFFFFFFULL);
		carry = sum >> 32;
	}
}

void Kib256AddMagnitudes(Kib256 res, Kib256 a, Kib256 b) 
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
void Kib512SubtractMagnitudes(Kib512 res, Kib512 a, Kib512 b) 
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
void Kib256SubtractMagnitudes(Kib256 res, Kib256 a, Kib256 b) 
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

void Kib512Add(Kib512 res, Kib512 a, Kib512 b) 
{
	if(a[0].sign == b[0].sign)
	{
		Kib512AddMagnitudes(res, a, b);
		res[0].sign = a[0].sign;
	}
	else
	{
		int cmp = Kib512CompareMagnitude(a, b);
		if(cmp == 0)
		{
			Kib512Zero(res);
		}
		else if
		(cmp > 0)
		{
			Kib512SubtractMagnitudes(res, a, b);
			res[0].sign = a[0].sign;
		}
		else
		{
			Kib512SubtractMagnitudes(res, b, a);
			res[0].sign = b[0].sign;
		}
	}
}


void Kib256Add(Kib256 res, Kib256 a, Kib256 b) 
{
	if(a[0].sign == b[0].sign)
	{
		Kib256AddMagnitudes(res, a, b);
		res[0].sign = a[0].sign;
	}
	else
	{
		int cmp = Kib256CompareMagnitude(a, b);
		if(cmp == 0)
		{
			Kib256Zero(res);
		}
		else if
		(cmp > 0)
		{
			Kib256SubtractMagnitudes(res, a, b);
			res[0].sign = a[0].sign;
		}
		else
		{
			Kib256SubtractMagnitudes(res, b, a);
			res[0].sign = b[0].sign;
		}
	}
}

void Kib256Mul(Kib512 res, Kib256 a, Kib256 b)
{
	Kib512Zero(res);
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




