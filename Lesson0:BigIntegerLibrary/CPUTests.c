#include "BigInteger.h"
#include <gmp.h>

//clear && gcc CPUTests.c -lm -lgmp -o m.o && ./m.o
void Kib256ToMPZ(mpz_t z,  Kib256 a)
{
	mpz_set_ui(z, 0);
	for(int i = MAX_LIMBS_256 - 1; i >= 0; i--)
	{
		mpz_mul_2exp(z, z, 32);
		mpz_add_ui(z, z, a[0].digits[i]);	
	}
	if(a[0].sign)mpz_neg(z, z);	
}

void Kib512ToMPZ(mpz_t z,  Kib512 a)
{
	mpz_set_ui(z, 0);
	for(int i = MAX_LIMBS_512 - 1; i >= 0; i--)
	{
		mpz_mul_2exp(z, z, 32);
		mpz_add_ui(z, z, a[0].digits[i]);	
	}
	if(a[0].sign)mpz_neg(z, z);	
}

int Kib256EqualToMPZ(Kib256 x, mpz_t z)
{
	mpz_t tmp; mpz_init(tmp);
	Kib256ToMPZ(tmp, x);
	int ok = (mpz_cmp(tmp, z) == 0);
	mpz_clear(tmp);
	return ok;
}

int Kib512EqualToMPZ(Kib512 x, mpz_t z)
{
	mpz_t tmp; mpz_init(tmp);
	Kib512ToMPZ(tmp, x);
	int ok = (mpz_cmp(tmp, z) == 0);
	mpz_clear(tmp);
	return ok;
}

void TestAddSubCorrectness()
{
	srand(65342);
	printf("=== Add/Sub Correctness checks ===\n");
	int testCount = 1000000;
	int pass = 1;
	Kib256 a, b, res;
	mpz_t za, zb, zr;
	mpz_inits(za, zb, zr, NULL);
	for(int i = 0; i < testCount; i++)
	{
		Kib256Random(a);
		Kib256Random(b);
		Kib256Add(res, a, b);
		
		Kib256ToMPZ(za, a);
		Kib256ToMPZ(zb, b);
		mpz_add(zr, za, zb);
		
		if(!Kib256EqualToMPZ(res, zr))
		{
			printf("Mismatch at test %d!\n", i);
			pass = 0;
			break;
		}
	}
	mpz_clears(za, zb, zr, NULL);
	printf(pass ? "All correctness tests passed!\n" : "Some tests failed!\n");
}

void TestAddSpeed()
{
	srand(65342);
	printf("=== Add/Sub Speed checks ===\n");
	int testCount = 1000000;
	int pass = 1;
	Kib256 a, b, res;
	mpz_t za, zb, zr;
	gmp_randstate_t mpzRandomState;
	gmp_randinit_default(mpzRandomState);
	
	mpz_inits(za, zb, zr, NULL);
	clock_t t1, t2;
	t1 = clock();
	for(int i = 0; i < testCount; i++)
	{
		Kib256Random(a);
		Kib256Random(b);
		Kib256Add(res, a, b);
	}
	t2 = clock();
	double kib256Speed = (double)(t2 - t1) / CLOCKS_PER_SEC;
	t1 = clock();
	for(int i = 0; i < testCount; i++)
	{
		mpz_urandomb(za,mpzRandomState, 256);
		mpz_urandomb(zb,mpzRandomState, 256);
		mpz_add(zr, za, zb);
	}
	t2 = clock();
	double mpzSpeed = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("Kib256 total time: %.6f sec\n", kib256Speed);
	printf("GMP total time:      %.6f sec\n", mpzSpeed);	
	mpz_clears(za, zb, zr, NULL);
	gmp_randclear(mpzRandomState);
}


void TestMulCorrectness()
{
	srand(65342);
	printf("=== Mul Correctness checks ===\n");
	int testCount = 1000000;
	int pass = 1;
	Kib256 a, b;
	Kib512 res;
	mpz_t za, zb, zr;
	mpz_inits(za, zb, zr, NULL);
	for(int i = 0; i < testCount; i++)
	{
		Kib256Random(a);
		Kib256Random(b);
		Kib256Mul(res, a, b);
		
		Kib256ToMPZ(za, a);
		Kib256ToMPZ(zb, b);
		mpz_mul(zr, za, zb);
		
		if(!Kib512EqualToMPZ(res, zr))
		{
			printf("Mismatch at test %d!\n", i);
			pass = 0;
			break;
		}
	}
	mpz_clears(za, zb, zr, NULL);
	printf(pass ? "All correctness tests passed!\n" : "Some tests failed!\n");
}

void TestMulSpeed()
{
	srand(65342);
	printf("=== Mul Speed checks ===\n");
	int testCount = 1000000;
	int pass = 1;
	Kib256 a, b;
	Kib512 res;
	mpz_t za, zb, zr;
	gmp_randstate_t mpzRandomState;
	gmp_randinit_default(mpzRandomState);
	
	mpz_inits(za, zb, zr, NULL);
	clock_t t1, t2;
	t1 = clock();
	for(int i = 0; i < testCount; i++)
	{
		Kib256Random(a);
		Kib256Random(b);
		Kib256Mul(res, a, b);
	}
	t2 = clock();
	double kib256Speed = (double)(t2 - t1) / CLOCKS_PER_SEC;
	t1 = clock();
	for(int i = 0; i < testCount; i++)
	{
		mpz_urandomb(za,mpzRandomState, 256);
		mpz_urandomb(zb,mpzRandomState, 256);
		mpz_mul(zr, za, zb);
	}
	t2 = clock();
	double mpzSpeed = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("Kib256 total time: %.6f sec\n", kib256Speed);
	printf("GMP total time:      %.6f sec\n", mpzSpeed);	
	mpz_clears(za, zb, zr, NULL);
	gmp_randclear(mpzRandomState);
}



int main()
{
	TestAddSubCorrectness();
	TestAddSpeed();
	TestMulCorrectness();
	TestMulSpeed();
	return 0;
}
