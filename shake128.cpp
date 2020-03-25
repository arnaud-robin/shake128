#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <bitset>
#include<iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

// Round constants for faster computation
const uint64_t keccak_RC[24] = {
        0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
        0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
        0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
        0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
        0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
        0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
        0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
        0x8000000000008080, 0x0000000080000001, 0x8000000080008008
};

int index(int x, int y)
{
	return 5*y +x;
}

uint64_t rot64L(uint64_t ac , int l)
{
	l = l%64;
	return (ac << l) | (ac >> (64-l));
}


uint64_t rot64R(uint64_t ac , int l)
{
	l = l%64;
	return (ac >> l) | (ac << (64-l));
}

void print_state(uint64_t* A)
{
	for(int i=0; i<5; i++)
	{
		for(int j=0; j<5; j++)
			printf("%#018" PRIx64 " ", A[index(j,i)]);
		printf("\n");
	}
}

void theta(uint64_t * A)
{
	uint64_t C[5];
	for(int x=0; x<5;x++)
	{
		C[x] = A[index(x,0)] ^  A[index(x,1)] ^ A[index(x,2)] ^ A[index(x,3)] ^ A[index(x,4)];
	}

	uint64_t D[5];
	for(int x=0; x<5;x++)
	{
		D[x] = C[(x+4)%5] ^ rot64L(C[(x+1)%5],1);
	}

	for(int x=0;x<5; x++)
		for(int y=0; y<5; y++)
			A[index(x,y)] ^= D[x];
}


void rho(uint64_t * A)
{
	uint64_t Ap[25] =  {0};
	Ap[0] = A[0];

	int x = 1; int y=0;
	for(int t =0; t<=23; t++)
	{
		Ap[index(x,y)] = rot64L(A[index(x,y)],(t+1)*(t+2)/2);
		
		int yp = (2*x + 3*y) %5;
		x = y;
		y = yp;
	}
	copy(Ap, Ap +25, A);
}


void pi(uint64_t * A)
{
	uint64_t Ap[25] =  {0};
	for(int x=0; x<5; x++)
		for(int y=0; y<5; y++)
			Ap[index(x,y)] = A[index((x+3*y)%5,x)];

	copy(Ap, Ap +25, A);
}

void chi(uint64_t * A)
{
	uint64_t Ap[25] =  {0};
	for(int x=0; x<5; x++)
		for(int y=0; y<5; y++)
			Ap[index(x,y)] = A[index(x,y)] ^ ((~A[index((x+1)%5,y)]) & A[index((x+2)%5,y)]) ;

	copy(Ap, Ap +25, A);
}

uint64_t rc(int t)
{
	t = t%255;
	if(t==0)
		return 1;

	bitset<9> R(1);

	for(int i=1 ; i<= t; i++)
	{
		R <<=1;
		R[0] = R[0] ^ R[8];
		R[4] = R[4] ^ R[8];
		R[5] = R[5] ^ R[8];
		R[6] = R[6] ^ R[8];
	}
	return R[0];
}

void iota(uint64_t * A, int ir)
{
	A[0] ^= keccak_RC[ir];
}

void rnd(uint64_t * A, int ir)
{
	theta(A);
	rho(A);
	pi(A);
	chi(A);
	iota(A,ir);
}

void Keccakp(int nr, vector<u_char> &S)
{
	int l =6;
	uint64_t A[25]  = {0};
	for(int i=0; i<25; i++)
	{
		for(int j=0; j<8; j++)
		{
			A[i] += ((u_int64_t)S[8*i+j])<<(8*j);
		}
	}

	for(int ir = 12 + 2*l - nr; ir< 12 +2*l; ir++)
	{
		rnd(A,ir);
	}
	
	u_char mask = 0xff;
	for(int i=0; i<25; i++)
	{
		for(int j=0; j<8; j++)
		{
			S[8*i + j] = (A[i]>>(8*j)) & mask;
		}
	}

}

void pad0star_1(int x, int m, vector<u_char> &out)
{
	m-=4; //We already added the prefix
	int j = ((-m-2)%x+x)%x;
	int l = (j-10)/8;
	for(int k=1; k<=l;k++)
	{
		out.push_back(0);
	}
	out.push_back(0x80);
	
} 

void sponge(vector<u_char>  &N, int d_in_bytes, int r, vector<u_char>  &out)
{
	pad0star_1(r, 8 * N.size(), N);

	int n = (8*N.size())/r;
	int c = 1600-r;

	vector<u_char>  S(200,0);

	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<r/8; j++)
			S[j] ^= N[i*r/8 + j] ;
		Keccakp(24,S);
	}

	vector<u_char> Z;

	while(Z.size() < d_in_bytes)
	{
		Z.insert(Z.end(), S.begin(), S.begin() + r/8);
		Keccakp(24,S);
	}

	vector<u_char> rep(Z.begin(), Z.begin()+ d_in_bytes);
	out.clear();
	out = rep;
}


void keccak(int c, vector<u_char> &S, int d, vector<u_char> &out)
{
	sponge(S, d, 1600-c, out);
}

void shake128(vector<u_char> in, int d, vector<u_char> &out)
{
	u_char pad = 0x1f;
	in.push_back(pad);
	keccak(256, in, d, out);
}

int main(int argc, char * argv[])
{
	if(argc != 2)
	{
		cout<<"Wrong parameters"<<endl;
		return 1;
	}
	
	int d_in_bytes = stoi(argv[1]);
	vector<u_char> in;
	vector<u_char> out;

	char c;
	while (cin.get(c))
	{
		in.push_back(c);
	}
	
	shake128(in, d_in_bytes, out);
	for(int i=0; i<d_in_bytes; ++i)
    	cout << std::setfill('0') << std::setw(2)<< std::hex << (int)out[i];
	cout<<endl;
	return 0;
}