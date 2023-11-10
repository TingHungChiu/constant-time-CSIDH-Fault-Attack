/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2017 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "fp.h"
#include "csidh.h"
#include "mont.h"
#include "uint.h"



int main(void)
{
    public_key shared_alice, shared_bob;
    (void)shared_alice;
    (void)shared_bob;
    
    uint8_t num_batches = 3;
    uint8_t my = 0;
    
    
    int8_t max[num_primes] = {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
                             5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8,
                             9, 9, 9, 10, 10, 10, 10, 9, 8, 8, 8, 7, 7, 7, 7, 7, 6, 5,
                             1, 2, 2};
    
    //int8_t max[num_primes] = {2,2};
    int8_t* x[num_primes];
    for (int i =0; i<num_primes;i++){x[i]=(int8_t*)malloc(max[i]*sizeof(int8_t));}
    for (int i =0; i<num_primes;i++){for (int j =0; j<num_primes;j++){x[i][j]=j;}}
    /*
    for (int i = 0; i < num_primes; i++) {
        // Loop through the columns
        for (int j = 0; j < num_primes; j++) {
            // Print each element
            printf("%d\t", x[i][j]);
        }
        // Move to the next row
        printf("\n");
    }
    */
    unsigned int num_isogenies = 4;

    public_key pub_bob = {{{0xa7071cf2062c5b28, 0x4ef6c4e374631ad5, 0x075a4dd6d3013833, 0xa3c0a67a26b9e943, 0x51601a8a437952f2, 0x4f45902681f6b516, 0xc8364e54abb888f4, 0x0266ebb102f36783}}};
    //public_key pub_bob ={{{0xF1183E6C4C9D2F5E,0x4964FF3249317265,0x86B3C76D00EF8204,0xD4C29C88FA51822C,0xAF4D884A13871716,0xA7B6AF476FEA943F,0x4E2ED542C116F131,0x53160D344B0DA451}}};
    private_key priv_alice = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    public_key shared_secret = {{{0xa7071cf2062c5b28, 0x4ef6c4e374631ad5, 0x075a4dd6d3013833, 0xa3c0a67a26b9e943, 0x51601a8a437952f2, 0x4f45902681f6b516, 0xc8364e54abb888f4, 0x0266ebb102f36783}}};
    
    //private_key priv_alice = {{1, 1}};
    //public_key shared_secret = {{{0x641f76a330197f23, 0x4efe042d8e2784a7, 0xd2dc3abc7ac942da, 0xeb8826cd70281773, 0x7b0d8f6688dcff83, 0x483b065d1c8e9f09, 0x9e42b5b3cb804006, 0x24889036edf112d1}}};

    clock_t start, end;
    start = clock();
    
    csidh(&shared_alice, &pub_bob, &priv_alice, num_batches, max, num_isogenies, my, x);
    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf ("%f\n",time_taken);
    
    for(int i = 0; i < 8; i++) {
        printf("%02llX\n",shared_alice.A.c[i]);
    }
    
}
