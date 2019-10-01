/* This file intends to read in the initial conditions from a text
file. The initial conditions should be in the following form:

x1	y1	z1	vx1	vy1	vz1	m1
.	.	.	.	.	.	.
.	.	.	.	.	.	.
.	.	.	.	.	.	.
xN	yN	zN	vxN	vyN	vzN	mN

The file will read in this data and then write a binary file in
the Gadget-2 default snapshot format described in Volker Springle's
Users Guide for the public release of Gadget-2.

-- User Guide --
https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=
rja&uact=8&ved=2ahUKEwj9tYOd4-biAhXJmq0KHVe9BtkQFjAAegQIABAC&url=
https%3A%2F%2Fwwwmpa.mpa-garching.mpg.de%2Fgadget%2Fusers-guide.
pdf&usg=AOvVaw01eHrnWA78FiW6TFMnmEXA
----------------

The code assume float precision is being used in the Gadget
Makefile.

Author: Michael Wentzel
Eamil: wentmich@umich.edu
Date: 12 June 2019												*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

/* GLOBALS
iFNAME:			txt file containing the initial conditions in the form above
oFNAME:			name of output snapshot file
N_PARTICLES:	        total number of particles to be writte to oFNAME
leftover:		number of bytes leftover in the header (from total of 256)
filler:			character to be used to fill up the header (can be arbitrary) */
//const char* iFNAME = "/nfs/turbo/bsafdi/wentmich/Axion_Structure_Sims/simulations/miniclustersICs/initial_data_in_txt_form.txt";
//const char* oFNAME = "/nfs/turbo/bsafdi/wentmich/Axion_Structure_Sims/simulations/miniclustersICs/initial_conditions_snapshot.dat";
const unsigned long N_FULL = 8000000000; // total number of type 1 particles in the simulation
const double cfact = 19637053.537620023; //13434500.510994378 * 1.4606070432329787;// * 205.40494734447185;
const int leftover = 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 8 - 2 * 4 - 4 * 8 - 3 * 4 - 6 * 4;
const int NFILES = 1000;
const int N_PARTICLES = 8000000;
const float MASS = 32051887299076.293;
const float BOX_LENGTH = 26262918072.827072;
const char* filler = "a";
char* x;
char* y;
char* z;
char* vx;
char* vy;
char* vz;
char* m;
/* DATA STRUCTURES
particle_data:		each line from iFNAME goes into a particle_data structure
io_header_1:		structure of the header written to oFNAME  				*/
struct particle_data
{
    unsigned int ID;
    float pos[3];
    float vel[3];
    float m, u, rho, hsml, pot, acce, endt, tstp;
} particles [N_PARTICLES];

struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned long npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int flag_entr_ics;
  char fill[leftover];	// fills to 256 Bytes
} header;

/* USER FUNCTIONS
initialize_header:			for this function, all of the parameters for the header 
							file must be specified. The function takes no inputs, so
							they all must be specified manually. All quantities should
							be expressed in cgs!!!               
load_particle_data:			load in the data from iFNAME into particle_data structures
write_unformatted_binary:	write the data from all N_PARTICLES particle_data 
							structures to oFNAME 									
make_oFNAME:				make the file oFNAME and then close it 					*/
void initialize_header()
{
    header.npart[0] = 0; // number of gas particles
    header.npart[1] = N_PARTICLES; // number of halo particles (axions in this scenario)
    header.npart[2] = 0; // number of disk particles
    header.npart[3] = 0; // number of bulge particles
    header.npart[4] = 0; // number of stars
    header.npart[5] = 0; // number of boundary particles

    header.mass[0]  = 0.0; // mass of gas
    header.mass[1]  = MASS; // mass of halo
    header.mass[2]  = 0.0; // mass of disk
    header.mass[3]  = 0.0; // mass of bulge
    header.mass[4]  = 0.0; // mass of star
    header.mass[5]  = 0.0; // mass of boundary particle

    header.time = 1.0 / 3412.0; // t_0 for Newtonian or a(t_0) for cosmological
    header.redshift = 1 / header.time - 1; // redshift at t_0 used for cosmological simulation
    
    header.flag_sfr = 0; // unused
    header.flag_feedback = 0; // unused

    header.npartTotal[0] = 0;//header.npart[0]; // total number of gas particles
    header.npartTotal[1] = N_FULL; // total number of halo particles
    header.npartTotal[2] = 0;//header.npart[2]; // total number of disk particles
    header.npartTotal[3] = 0;//header.npart[3]; // total number of bulge particles
    header.npartTotal[4] = 0;//header.npart[4]; // total number of stars
    header.npartTotal[5] = 0;//header.npart[5]; // total number of boundary particles

    header.flag_cooling = 0; // unused

    header.num_files = NFILES; // number of files per snapshot

    header.BoxSize = BOX_LENGTH * cfact; // size of the periodic box from [0, L]
    header.Omega0 = 0.31; // matter density at z = 0 in units of critical density
    header.OmegaLambda = 0.31 / 3412.0; // vacuum energy density at z = 0 in units of critical density
    header.HubbleParam = 3.2407789e-18 * 0.5665; // Hubble constant in units of km s^-1 Mpc^-1
    header.FlagAge = 0;  // unused
    header.FlagMetals = 0; // unused
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.NallHW[0] = 0; // 0 for less than 2^32 particles
    header.flag_entr_ics = 0; // 0 if thermal energy u instead of entropy
}

void load_particle_data(char * FileName)
{
    const int NLINE = 1024;
    char line_content[NLINE];
    float data[7];

    FILE * file;
    file = fopen(FileName, "r");
    for (unsigned int i = 0; i < N_PARTICLES; ++i)
    {
        fgets(line_content, NLINE, file);
        stringstream ss(line_content);
        for (int j = 0; j < 7; ++j)
            ss >> data[j];
        particles[i].pos[0] = data[0] * cfact; particles[i].pos[1] = data[1] * cfact; particles[i].pos[2] = data[2] * cfact;
        particles[i].vel[0] = data[3]; particles[i].vel[1] = data[4]; particles[i].vel[2] = data[5];
        particles[i].m = data[6]; particles[i].ID = i;
        //printf(line_content);
        //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\n", particles[i].pos[0], particles[i].pos[1], particles[i].pos[2],
        //        particles[i].vel[0], particles[i].vel[1], particles[i].vel[2], particles[i].m, particles[i].ID);
    }
}

void write_unformatted_binary(char * FileName)
{
	// open the file "oFNAME"
	ofstream binFile;
	binFile.open(FileName, ios::in | ios::binary);

    if (!binFile)
    {
        cout << "Cannot open file.\n";
    }
    
    // write the header info to the first block
    int block_size = 256; // size of block in bytes
    binFile.write((char*) &block_size, sizeof(int)); // write block size
    binFile.write((char*) &header, sizeof(io_header_1));
    binFile.write((char*) &block_size, sizeof(int)); // write block size

    // write the positions to the second block
    block_size = sizeof(double) * 3 * N_PARTICLES;
    binFile.write((char*) &block_size, sizeof(int)); // write block size
    for (int i = 0; i < N_PARTICLES; ++i)
    	binFile.write((char*) &particles[i].pos , sizeof(float) * 3);
    binFile.write((char*) &block_size, sizeof(int)); // write block size

    // write the velocities to the third block
    block_size = sizeof(double) * 3 * N_PARTICLES;
    binFile.write((char*) &block_size, sizeof(int)); // write block size
    for (int i = 0; i < N_PARTICLES; ++i)
    	binFile.write((char*) &particles[i].vel , sizeof(float) * 3);
    binFile.write((char*) &block_size, sizeof(int)); // write block size

    // write the particle IDs to the fourth block
    block_size = sizeof(int) * N_PARTICLES;
    binFile.write((char*) &block_size, sizeof(int)); // write block size
    for (int i = 0; i < N_PARTICLES; ++i)
    	binFile.write((char*) &particles[i].ID, sizeof(int));
    binFile.write((char*) &block_size, sizeof(int)); // write block size

    // in this case, masses are specified in the header
    // there are no gas particles, so there's no gas temperature
    // nothing else needs to be written for the initial conditions
}

void make_oFNAME(char* oFNAME)
{
	ofstream myfile (oFNAME);
	myfile.close();
}

//------------------------------ MAIN FUNCTION -----------------------------------
int main()
{
    char iFNAME[200];
    char oFNAME[200];
    int i;
    for (i = NFILES - 1; i >= 0; --i)
    {
        sprintf(iFNAME, "/scratch/bsafdi_root/bsafdi/wentmich/ics/ics.%i.txt", i);
        load_particle_data(iFNAME);
        initialize_header();
	printf("box size: %f \n", header.BoxSize);
        sprintf(oFNAME, "/scratch/bsafdi_root/bsafdi/wentmich/miniclusters-full/ics/ics.%i", i);
        make_oFNAME(oFNAME);
	write_unformatted_binary(oFNAME);
        //iFNAME.clear();
        //oFNAME.clear();
    }
}
