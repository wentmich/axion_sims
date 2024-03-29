#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>


int saveSnapshotInTxtFile(int suffix);
int load_snapshot(char *fname, int files);
int allocate_memory(void);

using namespace std;

struct io_header_1
{
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    int npartTotal[6];
    int flag_cooling;
    int num_files;
    double Boaxsize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];
} header1;


int NumPart, Ngas;

struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;
    
    float Rho, U, Temp, Ne;
} *P;

int *Id;

double Time, Redshift, L;

/* Here we load a snapshot file. It can be distributed onto
 * several files (for files > 1). The particles are brought
 * back into the order implied by their ID's.
 * A unit conversion routune is called to do unit \
 * conversion and to evaluate gas temperature. */
int main(int argc, char **argv)
{
    char path[200], input_fname[200], basename[200];
    int type, snapshot_number, files;
    
    sprintf(path, "./");
    sprintf(basename, "snapshot");
    snapshot_number = 0;
    files = 8;
    
    //printf("read snapshot total: %i\n\n", atoi(argv[1]));
    
    for(int i = 0; i < atoi(argv[1]); i++){
        int fnum = i * 5;
        //printf("Entering for loop\n"); 
        //printf("read snapshot%d\n", i);
        sprintf(input_fname, "%s/%s_%03d", path, basename, fnum);
        //printf("Filename: %s\n", input_fname);
        load_snapshot(input_fname, files);
        //reordering();
        //unit_conversion();
        saveSnapshotInTxtFile(i);
    }
}


/* Here the particle data is at your disposal */
int saveSnapshotInTxtFile(int suffix)
{
    	//printf("do_what_you_want.\n");

	double MaxPos[3]={-1e300,-1e300,-1e300},MinPos[3]={1e300,1e300,1e300};
	double MaxVel=0;
/*
	for(int i = 0; i < NumPart; i++)
	{

		for(int j = 0; j < 3; j++)
		{  		
			if( P[i].Pos[j] > MaxPos[j] ) MaxPos[j] = P[i].Pos[j];
			if( P[i].Pos[j] < MinPos[j] ) MinPos[j] = P[i].Pos[j];
		}
		double v = sqrt( P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2] );
		if( v > MaxVel ) MaxVel = v; 
	}

	printf("Boaxsize: %f\n",L);	
      	printf("MaxPosx: %f\n",MaxPos[0]);
      	printf("MaxPosy: %f\n",MaxPos[1]);
      	printf("MaxPosz: %f\n",MaxPos[2]);
      	printf("MinPosx: %f\n",MinPos[0]);
      	printf("MinPosy: %f\n",MinPos[1]);
      	printf("MinPosz: %f\n",MinPos[2]);
      	printf("MaxVel: %f\n",MaxVel);
*/
    
    	FILE *saved = stdout;
    	char filename[200];
    	sprintf(filename, "data_%i.inc", suffix);
    	//sprintf(filename, "data.inc");
    	//stdout = fopen(filename, "a");

	fstream fout;
	fout.open(filename,ios::out);

	for(int i = 0; i < NumPart; i++){
	//	double px = (P[i].Pos[0]-MinPos[0])/(MaxPos[0]-MinPos[0]);
 	//	double py = (P[i].Pos[1]-MinPos[1])/(MaxPos[1]-MinPos[1]);
	//	double pz = (P[i].Pos[2]-MinPos[2])/(MaxPos[2]-MinPos[2]);
		double px = P[i].Pos[0] / L;
		double py = P[i].Pos[1] / L;
		double pz = P[i].Pos[2] / L;
		//printf("sphere{<%f,%f,%f>,R}\n",px,py,pz);
		fout << "sphere{<"<<px<<","<<py<<","<<pz<<">,R}" << endl;
	}
    
    	//for(int i = 0; i < NumPart; i++){
    	//    printf("%i %f %f %f %f %f %f %f %f\n", Id[i], Time, Redshift, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]);	
    	//}

        //fflush(stdout);
    	//fclose(stdout);
    	//stdout = saved;
    	fout.close();
    	//printf("number particles %i; time %f and redshift %f\n", NumPart, Time, Redshift);
}


/* this template shows how one may convert from Gadget's units 
 * to cgs units. In this example, the temperature of the gas 
 * is computed (assuming that the electron density in units of
 * the hydrogen density was computed by the code. This is done
 * if cooling is enabled). */
int unit_conversion(void)
{
    double GRAVITY, BOLTZMANN, PROTONMASS;
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
    double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
    double G, Xh, HubbleParam;

    int i;
    double MeanWeight, u, gamma;
    
    /* physical constants in cgs units */
    GRAVITY = 6.672e-8;
    BOLTZMANN = 1.3806e-16;
    PROTONMASS = 1.6726e-24;

    /* internal unit system of the code */
    UnitLength_in_cm = 3.085678e21;    /*  code length unit in cm/h */
    UnitMass_in_g = 1.989e43;    /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s = 1.0e5;
    
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs = UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm * UnitLength_in_cm);
    UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / (UnitTime_in_s * UnitTime_in_s);
    UnitEnergy_in_cgs = UnitMass_in_g * (UnitLength_in_cm * UnitLength_in_cm) / (UnitTime_in_s * UnitTime_in_s);

    G = GRAVITY / (UnitLength_in_cm * UnitLength_in_cm * UnitLength_in_cm) * UnitMass_in_g * (UnitTime_in_s * UnitTime_in_s);
    
    
    Xh = 0.76;            /* mass fraction of hydrogen */
    HubbleParam = 0.65;

    for(i = 1; i <= NumPart; i++)
    {
        if(P[i].Type == 0)    /* gas particle */
        {
            MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * P[i].Ne) * PROTONMASS;
            
            /* convert internal energy to cgs units */
            
            u = P[i].U * UnitEnergy_in_cgs / UnitMass_in_g;
            
            gamma = 5.0 / 3;
            
            /* get temperature in Kelvin */
            
            P[i].Temp = MeanWeight / BOLTZMANN * (gamma - 1) * u;
        }
    }
}


/* this routine loads particle data from Gadget's default
 *  * binary file format. (A snapshot may be distributed
 *   * into multiple files.
 *    */
int load_snapshot(char *fname, int files)
{
    FILE *fd;
    char buf[200];
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
    
    for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
        if(files > 1)
            sprintf(buf, "%s.%d", fname, i);
        else
            sprintf(buf, "%s", fname);
        
        if(!(fd = fopen(buf, "r")))
        {
            printf("can't open file `%s`\n", buf);
            exit(0);
        }
        
        printf("reading `%s' ...\n", buf);
        fflush(stdout);
        
        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&header1, sizeof(header1), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);
        
        if(files == 1)
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header1.npart[k];
            Ngas = header1.npart[0];
        }
        else
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header1.npartTotal[k];
            Ngas = header1.npartTotal[0];
        }
        
        for(k = 0, ntot_withmasses = 0; k < 6; k++)
        {
            if(header1.mass[k] == 0)
                ntot_withmasses += header1.npart[k];
        }
        
        if(i == 0)
            allocate_memory();
        
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                //printf("read pos %i\n",pc_new);
                fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                //printf("read vel %i\n",pc_new);
                
                fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
        
        
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(&Id[pc_new], sizeof(int), 1, fd);
                pc_new++;
            }
        }
        SKIP;

        if(ntot_withmasses > 0)
            SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                P[pc_new].Type = k;
                
                if(header1.mass[k] == 0)
                    fread(&P[pc_new].Mass, sizeof(float), 1, fd);
                else
                    P[pc_new].Mass = header1.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses > 0)
            SKIP;
        
        
        if(header1.npart[1] > 0)
        {
            SKIP;
            for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
            {
                fread(&P[pc_sph].U, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            SKIP;
            for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
            {
                fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            if(header1.flag_cooling)
            {
                SKIP;
                for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
                {
                    fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
                    pc_sph++;
                }
                SKIP;
            }
            else
                for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
                {
                    P[pc_sph].Ne = 1.0;
                    pc_sph++;
                }
        }
        
        fclose(fd);
    }
    
    Time = header1.time;
    Redshift = header1.redshift;
    L = header1.Boaxsize;
    
}



/* this routine allocates the memory for the
 *  * particle data.
 *   */
int allocate_memory(void)
{
    //printf("allocating memory...\n");
    
    if(!(P = static_cast<particle_data*>(malloc(NumPart * sizeof(struct particle_data)))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    
    P--;                /* start with offset 1 */
    
    
    if(!(Id = static_cast<int*>(malloc(NumPart * sizeof(int)))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    
    Id--;                /* start with offset 1 */
    
    //printf("allocating memory...done\n");
}



/* This routine brings the particles back into
 *  * the order of their ID's.
 *   * NOTE: The routine only works if the ID's cover
 *    * the range from 1 to NumPart !
 *     * In other cases, one has to use more general
 *      * sorting routines.
 *       */
int reordering(void)
{
    int i, j;
    int idsource, idsave, dest;
    struct particle_data psave, psource;
    
    
    //printf("reordering....\n");
    
    for(i = 1; i <= NumPart; i++)
    {
        if(Id[i] != i)
        {
            psource = P[i];
            idsource = Id[i];
            dest = Id[i];
            
            do
            {
                psave = P[dest];
                idsave = Id[dest];
                
                P[dest] = psource;
                Id[dest] = idsource;
                
                if(dest == i)
                    break;
                
                psource = psave;
                idsource = idsave;
                
                dest = idsource;
            }
            while(1);
        }
    }
    
    //printf("done.\n");
    
    Id++;
    free(Id);
    
    //printf("space for particle ID freed\n");
}
