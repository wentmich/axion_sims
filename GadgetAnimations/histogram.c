#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <string>
using namespace std;

const short Nbins = 30;

int saveSnapshotInTxtFile(int suffix);
int load_snapshot(char *fname, int files, int indicator, int max);
int allocate_memory(int N);
int get_histogram_index(int N, float x, float boxlength);

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

struct histogram
{
    float grid[Nbins][Nbins][Nbins];
} histogram1;

float get_histogram_max(histogram hist);
float get_histogram_min(histogram hist);
float get_histogram_avg(histogram hist);
int saveSnapshotInBinaryFile(int suffix);
void make_oFNAME(char* oFNAME);
int NumPart, Ngas;

struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;
    
    float Rho, U, Temp, Ne;
} *P, single_particle;

int *Id;
float maximum, minimum;
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
    files = 1;
    
    //printf("read snapshot total: %i\n\n", atoi(argv[1]));
    
    for(int i = (atoi(argv[1])); i >= 0; i--){
        int fnum = i;
        sprintf(input_fname, "%s/%s_%03d", path, basename, fnum);
        load_snapshot(input_fname, files, i, atoi(argv[1]));
        saveSnapshotInBinaryFile(i);
    }
}

int saveSnapshotInBinaryFile(int suffix)
{
    int x, y, z;
    int value;
    char filename[200];
    FILE *my_file;

    sprintf(filename, "data_%i.df3", suffix);
    make_oFNAME(filename);

    if ((my_file = fopen(filename,"w")) == NULL)
      return 0;
    
    fputc(Nbins >> 8, my_file);
    fputc(Nbins & 0xff, my_file);
    fputc(Nbins >> 8, my_file);
    fputc(Nbins & 0xff, my_file);
    fputc(Nbins >> 8, my_file);
    fputc(Nbins & 0xff, my_file);

    for (z = 0; z < Nbins; z++) {
      for (y  = 0; y < Nbins; y++) {
         for (x = 0; x < Nbins; x++) {
            value = (int)histogram1.grid[x][y][z];
            fputc(value, my_file);
         }
      }
    }
    fclose(my_file);
}

/* this routine loads particle data from Gadget's default
 *  * binary file format. (A snapshot may be distributed
 *   * into multiple files.
 *    */
int load_snapshot(char *fname, int files, int indicator, int max)
{
    FILE *fd;
    char buf[200];
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;
    int x, y, z;
    int count = 0;
    for (x = 0; x < Nbins; x++)
    {
        for (y = 0; y < Nbins; y++)
        {
            for (z = 0; z < Nbins; z++)
            {
                histogram1.grid[x][y][z] = 0.0;
                count = count + 1;
            }
        }
    }
    printf("%i\n", count);

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
            allocate_memory(Nbins * Nbins * Nbins);
        
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(&single_particle.Pos[0], sizeof(float), 3, fd);
                x = get_histogram_index(Nbins, single_particle.Pos[0], header1.Boaxsize);
                y = get_histogram_index(Nbins, single_particle.Pos[1], header1.Boaxsize);
                z = get_histogram_index(Nbins, single_particle.Pos[2], header1.Boaxsize);
                histogram1.grid[x][y][z] = histogram1.grid[x][y][z] + 1.0;
                pc_new++;
            }
        }

        fclose(fd);
    }

    count = 0;
    n = 0;
    if (indicator == max)
    {
        maximum = get_histogram_max(histogram1);
        minimum = get_histogram_min(histogram1);
        printf("Maximum = %f \t Minimum = %f \n", maximum, minimum);
    }
    for (x = 0; x < Nbins; x++)
    {
        for (y = 0; y < Nbins; y++)
        {
            for (z = 0; z < Nbins; z++)
            {
                count  = count + histogram1.grid[x][y][z];
                histogram1.grid[x][y][z] = (histogram1.grid[x][y][z] - minimum) * 255.0 / (maximum - minimum);
                //histogram1.grid[x][y][z] = histogram1.grid[x][y][z] * 255.0 / maximum;
                //printf("%f\n", histogram1.grid[x][y][z]);
            }
        }
    }
    printf("Histogram Count = %i\n", count);
    Time = header1.time;
    Redshift = header1.redshift;
    L = header1.Boaxsize;
}



/* this routine allocates the memory for the
 *  * particle data.
 *   */
int allocate_memory(int N)
{
    if(!(P = static_cast<particle_data*>(malloc(N * sizeof(struct particle_data)))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    P--;                /* start with offset 1 */
    if(!(Id = static_cast<int*>(malloc(N * sizeof(int)))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    Id--;                /* start with offset 1 */
}


int get_histogram_index(int N, float x, float boxlength)
{
    int next_index;
    float integer_position;
    for (int i = 0; i < N; i++)
    {
        integer_position = x * N / boxlength;
        next_index  = i + 1;
        if (i <= integer_position && integer_position <= next_index)
            return i;
    }
    return -1;
}

float get_histogram_max(histogram hist)
{
    float max = 0.0;
    int i, j, k;
    for (i = 0; i < Nbins; i++)
    {
        for (j = 0; j < Nbins; j++)
        {
            for (k = 0; k < Nbins; k++)
            {
                if (hist.grid[i][j][k] > max)
                    max = hist.grid[i][j][k];
            }
        }
    }
    return max;
}

float get_histogram_min(histogram hist)
{
    float min = 10000000000.0;
    int i, j, k;
    for (i = 0; i < Nbins; i++)
    {
        for (j = 0; j < Nbins; j++)
        {
            for (k = 0; k < Nbins; k++)
            {
                if (hist.grid[i][j][k] < min)
                    min = hist.grid[i][j][k];
            }
        }
    }
    return min;
}

void make_oFNAME(char* oFNAME)
{
        ofstream myfile (oFNAME);
        myfile.close();
}

float get_histogram_avg(histogram hist)
{
    int i, j, k;
    float avg = 0.0;
    for (i = 0; i < Nbins; i++)
    {
        for (j = 0; j < Nbins; j++)
        {
            for (k = 0; k < Nbins; k++)
            {
                avg = avg + hist.grid[i][j][k];
            }
        }
    }
    avg = avg / (Nbins * Nbins * Nbins);
    return avg;
}
