/*
===============================================================================
         ___________.__            .____________.__  __          
     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
        \/    \/            \/      \/       \/               \/ 

                                                  
   eFindSite - ligand-binding site prediction from meta-threading

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org wfeinstein@lsu.edu

   Copyright 2013 Michal Brylinski, Wei Feinstein

   This file is part of eFindSite.

   eFindSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eFindSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#include "efindsite_omp.h"

using namespace std;

int main(int argc, char *argv[])
{
	printf("BEGIN PART2\n");
	time_t  t_beg, t_end;

	int     t_len1;
	char    *t_seq1;
	int     *t_res1;
	double  *t_xyz1;

	int     n_offset;
	int     s_offset;
	int     t_offset[MAXTPL];

	int     *t_len2;
	char    *t_seq2;
	int     *t_res2;
	double  *t_xyz2;

	int     *t_sco1;
	double  *t_sco2;
	double  *t_sco3;
	double  *t_sco4;
	int     *t_alig;
	double  *t_rmat;

	string output_name;
	bool output_opt = false;

	for ( int i = 0; i < argc; i++ ) 
	{
		if ( !strcmp(argv[i],"-j") && i < argc )
		{ output_name = string(argv[i+1]); output_opt = true; }
	}

	if ( !output_opt )
	{
		cout << "Provide job tag" << endl;
		exit(EXIT_FAILURE);
	}

	ifstream from_disk(output_name+".AtoB.out");
	printf("LOAD AtoB\n");

	string line;
	stringstream ss_line;
	string str_val;

	getline(from_disk, line);
	sscanf(line.c_str(), "%d", &t_len1);

	t_seq1 = new char[t_len1];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; ss_line.get(t_seq1[i]); i++ );

	t_res1 = new int[t_len1];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%d", &t_res1[i]);
	}

	t_xyz1 = new double[t_len1*3];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%lf", &t_xyz1[i]);
	}

	getline(from_disk, line);
	sscanf(line.c_str(), "%d", &n_offset);

	getline(from_disk, line);
	sscanf(line.c_str(), "%d", &s_offset);

	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%d", &t_offset[i]);
	} 

	t_len2 = new int[n_offset];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%d", &t_len2[i]);
	}

	t_seq2 = new char[s_offset];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; ss_line.get(t_seq2[i]); i++ );

	t_res2 = new int[s_offset];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%d", &t_res2[i]);
	}

	t_xyz2 = new double[s_offset*3];
	getline(from_disk, line);
	ss_line.str(line);
	ss_line.clear();
	for( int i=0; getline(ss_line, str_val, ' '); i++ )
	{
		sscanf(str_val.c_str(), "%lf", &t_xyz2[i]);
	}

	from_disk.close();


	t_sco1 = new int[n_offset];
	t_sco2 = new double[n_offset];
	t_sco3 = new double[n_offset];
	t_sco4 = new double[n_offset];
	t_alig = new int[t_len1 * n_offset];
	t_rmat = new double[n_offset *12];

	cout << "Calculating alignments ... " << flush;

 time(&t_beg);
 
 int tm_i;
 
 #pragma omp parallel for schedule(dynamic) private(tm_i)
 for ( tm_i = 0; tm_i < n_offset; tm_i++ )
 {
  int    o_sco1;
  double o_sco2, o_sco3, o_sco4;
  
  int o_len1 = t_len1;
  int o_len2 = t_len2[tm_i];
  
  char o_seq1[MAXPRO];
  char o_seq2[MAXPRO];
  
  int o_res1[MAXPRO];
  int o_res2[MAXPRO];
  
  double o_xyz1[MAXPRO][3];
  double o_xyz2[MAXPRO][3];
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
  {
   o_seq1[t_i] = t_seq1[t_i];
   
   o_res1[t_i] = t_res1[t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    o_xyz1[t_i][t_j] = t_xyz1[t_i*3+t_j];
  }
  
  for ( int t_i = 0; t_i < o_len2; t_i++ )
  {
   o_seq2[t_i] = t_seq2[t_offset[tm_i]+t_i];
   
   o_res2[t_i] = t_res2[t_offset[tm_i]+t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    o_xyz2[t_i][t_j] = t_xyz2[t_offset[tm_i]*3+t_i*3+t_j];
  }
  
  int o_alig[MAXPRO];
  double o_t[3];
  double o_u[3][3];
  
  frtmalign_( &o_sco1, &o_sco2, &o_sco3, &o_sco4, &o_len2, &o_len1, &o_seq2, &o_seq1, &o_alig, &o_res2, &o_res1, &o_xyz2, &o_xyz1, &o_t, &o_u, &o_len1 );
  
  t_sco1[tm_i] = o_sco1;
  t_sco2[tm_i] = o_sco2;
  t_sco3[tm_i] = o_sco3;
  t_sco4[tm_i] = 0.0;
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
   t_alig[tm_i*o_len1+t_i] = o_alig[t_i];
  
  for ( int t_i = 0; t_i < 3; t_i++ )
  {
   t_rmat[tm_i*12+t_i] = o_t[t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    t_rmat[tm_i*12+3+3*t_i+t_j] = o_u[t_i][t_j];
  }
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
   if ( o_alig[t_i] > -1 )
    if ( o_seq1[t_i] == o_seq2[o_alig[t_i]] )
     t_sco4[tm_i]++;
  
  t_sco4[tm_i] = t_sco4[tm_i] / ( (double) t_sco1[tm_i] );
 }

	time(&t_end);
	printf("\n\t******* alignment time: %.2f s *******\n", difftime(t_end, t_beg));
	FILE *to_disk;

	string file_BtoC(output_name+".BtoC.out");
	printf("SAVE BtoC\n");
	to_disk = fopen(file_BtoC.c_str(), "w");

	fprintf(to_disk, "%d\n", t_len1);

	fprintf(to_disk, "%d\n", n_offset);

	for ( int i=0; i < t_len1*n_offset; i++ )
	{
		fprintf(to_disk, "%d ", t_alig[i]);
	}
	fprintf(to_disk, "\n");

	for( int i=0; i<n_offset*12; i++)
	{
		fprintf(to_disk, "%f ", t_rmat[i]);
	}
	fprintf(to_disk, "\n");

	for( int i=0; i<n_offset; i++ )
	{
		fprintf(to_disk, "%d ", t_sco1[i]);
	}
	fprintf(to_disk, "\n");

	for( int i=0; i<n_offset; i++ )
	{
		fprintf(to_disk, "%f ", t_sco2[i]);
	}
	fprintf(to_disk, "\n");

	for( int i=0; i<n_offset; i++ )
	{
		fprintf(to_disk, "%f ", t_sco3[i]);
	}
	fprintf(to_disk, "\n");

	for( int i=0; i<n_offset; i++ )
	{
		fprintf(to_disk, "%f ", t_sco4[i]);
	}
	fprintf(to_disk, "\n");

	fclose(to_disk);

	delete [] t_seq1;
	delete [] t_res1;
	delete [] t_xyz1;
	delete [] t_len2;
	delete [] t_seq2;
	delete [] t_res2;
	delete [] t_xyz2;

	delete [] t_sco1;
	delete [] t_sco2;
	delete [] t_sco3;
	delete [] t_sco4;
	delete [] t_alig;
	delete [] t_rmat;

	printf("ENDING Part2\n");
	return 0;
}
