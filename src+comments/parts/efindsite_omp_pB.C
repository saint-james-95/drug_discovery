#include "efindsite_omp.h"
#define DEBUG_LEE 0

using namespace std;

int main(int argc, char *argv[])
{
    printf("BEGIN PART2\n");
    time_t  t_beg, t_end;
    /* declare variables for reading data from file */
    int     t_len1;
    char    *t_seq1;
    int     *t_res1;
    double  *t_xyz1;

    int     n_offset; // size of template_set
    int     s_offset; // total number of residues
    int     t_offset[MAXTPL]; // index of starting position for ith protein in template_set

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

    /* ***WARNGING***
     * input code expects a specific ordering of variables
     * it does no error checking
     */

    string line;
    stringstream ss_line;
    string str_val;

    /* load t_len1 */
    getline(from_disk, line);
    sscanf(line.c_str(), "%d", &t_len1);
    if(DEBUG_LEE) printf("from disk, t_len: %d\n", t_len1);

    /* load t_seq1 */
    t_seq1 = new char[t_len1];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; ss_line.get(t_seq1[i]); i++ );
    if(DEBUG_LEE) printf("t_seq1 redo this output\n"); //printf("from disk, t_seq1: %s\n", t_seq1);

    /* load t_res1 */
    t_res1 = new int[t_len1];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%d", &t_res1[i]);
    }
    if(DEBUG_LEE) printf("from disk, t_res1: ");
    if(DEBUG_LEE) for( int i=0; i<t_len1; i++ )
        printf("%d ", t_res1[i]);
    if(DEBUG_LEE) printf("\n");

    /* load t_xyz1 */
    t_xyz1 = new double[t_len1*3];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%lf", &t_xyz1[i]);
    }
    if(DEBUG_LEE) printf("from disk, t_xyz1: ");
    if(DEBUG_LEE) for( int i=0; i<t_len1*3; i++ )
        printf("%f ", t_xyz1[i]);
    if(DEBUG_LEE) printf("\n");

    /* load n_offset */
    getline(from_disk, line);
    sscanf(line.c_str(), "%d", &n_offset);
    if(DEBUG_LEE) printf("from disk, n_offset: %d\n", n_offset);

    /* load s_offset */
    getline(from_disk, line);
    sscanf(line.c_str(), "%d", &s_offset);
    if(DEBUG_LEE) printf("from disk, s_offset: %d\n", s_offset);

    /* load t_offset */
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%d", &t_offset[i]);
    } 
    if(DEBUG_LEE) printf("from disk, t_offset: "); 
    if(DEBUG_LEE) for( int i=0; i<n_offset; i++ )
        printf("%d ", t_offset[i]);
    if(DEBUG_LEE) printf("\n");
    
    /* load t_len2 */
    t_len2 = new int[n_offset];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%d", &t_len2[i]);
    }
    if(DEBUG_LEE) printf("from disk, t_len2: ");
    if(DEBUG_LEE) for( int i=0; i<n_offset; i++ )
    {
        printf("%d ", t_len2[i]);
    }
    if(DEBUG_LEE) printf("\n");

    /* load t_seq2 */
    t_seq2 = new char[s_offset];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; ss_line.get(t_seq2[i]); i++ );
    if(DEBUG_LEE) printf("redo this output line t_seq2\n"); //printf("%s\n", t_seq2);

    /* load t_res2 */
    t_res2 = new int[s_offset];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%d", &t_res2[i]);
    }
    if(DEBUG_LEE) printf("from disk, t_res2: ");
    if(DEBUG_LEE) for( int i=0; i<s_offset; i++ )
    {
        printf("%d ", t_res2[i]);
    }
    if(DEBUG_LEE) printf("\n");

    /* load t_xyz2 */
    t_xyz2 = new double[s_offset*3];
    getline(from_disk, line);
    ss_line.str(line);
    ss_line.clear();
    for( int i=0; getline(ss_line, str_val, ' '); i++ )
    {
        sscanf(str_val.c_str(), "%lf", &t_xyz2[i]);
    }
    if(DEBUG_LEE) printf("from disk, t_xyz2: ");
    if(DEBUG_LEE) for( int i=0; i<s_offset*3; i++ )
        printf("%f ", t_xyz2[i]);
    if(DEBUG_LEE) printf("\n");

    from_disk.close();


    t_sco1 = new int[n_offset];
    t_sco2 = new double[n_offset];
    t_sco3 = new double[n_offset];
    t_sco4 = new double[n_offset];
    t_alig = new int[t_len1 * n_offset];
    t_rmat = new double[n_offset *12];

    cout << "Calculating alignments ... " << flush;

    time(&t_beg); // start alignment calculations

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
    /* these files need to be written to disk for rest of program
     *      t_alig
     *      t_rmat
     *      t_sco1
     *      t_sco2
     *      t_sco3
     *      t_sco4
     *      t_len1
     */
    FILE *to_disk;
    string file_BtoC(output_name+".BtoC.out");
    printf("SAVE BtoC\n");
    to_disk = fopen(file_BtoC.c_str(), "w");

    // save t_len1 
    fprintf(to_disk, "%d\n", t_len1);

    // save n_offset
    fprintf(to_disk, "%d\n", n_offset);

    // save t_alig 
    for ( int i=0; i < t_len1*n_offset; i++ )
    {
        fprintf(to_disk, "%d ", t_alig[i]);
    }
    fprintf(to_disk, "\n");

    //save t_rmat 
    for( int i=0; i<n_offset*12; i++)
    {
        fprintf(to_disk, "%f ", t_rmat[i]);
    }
    fprintf(to_disk, "\n");

    // save t_sco1
    for( int i=0; i<n_offset; i++ )
    {
        fprintf(to_disk, "%d ", t_sco1[i]);
    }
    fprintf(to_disk, "\n");
    
    // save t_sco2
    for( int i=0; i<n_offset; i++ )
    {
        fprintf(to_disk, "%f ", t_sco2[i]);
    }
    fprintf(to_disk, "\n");
    
    // save t_sco3
    for( int i=0; i<n_offset; i++ )
    {
        fprintf(to_disk, "%f ", t_sco3[i]);
    }
    fprintf(to_disk, "\n");
    
    // save t_sco4
    for( int i=0; i<n_offset; i++ )
    {
        fprintf(to_disk, "%f ", t_sco4[i]);
    }
    fprintf(to_disk, "\n");

    fclose(to_disk);


    /* free memory */
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
