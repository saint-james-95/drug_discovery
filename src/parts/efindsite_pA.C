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
 printf("BEGIN PARTA\n");

 time_t t_start, t_end, t_bench1, t_bench2;
 //time(&t_start);
 
 double       cut_tmscore   = 0.40;
 double       cut_seqid     = 1.00;
 double       cut_clustdis  = 8.00;
 unsigned int cut_templates = MAXTPL;
 double       cut_binrest   = 0.25;
 int          cut_binresn   = 1;
 double       cut_clustlig  = 0.70;
 std::string  met_clustlig  = "T";
 
#ifdef __linux
 
 std::string  met_clustdis  = "P";
 
#else
 
 std::string  met_clustdis  = "L";
 
#endif
 
 cout << "------------------------------------------------------------" << endl
      << "                         efindsite" << endl
      << "                        version 1.2" << endl
      << "              ligand binding pocket prediction" << endl << endl
      << "       report bugs and issues to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 7 )
 {
  cout << " efindsite -s <target structure in PDB format>" << endl
       << "           -t <templates detected by eThread>" << endl
       << "           -i <secondary structure profile by psipred>" << endl
       << "           -e <sequence profile>" << endl
       << "           -o <output filename>" << endl << endl
       << " additional options:" << endl << endl
       << "           -l <auxiliary ligands in SDF format>" << endl
       << "           -b <sequence identity threshold (default 1.0)>" << endl
       << "           -m <TMscore threshold (default 0.4)>" << endl
       << "           -x <max number of templates (default " << MAXTPL << ")>" << endl
       << "           -r <binding residue threshold (default 0.25)>" << endl
       << "           -n <min # of binding residues (default 1)>" << endl
       << "           -d <pocket clustering cutoff (default 8.0)>" << endl

#ifdef __linux

       << "           -g <pocket clustering method (Linux only, default P)>" << endl
       << "               P - affinity propagation" << endl
       << "               L - average linkage" << endl

#endif

       << "           -f <fingerprint clustering cutoff (default 0.7)>" << endl
       << "           -c <fingerprint clustering method (default T)>" << endl
       << "               T - classical Tanimoto coeff" << endl
       << "               A - average Tanimoto coeff" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string target_name;
 bool target_opt = false;
 
 string psipred_name;
 bool psipred_opt = false;
 
 string sequence_name;
 bool sequence_opt = false;
 
 string templates_name;
 bool templates_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 string cmps_name;
 bool cmps_opt = false;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-s") && i < argc ) { target_name    = string(argv[i+1]); target_opt    = true; }
  if ( !strcmp(argv[i],"-i") && i < argc ) { psipred_name   = string(argv[i+1]); psipred_opt   = true; }
  if ( !strcmp(argv[i],"-e") && i < argc ) { sequence_name  = string(argv[i+1]); sequence_opt  = true; }
  if ( !strcmp(argv[i],"-t") && i < argc ) { templates_name = string(argv[i+1]); templates_opt = true; }
  if ( !strcmp(argv[i],"-o") && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
  if ( !strcmp(argv[i],"-l") && i < argc ) { cmps_name      = string(argv[i+1]); cmps_opt      = true; }
  if ( !strcmp(argv[i],"-b") && i < argc ) { cut_seqid      = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-m") && i < argc ) { cut_tmscore    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-x") && i < argc ) { cut_templates  = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-r") && i < argc ) { cut_binrest    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-n") && i < argc ) { cut_binresn    = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-f") && i < argc ) { cut_clustlig   = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-c") && i < argc ) { met_clustlig   = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-d") && i < argc ) { cut_clustdis   = atof(argv[i+1]);                         }
  
#ifdef __linux
  
  if ( !strcmp(argv[i],"-g") && i < argc ) { met_clustdis   = string(argv[i+1]);                       }
  
#endif
  
 }
 
 char * path1;
 
 path1 = getenv("EF_LIB"); if ( path1==NULL ) { cout << "EF_LIB is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EF_MAP"); if ( path1==NULL ) { cout << "EF_MAP is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EF_MOD"); if ( path1==NULL ) { cout << "EF_MOD is not set" << endl; exit(EXIT_FAILURE); }
 
#ifdef __linux
 
 path1 = getenv("EF_APL"); if ( path1==NULL ) { cout << "EF_APL is not set" << endl; exit(EXIT_FAILURE); }
 
#endif
 
 string lib_path;
 lib_path = getenv("EF_LIB");
 
 string cluster_def;
 cluster_def = getenv("EF_MAP");
 
 string model_path;
 model_path = getenv("EF_MOD");
 
#ifdef __linux
 
 string ap_lib;
 ap_lib = getenv("EF_APL");
 
#endif
 
 ifstream f01( (model_path+"/compositionSVM.model").c_str() );
 ifstream f02( (model_path+"/profileSVM.model").c_str() );
 ifstream f03( (model_path+"/residueSVM.model").c_str() );
 ifstream f04( (model_path+"/ligandSVM.model").c_str() );
 ifstream f05( (model_path+"/rank7SVM.model").c_str() );
 ifstream f06( (model_path+"/rank8SVM.model").c_str() );
 ifstream f07( (model_path+"/compositionSVM.scale").c_str() );
 ifstream f08( (model_path+"/profileSVM.scale").c_str() );
 ifstream f09( (model_path+"/residueSVM.scale").c_str() );
 ifstream f10( (model_path+"/ligandSVM.scale").c_str() );
 ifstream f11( (model_path+"/rank7SVM.scale").c_str() );
 ifstream f12( (model_path+"/rank8SVM.scale").c_str() );
 
 if ( !f01 || !f02 || !f03 || !f04 || !f05 || !f06 || !f07 || !f08 || !f09 || !f10 || !f11 || !f12 )
 {
  cout << "Could not find SVM models in " << model_path << endl;
  exit(EXIT_FAILURE);
 }
 
#ifdef __linux
 
 ifstream ap1( (ap_lib).c_str() );
 
 if ( !ap1 )
 {
  cout << "Could not find Affinity Propagation library: " << ap_lib << endl;
  exit(EXIT_FAILURE);
 }
 
#endif
 
 if ( !target_opt )
 {
  cout << "Provide target structure in PDB format" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !psipred_opt )
 {
  cout << "Provide secondary structure profile by psipred" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !sequence_opt )
 {
  cout << "Provide sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !templates_opt )
 {
  cout << "Provide templates detected by eThread" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outprot( (output_name+".templates.pdb").c_str() );
  outprot.close();
  
  ofstream outali( (output_name+".alignments.dat").c_str() );
  outali.close();
  
  ofstream outlig( (output_name+".ligands.sdf").c_str() );
  outlig.close();
  
  ofstream outpkt1( (output_name+".pockets.pdb").c_str() );
  outpkt1.close();
  
  ofstream outpkt2( (output_name+".pockets.dat").c_str() );
  outpkt2.close();
 }
 
 if ( cut_seqid < 1.0 )
  cout << "!!! Benchmarking mode activated with max sid of " << cut_seqid << " !!!" << endl << endl;
 
 if ( cut_tmscore < 0.4 )
  cout << "!!! TMscore of " << cut_tmscore << " is below the statistical significance threshold !!!" << endl << endl;
 
 if ( cut_binresn < 1 )
 {
  cout << "!!! Min # of binding residues must be >0, setting to 1 !!!" << endl << endl;
  
  cut_binresn = 1;
 }
 
 if ( cut_binrest < ( 1 / 1e6 ) )
 {
  cout << "!!! Threshold for binding residues must be >0, setting to 0.18 !!!" << endl << endl;
  
  cut_binrest = 0.18;
 }
 else if ( cut_binrest > 1 )
 {
  cout << "!!! Threshold for binding residues must be <=1, setting to 0.18 !!!" << endl << endl;
  
  cut_binrest = 0.18;
 }
 
 if ( cut_clustlig < ( 1 / 1e6 ) )
 {
  cout << "!!! Fingerprint clustering cutoff must be >0, setting to 0.1 !!!" << endl << endl;
  
  cut_clustlig = 0.1;
 }
 else if ( cut_clustlig > 1 )
 {
  cout << "!!! Fingerprint clustering cutoff must be <=1, setting to 1.0 !!!" << endl << endl;
  
  cut_clustlig = 1.0;
 }
 
 if ( met_clustlig != "T" && met_clustlig != "A" )
 {
  cout << "!!! Fingerprint clustering method must be either T or A, setting to T !!!" << endl << endl;
  
  met_clustlig = "T";
 }
 
#ifdef __linux
 
 if ( met_clustdis != "P" && met_clustdis != "L" )
 {
  cout << "!!! Pocket clustering method must be either P or L, setting to P !!!" << endl << endl;
  
  met_clustdis = "P";
 }
 
#endif
 
 if ( cut_templates > (int) MAXTPL )
 {
  cout << "!!! Max number of templates exceeded, setting to " << MAXTPL << " !!!" << endl << endl;
  
  cut_templates = MAXTPL;
 }
 
 if ( cut_templates < 1 )
 {
  cout << "!!! Max number of templates must be >0, setting to 1 !!!" << endl << endl;
  
  cut_templates = 1;
 }
 
#ifdef __linux
 
 cout << "Checking Affinity Propagation library ... " << flush;
 
 void *dlh = NULL;
 
 if ( !( dlh = dlopen( (ap_lib).c_str(), RTLD_LAZY ) ) )
 {
  cout << dlerror() << endl;
  exit(EXIT_FAILURE);
 }
 
 char *error;
 
 if ( ( error = dlerror() ) != NULL)
 {
  cout << error << endl;
  exit(EXIT_FAILURE);
 }
 
 dlclose(dlh);
 
 cout << "looks good" << endl << endl;
 
#endif
 
 Target * target;
 
 target = new Target( 0, 0 );
 
 if ( target->loadTarget(target_name) )
 {
  cout << "Cannot read target structure" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( target->loadPsipred(psipred_name) )
 {
  cout << "Cannot read psipred profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( target->loadSequence(sequence_name) )
 {
  cout << "Cannot read sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 Cmps * compounds;
 
 compounds = new Cmps( 0 );
 
 if ( cmps_opt )
 {
  if ( compounds->loadCompounds(cmps_name) )
  {
   cout << "Cannot read auxiliary compounds" << endl;
   exit(EXIT_FAILURE);
  }
 }

 string file_AtoB = output_name+".AtoB.out";
 string file_templates = output_name+".templates.out";
 string file_AtoC = output_name+".AtoC.out";
 printf("SAVE AtoC\n");
 to_disk = fopen(file_AtoC.c_str(), "w");
 fprintf(to_disk, "%f\n", cut_tmscore);
 fprintf(to_disk, "%f\n", cut_seqid);
 fprintf(to_disk, "%f\n", cut_clustdis);
 fprintf(to_disk, "%f\n", cut_binrest);
 fprintf(to_disk, "%d\n", cut_binresn);
 fprintf(to_disk, "%f\n", cut_clustlig);
 fprintf(to_disk, "%s\n", met_clustlig.c_str());
 fprintf(to_disk, "%s\n", met_clustdis.c_str());
 fprintf(to_disk, "%s\n", target_name.c_str());
 fprintf(to_disk, "%s\n", psipred_name.c_str());
 fprintf(to_disk, "%s\n", sequence_name.c_str());
 fprintf(to_disk, "%s\n", templates_name.c_str());

 fclose(to_disk);
 
 list<string> template_list;
 
 map<string,double> template_prob1;
 map<string,double> template_prob2;
 
 getList( templates_name, cluster_def, template_list, template_prob1, template_prob2 );
 
 cout << "eFindSite library: " << lib_path <<  endl
      << "eFindSite mapping: " << cluster_def << endl << endl
 
      << "Number of ligand-bound templates: " << setw(5) << template_list.size() << endl << endl;
 
 cout << "Template filtering ... " << flush;
 
 time(&t_bench1);
 
 multimap< int, Template *, greater<int> > template_set;
 
 list<string>::iterator it1;
 double template_load_sum = 0.;
 to_disk = fopen(file_templates.c_str(), "w");
 
 for ( it1 = template_list.begin(); it1 != template_list.end(); it1++ )
  if ( template_set.size() < cut_templates )
  {
   Template * template_tmp = new Template( 0, 0, 0, template_prob1.find(*it1)->second, template_prob2.find(*it1)->second );
   
   time(&t_start)
   bool load1 = template_tmp->loadTemplate( lib_path+"/data/"+(*it1).substr(1,2)+"/"+(*it1), template_list );
   time(&t_end)
   template_load_sum += difftime(t_end, t_start)
   
   if ( !load1 )
   {
    double sid1 = template_tmp->alignNW( target->getProteinSequence() );
    
    if ( sid1 <= cut_seqid )
    {
     fprintf(to_disk, "%s\n", (*it1).c_str());
     template_set.insert( std::pair< int,Template * >( template_tmp->getProteinResiduesTotal(), template_tmp ) );
    }
   }
  }
 
 fclose(to_disk);

 int t_offset[MAXTPL];
 int n_offset = 0;
 int s_offset = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator tpl1;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  t_offset[n_offset++] = s_offset;
  
  s_offset += (*tpl1).first;
 }
 
 time(&t_bench2);
 
 cout << ". " << template_set.size() << " templates survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;

 printf("\t*********\t total load time: %.0f\t*********\n\n", template_load_sum);

 cout << "Estimating data size ... " << flush;
 
 time(&t_bench1);
 
 int n_tar = target->getProteinResiduesTotal();
 int n_off = s_offset;
 int n_set = template_set.size();
 
 int  t_len1;
 int *t_len2 = new int [n_set];
 
 char *t_seq1 = new char [n_tar];
 char *t_seq2 = new char [n_off];
 
 int *t_res1 = new int [n_tar];
 int *t_res2 = new int [n_off];
 
 double *t_xyz1 = new double [n_tar * 3];
 double *t_xyz2 = new double [n_off * 3];
 
 int    *t_sco1 = new int    [n_set];
 double *t_sco2 = new double [n_set];
 double *t_sco3 = new double [n_set];
 double *t_sco4 = new double [n_set];
 
 int *t_alig = new int [n_tar * n_set];
 
 double *t_rmat = new double [n_set * 12];
 
 t_len1 = target->getProteinResiduesTotal();
 
 char t_seqs[MAXPRO];
 
 strcpy(t_seqs, (target->getProteinSequence()).c_str());
 
 for ( int t_i = 0; t_i < target->getProteinResiduesTotal(); t_i++ )
  t_seq1[t_i] = t_seqs[t_i];
  
 for ( int t_i = 0; t_i < t_len1; t_i++ )
  t_res1[t_i] = t_i + 1;
 
 target->getProteinCoords1D(t_xyz1);
 
 int i_offset = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  t_len2[i_offset] = ((*tpl1).second)->getProteinResiduesTotal();
  
  char t_seqt[MAXPRO];
  
  strcpy(t_seqt, (((*tpl1).second)->getProteinSequence()).c_str());
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal(); t_i++ )
   t_seq2[t_i+t_offset[i_offset]] = t_seqt[t_i];
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal(); t_i++ )
   t_res2[t_i+t_offset[i_offset]] = t_i + 1;
  
  double *t_xyzt = new double [((*tpl1).second)->getProteinResiduesTotal() * 3];
  
  ((*tpl1).second)->getProteinCoords1D(t_xyzt);
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal() * 3; t_i++ )
   t_xyz2[t_i+t_offset[i_offset]*3] = t_xyzt[t_i];
  
  delete [] t_xyzt;
  
  i_offset++;
 }
 
 int mem_len = sizeof(t_len1) + n_set * sizeof(t_len2);
 int mem_seq = n_tar * sizeof(t_seq1) + n_off * sizeof(t_seq2);
 int mem_res = n_tar * sizeof(t_res1) + n_off * sizeof(t_res2);
 int mem_xyz = n_tar * 3 * sizeof(t_xyz1) + n_off * 3 * sizeof(t_xyz2);
 int mem_sco = n_set * sizeof(t_sco1) + n_set * sizeof(t_sco2) + n_set * sizeof(t_sco3) + n_set * sizeof(t_sco4);
 int mem_ali = n_tar * n_set * sizeof(t_alig);
 int mem_mat = n_set * 12 * sizeof(t_rmat);
 
 time(&t_bench2);
 
 cout << setprecision(1) << ( (double) ( mem_len + mem_seq + mem_res + mem_xyz + mem_sco + mem_ali + mem_mat ) ) / 1e6 << "M (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;

 to_disk = fopen(file_AtoB.c_str(), "w");
 printf("SAVE AtoB\n");

 fprintf(to_disk, "%d\n", t_len1);

 for( int i=0; i<t_len1; i++)
  fprintf(to_disk, "%c", t_seq1[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<n_tar; i++ )
  fprintf(to_disk, "%d ", t_res1[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<n_tar*3; i++ )
  fprintf(to_disk, "%f ", t_xyz1[i]);
 fprintf(to_disk, "\n");

 fprintf(to_disk, "%d\n", n_offset);

 fprintf(to_disk, "%d\n", s_offset);

 for( int i=0; i<n_offset; i++ )
  fprintf(to_disk, "%d ", t_offset[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<n_offset; i++ )
  fprintf(to_disk, "%d ", t_len2[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<s_offset; i++ ) 
  fprintf(to_disk, "%c", t_seq2[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<s_offset; i++ )
  fprintf(to_disk, "%d ", t_res2[i]);
 fprintf(to_disk, "\n");

 for( int i=0; i<s_offset*3; i++ )
  fprintf(to_disk, "%f ", t_xyz2[i]);
 fprintf(to_disk, "\n");

 delete [] t_len2;
 delete [] t_seq1;
 delete [] t_seq2;
 delete [] t_res1;
 delete [] t_res2;
 delete [] t_xyz1;
 delete [] t_xyz2;
 delete [] t_sco1;
 delete [] t_sco2;
 delete [] t_sco3;
 delete [] t_sco4;
 delete [] t_alig;
 delete [] t_rmat;

 fclose(to_disk);

 printf("ENDING PART1\n");
 return EXIT_SUCCESS;
}

