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
 printf("BEGIN PARTC\n");

 Cmps * compounds;
 compounds = new Cmps( 0 );
 bool cmps_opt = false;
 string cmps_name;

 time_t t_bench1, t_bench2;

 double       cut_tmscore;
 double       cut_seqid;
 double       cut_clustdis;
 double       cut_binrest;
 int          cut_binresn;
 double       cut_clustlig;
 std::string  met_clustlig;
 std::string  met_clustdis;
 string       lib_path;
 string       cluster_def;
 string       model_path;
 string       ap_lib;
 string       target_name;
 string       psipred_name;
 string       sequence_name;
 string       templates_name;

 Target       *target;

 int     t_len1;
 int     n_offset;
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

 stringstream ss_line;
 string  str_val;

 lib_path = getenv("EF_LIB");
 cluster_def = getenv("EF_MAP");
 model_path = getenv("EF_MOD");
 ap_lib = getenv("EF_APL");

 printf("LOAD AtoC\n");
 ifstream from_disk(output_name+".AtoC.out");
 string line;

 getline(from_disk, line);
 sscanf(line.c_str(), "%lf", &cut_tmscore);

 getline(from_disk, line);
 sscanf(line.c_str(), "%lf", &cut_seqid);

 getline(from_disk, line);
 sscanf(line.c_str(), "%lf", &cut_clustdis);

 getline(from_disk, line);
 sscanf(line.c_str(), "%lf", &cut_binrest);

 getline(from_disk, line);
 sscanf(line.c_str(), "%d", &cut_binresn);

 getline(from_disk, line);
 sscanf(line.c_str(), "%lf", &cut_clustlig);

 getline(from_disk, met_clustlig);
 getline(from_disk, met_clustdis);
 getline(from_disk, target_name);
 getline(from_disk, psipred_name);
 getline(from_disk, sequence_name);
 getline(from_disk, templates_name);

 from_disk.close();

 printf("LOAD BtoC\n");
 from_disk.open(output_name+".BtoC.out");

 getline(from_disk, line);
 sscanf(line.c_str(), "%d", &t_len1);

 getline(from_disk, line);
 sscanf(line.c_str(), "%d", &n_offset);

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_alig = new int[t_len1*n_offset];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%d", &t_alig[i]);
 }

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_rmat = new double[n_offset*12];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%lf", &t_rmat[i]);
 }

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_sco1 = new int[n_offset];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%d", &t_sco1[i]);
 }

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_sco2 = new double[n_offset];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%lf", &t_sco2[i]);
 }

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_sco3 = new double[n_offset];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%lf", &t_sco3[i]);
 }

 getline(from_disk, line);
 ss_line.str(line);
 ss_line.clear();
 t_sco4 = new double[n_offset];
 for( int i=0; getline(ss_line, str_val, ' '); i++ )
 {
  sscanf(str_val.c_str(), "%lf", &t_sco4[i]);
 }

 from_disk.close();

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

 list<string> template_list;
 map<string,double> template_prob1;
 map<string,double> template_prob2;

 time(&t_bench1);
 getList( templates_name, cluster_def, template_list, template_prob1, template_prob2 );
 time(&t_bench2);
 printf("\tgetList time: %.2f\n", difftime(t_bench2, t_bench1) );


 multimap< int, Template *, greater<int> > template_set;
 list<string>::iterator it1;

 time(&t_bench1);
 string file_templates(output_name+".templates.out");
 printf("LOAD template_set\n");
 from_disk.open(file_templates);
 while(getline(from_disk, line))
 {
  Template * template_tmp = new Template( 0,  0, 0, template_prob1.find(line)->second, template_prob2.find(line)->second );
  template_tmp->loadTemplate( lib_path+"/data/"+(line).substr(1,2)+"/"+(line), template_list );
  template_set.insert( std::pair< int, Template * >( template_tmp->getProteinResiduesTotal(), template_tmp ) );
 }
 time(&t_bench2);
 from_disk.close();
 printf("\ttemplate_set loading time: %.2f (s)\n", difftime(t_bench2, t_bench1) );


 if ( cmps_opt )
 {
  if ( compounds->loadCompounds(cmps_name) )
  {
   cout << "Cannot read auxiliary compounds" << endl;
   exit(EXIT_FAILURE);
  }
 }

 int j_offset = 0;
 std::multimap< int, Template *, greater<int> >::iterator tpl1;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  ((*tpl1).second)->setProteinLengthTM( t_sco1[j_offset] );
  ((*tpl1).second)->setProteinRMSD( t_sco2[j_offset] );
  ((*tpl1).second)->setProteinTMscore( t_sco3[j_offset] );
  ((*tpl1).second)->setProteinSeqID2( t_sco4[j_offset] );
  
  int p_alig[MAXPRO];
  
  for ( int t_i = 0; t_i < t_len1; t_i++ )
   p_alig[t_i] = t_alig[j_offset*t_len1+t_i];
  
  ((*tpl1).second)->setTMalignment(p_alig, t_len1);
  
  double p_t[3];
  double p_u[3][3];
  
  for ( int t_i = 0; t_i < 3; t_i++ )
  {
   p_t[t_i] = t_rmat[j_offset*12+t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    p_u[t_i][t_j] = t_rmat[j_offset*12+3+3*t_i+t_j];
  }
  
  ((*tpl1).second)->setMatrix(p_t, p_u);
  
  j_offset++;
 }
 
 list<string> template_list_filtered;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); )
 {
  std::multimap< int, Template *, greater<int> >::iterator tpl6 = tpl1++;
  
  if ( ((*tpl6).second)->getProteinTMscore() < cut_tmscore )
   template_set.erase(tpl6);
  else
   template_list_filtered.push_back( ((*tpl6).second)->getProteinID() );
 }
 
 int ltot = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  ((*tpl1).second)->purgeAlignments(template_list_filtered);
  
  ((*tpl1).second)->calculateContacts();
  
  ltot += ((*tpl1).second)->getLigandsTotal();
 }
 
 if ( template_set.empty() )
 {
  cout << "no templates survived" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << template_set.size() << "/" << ltot << " templates/ligands survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Detecting pockets ... " << flush;
 
 time(&t_bench1);
 
 double *sim1 = new double[ltot*ltot];
 
 std::multimap< int, Template *, greater<int> >::iterator tpl2;
 std::multimap< int, Template *, greater<int> >::iterator tpl3;
 
 int nl1 = 0;
 int nl2 = 0;
 
 for ( tpl2 = template_set.begin(); tpl2 != template_set.end(); tpl2++ )
  for ( int il1 = 0; il1 < ((*tpl2).second)->getLigandsTotal(); il1++ )
  {
   for ( tpl3 = template_set.begin(); tpl3 != template_set.end(); tpl3++ )
    for ( int il2 = 0; il2 < ((*tpl3).second)->getLigandsTotal(); il2++ )
    {
     sim1[nl1*ltot+nl2] = getDistance( 1, (*tpl2).second, il1, (*tpl3).second, il2 );
     
     nl2++;
    }
   
   nl2 = 0;
   
   nl1++;
  }
 
 int * clu1 = new int [nl1];
 
 int clu2;
 
#ifdef __linux
 
 if ( met_clustdis == "P" && template_set.size() > 4 )
  clu2 = cluster_ap( sim1, clu1, nl1, cut_clustdis, ap_lib);
 
 else
  clu2 = cluster_avelink( sim1, clu1, nl1, cut_clustdis, "min" );
 
#else
 
 clu2 = cluster_avelink( sim1, clu1, nl1, cut_clustdis, "min" );
 
#endif
 
 delete [] sim1;
 
 if ( clu2 < 1 )
 {
  cout << "no pockets found" << endl << endl;
  
  template_set.clear();
  
  exit(EXIT_SUCCESS);
 }
 
 clu2 = refine_pockets(template_set, clu2, clu1, nl1, cut_clustdis);
 
 time(&t_bench2);
 
 cout << clu2 << " pockets found (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 int nl3 = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator tpl4;
 
 for ( tpl4 = template_set.begin(); tpl4 != template_set.end(); tpl4++ )
  for ( int il1 = 0; il1 < ((*tpl4).second)->getLigandsTotal(); il1++ )
   ((*tpl4).second)->setPocketNumber(il1, clu1[nl3++]);
 
 list< Pocket * > pocket_set;
 
 for ( int clu3 = 0; clu3 < clu2; clu3++ )
 {
  Pocket * pocket_tmp = new Pocket( clu3 );
  
  std::multimap< int, Template *, greater<int> >::iterator tpl5;
  
  for ( tpl5 = template_set.begin(); tpl5 != template_set.end(); tpl5++ )
   for ( int il1 = 0; il1 < ((*tpl5).second)->getLigandsTotal(); il1++ )
    if ( ((*tpl5).second)->getPocketNumber(il1) == clu3 )
     pocket_tmp->addTemplate((*tpl5).second);
  
  if ( pocket_tmp->getProteinsTotal() > 0 && pocket_tmp->getLigandsTotal() > 0 )
   pocket_set.push_back( pocket_tmp );
  else
   delete pocket_tmp;
 }
 
 delete [] clu1;
 
 cout << "Loading SVM models " << flush;
 
 ModelSVM * model_svm;
 
 model_svm = new ModelSVM( false, false, false, false, false, false, false, false, false );
 
 time(&t_bench1);
 
 model_svm->loadModel( 1, model_path+"/compositionSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 2, model_path+"/profileSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 3, model_path+"/residueSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 4, model_path+"/ligandSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 5, model_path+"/rank7SVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 6, model_path+"/rank8SVM.model" ); cout << '.' << flush;
 
 model_svm->loadScale( 1, model_path+"/compositionSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 2, model_path+"/profileSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 3, model_path+"/residueSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 4, model_path+"/ligandSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 5, model_path+"/rank7SVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 6, model_path+"/rank8SVM.scale" ); cout << '.' << flush;
 
 time(&t_bench2);
 
 if ( target->compositionSVM(model_svm) )
 {
  cout << "SVM for aa composition failed" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Predicting binding residues ... " << flush;
 
 time(&t_bench1);
 
 list< Pocket * > pocket_set_filtered;
 
 list< Pocket * >::iterator ipkt1;
 
 for ( ipkt1 = pocket_set.begin(); ipkt1 != pocket_set.end(); ipkt1++ )
 {
  (*ipkt1)->calculatePocketCenter();
  
  int bres1 = (*ipkt1)->calculateBindingResidues( target, model_svm, cut_binrest );
  
  if ( bres1 >= cut_binresn )
   pocket_set_filtered.push_back( *ipkt1 );
 }
 
 pocket_set.clear();
 
 if ( pocket_set_filtered.empty() )
 {
  cout << "no pockets survived" << endl << endl;
  
  template_set.clear();
  
  pocket_set_filtered.clear();
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << pocket_set_filtered.size() << " pockets survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 double fra1 = 0.0;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  fra1 += (double) (*ipkt1)->getLigandsTotal();
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double fra2 = 0.0;
  
  if ( fra1 > 0.0 )
   fra2 = ( (double) (*ipkt1)->getLigandsTotal() ) / fra1;
  
  (*ipkt1)->setPocketFraction( fra2 );
 }
 
 cout << "Constructing ligand fingerprints ... " << flush;
 
 time(&t_bench1);
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  (*ipkt1)->calculateFingerprintsSMILES( cut_clustlig, met_clustlig );
  (*ipkt1)->calculateFingerprintsMACCS( cut_clustlig, met_clustlig );
  
  if ( cmps_opt )
   (*ipkt1)->calculateCmpsScores( compounds, model_svm );
 }
 
 time(&t_bench2);
 
 cout << "done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Ranking pockets ... " << flush;
 
 time(&t_bench1);
 
 double rank1 = 0.0;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double rank2 = (*ipkt1)->calculateConfidence( cmps_opt, model_svm );
  
  if ( rank2 > rank1 )
   rank1 = rank2;
 }
 
 time(&t_bench2);
 
 cout << "top-ranked pocket has a confidence index of " << fixed << setprecision(1) << rank1 * 100 << "% (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 multimap<double,Pocket *> pocket_map_sorted;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  pocket_map_sorted.insert( pair<double,Pocket *>(-1.0*(*ipkt1)->getConfidence(),*ipkt1) );
 
 list< Pocket * > pocket_set_sorted;
 
 multimap<double,Pocket *>::iterator ipkt3;
 
 for ( ipkt3 = pocket_map_sorted.begin() ; ipkt3 != pocket_map_sorted.end(); ipkt3++ )
  pocket_set_sorted.push_back( (*ipkt3).second );
 
 pocket_set_filtered.clear();
 
 map<string,bool> chk1;
 map<string,bool> chk2;
 
 int ipkt2 = 1;
 
 for ( ipkt1 = pocket_set_sorted.begin(); ipkt1 != pocket_set_sorted.end(); ipkt1++ )
 {
  (*ipkt1)->setCenter( cut_binrest, cut_clustdis );
  
  (*ipkt1)->dumpProteinsAlignments( output_name, chk1, target );
  
  (*ipkt1)->dumpPocket( output_name, target, cut_binrest, ipkt2 );
  
  (*ipkt1)->dumpLigands( output_name, chk2, ipkt2 );
  
  ipkt2++;
 }
 
 template_set.clear();
 
 pocket_set_sorted.clear();
 printf("ENDING PARTC\n");
 return 0;
}
