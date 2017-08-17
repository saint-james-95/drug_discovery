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


#include "list.h"

using namespace std;

void getList( string templates_name, string cluster_def, list<string> &template_list, map<string,double> &template_prob1, map<string,double> &template_prob2 )
{
    multimap<string,string> cluster_map;

    /* these will hold probabilities of target protein ligand bounding
     * key is the template structure type
     * val is the probability
     */
    multimap<string,double> template_tmp1;
    multimap<string,double> template_tmp2;

    string line1;

    ifstream clusters_file( cluster_def.c_str() );

    if ( !clusters_file.is_open() )  { cout << "Cannot open " << cluster_def << endl; exit(EXIT_FAILURE); }

    /*
     * store contents of the map file (e.g. ethread-lib-jan2012.map)
     * into memory as cluster_map
     */
    while (getline(clusters_file,line1))
    {
        string tpl1;
        int    tpl3;

        tpl1 = line1.substr(0,15);

        string::size_type tpl6 = 0;

        bool tpl7 = true;

        /* remove whitespace */
        while( tpl7 )
        {
            tpl6 = tpl1.find(" ");
            if( tpl6 != string::npos )
                tpl1.erase( tpl6, 1 );
            else
                tpl7 = false;
        }

        tpl3 = atoi(line1.substr(15,5).c_str());

        /* for each protein to template insert into cluster_map i
         * e.g.: 
         *       // from ethread-lib-jan2012.map
         *       >          19hcA    2 1OFWB 1DUWA
         *       <19hcA, 10FWB> and <19hcA, 1DUWA> would be inserted into cluster_map
         */
        for ( int tpl4 = 0; tpl4 < tpl3; tpl4++ )
        {
            string tpl5;

            tpl5 = line1.substr(21+6*tpl4,5);

            cluster_map.insert( pair<string,string>(tpl1,tpl5) );
        }
    }

    clusters_file.close();

    /* 
     * a multimap<int,string> that orders by key value, largest to smallest
     * int will be the decimal part of a probability in etherad-fun
     * string will be the protein family name
     */
    multimap<int,string,std::greater<int> > template_fun;

    string line2;

    ifstream template_file( templates_name.c_str() );

    if ( !template_file.is_open() )  { cout << "Cannot open " << templates_name << endl; exit(EXIT_FAILURE); }

    /*
     * for each line in ethread-fun file
     */
    while (getline(template_file,line2))
        if ( line2.length() > 71 )
        {
            /* ignore lines beginning with # */
            if ( line2.compare(0, 1, "#") != 0 )
            {
                string tpl1;
                int    tpl2;
                double tpl3;
                double tpl4;

                /* template protein family (found to be similar to target unknown protein) */
                tpl1 = line2.substr(0,15);

                string::size_type tpl6 = 0;

                bool tpl7 = true;

                /* remove spaces */
                while( tpl7 )
                {
                    tpl6 = tpl1.find(" ");
                    if( tpl6 != string::npos )
                        tpl1.erase( tpl6, 1 );
                    else
                        tpl7 = false;
                }

                /* 
                 * tpl2 is the decimal part under ligand 
                 * according to ethread authors
                 * "a probability that the template binds ligands within 4 angstroms"
                 */ 
                tpl2 = atoi(line2.substr(57,6).c_str());

                /*
                 * same as tpl2 but includes 0 and decimal point
                 */
                tpl3 = atof(line2.substr(55,8).c_str());
                /*
                 * tpl4 is another probability
                 * according to ehtread authors
                 *
                 * "a probability that the Tanimoto coefficient between 
                 *  the template- and target-bound ligands is ≥0.5"
                 *
                 */
                tpl4 = atof(line2.substr(64,8).c_str());

                multimap<string,string>::iterator it1;
                pair<multimap<string,string>::iterator,multimap<string,string>::iterator> it2;

                /* lookup tpl1 aka the template protein family in cluster_map */
                it2 = cluster_map.equal_range(tpl1);

                /* for each protein and ligand of cluster_map that matches a protein 
                 * in ethread result
                 */
                for ( it1 = it2.first; it1 != it2.second; it1++ )
                {
                    /*
                     * insert to template_fun the probability of ligand bounding
                     * mapped to the protein family
                     */
                    template_fun.insert( pair<int,string>(tpl2,(*it1).second) );

                    /* store probabilities for ligand binding (tmp1) and 
                     * Tanimoto threshold surpassed (tmp2)
                     */
                    template_tmp1.insert( pair<string,double>((*it1).second, tpl3) );
                    template_tmp2.insert( pair<string,double>((*it1).second, tpl4) );
                }
            }
        }
        else if ( line2.length() ) // if line is shorter than 71 characters
        {                          // probably some different format unused
            if ( line2.compare(0, 1, "#") != 0 )
            {
                string tpl1;
                int    tpl2;
                double tpl3;
                double tpl4;

                tpl1 = line2;

                string::size_type tpl6 = 0;

                bool tpl7 = true;

                while( tpl7 )
                {
                    tpl6 = tpl1.find(" ");
                    if( tpl6 != string::npos )
                        tpl1.erase( tpl6, 1 );
                    else
                        tpl7 = false;
                }

                tpl2 = rand();

                tpl3 = 1.0;
                tpl4 = 1.0;

                multimap<string,string>::iterator it1;
                pair<multimap<string,string>::iterator,multimap<string,string>::iterator> it2;

                it2 = cluster_map.equal_range(tpl1);

                for ( it1 = it2.first; it1 != it2.second; it1++ )
                {
                    template_fun.insert( pair<int,string>(tpl2,(*it1).second) );

                    template_tmp1.insert( pair<string,double>((*it1).second, tpl3) );
                    template_tmp2.insert( pair<string,double>((*it1).second, tpl4) );
                }
            }
        }

    template_file.close();

    multimap<int,string>::iterator it3;

    /* 
     * this for loop creates a list of protein families and their relevant probabilities.
     * the loop makes sure that there are no duplicate entries of a protein family
     *
     * each entry is a line from the ethread-fun file
     *
     * note that template_fun is ordered based on greatest to least
     */
    for ( it3 = template_fun.begin() ; it3 != template_fun.end(); it3++ )
    {
        bool w1 = false; // signals no match found

        /* searches through template_list for a match with (*it3).second 
         * note (*it3).second would be the seq. of the highest probability 
         */
        if ( template_list.size() )
        {
            list<string>::iterator it4;

            it4 = find(template_list.begin(), template_list.end(), (*it3).second);

            if ( it4 == template_list.end() )
                w1 = true;
        }
        else
            w1 = true;

        /* if no match found in template_list */
        if ( w1 )
        {
            double prob1 = 0.0;
            double prob2 = 0.0;

            multimap<string,double>::iterator it5;
            pair<multimap<string,double>::iterator,multimap<string,double>::iterator> it6;

            /* 
             * returns range of indexes of all that match the seq 
             * note that this works because multipmap is ordered by seq name.
             */
            it6 = template_tmp1.equal_range( (*it3).second );

            /* for all matches in ligand binding probabilities
             * if its probability is larger than current prob1
             * then store it into prob1
             */
            for ( it5 = it6.first; it5 != it6.second; it5++ )
                if ( (*it5).second > prob1 )
                    prob1 = (*it5).second;

            /* repeate above steps for the tanimoto >= .5 probabilities
             */
            it6 = template_tmp2.equal_range( (*it3).second );

            for ( it5 = it6.first; it5 != it6.second; it5++ )
                if ( (*it5).second > prob2 )
                    prob2 = (*it5).second;

            template_prob1.insert( pair<string,double>((*it3).second, prob1) );
            template_prob2.insert( pair<string,double>((*it3).second, prob2) );

            template_list.push_back( (*it3).second );
        }
    }
}
